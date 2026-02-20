use std::hash::{Hash, Hasher};

use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::mol::{AtomId, Mol};
use crate::traits::{
    HasAromaticity, HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount,
    HasIsotope,
};

struct Fnv1aHasher(u64);

impl Fnv1aHasher {
    fn new() -> Self {
        Self(0xcbf29ce484222325)
    }
}

impl Hasher for Fnv1aHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= b as u64;
            self.0 = self.0.wrapping_mul(0x100000001b3);
        }
    }
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct AtomInvariant {
    atomic_num: u8,
    degree: u8,
    hydrogen_count: u8,
    formal_charge: i8,
    is_aromatic: bool,
    isotope: u16,
    singles: u8,
    doubles: u8,
    triples: u8,
    aromatic_bonds: u8,
}

fn atom_invariant<A, B>(mol: &Mol<A, B>, idx: NodeIndex) -> AtomInvariant
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + HasAromaticity + HasIsotope,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    let is_aromatic = atom.is_aromatic();
    let degree = mol.neighbors(idx).count() as u8;
    let mut singles: u8 = 0;
    let mut doubles: u8 = 0;
    let mut triples: u8 = 0;
    let mut aromatic_bonds: u8 = 0;
    for edge in mol.bonds_of(idx) {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        let neighbor = if a == idx { b } else { a };
        if is_aromatic && mol.atom(neighbor).is_aromatic() {
            aromatic_bonds += 1;
        } else {
            match mol.bond(edge).bond_order() {
                BondOrder::Single => singles += 1,
                BondOrder::Double => doubles += 1,
                BondOrder::Triple => triples += 1,
            }
        }
    }
    AtomInvariant {
        atomic_num: atom.atomic_num(),
        degree,
        hydrogen_count: atom.hydrogen_count(),
        formal_charge: atom.formal_charge(),
        is_aromatic,
        isotope: atom.isotope(),
        singles,
        doubles,
        triples,
        aromatic_bonds,
    }
}

fn hash_invariant(inv: &AtomInvariant) -> u64 {
    let mut h = Fnv1aHasher::new();
    inv.hash(&mut h);
    h.finish()
}

fn ranks_from_values(values: &[u64]) -> Vec<usize> {
    let n = values.len();
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by_key(|&i| values[i]);
    let mut ranks = vec![0usize; n];
    if n == 0 {
        return ranks;
    }
    ranks[indices[0]] = 0;
    for i in 1..n {
        ranks[indices[i]] = if values[indices[i]] == values[indices[i - 1]] {
            ranks[indices[i - 1]]
        } else {
            i
        };
    }
    ranks
}

fn count_distinct(ranks: &[usize]) -> usize {
    let mut sorted: Vec<usize> = ranks.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    sorted.len()
}

fn morgan_refine<A, B>(mol: &Mol<A, B>, ranks: &mut Vec<usize>)
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + HasAromaticity + HasIsotope,
    B: HasBondOrder,
{
    let n = mol.atom_count();
    let mut prev_distinct = count_distinct(ranks);

    loop {
        let mut new_values = vec![0u64; n];
        for node in mol.atoms() {
            let i = node.index();
            let mut neighbor_ranks: Vec<usize> =
                mol.neighbors(node).map(|nb| ranks[nb.index()]).collect();
            neighbor_ranks.sort_unstable();

            let mut h = Fnv1aHasher::new();
            ranks[i].hash(&mut h);
            neighbor_ranks.hash(&mut h);
            new_values[i] = h.finish();
        }
        let new_ranks = ranks_from_values(&new_values);
        let distinct = count_distinct(&new_ranks);
        if distinct <= prev_distinct {
            return;
        }
        *ranks = new_ranks;
        prev_distinct = distinct;
    }
}

fn chirality_refine<A, B>(mol: &Mol<A, B>, ranks: &mut Vec<usize>)
where
    A: HasHydrogenCount,
{
    let n = ranks.len();
    let mut values = vec![0u64; n];
    let mut any_stereo = false;

    for stereo in mol.tetrahedral_stereo() {
        let center = stereo.center;
        if center.index() >= n {
            continue;
        }

        let full = stereo.above;

        let rank_of = |id: &AtomId| -> usize {
            match id {
                AtomId::Node(idx) => ranks[idx.index()],
                AtomId::VirtualH(_, _) => n,
            }
        };

        let ranks4 = [rank_of(&full[0]), rank_of(&full[1]), rank_of(&full[2]), rank_of(&full[3])];

        // Need all 4 neighbor ranks distinct for unambiguous parity.
        let mut sorted4 = ranks4;
        sorted4.sort_unstable();
        if sorted4.windows(2).any(|w| w[0] == w[1]) {
            continue;
        }

        let perm: [usize; 4] = std::array::from_fn(|i| {
            sorted4.iter().position(|&s| s == ranks4[i]).unwrap()
        });
        let mut visited = [false; 4];
        let mut swaps = 0usize;
        for i in 0..4 {
            if visited[i] {
                continue;
            }
            let mut cycle_len = 0;
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
                cycle_len += 1;
            }
            swaps += cycle_len - 1;
        }

        let parity: u64 = if swaps.is_multiple_of(2) { 1 } else { 2 };

        let mut h = Fnv1aHasher::new();
        ranks[center.index()].hash(&mut h);
        parity.hash(&mut h);
        values[center.index()] = h.finish();
        any_stereo = true;
    }

    if !any_stereo {
        return;
    }

    for i in 0..n {
        if values[i] == 0 {
            let mut h = Fnv1aHasher::new();
            ranks[i].hash(&mut h);
            values[i] = h.finish();
        }
    }

    *ranks = ranks_from_values(&values);
}

pub fn canonical_ordering<A, B>(mol: &Mol<A, B>) -> Vec<usize>
where
    A: HasAtomicNum
        + HasHydrogenCount
        + HasFormalCharge
        + HasAromaticity
        + HasIsotope,
    B: HasBondOrder,
{
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    let invariants: Vec<AtomInvariant> = (0..n)
        .map(|i| atom_invariant(mol, NodeIndex::new(i)))
        .collect();

    let initial_values: Vec<u64> = invariants.iter().map(hash_invariant).collect();
    let mut ranks = ranks_from_values(&initial_values);

    morgan_refine(mol, &mut ranks);

    if count_distinct(&ranks) < n {
        break_ties(mol, &mut ranks, &invariants);
    }

    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by_key(|&i| ranks[i]);
    let mut final_ranks = vec![0usize; n];
    for (rank, &atom_idx) in indices.iter().enumerate() {
        final_ranks[atom_idx] = rank;
    }
    final_ranks
}

fn break_ties<A, B>(mol: &Mol<A, B>, ranks: &mut Vec<usize>, invariants: &[AtomInvariant])
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + HasAromaticity + HasIsotope,
    B: HasBondOrder,
{
    let n = ranks.len();

    loop {
        if count_distinct(ranks) == n {
            return;
        }

        let min_tied_rank = find_min_tied_rank(ranks);

        let tied_atoms: Vec<usize> = (0..n)
            .filter(|&i| ranks[i] == min_tied_rank)
            .collect();

        // Try promoting each tied atom and pick the one that yields the
        // lexicographically smallest invariant trace (atom invariants sorted
        // by resulting rank). This is independent of atom numbering.
        let max_rank = *ranks.iter().max().unwrap();
        let mut best_trace: Option<Vec<u64>> = None;
        let mut best_ranks: Option<Vec<usize>> = None;

        for &candidate in &tied_atoms {
            let mut trial = ranks.clone();
            trial[candidate] = max_rank + 1;
            morgan_refine(mol, &mut trial);
            loop {
                let prev_distinct = count_distinct(&trial);
                chirality_refine(mol, &mut trial);
                morgan_refine(mol, &mut trial);
                if count_distinct(&trial) <= prev_distinct {
                    break;
                }
            }
            // Build trace: sort atoms by their trial rank, then hash each
            // atom's invariant + neighbor-rank signature.
            let mut indexed: Vec<(usize, usize)> = trial.iter().copied().enumerate().collect();
            indexed.sort_by_key(|&(_, r)| r);
            let trace: Vec<u64> = indexed
                .iter()
                .map(|&(atom_i, _)| {
                    let mut h = Fnv1aHasher::new();
                    invariants[atom_i].hash(&mut h);
                    let mut nb_ranks: Vec<usize> = mol
                        .neighbors(NodeIndex::new(atom_i))
                        .map(|nb| trial[nb.index()])
                        .collect();
                    nb_ranks.sort_unstable();
                    nb_ranks.hash(&mut h);
                    h.finish()
                })
                .collect();
            if best_trace.as_ref().is_none_or(|best| trace < *best) {
                best_trace = Some(trace);
                best_ranks = Some(trial);
            }
        }

        *ranks = best_ranks.unwrap();
    }
}

fn find_min_tied_rank(ranks: &[usize]) -> usize {
    let mut counts = std::collections::HashMap::new();
    for &r in ranks {
        *counts.entry(r).or_insert(0usize) += 1;
    }
    *counts
        .iter()
        .filter(|(_, &count)| count > 1)
        .map(|(rank, _)| rank)
        .min()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    #[test]
    fn empty_mol() {
        let mol = Mol::<crate::atom::Atom, crate::bond::Bond>::new();
        assert!(canonical_ordering(&mol).is_empty());
    }

    #[test]
    fn single_atom() {
        let mol = from_smiles("C").unwrap();
        let ranks = canonical_ordering(&mol);
        assert_eq!(ranks.len(), 1);
        assert_eq!(ranks[0], 0);
    }

    #[test]
    fn ethanol_all_distinct() {
        let mol = from_smiles("CCO").unwrap();
        let ranks = canonical_ordering(&mol);
        assert_eq!(ranks.len(), 3);
        let mut sorted = ranks.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    #[test]
    fn benzene_total_ordering() {
        let mol = from_smiles("c1ccccc1").unwrap();
        let ranks = canonical_ordering(&mol);
        assert_eq!(ranks.len(), 6);
        let mut sorted = ranks.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4, 5]);
    }
}
