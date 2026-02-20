use std::hash::{Hash, Hasher};

use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::mol::{AtomId, Mol};
use crate::traits::{
    HasAromaticity, HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount, HasIsotope,
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

fn above_rank_seq<A, B>(
    mol: &Mol<A, B>,
    stereo: &crate::mol::TetrahedralStereo,
    ranks: &[usize],
) -> [usize; 4]
where
    A: HasHydrogenCount,
{
    let n = mol.atom_count();
    std::array::from_fn(|i| match stereo.above[i] {
        AtomId::Node(idx) => ranks[idx.index()],
        AtomId::VirtualH(_, _) => n,
    })
}

fn hash_stereo_for_ranks(ranks4: &[usize; 4], center_rank: usize) -> u64 {
    let mut sorted4 = *ranks4;
    sorted4.sort_unstable();
    let has_ties = sorted4.windows(2).any(|w| w[0] == w[1]);

    let mut h = Fnv1aHasher::new();
    center_rank.hash(&mut h);

    if has_ties {
        sorted4.hash(&mut h);
        3u64.hash(&mut h);
    } else {
        let perm: [usize; 4] =
            std::array::from_fn(|i| sorted4.iter().position(|&s| s == ranks4[i]).unwrap());
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
        parity.hash(&mut h);
    }

    h.finish()
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

        let ranks4 = above_rank_seq(mol, stereo, ranks);
        values[center.index()] = hash_stereo_for_ranks(&ranks4, ranks[center.index()]);
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

fn ez_refine<A, B>(mol: &Mol<A, B>, ranks: &mut Vec<usize>)
where
    A: HasHydrogenCount,
{
    let n = ranks.len();
    let mut values = vec![0u64; n];
    let mut any_stereo = false;

    for ez in mol.ez_stereo() {
        let a = ez.bond.0;
        let b = ez.bond.1;
        if a.index() >= n || b.index() >= n {
            continue;
        }

        let highest_ref_a = mol
            .neighbors(a)
            .filter(|&nb| nb != b)
            .max_by_key(|nb| ranks[nb.index()]);
        let highest_ref_b = mol
            .neighbors(b)
            .filter(|&nb| nb != a)
            .max_by_key(|nb| ranks[nb.index()]);

        let (canon_ref_a, canon_ref_b) = match (highest_ref_a, highest_ref_b) {
            (Some(ra), Some(rb)) => (ra, rb),
            _ => continue,
        };

        let same_side_a = match &ez.refs[0] {
            AtomId::Node(idx) => *idx == canon_ref_a,
            AtomId::VirtualH(_, _) => false,
        };
        let same_side_b = match &ez.refs[1] {
            AtomId::Node(idx) => *idx == canon_ref_b,
            AtomId::VirtualH(_, _) => false,
        };

        let canon_refs_cis = same_side_a == same_side_b;
        let parity: u64 = if canon_refs_cis { 1 } else { 2 };

        for &atom in &[a, b] {
            let mut h = Fnv1aHasher::new();
            ranks[atom.index()].hash(&mut h);
            parity.hash(&mut h);
            values[atom.index()] = h.finish();
            any_stereo = true;
        }
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
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + HasAromaticity + HasIsotope,
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

    loop {
        let prev = count_distinct(&ranks);
        chirality_refine(mol, &mut ranks);
        morgan_refine(mol, &mut ranks);
        if count_distinct(&ranks) <= prev {
            break;
        }
    }

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

        let min_tied_rank = find_best_tied_rank(mol, ranks);

        let tied_atoms: Vec<usize> = (0..n).filter(|&i| ranks[i] == min_tied_rank).collect();

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
                ez_refine(mol, &mut trial);
                morgan_refine(mol, &mut trial);
                if count_distinct(&trial) <= prev_distinct {
                    break;
                }
            }
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
                    if let Some(stereo) = mol.tetrahedral_stereo_for(NodeIndex::new(atom_i)) {
                        let ranks4 = above_rank_seq(mol, stereo, &trial);
                        let stereo_hash = hash_stereo_for_ranks(&ranks4, trial[atom_i]);
                        stereo_hash.hash(&mut h);
                    }
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

fn find_best_tied_rank<A, B>(mol: &Mol<A, B>, ranks: &[usize]) -> usize
where
    A: HasHydrogenCount,
{
    let mut counts = std::collections::HashMap::new();
    for &r in ranks {
        *counts.entry(r).or_insert(0usize) += 1;
    }
    let tied_ranks: Vec<usize> = counts
        .into_iter()
        .filter(|&(_, count)| count > 1)
        .map(|(rank, _)| rank)
        .collect();

    let has_stereo_at_rank = |r: usize| -> bool {
        mol.tetrahedral_stereo()
            .iter()
            .any(|s| s.center.index() < ranks.len() && ranks[s.center.index()] == r)
    };

    let non_stereo: Option<usize> = tied_ranks
        .iter()
        .copied()
        .filter(|&r| !has_stereo_at_rank(r))
        .min();

    non_stereo.unwrap_or_else(|| *tied_ranks.iter().min().unwrap())
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

    #[test]
    fn canonical_135_cyclohexane_triol_idempotent() {
        let cases = [
            "[C@@H]1(O)C[C@H](O)C[C@@H](O)C1",
            "[C@H]1(O)C[C@@H](O)C[C@H](O)C1",
            "[C@@H]1(O)C[C@@H](O)C[C@@H](O)C1",
            "[C@H]1(O)C[C@H](O)C[C@H](O)C1",
            "[C@@H]1(F)C[C@H](F)C[C@@H](F)C1",
            "[C@@H]1(Cl)C[C@H](Cl)C[C@@H](Cl)C1",
        ];
        for smi in &cases {
            let mol = from_smiles(smi).unwrap();
            let c1 = crate::smiles::to_canonical_smiles(&mol);
            let reparsed = from_smiles(&c1).unwrap();
            let c2 = crate::smiles::to_canonical_smiles(&reparsed);
            assert_eq!(c1, c2, "round-trip failed for {smi}: '{c1}' vs '{c2}'");
        }
    }

    #[test]
    fn canonical_spiro_stereo_idempotent() {
        let cases = [
            "[C@@H]1(CC1)[C@H]2CC2",
            "[C@H]1(CC1)[C@H]2CC2",
            "[C@@H]1(CCC1)[C@H]2CCC2",
            "[C@@H]1(CC1)[C@H]2CCC2",
            "[C@@H]1(CCC1)[C@H]2CC2",
            "[C@@H]1(CC1)[C@@H](F)CC",
            "[C@@H]1(CCCC1)[C@H]2CCCC2",
        ];
        for smi in &cases {
            let mol = from_smiles(smi).unwrap();
            let c1 = crate::smiles::to_canonical_smiles(&mol);
            let reparsed = from_smiles(&c1).unwrap();
            let c2 = crate::smiles::to_canonical_smiles(&reparsed);
            assert_eq!(c1, c2, "round-trip failed for {smi}: '{c1}' vs '{c2}'");
        }
    }

    #[test]
    fn canonical_symmetric_chirality_permutation_invariant() {
        use crate::graph_ops::renumber_atoms;
        let cases = [
            "[C@@H]1(O)C[C@H](O)C[C@@H](O)C1",
            "[C@@H]1(CC1)[C@H]2CC2",
            "[C@@H]1(CCC1)[C@H]2CCC2",
        ];
        for smi in &cases {
            let mol = from_smiles(smi).unwrap();
            let c1 = crate::smiles::to_canonical_smiles(&mol);
            let n = mol.atom_count();
            for offset in 1..n {
                let perm: Vec<usize> = (0..n).map(|i| (i + offset) % n).collect();
                let renum = renumber_atoms(&mol, &perm).unwrap();
                let c2 = crate::smiles::to_canonical_smiles(&renum);
                assert_eq!(c1, c2, "permutation invariance failed for {smi} with offset {offset}: '{c1}' vs '{c2}'");
            }
        }
    }
}
