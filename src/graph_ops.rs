use std::collections::VecDeque;

use petgraph::graph::NodeIndex;

use crate::bond::BondStereo;
use crate::canonical::canonical_ordering;
use crate::mol::{AtomId, Mol};
use crate::traits::{
    HasAromaticity, HasAtomicNum, HasBondOrder, HasBondStereoMut,
    HasFormalCharge, HasHydrogenCount, HasIsotope,
};

pub fn adjacency_matrix<A, B>(mol: &Mol<A, B>) -> Vec<Vec<bool>> {
    let n = mol.atom_count();
    let mut matrix = vec![vec![false; n]; n];
    for edge in mol.bonds() {
        if let Some((a, b)) = mol.bond_endpoints(edge) {
            matrix[a.index()][b.index()] = true;
            matrix[b.index()][a.index()] = true;
        }
    }
    matrix
}

pub fn distance_matrix<A, B>(mol: &Mol<A, B>) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    let mut dist = vec![vec![usize::MAX; n]; n];
    for start in mol.atoms() {
        let si = start.index();
        dist[si][si] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(current) = queue.pop_front() {
            let d = dist[si][current.index()];
            for neighbor in mol.neighbors(current) {
                if dist[si][neighbor.index()] == usize::MAX {
                    dist[si][neighbor.index()] = d + 1;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    dist
}

pub fn shortest_path<A, B>(
    mol: &Mol<A, B>,
    from: NodeIndex,
    to: NodeIndex,
) -> Option<Vec<NodeIndex>> {
    if from == to {
        return Some(vec![from]);
    }
    let n = mol.atom_count();
    let mut pred = vec![None; n];
    let mut visited = vec![false; n];
    visited[from.index()] = true;
    let mut queue = VecDeque::new();
    queue.push_back(from);
    while let Some(current) = queue.pop_front() {
        for neighbor in mol.neighbors(current) {
            if !visited[neighbor.index()] {
                visited[neighbor.index()] = true;
                pred[neighbor.index()] = Some(current);
                if neighbor == to {
                    let mut path = vec![to];
                    let mut node = to;
                    while let Some(p) = pred[node.index()] {
                        path.push(p);
                        node = p;
                    }
                    path.reverse();
                    return Some(path);
                }
                queue.push_back(neighbor);
            }
        }
    }
    None
}

pub fn connected_components<A, B>(mol: &Mol<A, B>) -> Vec<Vec<NodeIndex>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components = Vec::new();
    for node in mol.atoms() {
        if visited[node.index()] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![node];
        while let Some(current) = stack.pop() {
            if visited[current.index()] {
                continue;
            }
            visited[current.index()] = true;
            component.push(current);
            for neighbor in mol.neighbors(current) {
                if !visited[neighbor.index()] {
                    stack.push(neighbor);
                }
            }
        }
        component.sort();
        components.push(component);
    }
    components
}

pub fn num_components<A, B>(mol: &Mol<A, B>) -> usize {
    connected_components(mol).len()
}

pub fn get_fragments<A: Clone, B: Clone>(mol: &Mol<A, B>) -> Vec<Mol<A, B>> {
    let components = connected_components(mol);
    let mut fragments = Vec::with_capacity(components.len());
    for component in &components {
        let mut frag = Mol::new();
        let mut index_map = vec![NodeIndex::new(0); mol.atom_count()];
        for &old_idx in component {
            let new_idx = frag.add_atom(mol.atom(old_idx).clone());
            index_map[old_idx.index()] = new_idx;
        }
        for &old_idx in component {
            for edge in mol.bonds_of(old_idx) {
                if let Some((a, b)) = mol.bond_endpoints(edge) {
                    if a == old_idx && a.index() < b.index() {
                        frag.add_bond(
                            index_map[a.index()],
                            index_map[b.index()],
                            mol.bond(edge).clone(),
                        );
                    }
                }
            }
        }
        fragments.push(frag);
    }
    fragments
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RenumberError {
    LengthMismatch { expected: usize, got: usize },
    InvalidPermutation,
}

impl std::fmt::Display for RenumberError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LengthMismatch { expected, got } => {
                write!(f, "new_order length {got} != atom count {expected}")
            }
            Self::InvalidPermutation => write!(f, "new_order is not a valid permutation"),
        }
    }
}

impl std::error::Error for RenumberError {}

fn validate_permutation(new_order: &[usize], n: usize) -> Result<(), RenumberError> {
    if new_order.len() != n {
        return Err(RenumberError::LengthMismatch {
            expected: n,
            got: new_order.len(),
        });
    }
    let mut seen = vec![false; n];
    for &idx in new_order {
        if idx >= n || seen[idx] {
            return Err(RenumberError::InvalidPermutation);
        }
        seen[idx] = true;
    }
    Ok(())
}

pub fn renumber_atoms<A: Clone, B: HasBondStereoMut + Clone>(
    mol: &Mol<A, B>,
    new_order: &[usize],
) -> Result<Mol<A, B>, RenumberError> {
    let n = mol.atom_count();
    validate_permutation(new_order, n)?;

    let mut new_mol = Mol::new();

    // new_order[new_idx] = old_idx
    for &old_idx in new_order {
        new_mol.add_atom(mol.atom(NodeIndex::new(old_idx)).clone());
    }

    // old_to_new[old_idx] = new_idx
    let mut old_to_new = vec![0usize; n];
    for (new_idx, &old_idx) in new_order.iter().enumerate() {
        old_to_new[old_idx] = new_idx;
    }

    for edge in mol.bonds() {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        let new_a = NodeIndex::new(old_to_new[a.index()]);
        let new_b = NodeIndex::new(old_to_new[b.index()]);
        new_mol.add_bond(new_a, new_b, mol.bond(edge).clone());
    }

    let old_to_new_nodes: Vec<NodeIndex> = old_to_new.iter().map(|&i| NodeIndex::new(i)).collect();
    remap_bond_stereo(&mut new_mol, &old_to_new_nodes);

    let remapped_stereo = mol
        .tetrahedral_stereo()
        .iter()
        .map(|s| s.map(|aid| remap_atom_id(aid, &old_to_new_nodes)))
        .collect();
    new_mol.set_tetrahedral_stereo(remapped_stereo);

    Ok(new_mol)
}

fn remap_atom_id(aid: AtomId, old_to_new: &[NodeIndex]) -> AtomId {
    match aid {
        AtomId::Node(idx) => AtomId::Node(old_to_new[idx.index()]),
        AtomId::VirtualH(parent, n) => AtomId::VirtualH(old_to_new[parent.index()], n),
    }
}

fn remap_bond_stereo<A, B>(mol: &mut Mol<A, B>, old_to_new: &[NodeIndex])
where
    B: HasBondStereoMut,
{
    for edge in mol.bonds().collect::<Vec<_>>() {
        let stereo = mol.bond(edge).bond_stereo();
        let remapped = match stereo {
            BondStereo::None => BondStereo::None,
            BondStereo::Cis(a, b) => BondStereo::Cis(old_to_new[a.index()], old_to_new[b.index()]),
            BondStereo::Trans(a, b) => {
                BondStereo::Trans(old_to_new[a.index()], old_to_new[b.index()])
            }
        };
        *mol.bond_mut(edge).bond_stereo_mut() = remapped;
    }
}

pub fn renumber_atoms_canonical<A, B>(mol: &Mol<A, B>) -> Mol<A, B>
where
    A: HasAtomicNum
        + HasFormalCharge
        + HasHydrogenCount
        + HasIsotope
        + HasAromaticity
        + Clone,
    B: HasBondOrder + HasBondStereoMut + Clone,
{
    let n = mol.atom_count();
    if n == 0 {
        return Mol::new();
    }
    let ranks = canonical_ordering(mol);
    let mut new_order = vec![0usize; n];
    for (old_idx, &rank) in ranks.iter().enumerate() {
        new_order[rank] = old_idx;
    }
    renumber_atoms(mol, &new_order).expect("canonical ordering is a valid permutation")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder, BondStereo};
    use crate::mol::AtomId;
    use crate::smiles::{from_smiles, to_canonical_smiles};

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    #[test]
    fn renumber_identity() {
        let mol = from_smiles("CCO").unwrap();
        let identity: Vec<usize> = (0..mol.atom_count()).collect();
        let renum = renumber_atoms(&mol, &identity).unwrap();
        assert_eq!(renum.atom_count(), mol.atom_count());
        assert_eq!(renum.bond_count(), mol.bond_count());
        for i in 0..mol.atom_count() {
            assert_eq!(
                renum.atom(n(i)).atomic_num,
                mol.atom(n(i)).atomic_num
            );
        }
    }

    #[test]
    fn renumber_reversed() {
        let mol = from_smiles("CCO").unwrap();
        let n_atoms = mol.atom_count();
        let reversed: Vec<usize> = (0..n_atoms).rev().collect();
        let renum = renumber_atoms(&mol, &reversed).unwrap();
        assert_eq!(renum.atom_count(), n_atoms);
        assert_eq!(renum.bond_count(), mol.bond_count());
        // new[0] should be old[2] (oxygen)
        assert_eq!(renum.atom(n(0)).atomic_num, 8);
        // new[2] should be old[0] (carbon)
        assert_eq!(renum.atom(n(2)).atomic_num, 6);
    }

    #[test]
    fn renumber_preserves_bond_connectivity() {
        let mol = from_smiles("C-C=O").unwrap();
        // reverse: new_order = [2, 1, 0] → O, C, C
        let renum = renumber_atoms(&mol, &[2, 1, 0]).unwrap();
        // old bond 0-1 (single) → new 2-1
        // old bond 1-2 (double) → new 1-0
        let adj = adjacency_matrix(&renum);
        assert!(adj[0][1]); // O=C
        assert!(adj[1][2]); // C-C
        assert!(!adj[0][2]);
    }

    #[test]
    fn renumber_preserves_bond_stereo() {
        // Build a molecule with E/Z stereo manually
        let mut mol = Mol::new();
        let c0 = mol.add_atom(Atom { atomic_num: 6, hydrogen_count: 1, ..Atom::default() });
        let c1 = mol.add_atom(Atom { atomic_num: 6, hydrogen_count: 1, ..Atom::default() });
        let f2 = mol.add_atom(Atom { atomic_num: 9, ..Atom::default() });
        let cl3 = mol.add_atom(Atom { atomic_num: 17, ..Atom::default() });
        mol.add_bond(c0, c1, Bond { order: BondOrder::Double, stereo: BondStereo::Trans(n(2), n(3)) });
        mol.add_bond(c0, f2, Bond { order: BondOrder::Single, stereo: BondStereo::None });
        mol.add_bond(c1, cl3, Bond { order: BondOrder::Single, stereo: BondStereo::None });

        // Renumber: reverse [3, 2, 1, 0]
        // old_to_new: 0→3, 1→2, 2→1, 3→0
        let new_order = vec![3, 2, 1, 0];
        let renum_basic = renumber_atoms(&mol, &new_order).unwrap();

        // The basic renumber clones bonds without remapping stereo refs
        // Use renumber_atoms_canonical which does the full remap
        let canonical = renumber_atoms_canonical(&mol);
        // Verify stereo still references valid atoms
        for edge in canonical.bonds() {
            let bond = canonical.bond(edge);
            match bond.stereo {
                BondStereo::Cis(a, b) | BondStereo::Trans(a, b) => {
                    assert!(a.index() < canonical.atom_count());
                    assert!(b.index() < canonical.atom_count());
                }
                BondStereo::None => {}
            }
        }
        let double_edge = renum_basic.bonds().find(|&e| renum_basic.bond(e).order == BondOrder::Double).unwrap();
        let stereo = renum_basic.bond(double_edge).stereo;
        // old_to_new: 0→3, 1→2, 2→1, 3→0. Trans(2,3) → Trans(1,0)
        assert_eq!(stereo, BondStereo::Trans(n(1), n(0)));
    }

    #[test]
    fn renumber_preserves_chirality() {
        let mol = from_smiles("[C@@](F)(Cl)(Br)I").unwrap();
        let smiles_orig = to_canonical_smiles(&mol);

        let canonical = renumber_atoms_canonical(&mol);
        let smiles_renum = to_canonical_smiles(&canonical);
        assert_eq!(smiles_orig, smiles_renum);
    }

    #[test]
    fn renumber_chirality_manual() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let f = mol.add_atom(Atom { atomic_num: 9, ..Atom::default() });
        let cl = mol.add_atom(Atom { atomic_num: 17, ..Atom::default() });
        let br = mol.add_atom(Atom { atomic_num: 35, ..Atom::default() });
        let iodine = mol.add_atom(Atom { atomic_num: 53, ..Atom::default() });
        mol.add_bond(c, f, Bond::default());
        mol.add_bond(c, cl, Bond::default());
        mol.add_bond(c, br, Bond::default());
        mol.add_bond(c, iodine, Bond::default());
        mol.add_tetrahedral_stereo([
            AtomId::Node(c),
            AtomId::Node(f),
            AtomId::Node(cl),
            AtomId::Node(br),
        ]);

        // Identity renumber preserves chirality
        let id_order: Vec<usize> = (0..5).collect();
        let renum = renumber_atoms(&mol, &id_order).unwrap();
        assert!(renum.tetrahedral_stereo_for(n(0)).is_some());

        // Swap two neighbors: indices get remapped through the permutation
        let swap_order = vec![0, 2, 1, 3, 4];
        let renum2 = renumber_atoms(&mol, &swap_order).unwrap();
        assert!(renum2.tetrahedral_stereo_for(n(0)).is_some());

        // Canonical SMILES should be identical regardless of renumbering
        let smiles_orig = to_canonical_smiles(&mol);
        let smiles_renum = to_canonical_smiles(&renum2);
        assert_eq!(smiles_orig, smiles_renum);
    }

    #[test]
    fn renumber_invalid_length() {
        let mol = from_smiles("CC").unwrap();
        let result = renumber_atoms(&mol, &[0, 1, 2]);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            RenumberError::LengthMismatch { expected: 2, got: 3 }
        ));
    }

    #[test]
    fn renumber_invalid_duplicate() {
        let mol = from_smiles("CC").unwrap();
        let result = renumber_atoms(&mol, &[0, 0]);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            RenumberError::InvalidPermutation
        ));
    }

    #[test]
    fn renumber_invalid_out_of_range() {
        let mol = from_smiles("CC").unwrap();
        let result = renumber_atoms(&mol, &[0, 5]);
        assert!(result.is_err());
    }

    #[test]
    fn canonical_renumber_deterministic() {
        let mol1 = from_smiles("OCC").unwrap();
        let mol2 = from_smiles("CCO").unwrap();
        let can1 = renumber_atoms_canonical(&mol1);
        let can2 = renumber_atoms_canonical(&mol2);
        assert_eq!(can1.atom_count(), can2.atom_count());
        assert_eq!(can1.bond_count(), can2.bond_count());
        for i in 0..can1.atom_count() {
            assert_eq!(can1.atom(n(i)).atomic_num, can2.atom(n(i)).atomic_num);
        }
    }

    #[test]
    fn canonical_renumber_empty() {
        let mol = Mol::<Atom, Bond>::new();
        let renum = renumber_atoms_canonical(&mol);
        assert_eq!(renum.atom_count(), 0);
    }

    #[test]
    fn renumber_empty_mol() {
        let mol = Mol::<Atom, Bond>::new();
        let renum = renumber_atoms(&mol, &[]).unwrap();
        assert_eq!(renum.atom_count(), 0);
    }

    #[test]
    fn adjacency_ethane() {
        let mol = from_smiles("CC").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 2);
        assert!(!adj[0][0]);
        assert!(adj[0][1]);
        assert!(adj[1][0]);
        assert!(!adj[1][1]);
    }

    #[test]
    fn adjacency_methane() {
        let mol = from_smiles("C").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 1);
        assert!(!adj[0][0]);
    }

    #[test]
    fn adjacency_triangle() {
        let mol = from_smiles("C1CC1").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 3);
        assert!(adj[0][1]);
        assert!(adj[1][2]);
        assert!(adj[0][2]);
        assert!(adj[1][0]);
        assert!(adj[2][1]);
        assert!(adj[2][0]);
        for (i, row) in adj.iter().enumerate().take(3) {
            assert!(!row[i]);
        }
    }

    #[test]
    fn distance_linear_chain() {
        let mol = from_smiles("CCCC").unwrap();
        let dist = distance_matrix(&mol);
        assert_eq!(dist[0][3], 3);
        assert_eq!(dist[0][1], 1);
        assert_eq!(dist[1][3], 2);
        assert_eq!(dist[0][0], 0);
    }

    #[test]
    fn distance_cyclohexane() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        let dist = distance_matrix(&mol);
        assert_eq!(dist.len(), 6);
        assert_eq!(dist[0][3], 3);
        assert_eq!(dist[0][1], 1);
        assert_eq!(dist[0][2], 2);
    }

    #[test]
    fn shortest_path_linear() {
        let mol = from_smiles("CCCCC").unwrap();
        let path = shortest_path(&mol, NodeIndex::new(0), NodeIndex::new(4));
        assert_eq!(
            path,
            Some(vec![
                NodeIndex::new(0),
                NodeIndex::new(1),
                NodeIndex::new(2),
                NodeIndex::new(3),
                NodeIndex::new(4),
            ])
        );
    }

    #[test]
    fn shortest_path_no_path() {
        let mol = from_smiles("[Na+].[Cl-]").unwrap();
        let path = shortest_path(&mol, NodeIndex::new(0), NodeIndex::new(1));
        assert_eq!(path, None);
    }

    #[test]
    fn components_nacl() {
        let mol = from_smiles("[Na+].[Cl-]").unwrap();
        let comps = connected_components(&mol);
        assert_eq!(comps.len(), 2);
    }

    #[test]
    fn components_single() {
        let mol = from_smiles("CCO").unwrap();
        assert_eq!(num_components(&mol), 1);
    }

    #[test]
    fn components_empty() {
        let mol: Mol<(), ()> = Mol::new();
        assert_eq!(num_components(&mol), 0);
    }

    #[test]
    fn fragments_three() {
        let mol = from_smiles("[Na+].[Cl-].O").unwrap();
        let frags = get_fragments(&mol);
        assert_eq!(frags.len(), 3);
        let mut counts: Vec<usize> = frags.iter().map(|f| f.atom_count()).collect();
        counts.sort();
        assert_eq!(counts, vec![1, 1, 1]);
    }

    #[test]
    fn fragments_single() {
        let mol = from_smiles("CCO").unwrap();
        let frags = get_fragments(&mol);
        assert_eq!(frags.len(), 1);
        assert_eq!(frags[0].atom_count(), mol.atom_count());
        assert_eq!(frags[0].bond_count(), mol.bond_count());
    }
}
