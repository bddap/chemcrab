use std::collections::HashSet;

use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::mol::Mol;
use crate::rings::RingInfo;
use crate::traits::{HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    Huckel,
}

const SP2_CAPABLE: [u8; 9] = [
    5,  // B
    6,  // C
    7,  // N
    8,  // O
    15, // P
    16, // S
    33, // As
    34, // Se
    52, // Te
];

pub fn find_aromatic_atoms<A, B>(mol: &Mol<A, B>) -> Vec<bool>
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    find_aromatic_atoms_huckel(mol)
}

pub fn set_aromaticity(
    mol: &mut Mol<crate::atom::Atom, crate::bond::Bond>,
    model: AromaticityModel,
) {
    match model {
        AromaticityModel::Huckel => {
            let aromatic = find_aromatic_atoms_huckel(mol);
            let indices: Vec<_> = mol.atoms().collect();
            for idx in indices {
                mol.atom_mut(idx).is_aromatic = aromatic[idx.index()];
            }
        }
    }
}

fn find_aromatic_atoms_huckel<A, B>(mol: &Mol<A, B>) -> Vec<bool>
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    let n = mol.atom_count();
    let mut aromatic = vec![false; n];

    let ring_info = RingInfo::sssr(mol);
    let rings = ring_info.rings();

    for ring in rings {
        if is_aromatic_ring(mol, ring) {
            for &atom_idx in ring {
                aromatic[atom_idx.index()] = true;
            }
        }
    }

    for system in fused_ring_systems(rings) {
        if system.len() < 2 {
            continue;
        }
        mark_fused_system_aromatic(mol, &system, rings, &mut aromatic);
    }

    aromatic
}

fn fused_ring_systems(rings: &[Vec<NodeIndex>]) -> Vec<Vec<usize>> {
    let n = rings.len();
    let mut adj = vec![vec![false; n]; n];
    for i in 0..n {
        let set_i: HashSet<NodeIndex> = rings[i].iter().copied().collect();
        for j in (i + 1)..n {
            let shared = rings[j].iter().filter(|a| set_i.contains(a)).count();
            if shared >= 2 {
                adj[i][j] = true;
                adj[j][i] = true;
            }
        }
    }

    let mut visited = vec![false; n];
    let mut components = Vec::new();
    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![i];
        while let Some(cur) = stack.pop() {
            if visited[cur] {
                continue;
            }
            visited[cur] = true;
            component.push(cur);
            for j in 0..n {
                if adj[cur][j] && !visited[j] {
                    stack.push(j);
                }
            }
        }
        components.push(component);
    }
    components
}

fn mark_fused_system_aromatic<A, B>(
    mol: &Mol<A, B>,
    system: &[usize],
    rings: &[Vec<NodeIndex>],
    aromatic: &mut [bool],
) where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    loop {
        let mut changed = false;
        for &ring_idx in system {
            let ring = &rings[ring_idx];
            if ring.iter().all(|&a| aromatic[a.index()]) {
                continue;
            }
            if is_aromatic_ring_in_fused_system(mol, ring, aromatic) {
                for &atom_idx in ring {
                    if !aromatic[atom_idx.index()] {
                        aromatic[atom_idx.index()] = true;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
}

fn is_aromatic_ring_in_fused_system<A, B>(
    mol: &Mol<A, B>,
    ring: &[NodeIndex],
    aromatic: &[bool],
) -> bool
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    if ring.len() < 3 {
        return false;
    }

    for &atom_idx in ring {
        if !SP2_CAPABLE.contains(&mol.atom(atom_idx).atomic_num()) {
            return false;
        }
    }

    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some(edge) = mol.bond_between(a, b) {
            if mol.bond(edge).bond_order() == BondOrder::Triple {
                return false;
            }
        }
    }

    ring.iter().all(|&atom_idx| is_sp2_in_fused_system(mol, atom_idx, aromatic))
}

fn is_sp2_in_fused_system<A, B>(
    mol: &Mol<A, B>,
    atom_idx: NodeIndex,
    aromatic: &[bool],
) -> bool
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    if aromatic[atom_idx.index()] {
        return true;
    }

    let atom = mol.atom(atom_idx);
    let anum = atom.atomic_num();
    let charge = atom.formal_charge();
    let has_double = has_any_double_bond(mol, atom_idx);
    let total_degree = mol.neighbors(atom_idx).count() as u8 + atom.hydrogen_count();

    match anum {
        6 => match charge {
            0 => has_double,
            1 => true,
            -1 => true,
            _ => false,
        },
        7 => match charge {
            0 => has_double || total_degree <= 3,
            1 => has_double,
            _ => false,
        },
        8 | 16 | 34 | 52 => true,
        5 => has_double,
        15 | 33 => has_double || total_degree <= 3,
        _ => false,
    }
}

fn is_aromatic_ring<A, B>(mol: &Mol<A, B>, ring: &[NodeIndex]) -> bool
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    if ring.len() < 3 {
        return false;
    }

    for &atom_idx in ring {
        if !SP2_CAPABLE.contains(&mol.atom(atom_idx).atomic_num()) {
            return false;
        }
    }

    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some(edge) = mol.bond_between(a, b) {
            if mol.bond(edge).bond_order() == BondOrder::Triple {
                return false;
            }
        }
    }

    let mut pi_total: u8 = 0;
    for (i, &atom_idx) in ring.iter().enumerate() {
        match pi_electrons(mol, atom_idx, ring, i) {
            Some(e) => pi_total = pi_total.saturating_add(e),
            None => return false,
        }
    }

    is_huckel(pi_total)
}

fn pi_electrons<A, B>(
    mol: &Mol<A, B>,
    atom_idx: NodeIndex,
    ring: &[NodeIndex],
    pos_in_ring: usize,
) -> Option<u8>
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    let atom = mol.atom(atom_idx);
    let anum = atom.atomic_num();
    let charge = atom.formal_charge();

    let has_double = has_any_double_bond(mol, atom_idx);
    let has_double_in_ring = has_double_to_ring_neighbor(mol, atom_idx, ring, pos_in_ring);

    let total_degree = mol.neighbors(atom_idx).count() as u8 + atom.hydrogen_count();
    let ring_degree = ring_neighbor_count(ring);

    match anum {
        6 => match charge {
            0 => {
                if has_double_in_ring {
                    Some(1)
                } else {
                    None
                }
            }
            -1 => Some(2),
            1 => {
                if has_double {
                    Some(1)
                } else {
                    Some(0)
                }
            }
            _ => None,
        },
        7 => match charge {
            0 => {
                if has_double {
                    Some(1)
                } else if ring_degree == 2 && total_degree <= 3 {
                    Some(2)
                } else {
                    None
                }
            }
            1 => {
                if has_double_in_ring {
                    Some(1)
                } else {
                    None
                }
            }
            _ => None,
        },
        8 | 16 | 34 | 52 => {
            if has_double_in_ring {
                Some(1)
            } else if ring_degree == 2 {
                Some(2)
            } else {
                None
            }
        }
        5 => {
            if has_double {
                Some(1)
            } else {
                None
            }
        }
        15 | 33 => {
            if has_double {
                Some(1)
            } else if ring_degree == 2 && total_degree <= 3 {
                Some(2)
            } else {
                None
            }
        }
        _ => None,
    }
}

fn has_any_double_bond<A, B>(mol: &Mol<A, B>, atom_idx: NodeIndex) -> bool
where
    B: HasBondOrder,
{
    mol.bonds_of(atom_idx)
        .any(|e| mol.bond(e).bond_order() == BondOrder::Double)
}

fn has_double_to_ring_neighbor<A, B>(
    mol: &Mol<A, B>,
    atom_idx: NodeIndex,
    ring: &[NodeIndex],
    pos_in_ring: usize,
) -> bool
where
    B: HasBondOrder,
{
    let len = ring.len();
    let prev = ring[(pos_in_ring + len - 1) % len];
    let next = ring[(pos_in_ring + 1) % len];

    for neighbor in [prev, next] {
        if let Some(edge) = mol.bond_between(atom_idx, neighbor) {
            if mol.bond(edge).bond_order() == BondOrder::Double {
                return true;
            }
        }
    }
    false
}

fn ring_neighbor_count(ring: &[NodeIndex]) -> u8 {
    if ring.len() > 1 { 2 } else { 0 }
}

fn is_huckel(pi_electrons: u8) -> bool {
    if pi_electrons < 2 {
        return false;
    }
    (pi_electrons - 2).is_multiple_of(4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};
    use crate::smiles::from_smiles;

    #[test]
    fn benzene_all_aromatic() {
        let mol = from_smiles("c1ccccc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 6);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn cyclohexane_none_aromatic() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert!(arom.iter().all(|&a| !a));
    }

    #[test]
    fn pyridine_all_aromatic() {
        let mol = from_smiles("c1ccncc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 6);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn pyrrole_all_aromatic() {
        let mol = from_smiles("[nH]1cccc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 5);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn furan_all_aromatic() {
        let mol = from_smiles("o1cccc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 5);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn thiophene_all_aromatic() {
        let mol = from_smiles("s1cccc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 5);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn naphthalene_all_aromatic() {
        let mol = from_smiles("c1ccc2ccccc2c1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 10);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn cyclopentadienyl_anion_aromatic() {
        let mol = from_smiles("[C-]1=CC=CC=1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 5);
        assert!(arom.iter().all(|&a| a));
    }

    #[test]
    fn phenol_ring_aromatic_oxygen_not() {
        let mol = from_smiles("Oc1ccccc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 7);
        assert!(!arom[0]);
        for (i, &is_arom) in arom.iter().enumerate().take(7).skip(1) {
            assert!(is_arom, "ring atom {} should be aromatic", i);
        }
    }

    #[test]
    fn cyclopentadiene_not_aromatic() {
        let mol = from_smiles("C1=CCC=C1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert!(arom.iter().all(|&a| !a));
    }

    #[test]
    fn cyclooctatetraene_not_aromatic() {
        let mol = from_smiles("C1=CC=CC=CC=C1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert!(arom.iter().all(|&a| !a));
    }

    #[test]
    fn set_aromaticity_on_kekulized_benzene() {
        let mut mol = Mol::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| {
                mol.add_atom(Atom {
                    atomic_num: 6,
                    hydrogen_count: 1,
                    ..Atom::default()
                })
            })
            .collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            mol.add_bond(
                atoms[i],
                atoms[(i + 1) % 6],
                Bond { order },
            );
        }

        for idx in mol.atoms() {
            assert!(!mol.atom(idx).is_aromatic);
        }

        set_aromaticity(&mut mol, AromaticityModel::Huckel);

        for idx in mol.atoms() {
            assert!(
                mol.atom(idx).is_aromatic,
                "atom {} should be aromatic after set_aromaticity",
                idx.index()
            );
        }
    }

    #[test]
    fn caffeine_fused_system_aromatic() {
        let mol = from_smiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O").unwrap();
        let arom = find_aromatic_atoms(&mol);
        let aromatic_count = arom.iter().filter(|&&a| a).count();
        assert_eq!(
            aromatic_count, 9,
            "caffeine purine-like fused system should have 9 aromatic atoms, got {}",
            aromatic_count
        );
    }

    #[test]
    fn boronic_ester_phenyl_rings_only() {
        let mol = from_smiles("c1ccc(B2Oc3ccccc3O2)cc1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        let aromatic_count = arom.iter().filter(|&&a| a).count();
        assert_eq!(
            aromatic_count, 12,
            "boronic ester should have 12 aromatic atoms (two phenyl rings), got {}",
            aromatic_count
        );
    }


    #[test]
    fn pyrene_all_aromatic() {
        let mol = from_smiles("c1cc2ccc3cccc4ccc(c1)c2c34").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 16);
        assert!(
            arom.iter().all(|&a| a),
            "all 16 pyrene atoms should be aromatic"
        );
    }

    #[test]
    fn anthracene_all_aromatic() {
        let mol = from_smiles("c1ccc2cc3ccccc3cc2c1").unwrap();
        let arom = find_aromatic_atoms(&mol);
        assert_eq!(arom.len(), 14);
        assert!(
            arom.iter().all(|&a| a),
            "all 14 anthracene atoms should be aromatic"
        );
    }

    #[test]
    fn pyrene_canonical_round_trip() {
        let smiles = "c1cc2ccc3cccc4ccc(c1)c2c34";
        let mol = from_smiles(smiles).unwrap();
        let canonical = crate::to_canonical_smiles(&mol);
        let mol2 = from_smiles(&canonical).unwrap();
        let canonical2 = crate::to_canonical_smiles(&mol2);
        assert_eq!(canonical, canonical2, "pyrene canonical SMILES not idempotent");
    }

    #[test]
    fn huckel_rule() {
        assert!(!is_huckel(0));
        assert!(!is_huckel(1));
        assert!(is_huckel(2));
        assert!(!is_huckel(4));
        assert!(is_huckel(6));
        assert!(!is_huckel(8));
        assert!(is_huckel(10));
        assert!(is_huckel(14));
        assert!(is_huckel(18));
    }
}
