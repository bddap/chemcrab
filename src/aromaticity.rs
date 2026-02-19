use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::mol::Mol;
use crate::rings::RingInfo;
use crate::traits::{HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    Rdkit,
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
    find_aromatic_atoms_rdkit(mol)
}

pub fn set_aromaticity(
    mol: &mut Mol<crate::atom::Atom, crate::bond::Bond>,
    model: AromaticityModel,
) {
    match model {
        AromaticityModel::Rdkit => {
            let aromatic = find_aromatic_atoms_rdkit(mol);
            let indices: Vec<_> = mol.atoms().collect();
            for idx in indices {
                mol.atom_mut(idx).is_aromatic = aromatic[idx.index()];
            }
        }
    }
}

fn find_aromatic_atoms_rdkit<A, B>(mol: &Mol<A, B>) -> Vec<bool>
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    let n = mol.atom_count();
    let mut aromatic = vec![false; n];

    let ring_info = RingInfo::sssr(mol);

    for ring in ring_info.rings() {
        if is_aromatic_ring(mol, ring) {
            for &atom_idx in ring {
                aromatic[atom_idx.index()] = true;
            }
        }
    }

    aromatic
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
                if has_double {
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
                Some(0)
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
        for i in 1..7 {
            assert!(arom[i], "ring atom {} should be aromatic", i);
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
                Bond {
                    order,
                    ..Bond::default()
                },
            );
        }

        for idx in mol.atoms() {
            assert!(!mol.atom(idx).is_aromatic);
        }

        set_aromaticity(&mut mol, AromaticityModel::Rdkit);

        for idx in mol.atoms() {
            assert!(
                mol.atom(idx).is_aromatic,
                "atom {} should be aromatic after set_aromaticity",
                idx.index()
            );
        }
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
