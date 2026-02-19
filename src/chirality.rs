use crate::atom::{Atom, Chirality};
use crate::mol::Mol;
use crate::traits::HasBondOrder;

pub fn cleanup_chirality<B>(mol: &mut Mol<Atom, B>)
where
    B: HasBondOrder,
{
    let indices: Vec<_> = mol.atoms().collect();
    for idx in indices {
        let atom = mol.atom(idx);
        if atom.chirality == Chirality::None {
            continue;
        }

        let total_neighbors =
            mol.neighbors(idx).count() as u8 + atom.hydrogen_count;

        if total_neighbors < 3 {
            mol.atom_mut(idx).chirality = Chirality::None;
            continue;
        }

        if total_neighbors < 4 && !has_lone_pair_stereo(atom.atomic_num) {
            mol.atom_mut(idx).chirality = Chirality::None;
            continue;
        }

        if all_neighbors_identical(mol, idx) {
            mol.atom_mut(idx).chirality = Chirality::None;
        }
    }
}

fn has_lone_pair_stereo(atomic_num: u8) -> bool {
    matches!(atomic_num, 7 | 15 | 33 | 16 | 34 | 52)
}

fn all_neighbors_identical<B: HasBondOrder>(mol: &Mol<Atom, B>, idx: petgraph::graph::NodeIndex) -> bool {
    let atom = mol.atom(idx);
    let mut neighbor_keys: Vec<(u8, u8)> = mol
        .neighbors(idx)
        .map(|nb| {
            let nb_atom = mol.atom(nb);
            let bond_edge = mol.bond_between(idx, nb).unwrap();
            let bo = match mol.bond(bond_edge).bond_order() {
                crate::bond::BondOrder::Single => 1u8,
                crate::bond::BondOrder::Double => 2,
                crate::bond::BondOrder::Triple => 3,
            };
            (nb_atom.atomic_num, bo)
        })
        .collect();

    for _ in 0..atom.hydrogen_count {
        neighbor_keys.push((1, 1));
    }

    if neighbor_keys.len() < 2 {
        return true;
    }

    let first = neighbor_keys[0];
    neighbor_keys.iter().all(|k| *k == first)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;
    use crate::mol::Mol;
    use crate::smiles::from_smiles;

    #[test]
    fn chiral_center_preserved() {
        let mut mol = from_smiles("F[C@H](Cl)Br").unwrap();
        cleanup_chirality(&mut mol);
        let c_idx = mol
            .atoms()
            .find(|&idx| mol.atom(idx).atomic_num == 6)
            .unwrap();
        assert_ne!(mol.atom(c_idx).chirality, Chirality::None);
    }

    #[test]
    fn chirality_removed_too_few_neighbors() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 1,
            ..Atom::default()
        });
        let f = mol.add_atom(Atom {
            atomic_num: 9,
            ..Atom::default()
        });
        mol.add_bond(c, f, Bond::default());
        // C has 1 explicit + 1 implicit H = 2 neighbors. < 3 → remove chirality.
        cleanup_chirality(&mut mol);
        assert_eq!(mol.atom(c).chirality, Chirality::None);
    }

    #[test]
    fn chirality_removed_all_identical() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 0,
            ..Atom::default()
        });
        for _ in 0..4 {
            let f = mol.add_atom(Atom {
                atomic_num: 9,
                ..Atom::default()
            });
            mol.add_bond(c, f, Bond::default());
        }
        cleanup_chirality(&mut mol);
        assert_eq!(mol.atom(c).chirality, Chirality::None);
    }

    #[test]
    fn chirality_preserved_distinct_neighbors() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 1,
            ..Atom::default()
        });
        let f = mol.add_atom(Atom {
            atomic_num: 9,
            ..Atom::default()
        });
        let cl = mol.add_atom(Atom {
            atomic_num: 17,
            ..Atom::default()
        });
        let br = mol.add_atom(Atom {
            atomic_num: 35,
            ..Atom::default()
        });
        mol.add_bond(c, f, Bond::default());
        mol.add_bond(c, cl, Bond::default());
        mol.add_bond(c, br, Bond::default());
        cleanup_chirality(&mut mol);
        assert_eq!(mol.atom(c).chirality, Chirality::Cw);
    }

    #[test]
    fn nitrogen_with_three_neighbors_preserved() {
        let mut mol = Mol::<Atom, Bond>::new();
        let n = mol.add_atom(Atom {
            atomic_num: 7,
            chirality: Chirality::Cw,
            hydrogen_count: 0,
            ..Atom::default()
        });
        for anum in [6, 9, 17] {
            let a = mol.add_atom(Atom {
                atomic_num: anum,
                ..Atom::default()
            });
            mol.add_bond(n, a, Bond::default());
        }
        cleanup_chirality(&mut mol);
        assert_eq!(mol.atom(n).chirality, Chirality::Cw);
    }

    #[test]
    fn non_chiral_atom_untouched() {
        let mut mol = from_smiles("CC").unwrap();
        cleanup_chirality(&mut mol);
        for idx in mol.atoms() {
            assert_eq!(mol.atom(idx).chirality, Chirality::None);
        }
    }
}
