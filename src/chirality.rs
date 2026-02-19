use crate::atom::Atom;
use crate::mol::{AtomId, Mol};
use crate::traits::HasBondOrder;

pub fn cleanup_chirality<B>(mol: &mut Mol<Atom, B>)
where
    B: HasBondOrder,
{
    let to_remove: Vec<petgraph::graph::NodeIndex> = mol
        .tetrahedral_stereo()
        .iter()
        .filter_map(|s| {
            let center = match s[0] {
                AtomId::Node(idx) => idx,
                _ => return Some(petgraph::graph::NodeIndex::new(usize::MAX)),
            };
            let atom = mol.atom(center);
            let total_neighbors = mol.neighbors(center).count() as u8 + atom.hydrogen_count;

            if total_neighbors < 3 {
                return Some(center);
            }
            if total_neighbors < 4 && !has_lone_pair_stereo(atom.atomic_num) {
                return Some(center);
            }
            if all_neighbors_identical(mol, center) {
                return Some(center);
            }
            None
        })
        .collect();

    for center in to_remove {
        mol.remove_tetrahedral_stereo(center);
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
    use crate::mol::{AtomId, Mol};
    use crate::smiles::from_smiles;

    #[test]
    fn chiral_center_preserved() {
        let mut mol = from_smiles("F[C@H](Cl)Br").unwrap();
        cleanup_chirality(&mut mol);
        let c_idx = mol
            .atoms()
            .find(|&idx| mol.atom(idx).atomic_num == 6)
            .unwrap();
        assert!(mol.tetrahedral_stereo_for(c_idx).is_some());
    }

    #[test]
    fn chirality_removed_too_few_neighbors() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 1,
            ..Atom::default()
        });
        let f = mol.add_atom(Atom {
            atomic_num: 9,
            ..Atom::default()
        });
        mol.add_bond(c, f, Bond::default());
        mol.add_tetrahedral_stereo([
            AtomId::Node(c),
            AtomId::Node(f),
            AtomId::VirtualH(c, 0),
            AtomId::VirtualH(c, 0),
        ]);
        // C has 1 explicit + 1 implicit H = 2 neighbors. < 3 → remove chirality.
        cleanup_chirality(&mut mol);
        assert!(mol.tetrahedral_stereo_for(c).is_none());
    }

    #[test]
    fn chirality_removed_all_identical() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let mut fs = Vec::new();
        for _ in 0..4 {
            let f = mol.add_atom(Atom {
                atomic_num: 9,
                ..Atom::default()
            });
            mol.add_bond(c, f, Bond::default());
            fs.push(f);
        }
        mol.add_tetrahedral_stereo([
            AtomId::Node(c),
            AtomId::Node(fs[0]),
            AtomId::Node(fs[1]),
            AtomId::Node(fs[2]),
        ]);
        cleanup_chirality(&mut mol);
        assert!(mol.tetrahedral_stereo_for(c).is_none());
    }

    #[test]
    fn chirality_preserved_distinct_neighbors() {
        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
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
        mol.add_tetrahedral_stereo([
            AtomId::Node(c),
            AtomId::Node(f),
            AtomId::Node(cl),
            AtomId::Node(br),
        ]);
        cleanup_chirality(&mut mol);
        assert!(mol.tetrahedral_stereo_for(c).is_some());
    }

    #[test]
    fn nitrogen_with_three_neighbors_preserved() {
        let mut mol = Mol::<Atom, Bond>::new();
        let n = mol.add_atom(Atom {
            atomic_num: 7,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let mut neighbors = Vec::new();
        for anum in [6, 9, 17] {
            let a = mol.add_atom(Atom {
                atomic_num: anum,
                ..Atom::default()
            });
            mol.add_bond(n, a, Bond::default());
            neighbors.push(a);
        }
        mol.add_tetrahedral_stereo([
            AtomId::Node(n),
            AtomId::Node(neighbors[0]),
            AtomId::Node(neighbors[1]),
            AtomId::Node(neighbors[2]),
        ]);
        cleanup_chirality(&mut mol);
        assert!(mol.tetrahedral_stereo_for(n).is_some());
    }

    #[test]
    fn non_chiral_atom_untouched() {
        let mut mol = from_smiles("CC").unwrap();
        cleanup_chirality(&mut mol);
        assert!(mol.tetrahedral_stereo().is_empty());
    }
}
