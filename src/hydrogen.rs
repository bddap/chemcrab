use crate::atom::{Atom, Chirality};
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::mol::Mol;

fn flip_chirality(c: Chirality) -> Chirality {
    match c {
        Chirality::Cw => Chirality::Ccw,
        Chirality::Ccw => Chirality::Cw,
        Chirality::None => Chirality::None,
    }
}

pub fn add_hs(mol: &Mol<Atom, Bond>) -> Mol<Atom, Bond> {
    let mut result = Mol::new();
    let mut index_map = Vec::new();

    for idx in mol.atoms() {
        let atom = mol.atom(idx);
        let new_idx = result.add_atom(Atom {
            hydrogen_count: 0,
            ..*atom
        });
        index_map.push(new_idx);
    }

    for edge in mol.bonds() {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        let bond = mol.bond(edge);
        result.add_bond(index_map[a.index()], index_map[b.index()], bond.clone());
    }

    for (idx, &parent) in index_map.iter().enumerate() {
        let orig_idx = petgraph::graph::NodeIndex::new(idx);
        let atom = mol.atom(orig_idx);
        let h_count = atom.hydrogen_count as usize;
        for _ in 0..h_count {
            let h = result.add_atom(Atom {
                atomic_num: 1,
                ..Atom::default()
            });
            result.add_bond(
                parent,
                h,
                Bond {
                    order: BondOrder::Single,
                    stereo: BondStereo::None,
                },
            );
        }
        if matches!(atom.chirality, Chirality::Cw | Chirality::Ccw) && h_count > 0 {
            let graph_neighbor_count = mol.neighbors(orig_idx).count();
            if (h_count * graph_neighbor_count) % 2 == 1 {
                result.atom_mut(parent).chirality = flip_chirality(atom.chirality);
            }
        }
    }

    result
}

pub fn remove_hs(mol: &Mol<Atom, Bond>) -> Mol<Atom, Bond> {
    let node_count = mol.atom_count();
    let mut removable = vec![false; node_count];
    let mut extra_h: Vec<u8> = vec![0; node_count];

    for idx in mol.atoms() {
        let atom = mol.atom(idx);
        if atom.atomic_num == 1
            && atom.isotope == 0
            && atom.formal_charge == 0
        {
            let neighbors: Vec<_> = mol.neighbors(idx).collect();
            if neighbors.len() == 1 {
                removable[idx.index()] = true;
                extra_h[neighbors[0].index()] += 1;
            }
        }
    }

    let mut result = Mol::new();
    let mut index_map = vec![None; node_count];

    for idx in mol.atoms() {
        if removable[idx.index()] {
            continue;
        }
        let atom = mol.atom(idx);
        let chirality = match atom.chirality {
            Chirality::Cw | Chirality::Ccw => {
                let neighbors: Vec<_> = mol.neighbors(idx).collect();
                let removable_after_count: usize = neighbors
                    .iter()
                    .enumerate()
                    .filter(|(_, n)| removable[n.index()])
                    .map(|(pos, _)| {
                        neighbors[pos + 1..]
                            .iter()
                            .filter(|n| !removable[n.index()])
                            .count()
                    })
                    .sum();
                if removable_after_count % 2 == 1 {
                    flip_chirality(atom.chirality)
                } else {
                    atom.chirality
                }
            }
            Chirality::None => Chirality::None,
        };
        let new_idx = result.add_atom(Atom {
            chirality,
            hydrogen_count: atom.hydrogen_count + extra_h[idx.index()],
            ..*atom
        });
        index_map[idx.index()] = Some(new_idx);
    }

    for edge in mol.bonds() {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        if let (Some(new_a), Some(new_b)) = (index_map[a.index()], index_map[b.index()]) {
            result.add_bond(new_a, new_b, mol.bond(edge).clone());
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bond::SmilesBondOrder;
    use crate::smiles::parse_smiles;
    use crate::SmilesBond;
    use petgraph::graph::NodeIndex;

    fn smiles_to_mol(s: &str) -> Mol<Atom, Bond> {
        let smol = parse_smiles(s).unwrap();
        let mut mol = Mol::new();
        let mut map = Vec::new();
        for idx in smol.atoms() {
            map.push(mol.add_atom(smol.atom(idx).clone()));
        }
        for edge in smol.bonds() {
            let (a, b) = smol.bond_endpoints(edge).unwrap();
            let sb: &SmilesBond = smol.bond(edge);
            let order = match sb.order {
                SmilesBondOrder::Single | SmilesBondOrder::Implicit => BondOrder::Single,
                SmilesBondOrder::Double => BondOrder::Double,
                SmilesBondOrder::Triple => BondOrder::Triple,
                SmilesBondOrder::Aromatic => BondOrder::Single,
            };
            mol.add_bond(
                map[a.index()],
                map[b.index()],
                Bond {
                    order,
                    stereo: BondStereo::None,
                },
            );
        }
        mol
    }

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    #[test]
    fn add_hs_methane() {
        let mol = smiles_to_mol("C");
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.atom(n(0)).hydrogen_count, 4);

        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom_count(), 5);
        assert_eq!(explicit.bond_count(), 4);
        assert_eq!(explicit.atom(n(0)).hydrogen_count, 0);
        assert_eq!(explicit.atom(n(0)).atomic_num, 6);
        for i in 1..5 {
            assert_eq!(explicit.atom(n(i)).atomic_num, 1);
        }
    }

    #[test]
    fn add_hs_ethane() {
        let mol = smiles_to_mol("CC");
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom_count(), 8);
        assert_eq!(explicit.bond_count(), 7);
    }

    #[test]
    fn add_hs_water() {
        let mol = smiles_to_mol("O");
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom_count(), 3);
        assert_eq!(explicit.bond_count(), 2);
        assert_eq!(explicit.atom(n(0)).atomic_num, 8);
    }

    #[test]
    fn add_hs_iron_no_virtual() {
        let mol = smiles_to_mol("[Fe]");
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom_count(), 1);
        assert_eq!(explicit.bond_count(), 0);
    }

    #[test]
    fn add_hs_preserves_chirality_no_virtual_h() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 0,
            ..Atom::default()
        });
        for _ in 0..4 {
            let n = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, n, Bond::default());
        }
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom(n(0)).chirality, Chirality::Cw);
    }

    #[test]
    fn add_hs_remaps_chirality_with_virtual_h() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Ccw,
            hydrogen_count: 1,
            ..Atom::default()
        });
        for _ in 0..3 {
            let n = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, n, Bond::default());
        }
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom(n(0)).chirality, Chirality::Cw);
    }

    #[test]
    fn remove_hs_preserves_chirality_no_removable_h_neighbor() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 0,
            ..Atom::default()
        });
        for _ in 0..4 {
            let n = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, n, Bond::default());
        }
        let result = remove_hs(&mol);
        assert_eq!(result.atom(n(0)).chirality, Chirality::Cw);
    }

    #[test]
    fn chirality_round_trip_no_h_on_center() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Ccw,
            hydrogen_count: 0,
            ..Atom::default()
        });
        for _ in 0..4 {
            let n = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, n, Bond::default());
        }
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom(n(0)).chirality, Chirality::Ccw);
        let collapsed = remove_hs(&explicit);
        assert_eq!(collapsed.atom(n(0)).chirality, Chirality::Ccw);
    }

    #[test]
    fn chirality_round_trip_with_h_on_center() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            chirality: Chirality::Cw,
            hydrogen_count: 1,
            ..Atom::default()
        });
        for _ in 0..3 {
            let nn = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nn, Bond::default());
        }
        let explicit = add_hs(&mol);
        assert_eq!(explicit.atom(n(0)).chirality, Chirality::Ccw);
        let collapsed = remove_hs(&explicit);
        assert_eq!(collapsed.atom(n(0)).chirality, Chirality::Cw);
        assert_eq!(collapsed.atom(n(0)).hydrogen_count, 1);
    }

    #[test]
    fn remove_hs_round_trip() {
        let mol = smiles_to_mol("CC");
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);

        assert_eq!(collapsed.atom_count(), 2);
        assert_eq!(collapsed.bond_count(), 1);
        assert_eq!(collapsed.atom(n(0)).atomic_num, 6);
        assert_eq!(collapsed.atom(n(0)).hydrogen_count, 3);
        assert_eq!(collapsed.atom(n(1)).atomic_num, 6);
        assert_eq!(collapsed.atom(n(1)).hydrogen_count, 3);
    }

    #[test]
    fn remove_hs_keeps_deuterium() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let d = mol.add_atom(Atom {
            atomic_num: 1,
            isotope: 2,
            ..Atom::default()
        });
        let h1 = mol.add_atom(Atom {
            atomic_num: 1,
            ..Atom::default()
        });
        let h2 = mol.add_atom(Atom {
            atomic_num: 1,
            ..Atom::default()
        });
        let h3 = mol.add_atom(Atom {
            atomic_num: 1,
            ..Atom::default()
        });
        let single = Bond::default();
        mol.add_bond(c, d, single.clone());
        mol.add_bond(c, h1, single.clone());
        mol.add_bond(c, h2, single.clone());
        mol.add_bond(c, h3, single);

        let result = remove_hs(&mol);
        assert_eq!(result.atom_count(), 2);
        assert_eq!(result.atom(n(0)).atomic_num, 6);
        assert_eq!(result.atom(n(0)).hydrogen_count, 3);
        assert_eq!(result.atom(n(1)).atomic_num, 1);
        assert_eq!(result.atom(n(1)).isotope, 2);
    }

    #[test]
    fn remove_hs_keeps_charged_h() {
        let mut mol = Mol::new();
        mol.add_atom(Atom {
            atomic_num: 1,
            formal_charge: 1,
            ..Atom::default()
        });
        let result = remove_hs(&mol);
        assert_eq!(result.atom_count(), 1);
        assert_eq!(result.atom(n(0)).atomic_num, 1);
        assert_eq!(result.atom(n(0)).formal_charge, 1);
    }

    #[test]
    fn remove_hs_no_explicit_hs() {
        let mol = smiles_to_mol("C");
        let result = remove_hs(&mol);
        assert_eq!(result.atom_count(), 1);
        assert_eq!(result.atom(n(0)).hydrogen_count, 4);
    }

    #[test]
    fn empty_molecule_add_hs() {
        let mol = Mol::<Atom, Bond>::new();
        let result = add_hs(&mol);
        assert_eq!(result.atom_count(), 0);
        assert_eq!(result.bond_count(), 0);
    }

    #[test]
    fn empty_molecule_remove_hs() {
        let mol = Mol::<Atom, Bond>::new();
        let result = remove_hs(&mol);
        assert_eq!(result.atom_count(), 0);
        assert_eq!(result.bond_count(), 0);
    }
}
