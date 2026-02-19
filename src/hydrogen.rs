use petgraph::graph::NodeIndex;

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

fn remap_stereo(stereo: BondStereo, map: &dyn Fn(NodeIndex) -> Option<NodeIndex>) -> BondStereo {
    match stereo {
        BondStereo::Cis(a, b) => match (map(a), map(b)) {
            (Some(na), Some(nb)) => BondStereo::Cis(na, nb),
            _ => BondStereo::None,
        },
        BondStereo::Trans(a, b) => match (map(a), map(b)) {
            (Some(na), Some(nb)) => BondStereo::Trans(na, nb),
            _ => BondStereo::None,
        },
        BondStereo::None => BondStereo::None,
    }
}

fn stereo_referenced_atoms(mol: &Mol<Atom, Bond>) -> Vec<bool> {
    let mut referenced = vec![false; mol.atom_count()];
    for edge in mol.bonds() {
        match mol.bond(edge).stereo {
            BondStereo::Cis(a, b) | BondStereo::Trans(a, b) => {
                referenced[a.index()] = true;
                referenced[b.index()] = true;
            }
            BondStereo::None => {}
        }
    }
    referenced
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
        let stereo = remap_stereo(bond.stereo, &|n| Some(index_map[n.index()]));
        result.add_bond(
            index_map[a.index()],
            index_map[b.index()],
            Bond { stereo, ..bond.clone() },
        );
    }

    for (idx, &parent) in index_map.iter().enumerate() {
        let orig_idx = NodeIndex::new(idx);
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

#[derive(Debug, Clone, Default)]
pub struct RemoveHsOptions {
    pub keep_stereo_hs: bool,
}

pub fn remove_hs(mol: &Mol<Atom, Bond>) -> Mol<Atom, Bond> {
    remove_hs_with(mol, &RemoveHsOptions::default())
}

pub fn remove_hs_with(mol: &Mol<Atom, Bond>, opts: &RemoveHsOptions) -> Mol<Atom, Bond> {
    let node_count = mol.atom_count();
    let stereo_refs = if opts.keep_stereo_hs {
        stereo_referenced_atoms(mol)
    } else {
        vec![false; node_count]
    };

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
                if opts.keep_stereo_hs && stereo_refs[idx.index()] {
                    continue;
                }
                removable[idx.index()] = true;
                extra_h[neighbors[0].index()] += 1;
            }
        }
    }

    let mut result = Mol::new();
    let mut index_map: Vec<Option<NodeIndex>> = vec![None; node_count];

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
            let bond = mol.bond(edge);
            let stereo = remap_stereo(bond.stereo, &|n| index_map[n.index()]);
            result.add_bond(new_a, new_b, Bond { stereo, ..bond.clone() });
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bond::SmilesBondOrder;
    use crate::smiles::{from_smiles, parse_smiles, to_smiles};
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

    fn find_double_bond_stereo(mol: &Mol<Atom, Bond>) -> Option<BondStereo> {
        mol.bonds().find_map(|e| match mol.bond(e).stereo {
            BondStereo::None => None,
            s => Some(s),
        })
    }

    fn stereo_is_trans(s: BondStereo) -> bool {
        matches!(s, BondStereo::Trans(_, _))
    }

    fn stereo_is_cis(s: BondStereo) -> bool {
        matches!(s, BondStereo::Cis(_, _))
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

    // ---- E/Z stereo preservation ----

    #[test]
    fn add_hs_preserves_ez_trans() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let stereo = find_double_bond_stereo(&mol).unwrap();
        assert!(stereo_is_trans(stereo));

        let explicit = add_hs(&mol);
        let stereo2 = find_double_bond_stereo(&explicit)
            .expect("E/Z stereo lost after add_hs");
        assert!(stereo_is_trans(stereo2));

        if let BondStereo::Trans(a, b) = stereo2 {
            assert!(a.index() < explicit.atom_count());
            assert!(b.index() < explicit.atom_count());
        }
    }

    #[test]
    fn add_hs_preserves_ez_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        let stereo = find_double_bond_stereo(&mol).unwrap();
        assert!(stereo_is_cis(stereo));

        let explicit = add_hs(&mol);
        let stereo2 = find_double_bond_stereo(&explicit)
            .expect("E/Z stereo lost after add_hs");
        assert!(stereo_is_cis(stereo2));
    }

    #[test]
    fn remove_hs_preserves_ez_trans() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let stereo = find_double_bond_stereo(&collapsed)
            .expect("E/Z stereo lost after remove_hs");
        assert!(stereo_is_trans(stereo));
    }

    #[test]
    fn remove_hs_preserves_ez_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let stereo = find_double_bond_stereo(&collapsed)
            .expect("E/Z stereo lost after remove_hs");
        assert!(stereo_is_cis(stereo));
    }

    #[test]
    fn ez_round_trip_smiles_f_trans() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let smiles = to_smiles(&collapsed);
        let reparsed = from_smiles(&smiles).unwrap();
        let stereo = find_double_bond_stereo(&reparsed)
            .expect("E/Z stereo lost in round trip");
        assert!(stereo_is_trans(stereo));
    }

    #[test]
    fn ez_round_trip_smiles_f_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let smiles = to_smiles(&collapsed);
        let reparsed = from_smiles(&smiles).unwrap();
        let stereo = find_double_bond_stereo(&reparsed)
            .expect("E/Z stereo lost in round trip");
        assert!(stereo_is_cis(stereo));
    }

    #[test]
    fn ez_round_trip_cl_trans() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let stereo = find_double_bond_stereo(&collapsed)
            .expect("E/Z stereo lost after round trip");
        assert!(stereo_is_trans(stereo));
    }

    #[test]
    fn ez_round_trip_alkyl() {
        let mol = from_smiles("CC/C=C/CC").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let stereo = find_double_bond_stereo(&collapsed)
            .expect("E/Z stereo lost after round trip");
        assert!(stereo_is_trans(stereo));
    }

    #[test]
    fn ez_stereo_ref_atoms_valid_after_add_hs() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        if let Some(BondStereo::Trans(a, b) | BondStereo::Cis(a, b)) =
            find_double_bond_stereo(&explicit)
        {
            assert!(
                a.index() < explicit.atom_count(),
                "ref atom {a:?} out of bounds"
            );
            assert!(
                b.index() < explicit.atom_count(),
                "ref atom {b:?} out of bounds"
            );
        }
    }

    #[test]
    fn ez_stereo_ref_atoms_valid_after_remove_hs() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        if let Some(BondStereo::Trans(a, b) | BondStereo::Cis(a, b)) =
            find_double_bond_stereo(&collapsed)
        {
            assert!(
                a.index() < collapsed.atom_count(),
                "ref atom {a:?} out of bounds"
            );
            assert!(
                b.index() < collapsed.atom_count(),
                "ref atom {b:?} out of bounds"
            );
        }
    }

    // ---- keep_stereo_hs option ----

    #[test]
    fn remove_hs_keep_stereo_hs_preserves_ref_atoms() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);

        let stereo = find_double_bond_stereo(&result)
            .expect("E/Z stereo lost with keep_stereo_hs");
        assert!(stereo_is_trans(stereo));

        if let BondStereo::Trans(a, b) = stereo {
            assert_eq!(result.atom(a).atomic_num, 9, "left ref should be F");
            assert_eq!(result.atom(b).atomic_num, 9, "right ref should be F");
        }
    }

    #[test]
    fn keep_stereo_hs_retains_h_on_double_bond_carbon() {
        // In CC/C=C/CC, the stereo refs are the C neighbors, not H.
        // After add_hs, H atoms bonded to the double-bond carbons may become
        // stereo references. keep_stereo_hs should retain those.
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);

        let stereo_before = find_double_bond_stereo(&explicit).unwrap();

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);
        let stereo_after = find_double_bond_stereo(&result).unwrap();

        // Both should be trans
        assert!(stereo_is_trans(stereo_before));
        assert!(stereo_is_trans(stereo_after));
    }

    #[test]
    fn keep_stereo_hs_still_removes_non_stereo_hs() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        // explicit has F, C(H), =, C(H), F => 4 heavy + some H atoms
        let total_explicit = explicit.atom_count();
        assert!(total_explicit > 4);

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);
        // Should have fewer atoms than fully explicit but may keep stereo H
        assert!(result.atom_count() <= total_explicit);
    }

    #[test]
    fn remove_hs_default_does_not_keep_stereo_hs() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        // Default: all plain H removed, back to implicit
        assert_eq!(collapsed.atom_count(), 4);
    }
}
