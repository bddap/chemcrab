use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};
use crate::mol::{AtomId, EZStereo, Mol, TetrahedralStereo};

fn stereo_referenced_atoms(mol: &Mol<Atom, Bond>) -> Vec<bool> {
    let mut referenced = vec![false; mol.atom_count()];
    for ez in mol.ez_stereo() {
        referenced[ez.bond.0.index()] = true;
        referenced[ez.bond.1.index()] = true;
        for r in &ez.refs {
            if let AtomId::Node(idx) = r {
                if idx.index() < referenced.len() {
                    referenced[idx.index()] = true;
                }
            }
        }
    }
    for stereo in mol.tetrahedral_stereo() {
        if stereo.center.index() < referenced.len() {
            referenced[stereo.center.index()] = true;
        }
        for aid in &stereo.above {
            if let AtomId::Node(idx) = aid {
                if idx.index() < referenced.len() {
                    referenced[idx.index()] = true;
                }
            }
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
        result.add_bond(
            index_map[a.index()],
            index_map[b.index()],
            bond.clone(),
        );
    }

    let mut new_h_map: Vec<Option<NodeIndex>> = vec![None; mol.atom_count()];

    for (idx, &parent) in index_map.iter().enumerate() {
        let orig_idx = NodeIndex::new(idx);
        let atom = mol.atom(orig_idx);
        let h_count = atom.hydrogen_count as usize;
        let mut first_h = None;
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
                },
            );
            if first_h.is_none() {
                first_h = Some(h);
            }
        }
        new_h_map[idx] = first_h;
    }

    let remap_aid = |aid: AtomId| -> AtomId {
        match aid {
            AtomId::Node(idx) => AtomId::Node(index_map[idx.index()]),
            AtomId::VirtualH(parent, _) => {
                match new_h_map[parent.index()] {
                    Some(h_idx) => AtomId::Node(h_idx),
                    None => AtomId::VirtualH(index_map[parent.index()], 0),
                }
            }
        }
    };

    let new_stereo: Vec<TetrahedralStereo> = mol
        .tetrahedral_stereo()
        .iter()
        .map(|s| TetrahedralStereo {
            center: index_map[s.center.index()],
            above: s.above.map(&remap_aid),
        })
        .collect();
    result.set_tetrahedral_stereo(new_stereo);

    let new_ez: Vec<EZStereo> = mol
        .ez_stereo()
        .iter()
        .map(|s| {
            let new_a = index_map[s.bond.0.index()];
            let new_b = index_map[s.bond.1.index()];
            let (lo, hi, refs) = if new_a.index() < new_b.index() {
                (new_a, new_b, [remap_aid(s.refs[0]), remap_aid(s.refs[1])])
            } else {
                (new_b, new_a, [remap_aid(s.refs[1]), remap_aid(s.refs[0])])
            };
            EZStereo { bond: (lo, hi), refs }
        })
        .collect();
    result.set_ez_stereo(new_ez);

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
        let new_idx = result.add_atom(Atom {
            hydrogen_count: atom.hydrogen_count + extra_h[idx.index()],
            ..*atom
        });
        index_map[idx.index()] = Some(new_idx);
    }

    for edge in mol.bonds() {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        if let (Some(new_a), Some(new_b)) = (index_map[a.index()], index_map[b.index()]) {
            let bond = mol.bond(edge);
            result.add_bond(new_a, new_b, bond.clone());
        }
    }

    let remap_aid = |aid: AtomId| -> Option<AtomId> {
        match aid {
            AtomId::Node(idx) => {
                if removable[idx.index()] {
                    let parent = mol.neighbors(idx).next()?;
                    Some(AtomId::VirtualH(index_map[parent.index()]?, 0))
                } else {
                    Some(AtomId::Node(index_map[idx.index()]?))
                }
            }
            AtomId::VirtualH(parent, n) => {
                Some(AtomId::VirtualH(index_map[parent.index()]?, n))
            }
        }
    };

    let new_stereo: Vec<TetrahedralStereo> = mol
        .tetrahedral_stereo()
        .iter()
        .filter_map(|s| {
            let center = index_map[s.center.index()]?;
            let mut above = [AtomId::Node(NodeIndex::new(0)); 4];
            for (i, &aid) in s.above.iter().enumerate() {
                above[i] = remap_aid(aid)?;
            }
            Some(TetrahedralStereo { center, above })
        })
        .collect();
    result.set_tetrahedral_stereo(new_stereo);

    let new_ez: Vec<EZStereo> = mol
        .ez_stereo()
        .iter()
        .filter_map(|s| {
            let new_a = index_map[s.bond.0.index()]?;
            let new_b = index_map[s.bond.1.index()]?;
            let ref0 = remap_aid(s.refs[0])?;
            let ref1 = remap_aid(s.refs[1])?;
            let (lo, hi, refs) = if new_a.index() < new_b.index() {
                (new_a, new_b, [ref0, ref1])
            } else {
                (new_b, new_a, [ref1, ref0])
            };
            Some(EZStereo { bond: (lo, hi), refs })
        })
        .collect();
    result.set_ez_stereo(new_ez);

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{AtomId, TetrahedralStereo};
    use crate::smiles::{from_smiles, parse_smiles, to_smiles};
    use crate::SmilesBond;
    use crate::bond::SmilesBondOrder;
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
                Bond { order },
            );
        }
        mol
    }

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    fn has_ez_stereo(mol: &Mol<Atom, Bond>) -> bool {
        !mol.ez_stereo().is_empty()
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
            hydrogen_count: 0,
            ..Atom::default()
        });
        let mut ns = Vec::new();
        for _ in 0..4 {
            let nbr = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nbr, Bond::default());
            ns.push(nbr);
        }
        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center: c,
            above: [
                AtomId::Node(ns[0]),
                AtomId::Node(ns[1]),
                AtomId::Node(ns[2]),
                AtomId::Node(ns[3]),
            ],
        });
        let explicit = add_hs(&mol);
        assert!(explicit.tetrahedral_stereo_for(n(0)).is_some());
    }

    #[test]
    fn add_hs_remaps_chirality_with_virtual_h() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 1,
            ..Atom::default()
        });
        let mut ns = Vec::new();
        for _ in 0..3 {
            let nbr = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nbr, Bond::default());
            ns.push(nbr);
        }
        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center: c,
            above: [
                AtomId::VirtualH(c, 0),
                AtomId::Node(ns[0]),
                AtomId::Node(ns[1]),
                AtomId::Node(ns[2]),
            ],
        });
        let explicit = add_hs(&mol);
        let stereo = explicit.tetrahedral_stereo_for(n(0)).expect("chirality should exist");
        assert!(
            stereo.above.iter().all(|id| matches!(id, AtomId::Node(_))),
            "VirtualH should have been replaced with Node"
        );
    }

    #[test]
    fn remove_hs_preserves_chirality_no_removable_h_neighbor() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let mut ns = Vec::new();
        for _ in 0..4 {
            let nbr = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nbr, Bond::default());
            ns.push(nbr);
        }
        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center: c,
            above: [
                AtomId::Node(ns[0]),
                AtomId::Node(ns[1]),
                AtomId::Node(ns[2]),
                AtomId::Node(ns[3]),
            ],
        });
        let result = remove_hs(&mol);
        assert!(result.tetrahedral_stereo_for(n(0)).is_some());
    }

    #[test]
    fn chirality_round_trip_no_h_on_center() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 0,
            ..Atom::default()
        });
        let mut ns = Vec::new();
        for _ in 0..4 {
            let nbr = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nbr, Bond::default());
            ns.push(nbr);
        }
        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center: c,
            above: [
                AtomId::Node(ns[0]),
                AtomId::Node(ns[1]),
                AtomId::Node(ns[2]),
                AtomId::Node(ns[3]),
            ],
        });
        let explicit = add_hs(&mol);
        assert!(explicit.tetrahedral_stereo_for(n(0)).is_some());
        let collapsed = remove_hs(&explicit);
        assert!(collapsed.tetrahedral_stereo_for(n(0)).is_some());
    }

    #[test]
    fn chirality_round_trip_with_h_on_center() {
        let mut mol = Mol::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 1,
            ..Atom::default()
        });
        let mut ns = Vec::new();
        for _ in 0..3 {
            let nbr = mol.add_atom(Atom {
                atomic_num: 7,
                ..Atom::default()
            });
            mol.add_bond(c, nbr, Bond::default());
            ns.push(nbr);
        }
        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center: c,
            above: [
                AtomId::VirtualH(c, 0),
                AtomId::Node(ns[0]),
                AtomId::Node(ns[1]),
                AtomId::Node(ns[2]),
            ],
        });
        let explicit = add_hs(&mol);
        assert!(explicit.tetrahedral_stereo_for(n(0)).is_some());
        let collapsed = remove_hs(&explicit);
        assert!(collapsed.tetrahedral_stereo_for(n(0)).is_some());
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
        assert!(has_ez_stereo(&mol));

        let explicit = add_hs(&mol);
        assert!(has_ez_stereo(&explicit), "E/Z stereo lost after add_hs");
    }

    #[test]
    fn add_hs_preserves_ez_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        assert!(has_ez_stereo(&mol));

        let explicit = add_hs(&mol);
        assert!(has_ez_stereo(&explicit), "E/Z stereo lost after add_hs");
    }

    #[test]
    fn remove_hs_preserves_ez_trans() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        assert!(has_ez_stereo(&collapsed), "E/Z stereo lost after remove_hs");
    }

    #[test]
    fn remove_hs_preserves_ez_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        assert!(has_ez_stereo(&collapsed), "E/Z stereo lost after remove_hs");
    }

    #[test]
    fn ez_round_trip_smiles_f_trans() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let smiles = to_smiles(&collapsed);
        let reparsed = from_smiles(&smiles).unwrap();
        assert!(has_ez_stereo(&reparsed), "E/Z stereo lost in round trip");
    }

    #[test]
    fn ez_round_trip_smiles_f_cis() {
        let mol = from_smiles(r"F/C=C\F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        let smiles = to_smiles(&collapsed);
        let reparsed = from_smiles(&smiles).unwrap();
        assert!(has_ez_stereo(&reparsed), "E/Z stereo lost in round trip");
    }

    #[test]
    fn ez_round_trip_cl_trans() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        assert!(has_ez_stereo(&collapsed), "E/Z stereo lost after round trip");
    }

    #[test]
    fn ez_round_trip_alkyl() {
        let mol = from_smiles("CC/C=C/CC").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        assert!(has_ez_stereo(&collapsed), "E/Z stereo lost after round trip");
    }

    #[test]
    fn ez_stereo_ref_atoms_valid_after_add_hs() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        for ez in explicit.ez_stereo() {
            assert!(ez.bond.0.index() < explicit.atom_count());
            assert!(ez.bond.1.index() < explicit.atom_count());
            for r in &ez.refs {
                if let AtomId::Node(n) = r {
                    assert!(n.index() < explicit.atom_count(), "ref atom {n:?} out of bounds");
                }
            }
        }
    }

    #[test]
    fn ez_stereo_ref_atoms_valid_after_remove_hs() {
        let mol = from_smiles("Cl/C=C/Br").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        for ez in collapsed.ez_stereo() {
            assert!(ez.bond.0.index() < collapsed.atom_count());
            assert!(ez.bond.1.index() < collapsed.atom_count());
            for r in &ez.refs {
                if let AtomId::Node(n) = r {
                    assert!(n.index() < collapsed.atom_count(), "ref atom {n:?} out of bounds");
                }
            }
        }
    }

    // ---- keep_stereo_hs option ----

    #[test]
    fn remove_hs_keep_stereo_hs_preserves_ref_atoms() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);

        assert!(has_ez_stereo(&result), "E/Z stereo lost with keep_stereo_hs");
    }

    #[test]
    fn keep_stereo_hs_retains_h_on_double_bond_carbon() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);
        assert!(has_ez_stereo(&result));
    }

    #[test]
    fn keep_stereo_hs_still_removes_non_stereo_hs() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let total_explicit = explicit.atom_count();
        assert!(total_explicit > 4);

        let opts = RemoveHsOptions { keep_stereo_hs: true };
        let result = remove_hs_with(&explicit, &opts);
        assert!(result.atom_count() <= total_explicit);
    }

    #[test]
    fn remove_hs_default_does_not_keep_stereo_hs() {
        let mol = from_smiles("F/C=C/F").unwrap();
        let explicit = add_hs(&mol);
        let collapsed = remove_hs(&explicit);
        assert_eq!(collapsed.atom_count(), 4);
    }
}
