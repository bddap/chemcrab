use crate::*;

#[test]
fn mol_add_atoms_and_bonds() {
    let mut mol = Mol::<Atom, Bond>::new();
    let c = mol.add_atom(Atom {
        atomic_num: 6,
        ..Atom::default()
    });
    let o = mol.add_atom(Atom {
        atomic_num: 8,
        ..Atom::default()
    });
    let bond_idx = mol.add_bond(
        c,
        o,
        Bond {
            order: BondOrder::Double,
            ..Bond::default()
        },
    );

    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.atom(c).atomic_num, 6);
    assert_eq!(mol.atom(o).atomic_num, 8);
    assert_eq!(mol.bond(bond_idx).order, BondOrder::Double);
}

#[test]
fn mol_neighbors_and_bonds_of() {
    let mut mol = Mol::<Atom, Bond>::new();
    let a = mol.add_atom(Atom::default());
    let b = mol.add_atom(Atom::default());
    let c = mol.add_atom(Atom::default());
    mol.add_bond(a, b, Bond::default());
    mol.add_bond(a, c, Bond::default());

    let neighbors: Vec<_> = mol.neighbors(a).collect();
    assert_eq!(neighbors.len(), 2);

    let incident: Vec<_> = mol.bonds_of(a).collect();
    assert_eq!(incident.len(), 2);
}

#[test]
fn mol_bond_between_and_endpoints() {
    let mut mol = Mol::<Atom, Bond>::new();
    let a = mol.add_atom(Atom::default());
    let b = mol.add_atom(Atom::default());
    let c = mol.add_atom(Atom::default());
    let e = mol.add_bond(a, b, Bond::default());

    assert_eq!(mol.bond_between(a, b), Some(e));
    assert_eq!(mol.bond_between(a, c), None);

    let (src, dst) = mol.bond_endpoints(e).unwrap();
    assert!((src == a && dst == b) || (src == b && dst == a));
}

#[test]
fn mol_iterators() {
    let mut mol = Mol::<Atom, Bond>::new();
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());

    assert_eq!(mol.atoms().count(), 2);
    assert_eq!(mol.bonds().count(), 0);
}

#[test]
fn mol_atom_mut() {
    let mut mol = Mol::<Atom, Bond>::new();
    let idx = mol.add_atom(Atom::default());
    mol.atom_mut(idx).atomic_num = 7;
    assert_eq!(mol.atom(idx).atomic_num, 7);
}

#[test]
fn atom_trait_impls() {
    let atom = Atom {
        atomic_num: 6,
        formal_charge: -1,
        isotope: 13,
        chirality: Chirality::Cw,
        hydrogen_count: 3,
        is_aromatic: true,
    };

    assert_eq!(HasAtomicNum::atomic_num(&atom), 6);
    assert_eq!(HasFormalCharge::formal_charge(&atom), -1);
    assert_eq!(HasIsotope::isotope(&atom), 13);
    assert_eq!(HasChirality::chirality(&atom), Chirality::Cw);
    assert_eq!(HasHydrogenCount::hydrogen_count(&atom), 3);
    assert!(HasAromaticity::is_aromatic(&atom));
}

#[test]
fn bond_trait_impls() {
    let bond = Bond {
        order: BondOrder::Triple,
        stereo: BondStereo::None,
    };

    assert_eq!(HasBondOrder::bond_order(&bond), BondOrder::Triple);
    assert_eq!(HasBondStereo::bond_stereo(&bond), BondStereo::None);
}

#[test]
fn smiles_bond_has_bond_stereo_but_not_order() {
    use petgraph::graph::NodeIndex;

    let sb = SmilesBond {
        order: SmilesBondOrder::Aromatic,
        stereo: BondStereo::Cis(NodeIndex::new(0), NodeIndex::new(1)),
    };

    assert_eq!(
        HasBondStereo::bond_stereo(&sb),
        BondStereo::Cis(NodeIndex::new(0), NodeIndex::new(1))
    );
}

#[test]
fn with_valence_wraps_atom() {
    let atom = Atom {
        atomic_num: 6,
        hydrogen_count: 2,
        ..Atom::default()
    };
    let enriched = WithValence {
        inner: atom,
        valence: 4,
    };

    assert_eq!(HasValence::valence(&enriched), 4);
    assert_eq!(HasAtomicNum::atomic_num(&enriched), 6);
    assert_eq!(HasHydrogenCount::hydrogen_count(&enriched), 2);
}

#[test]
fn with_hybridization_wraps_atom() {
    let atom = Atom {
        atomic_num: 7,
        ..Atom::default()
    };
    let enriched = WithHybridization {
        inner: atom,
        hybridization: Hybridization::SP2,
    };

    assert_eq!(HasHybridization::hybridization(&enriched), Hybridization::SP2);
    assert_eq!(HasAtomicNum::atomic_num(&enriched), 7);
}

#[test]
fn nested_wrappers() {
    let atom = Atom {
        atomic_num: 8,
        ..Atom::default()
    };
    let enriched = WithHybridization {
        inner: WithValence {
            inner: atom,
            valence: 2,
        },
        hybridization: Hybridization::SP3,
    };

    assert_eq!(HasHybridization::hybridization(&enriched), Hybridization::SP3);
    assert_eq!(HasValence::valence(&enriched), 2);
    assert_eq!(HasAtomicNum::atomic_num(&enriched), 8);
}

#[test]
fn with_position_2d() {
    let atom = Atom {
        atomic_num: 6,
        ..Atom::default()
    };
    let mut enriched = WithPosition2D {
        inner: atom,
        position_2d: Some([1.0, 2.0]),
    };

    assert_eq!(HasPosition2D::position_2d(&enriched), Some([1.0, 2.0]));
    HasPosition2D::set_position_2d(&mut enriched, None);
    assert_eq!(HasPosition2D::position_2d(&enriched), None);
    assert_eq!(HasAtomicNum::atomic_num(&enriched), 6);
}

#[test]
fn with_position_3d() {
    let atom = Atom {
        atomic_num: 6,
        ..Atom::default()
    };
    let mut enriched = WithPosition3D {
        inner: atom,
        position_3d: Some([1.0, 2.0, 3.0]),
    };

    assert_eq!(
        HasPosition3D::position_3d(&enriched),
        Some([1.0, 2.0, 3.0])
    );
    HasPosition3D::set_position_3d(&mut enriched, None);
    assert_eq!(HasPosition3D::position_3d(&enriched), None);
}

#[test]
fn chirality_default_is_none() {
    assert_eq!(Chirality::default(), Chirality::None);
}

#[test]
fn bond_order_default_is_single() {
    assert_eq!(BondOrder::default(), BondOrder::Single);
}

#[test]
fn atom_default() {
    let atom = Atom::default();
    assert_eq!(atom.atomic_num, 0);
    assert_eq!(atom.formal_charge, 0);
    assert_eq!(atom.isotope, 0);
    assert_eq!(atom.chirality, Chirality::None);
    assert_eq!(atom.hydrogen_count, 0);
    assert!(!atom.is_aromatic);
}

#[test]
fn mol_default() {
    let mol = Mol::<Atom, Bond>::default();
    assert_eq!(mol.atom_count(), 0);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn mol_graph_access() {
    let mut mol = Mol::<Atom, Bond>::new();
    mol.add_atom(Atom::default());
    assert_eq!(mol.graph().node_count(), 1);
}
