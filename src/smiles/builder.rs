use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{SmilesBond, SmilesBondOrder};
use crate::element::Element;
use crate::mol::{AtomId, EZStereo, Mol, TetrahedralStereo};
use crate::smiles::parse_tree::{ParseAtom, ParseTree};
use crate::smiles::tokenizer::{BondToken, ChiralityToken};

pub fn build_mol(tree: &ParseTree) -> Mol<Atom, SmilesBond> {
    let mut mol = Mol::new();
    let mut node_indices: Vec<NodeIndex> = Vec::with_capacity(tree.atoms.len());

    for parse_atom in &tree.atoms {
        let atom = Atom {
            atomic_num: parse_atom.element.atomic_num(),
            formal_charge: parse_atom.charge,
            isotope: parse_atom.isotope,
            hydrogen_count: 0,
            is_aromatic: parse_atom.is_aromatic,
        };
        node_indices.push(mol.add_atom(atom));
    }

    let mut added_edges: Vec<Vec<usize>> = vec![Vec::new(); tree.atoms.len()];

    for (i, parse_atom) in tree.atoms.iter().enumerate() {
        for neighbor in &parse_atom.neighbors {
            let j = neighbor.atom_idx;
            if !added_edges[i].contains(&j) {
                let order = resolve_bond_order(
                    &neighbor.bond,
                    parse_atom.is_aromatic,
                    tree.atoms[j].is_aromatic,
                );
                let bond = SmilesBond { order };
                mol.add_bond(node_indices[i], node_indices[j], bond);
                added_edges[i].push(j);
                added_edges[j].push(i);
            }
        }
    }

    resolve_chirality(&mut mol, tree, &node_indices);
    resolve_ez_stereo(&mut mol, tree, &node_indices);
    resolve_hydrogen_counts(&mut mol, tree, &node_indices);

    mol
}

fn resolve_bond_order(
    bond_tok: &Option<BondToken>,
    from_aromatic: bool,
    to_aromatic: bool,
) -> SmilesBondOrder {
    match bond_tok {
        Some(BondToken::Single) => SmilesBondOrder::Single,
        Some(BondToken::Double) => SmilesBondOrder::Double,
        Some(BondToken::Triple) => SmilesBondOrder::Triple,
        Some(BondToken::Aromatic) => SmilesBondOrder::Aromatic,
        Some(BondToken::Up) | Some(BondToken::Down) => SmilesBondOrder::Single,
        None => {
            if from_aromatic && to_aromatic {
                SmilesBondOrder::Aromatic
            } else {
                SmilesBondOrder::Implicit
            }
        }
    }
}

fn resolve_chirality(mol: &mut Mol<Atom, SmilesBond>, tree: &ParseTree, indices: &[NodeIndex]) {
    for (i, parse_atom) in tree.atoms.iter().enumerate() {
        if parse_atom.chirality == ChiralityToken::None {
            continue;
        }

        let center = indices[i];
        let has_bracket_h = parse_atom.is_bracket && parse_atom.hcount.unwrap_or(0) > 0;
        let has_preceding =
            !parse_atom.neighbors.is_empty() && parse_atom.neighbors[0].atom_idx < i;

        let explicit_smiles: Vec<AtomId> = parse_atom
            .neighbors
            .iter()
            .map(|n| AtomId::Node(indices[n.atom_idx]))
            .collect();

        let mut smiles_order: Vec<AtomId> = Vec::with_capacity(explicit_smiles.len() + 1);
        if has_bracket_h {
            if has_preceding {
                smiles_order.push(explicit_smiles[0]);
                smiles_order.push(AtomId::VirtualH(center, 0));
                smiles_order.extend_from_slice(&explicit_smiles[1..]);
            } else {
                smiles_order.push(AtomId::VirtualH(center, 0));
                smiles_order.extend_from_slice(&explicit_smiles);
            }
        } else {
            smiles_order.extend_from_slice(&explicit_smiles);
        }

        if smiles_order.len() < 3 {
            continue;
        }

        let mut above = if smiles_order.len() >= 4 {
            [
                smiles_order[0],
                smiles_order[1],
                smiles_order[2],
                smiles_order[3],
            ]
        } else {
            [
                AtomId::VirtualH(center, 0),
                smiles_order[0],
                smiles_order[1],
                smiles_order[2],
            ]
        };

        if parse_atom.chirality == ChiralityToken::Clockwise {
            above.swap(2, 3);
        }

        mol.add_tetrahedral_stereo(TetrahedralStereo {
            center,
            above,
        });
    }
}

fn resolve_ez_stereo(mol: &mut Mol<Atom, SmilesBond>, tree: &ParseTree, indices: &[NodeIndex]) {
    for (i, parse_atom) in tree.atoms.iter().enumerate() {
        for neighbor in &parse_atom.neighbors {
            let j = neighbor.atom_idx;
            if i >= j {
                continue;
            }

            let edge_idx = match mol.bond_between(indices[i], indices[j]) {
                Some(e) => e,
                None => continue,
            };

            if mol.bond(edge_idx).order != SmilesBondOrder::Double {
                continue;
            }

            let left_dir = find_directional_neighbor(tree, i, j);
            let right_dir = find_directional_neighbor(tree, j, i);

            if let (Some((left_atom, left_bond)), Some((right_atom, right_bond))) =
                (left_dir, right_dir)
            {
                let same_direction = match (left_bond, right_bond) {
                    (BondToken::Up, BondToken::Up) | (BondToken::Down, BondToken::Down) => true,
                    (BondToken::Up, BondToken::Down) | (BondToken::Down, BondToken::Up) => false,
                    _ => continue,
                };

                let (cis_left, cis_right) = if same_direction {
                    let other = other_substituent(mol, indices[j], indices[i], indices[right_atom]);
                    (AtomId::Node(indices[left_atom]), other)
                } else {
                    (AtomId::Node(indices[left_atom]), AtomId::Node(indices[right_atom]))
                };

                let (lo, hi) = if indices[i].index() < indices[j].index() {
                    (indices[i], indices[j])
                } else {
                    (indices[j], indices[i])
                };

                let refs = if lo == indices[i] {
                    [cis_left, cis_right]
                } else {
                    [cis_right, cis_left]
                };

                mol.add_ez_stereo(EZStereo { bond: (lo, hi), refs });
            }
        }
    }
}

fn other_substituent(
    mol: &Mol<Atom, SmilesBond>,
    db_atom: NodeIndex,
    db_partner: NodeIndex,
    known_ref: NodeIndex,
) -> AtomId {
    for neighbor in mol.neighbors(db_atom) {
        if neighbor != db_partner && neighbor != known_ref {
            return AtomId::Node(neighbor);
        }
    }
    AtomId::VirtualH(db_atom, 0)
}

fn find_directional_neighbor(
    tree: &ParseTree,
    db_atom: usize,
    other_db_atom: usize,
) -> Option<(usize, BondToken)> {
    for neighbor in &tree.atoms[db_atom].neighbors {
        if neighbor.atom_idx == other_db_atom {
            continue;
        }
        if let Some(bond) = &neighbor.bond {
            match bond {
                BondToken::Up | BondToken::Down => {
                    return Some((neighbor.atom_idx, *bond));
                }
                _ => {}
            }
        }
    }

    for neighbor in &tree.atoms[db_atom].neighbors {
        if neighbor.atom_idx == other_db_atom {
            continue;
        }
        let other_atom = &tree.atoms[neighbor.atom_idx];
        for n2 in &other_atom.neighbors {
            if n2.atom_idx == db_atom {
                if let Some(bond) = &n2.bond {
                    match bond {
                        BondToken::Up => return Some((neighbor.atom_idx, BondToken::Down)),
                        BondToken::Down => return Some((neighbor.atom_idx, BondToken::Up)),
                        _ => {}
                    }
                }
            }
        }
    }

    None
}

fn resolve_hydrogen_counts(
    mol: &mut Mol<Atom, SmilesBond>,
    tree: &ParseTree,
    indices: &[NodeIndex],
) {
    for (i, parse_atom) in tree.atoms.iter().enumerate() {
        let h_count = if parse_atom.is_bracket {
            parse_atom.hcount.unwrap_or(0)
        } else {
            compute_implicit_h(mol, indices[i], parse_atom)
        };
        mol.atom_mut(indices[i]).hydrogen_count = h_count;
    }
}

fn compute_implicit_h(
    mol: &Mol<Atom, SmilesBond>,
    node: NodeIndex,
    parse_atom: &ParseAtom,
) -> u8 {
    let valences = parse_atom.element.default_valences();
    if valences.is_empty() {
        return 0;
    }

    let bond_order_sum = bond_order_sum(mol, node);

    let charge = parse_atom.charge;
    let effective_valences = adjust_valences_for_charge(valences, parse_atom.element, charge);

    let target = effective_valences
        .iter()
        .find(|&&v| v >= bond_order_sum)
        .copied()
        .unwrap_or(0);

    if target < bond_order_sum {
        return 0;
    }

    let mut h = target - bond_order_sum;

    if parse_atom.is_aromatic && h > 0 {
        h -= 1;
    }

    h
}

fn bond_order_sum(mol: &Mol<Atom, SmilesBond>, node: NodeIndex) -> u8 {
    let mut sum: u8 = 0;
    for edge_idx in mol.bonds_of(node) {
        let order = match mol.bond(edge_idx).order {
            SmilesBondOrder::Single => 1,
            SmilesBondOrder::Double => 2,
            SmilesBondOrder::Triple => 3,
            SmilesBondOrder::Aromatic => 1,
            SmilesBondOrder::Implicit => 1,
        };
        sum = sum.saturating_add(order);
    }
    sum
}

fn adjust_valences_for_charge(
    valences: &[u8],
    _element: Element,
    charge: i8,
) -> Vec<u8> {
    debug_assert_eq!(charge, 0, "bare atoms in SMILES never carry charge");
    valences.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_tree::build_parse_tree;
    use crate::smiles::tokenizer::tokenize;

    fn parse(s: &str) -> Mol<Atom, SmilesBond> {
        let tokens = tokenize(s).unwrap();
        let tree = build_parse_tree(&tokens).unwrap();
        build_mol(&tree)
    }

    #[test]
    fn methane_h_count() {
        let mol = parse("C");
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.atom(NodeIndex::new(0)).hydrogen_count, 4);
    }

    #[test]
    fn ethane_h_counts() {
        let mol = parse("CC");
        assert_eq!(mol.atom(NodeIndex::new(0)).hydrogen_count, 3);
        assert_eq!(mol.atom(NodeIndex::new(1)).hydrogen_count, 3);
    }

    #[test]
    fn ethene_h_counts() {
        let mol = parse("C=C");
        assert_eq!(mol.atom(NodeIndex::new(0)).hydrogen_count, 2);
        assert_eq!(mol.atom(NodeIndex::new(1)).hydrogen_count, 2);
    }

    #[test]
    fn bracket_atom_h() {
        let mol = parse("[CH4]");
        assert_eq!(mol.atom(NodeIndex::new(0)).hydrogen_count, 4);
    }

    #[test]
    fn bracket_no_h() {
        let mol = parse("[C]");
        assert_eq!(mol.atom(NodeIndex::new(0)).hydrogen_count, 0);
    }

    #[test]
    fn aromatic_carbon_benzene() {
        let mol = parse("c1ccccc1");
        for i in 0..6 {
            let atom = mol.atom(NodeIndex::new(i));
            assert!(atom.is_aromatic);
            assert_eq!(atom.hydrogen_count, 1, "atom {} should have 1 H", i);
        }
    }
}
