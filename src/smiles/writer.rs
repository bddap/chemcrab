use petgraph::graph::NodeIndex;

use crate::atom::{Atom, Chirality};
use crate::bond::{Bond, BondOrder};
use crate::element::Element;
use crate::graph_ops::connected_components;
use crate::mol::Mol;

pub fn to_smiles(mol: &Mol<Atom, Bond>) -> String {
    let components = connected_components(mol);
    let mut parts = Vec::with_capacity(components.len());
    for component in &components {
        parts.push(write_fragment(mol, component));
    }
    parts.join(".")
}

struct RingClosure {
    ring_id: usize,
    order: BondOrder,
    other: NodeIndex,
}

fn write_fragment(mol: &Mol<Atom, Bond>, component: &[NodeIndex]) -> String {
    let n = mol.atom_count();
    let start = component[0];

    let mut visited = vec![false; n];
    let mut parent = vec![None::<NodeIndex>; n];
    let mut ring_opens: Vec<Vec<RingClosure>> = (0..n).map(|_| Vec::new()).collect();
    let mut ring_closes: Vec<Vec<RingClosure>> = (0..n).map(|_| Vec::new()).collect();
    let mut next_ring_id: usize = 1;
    let mut children: Vec<Vec<NodeIndex>> = (0..n).map(|_| Vec::new()).collect();

    let neighbor_lists: Vec<Vec<NodeIndex>> =
        (0..n).map(|i| mol.neighbors(NodeIndex::new(i)).collect()).collect();

    let mut stack: Vec<(NodeIndex, usize)> = Vec::new();
    visited[start.index()] = true;
    stack.push((start, 0));

    loop {
        let Some(&mut (node, ref mut ni)) = stack.last_mut() else {
            break;
        };
        let neighbors = &neighbor_lists[node.index()];
        if *ni >= neighbors.len() {
            stack.pop();
            continue;
        }
        let neighbor = neighbors[*ni];
        *ni += 1;

        if !visited[neighbor.index()] {
            visited[neighbor.index()] = true;
            parent[neighbor.index()] = Some(node);
            children[node.index()].push(neighbor);
            stack.push((neighbor, 0));
        } else if parent[node.index()] != Some(neighbor) {
            let already = ring_opens[neighbor.index()]
                .iter()
                .any(|rc| ring_closes[node.index()].iter().any(|rc2| rc2.ring_id == rc.ring_id))
                || ring_opens[node.index()]
                    .iter()
                    .any(|rc| {
                        ring_closes[neighbor.index()]
                            .iter()
                            .any(|rc2| rc2.ring_id == rc.ring_id)
                    });
            if !already {
                let edge = mol.bond_between(node, neighbor).unwrap();
                let order = mol.bond(edge).order;
                let ring_id = next_ring_id;
                next_ring_id += 1;
                ring_opens[neighbor.index()].push(RingClosure {
                    ring_id,
                    order,
                    other: node,
                });
                ring_closes[node.index()].push(RingClosure {
                    ring_id,
                    order,
                    other: neighbor,
                });
            }
        }
    }

    let mut out = String::new();
    write_node(mol, start, &children, &ring_opens, &ring_closes, &mut out);
    out
}

fn write_node(
    mol: &Mol<Atom, Bond>,
    node: NodeIndex,
    children: &[Vec<NodeIndex>],
    ring_opens: &[Vec<RingClosure>],
    ring_closes: &[Vec<RingClosure>],
    out: &mut String,
) {
    let atom = mol.atom(node);
    write_atom_symbol(mol, node, out);

    for rc in &ring_opens[node.index()] {
        write_bond_symbol(rc.order, atom, mol.atom(rc.other), out);
        write_ring_digit(rc.ring_id, out);
    }

    for rc in &ring_closes[node.index()] {
        write_bond_symbol(rc.order, atom, mol.atom(rc.other), out);
        write_ring_digit(rc.ring_id, out);
    }

    let kids = &children[node.index()];
    if kids.is_empty() {
        return;
    }

    let last = kids.len() - 1;
    for (i, &child) in kids.iter().enumerate() {
        let is_branch = i < last;
        if is_branch {
            out.push('(');
        }
        let edge = mol.bond_between(node, child).unwrap();
        write_bond_symbol(mol.bond(edge).order, atom, mol.atom(child), out);
        write_node(mol, child, children, ring_opens, ring_closes, out);
        if is_branch {
            out.push(')');
        }
    }
}

fn write_bond_symbol(order: BondOrder, from: &Atom, to: &Atom, out: &mut String) {
    if from.is_aromatic && to.is_aromatic {
        return;
    }
    match order {
        BondOrder::Single => {}
        BondOrder::Double => out.push('='),
        BondOrder::Triple => out.push('#'),
    }
}

fn write_ring_digit(id: usize, out: &mut String) {
    assert!(id <= 99, "ring id {id} exceeds SMILES maximum of 99");
    if id <= 9 {
        out.push(char::from(b'0' + id as u8));
    } else {
        out.push('%');
        out.push(char::from(b'0' + (id / 10) as u8));
        out.push(char::from(b'0' + (id % 10) as u8));
    }
}

fn write_atom_symbol(mol: &Mol<Atom, Bond>, node: NodeIndex, out: &mut String) {
    let atom = mol.atom(node);
    let elem = Element::from_atomic_num(atom.atomic_num);

    if can_write_bare(mol, node) {
        let symbol = elem.unwrap().symbol();
        if atom.is_aromatic {
            for c in symbol.chars() {
                out.push(c.to_ascii_lowercase());
            }
        } else {
            out.push_str(symbol);
        }
    } else {
        write_bracket_atom(atom, elem, out);
    }
}

fn can_write_bare(mol: &Mol<Atom, Bond>, node: NodeIndex) -> bool {
    let atom = mol.atom(node);

    let elem = match Element::from_atomic_num(atom.atomic_num) {
        Some(e) => e,
        None => return false,
    };

    if !elem.is_organic_subset() {
        return false;
    }
    if atom.isotope != 0 || atom.formal_charge != 0 || atom.chirality != Chirality::None {
        return false;
    }

    let expected_h =
        implicit_h_for_bare_atom(elem, atom.is_aromatic, reader_bond_order_sum(mol, node));
    atom.hydrogen_count == expected_h
}

/// Mirrors `compute_implicit_h` in the parser's builder: given what the reader
/// would see for a bare atom, compute the implicit hydrogen count.
fn implicit_h_for_bare_atom(elem: Element, is_aromatic: bool, bos: u8) -> u8 {
    let valences = elem.default_valences();
    if valences.is_empty() {
        return 0;
    }
    let target = valences.iter().find(|&&v| v >= bos).copied().unwrap_or(0);
    if target < bos {
        return 0;
    }
    let mut h = target - bos;
    if is_aromatic && h > 0 {
        h -= 1;
    }
    h
}

/// Bond order sum as the reader would compute it for a bare atom: bonds between
/// two aromatic atoms count as 1 (aromatic/implicit), others use actual order.
fn reader_bond_order_sum(mol: &Mol<Atom, Bond>, node: NodeIndex) -> u8 {
    let atom = mol.atom(node);
    let mut sum: u8 = 0;
    for edge_idx in mol.bonds_of(node) {
        let (a, b) = mol.bond_endpoints(edge_idx).unwrap();
        let neighbor = mol.atom(if a == node { b } else { a });

        let contribution = if atom.is_aromatic && neighbor.is_aromatic {
            1
        } else {
            match mol.bond(edge_idx).order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
            }
        };
        sum = sum.saturating_add(contribution);
    }
    sum
}

fn write_bracket_atom(atom: &Atom, elem: Option<Element>, out: &mut String) {
    out.push('[');

    if atom.isotope != 0 {
        out.push_str(&atom.isotope.to_string());
    }

    match elem {
        Some(e) => {
            let symbol = e.symbol();
            if atom.is_aromatic {
                for c in symbol.chars() {
                    out.push(c.to_ascii_lowercase());
                }
            } else {
                out.push_str(symbol);
            }
        }
        None => out.push('*'),
    }

    if atom.hydrogen_count > 0 {
        out.push('H');
        if atom.hydrogen_count > 1 {
            out.push_str(&atom.hydrogen_count.to_string());
        }
    }

    if atom.formal_charge > 0 {
        out.push('+');
        if atom.formal_charge > 1 {
            out.push_str(&atom.formal_charge.to_string());
        }
    } else if atom.formal_charge < 0 {
        out.push('-');
        if atom.formal_charge < -1 {
            out.push_str(&atom.formal_charge.abs().to_string());
        }
    }

    out.push(']');
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    fn round_trip(smiles: &str) -> (Mol<Atom, Bond>, Mol<Atom, Bond>, String) {
        let mol1 = from_smiles(smiles).unwrap();
        let written = to_smiles(&mol1);
        let mol2 = from_smiles(&written).unwrap_or_else(|e| {
            panic!("Failed to re-parse '{written}' (from '{smiles}'): {e}");
        });
        (mol1, mol2, written)
    }

    fn assert_same_structure(mol1: &Mol<Atom, Bond>, mol2: &Mol<Atom, Bond>, ctx: &str) {
        assert_eq!(mol1.atom_count(), mol2.atom_count(), "{ctx}: atom count");
        assert_eq!(mol1.bond_count(), mol2.bond_count(), "{ctx}: bond count");

        let mut e1: Vec<u8> = mol1.atoms().map(|n| mol1.atom(n).atomic_num).collect();
        let mut e2: Vec<u8> = mol2.atoms().map(|n| mol2.atom(n).atomic_num).collect();
        e1.sort();
        e2.sort();
        assert_eq!(e1, e2, "{ctx}: elements");
    }

    #[test]
    fn methane() {
        let (m1, m2, s) = round_trip("C");
        assert_eq!(s, "C");
        assert_same_structure(&m1, &m2, "methane");
    }

    #[test]
    fn ethane() {
        let (m1, m2, s) = round_trip("CC");
        assert_eq!(s, "CC");
        assert_same_structure(&m1, &m2, "ethane");
    }

    #[test]
    fn ethene() {
        let (m1, m2, s) = round_trip("C=C");
        assert!(s.contains('='));
        assert_same_structure(&m1, &m2, "ethene");
    }

    #[test]
    fn cyclohexane() {
        let (m1, m2, s) = round_trip("C1CCCCC1");
        assert!(s.contains('1'));
        assert_same_structure(&m1, &m2, "cyclohexane");
    }

    #[test]
    fn benzene() {
        let (m1, m2, s) = round_trip("c1ccccc1");
        assert_same_structure(&m1, &m2, "benzene");
        assert!(s.contains('c'), "expected lowercase aromatic: {s}");
    }

    #[test]
    fn water() {
        let (m1, m2, s) = round_trip("O");
        assert_eq!(s, "O");
        assert_same_structure(&m1, &m2, "water");
    }

    #[test]
    fn methanol() {
        let (m1, m2, _) = round_trip("CO");
        assert_same_structure(&m1, &m2, "methanol");
    }

    #[test]
    fn acetic_acid() {
        let (m1, m2, s) = round_trip("CC(=O)O");
        assert_same_structure(&m1, &m2, "acetic acid");
        assert!(s.contains('='));
    }

    #[test]
    fn sodium_chloride() {
        let (m1, m2, s) = round_trip("[Na+].[Cl-]");
        assert_same_structure(&m1, &m2, "NaCl");
        assert!(s.contains('.'));
    }

    #[test]
    fn iron() {
        let (m1, m2, s) = round_trip("[Fe]");
        assert_eq!(s, "[Fe]");
        assert_same_structure(&m1, &m2, "iron");
    }

    #[test]
    fn carbon_13() {
        let (m1, m2, s) = round_trip("[13C]");
        assert!(s.contains("13"));
        assert_same_structure(&m1, &m2, "13C");
    }

    #[test]
    fn ammonium() {
        let (m1, m2, s) = round_trip("[NH4+]");
        assert!(s.contains('+'));
        assert_same_structure(&m1, &m2, "NH4+");
    }

    #[test]
    fn pyridine() {
        let (m1, m2, _) = round_trip("c1ccncc1");
        assert_same_structure(&m1, &m2, "pyridine");
    }

    #[test]
    fn naphthalene() {
        let (m1, m2, _) = round_trip("c1ccc2ccccc2c1");
        assert_same_structure(&m1, &m2, "naphthalene");
    }

    #[test]
    fn empty_mol() {
        let mol = Mol::<Atom, Bond>::new();
        assert_eq!(to_smiles(&mol), "");
    }

    #[test]
    fn triple_bond() {
        let (m1, m2, s) = round_trip("C#C");
        assert!(s.contains('#'));
        assert_same_structure(&m1, &m2, "ethyne");
    }

    #[test]
    fn phenol() {
        let (m1, m2, _) = round_trip("Oc1ccccc1");
        assert_same_structure(&m1, &m2, "phenol");
    }

    #[test]
    fn three_fragments() {
        let (m1, m2, s) = round_trip("[Na+].[Cl-].O");
        assert_same_structure(&m1, &m2, "three fragments");
        assert_eq!(s.matches('.').count(), 2);
    }
}
