use petgraph::graph::NodeIndex;

use crate::element::Element;
use crate::mol::Mol;

use super::query::{AtomExpr, BondExpr};

pub fn to_smarts(mol: &Mol<AtomExpr, BondExpr>) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }

    let components = find_components(mol);
    let mut parts = Vec::new();
    for component in &components {
        parts.push(write_component(mol, component));
    }
    parts.join(".")
}

fn find_components(mol: &Mol<AtomExpr, BondExpr>) -> Vec<Vec<NodeIndex>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components = Vec::new();

    for start in mol.atoms() {
        if visited[start.index()] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        visited[start.index()] = true;
        while let Some(node) = stack.pop() {
            component.push(node);
            for nb in mol.neighbors(node) {
                if !visited[nb.index()] {
                    visited[nb.index()] = true;
                    stack.push(nb);
                }
            }
        }
        components.push(component);
    }
    components
}

fn write_component(mol: &Mol<AtomExpr, BondExpr>, component: &[NodeIndex]) -> String {
    let n = mol.atom_count();
    let start = component[0];

    let mut visited = vec![false; n];
    let mut parent = vec![None::<NodeIndex>; n];
    let mut children: Vec<Vec<NodeIndex>> = (0..n).map(|_| Vec::new()).collect();
    let mut ring_opens: Vec<Vec<(usize, NodeIndex)>> = (0..n).map(|_| Vec::new()).collect();
    let mut ring_closes: Vec<Vec<(usize, NodeIndex)>> = (0..n).map(|_| Vec::new()).collect();
    let mut next_ring_id: usize = 1;

    let mut stack: Vec<(NodeIndex, usize)> = Vec::new();
    let neighbor_lists: Vec<Vec<NodeIndex>> = (0..n)
        .map(|i| mol.neighbors(NodeIndex::new(i)).collect())
        .collect();

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
                .any(|(rid, _)| ring_closes[node.index()].iter().any(|(rid2, _)| rid2 == rid))
                || ring_opens[node.index()]
                    .iter()
                    .any(|(rid, _)| ring_closes[neighbor.index()].iter().any(|(rid2, _)| rid2 == rid));
            if !already {
                let ring_id = next_ring_id;
                next_ring_id += 1;
                ring_opens[neighbor.index()].push((ring_id, node));
                ring_closes[node.index()].push((ring_id, neighbor));
            }
        }
    }

    let mut out = String::new();
    write_node(mol, start, &children, &ring_opens, &ring_closes, &mut out);
    out
}

fn write_node(
    mol: &Mol<AtomExpr, BondExpr>,
    node: NodeIndex,
    children: &[Vec<NodeIndex>],
    ring_opens: &[Vec<(usize, NodeIndex)>],
    ring_closes: &[Vec<(usize, NodeIndex)>],
    out: &mut String,
) {
    write_atom_expr(mol.atom(node), out);

    for &(ring_id, other) in &ring_opens[node.index()] {
        if let Some(edge) = mol.bond_between(node, other) {
            write_bond_symbol(mol.bond(edge), out);
        }
        write_ring_digit(ring_id, out);
    }

    for &(ring_id, other) in &ring_closes[node.index()] {
        if let Some(edge) = mol.bond_between(node, other) {
            write_bond_symbol(mol.bond(edge), out);
        }
        write_ring_digit(ring_id, out);
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
        if let Some(edge) = mol.bond_between(node, child) {
            write_bond_symbol(mol.bond(edge), out);
        }
        write_node(mol, child, children, ring_opens, ring_closes, out);
        if is_branch {
            out.push(')');
        }
    }
}

fn write_bond_symbol(bond: &BondExpr, out: &mut String) {
    match bond {
        BondExpr::SingleOrAromatic => {}
        BondExpr::True => out.push('~'),
        BondExpr::Single => out.push('-'),
        BondExpr::Double => out.push('='),
        BondExpr::Triple => out.push('#'),
        BondExpr::Aromatic => out.push(':'),
        BondExpr::Ring => out.push('@'),
        BondExpr::Up => out.push('/'),
        BondExpr::Down => out.push('\\'),
        BondExpr::Not(inner) => {
            out.push('!');
            write_bond_symbol(inner, out);
        }
        BondExpr::And(exprs) => {
            for (i, e) in exprs.iter().enumerate() {
                if i > 0 {
                    out.push('&');
                }
                write_bond_symbol(e, out);
            }
        }
        BondExpr::Or(exprs) => {
            for (i, e) in exprs.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_bond_symbol(e, out);
            }
        }
    }
}

fn write_ring_digit(id: usize, out: &mut String) {
    if id <= 9 {
        out.push(char::from(b'0' + id as u8));
    } else {
        out.push('%');
        out.push(char::from(b'0' + (id / 10) as u8));
        out.push(char::from(b'0' + (id % 10) as u8));
    }
}

fn write_atom_expr(expr: &AtomExpr, out: &mut String) {
    match expr {
        AtomExpr::True => out.push('*'),
        AtomExpr::Aromatic => {
            out.push('[');
            out.push('a');
            out.push(']');
        }
        AtomExpr::Aliphatic => {
            out.push('[');
            out.push('A');
            out.push(']');
        }
        AtomExpr::Element { atomic_num, aromatic } => {
            if can_write_bare(*atomic_num, *aromatic) {
                write_bare_element(*atomic_num, *aromatic, out);
            } else {
                out.push('[');
                write_atom_expr_inner(expr, out);
                out.push(']');
            }
        }
        _ => {
            out.push('[');
            write_atom_expr_inner(expr, out);
            out.push(']');
        }
    }
}

fn can_write_bare(atomic_num: u8, aromatic: Option<bool>) -> bool {
    match aromatic {
        Some(is_arom) => {
            if let Some(elem) = Element::from_atomic_num(atomic_num) {
                if is_arom {
                    matches!(atomic_num, 6 | 7 | 8 | 15 | 16)
                } else {
                    elem.is_organic_subset() || atomic_num == 1
                }
            } else {
                false
            }
        }
        None => false,
    }
}

fn write_bare_element(atomic_num: u8, aromatic: Option<bool>, out: &mut String) {
    if let Some(elem) = Element::from_atomic_num(atomic_num) {
        let sym = elem.symbol();
        if aromatic == Some(true) {
            for c in sym.chars() {
                out.push(c.to_ascii_lowercase());
            }
        } else {
            out.push_str(sym);
        }
    }
}

fn write_element_symbol(atomic_num: u8, aromatic: Option<bool>, out: &mut String) {
    if let Some(elem) = Element::from_atomic_num(atomic_num) {
        let sym = elem.symbol();
        if aromatic == Some(true) {
            for c in sym.chars() {
                out.push(c.to_ascii_lowercase());
            }
        } else {
            out.push_str(sym);
        }
    }
}

fn write_atom_expr_inner(expr: &AtomExpr, out: &mut String) {
    match expr {
        AtomExpr::True => out.push('*'),
        AtomExpr::Element { atomic_num, aromatic } => match aromatic {
            Some(is_arom) => write_element_symbol(*atomic_num, Some(*is_arom), out),
            None => {
                out.push('#');
                out.push_str(&atomic_num.to_string());
            }
        },
        AtomExpr::Aromatic => out.push('a'),
        AtomExpr::Aliphatic => out.push('A'),
        AtomExpr::Isotope(iso) => {
            out.push_str(&iso.to_string());
        }
        AtomExpr::Degree(d) => {
            out.push('D');
            out.push_str(&d.to_string());
        }
        AtomExpr::Valence(v) => {
            out.push('v');
            out.push_str(&v.to_string());
        }
        AtomExpr::Connectivity(x) => {
            out.push('X');
            out.push_str(&x.to_string());
        }
        AtomExpr::TotalHCount(h) => {
            out.push('H');
            out.push_str(&h.to_string());
        }
        AtomExpr::ImplicitHCount(h) => {
            out.push('h');
            out.push_str(&h.to_string());
        }
        AtomExpr::RingMembership(n) => {
            out.push('R');
            out.push_str(&n.to_string());
        }
        AtomExpr::SmallestRingSize(r) => {
            out.push('r');
            out.push_str(&r.to_string());
        }
        AtomExpr::RingBondCount(x) => {
            out.push('x');
            out.push_str(&x.to_string());
        }
        AtomExpr::Chirality(chiral) => {
            use crate::atom::Chirality;
            match chiral {
                Chirality::Ccw => out.push('@'),
                Chirality::Cw => out.push_str("@@"),
                Chirality::None => {}
            }
        }
        AtomExpr::Charge(c) => {
            if *c >= 0 {
                out.push('+');
                out.push_str(&c.to_string());
            } else {
                out.push('-');
                out.push_str(&c.abs().to_string());
            }
        }
        AtomExpr::InRing => out.push('R'),
        AtomExpr::NotInRing => {
            out.push('R');
            out.push('0');
        }
        AtomExpr::Recursive(inner) => {
            out.push_str("$(");
            out.push_str(&super::writer::to_smarts(inner));
            out.push(')');
        }
        AtomExpr::And(exprs) => {
            for (i, e) in exprs.iter().enumerate() {
                if i > 0 {
                    out.push('&');
                }
                write_atom_expr_inner(e, out);
            }
        }
        AtomExpr::Or(exprs) => {
            for (i, e) in exprs.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_atom_expr_inner(e, out);
            }
        }
        AtomExpr::Not(inner) => {
            out.push('!');
            write_atom_expr_inner(inner, out);
        }
    }
}
