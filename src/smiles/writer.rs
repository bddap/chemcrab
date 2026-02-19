use std::collections::HashMap;

use petgraph::graph::NodeIndex;

use crate::atom::{Atom, Chirality};
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::canonical::canonical_ordering;
use crate::element::Element;
use crate::graph_ops::connected_components;
use crate::mol::Mol;

pub fn to_smiles(mol: &Mol<Atom, Bond>) -> String {
    let components = connected_components(mol);
    let mut parts = Vec::with_capacity(components.len());
    for component in &components {
        parts.push(write_fragment(mol, component, None));
    }
    parts.join(".")
}

pub fn to_canonical_smiles(mol: &Mol<Atom, Bond>) -> String {
    let ranks = canonical_ordering(mol);
    let components = connected_components(mol);
    let mut parts = Vec::with_capacity(components.len());
    for component in &components {
        parts.push(write_fragment(mol, component, Some(&ranks)));
    }
    parts.join(".")
}

struct RingClosure {
    ring_id: usize,
    order: BondOrder,
    other: NodeIndex,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Up,
    Down,
}

impl Direction {
    fn flip(self) -> Self {
        match self {
            Direction::Up => Direction::Down,
            Direction::Down => Direction::Up,
        }
    }

    fn as_char(self) -> char {
        match self {
            Direction::Up => '/',
            Direction::Down => '\\',
        }
    }
}

fn compute_bond_directions(
    mol: &Mol<Atom, Bond>,
    parent: &[Option<NodeIndex>],
    children: &[Vec<NodeIndex>],
    ring_opens: &[Vec<RingClosure>],
    ring_closes: &[Vec<RingClosure>],
) -> HashMap<(NodeIndex, NodeIndex), Direction> {
    let mut dirs: HashMap<(NodeIndex, NodeIndex), Direction> = HashMap::new();

    for edge in mol.bonds() {
        let bond = mol.bond(edge);
        if bond.order != BondOrder::Double {
            continue;
        }
        let (ref_left, ref_right, is_trans) = match bond.stereo {
            BondStereo::Trans(l, r) => (l, r, true),
            BondStereo::Cis(l, r) => (l, r, false),
            BondStereo::None => continue,
        };

        let (left, right) = mol.bond_endpoints(edge).unwrap();

        let left_write_dir =
            write_direction(left, ref_left, parent, children, ring_opens, ring_closes);
        let right_write_dir =
            write_direction(right, ref_right, parent, children, ring_opens, ring_closes);

        let left_char = Direction::Up;
        let right_char = if is_trans { left_char } else { left_char.flip() };

        set_bond_dir(&mut dirs, left, ref_left, left_char, left_write_dir);
        set_bond_dir(&mut dirs, right, ref_right, right_char, right_write_dir);
    }

    dirs
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum WriteDir {
    ParentToChild,
    ChildToParent,
    RingOpen,
    RingClose,
}

fn write_direction(
    db_atom: NodeIndex,
    ref_atom: NodeIndex,
    parent: &[Option<NodeIndex>],
    children: &[Vec<NodeIndex>],
    ring_opens: &[Vec<RingClosure>],
    ring_closes: &[Vec<RingClosure>],
) -> WriteDir {
    if children[db_atom.index()].contains(&ref_atom) {
        return WriteDir::ParentToChild;
    }
    if parent[db_atom.index()] == Some(ref_atom) {
        return WriteDir::ChildToParent;
    }
    if ring_opens[db_atom.index()].iter().any(|rc| rc.other == ref_atom) {
        return WriteDir::RingOpen;
    }
    if ring_closes[db_atom.index()].iter().any(|rc| rc.other == ref_atom) {
        return WriteDir::RingClose;
    }
    WriteDir::ParentToChild
}

fn set_bond_dir(
    dirs: &mut HashMap<(NodeIndex, NodeIndex), Direction>,
    db_atom: NodeIndex,
    ref_atom: NodeIndex,
    smiles_char: Direction,
    write_dir: WriteDir,
) {
    match write_dir {
        WriteDir::ParentToChild => {
            dirs.insert((db_atom, ref_atom), smiles_char);
            dirs.insert((ref_atom, db_atom), smiles_char.flip());
        }
        WriteDir::ChildToParent => {
            dirs.insert((ref_atom, db_atom), smiles_char);
            dirs.insert((db_atom, ref_atom), smiles_char.flip());
        }
        WriteDir::RingOpen | WriteDir::RingClose => {
            dirs.insert((db_atom, ref_atom), smiles_char);
            dirs.insert((ref_atom, db_atom), smiles_char.flip());
        }
    }
}

fn write_fragment(mol: &Mol<Atom, Bond>, component: &[NodeIndex], ranks: Option<&[usize]>) -> String {
    let n = mol.atom_count();
    let start = match ranks {
        Some(r) => *component.iter().min_by_key(|&&node| r[node.index()]).unwrap(),
        None => component[0],
    };

    let mut visited = vec![false; n];
    let mut parent = vec![None::<NodeIndex>; n];
    let mut ring_opens: Vec<Vec<RingClosure>> = (0..n).map(|_| Vec::new()).collect();
    let mut ring_closes: Vec<Vec<RingClosure>> = (0..n).map(|_| Vec::new()).collect();
    let mut next_ring_id: usize = 1;
    let mut children: Vec<Vec<NodeIndex>> = (0..n).map(|_| Vec::new()).collect();

    let neighbor_lists: Vec<Vec<NodeIndex>> = (0..n)
        .map(|i| {
            let mut neighbors: Vec<NodeIndex> = mol.neighbors(NodeIndex::new(i)).collect();
            if let Some(r) = ranks {
                neighbors.sort_by_key(|nb| r[nb.index()]);
            }
            neighbors
        })
        .collect();

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

    let bond_dirs =
        compute_bond_directions(mol, &parent, &children, &ring_opens, &ring_closes);

    let ctx = DfsContext {
        parent,
        children,
        ring_opens,
        ring_closes,
        bond_dirs,
    };

    let mut out = String::new();
    write_node(mol, start, &ctx, &mut out);
    out
}

struct DfsContext {
    parent: Vec<Option<NodeIndex>>,
    children: Vec<Vec<NodeIndex>>,
    ring_opens: Vec<Vec<RingClosure>>,
    ring_closes: Vec<Vec<RingClosure>>,
    bond_dirs: HashMap<(NodeIndex, NodeIndex), Direction>,
}

impl DfsContext {
    fn smiles_neighbor_order(&self, node: NodeIndex) -> Vec<NodeIndex> {
        let mut order = Vec::new();
        if let Some(p) = self.parent[node.index()] {
            order.push(p);
        }
        for rc in &self.ring_opens[node.index()] {
            order.push(rc.other);
        }
        for rc in &self.ring_closes[node.index()] {
            order.push(rc.other);
        }
        for &child in &self.children[node.index()] {
            order.push(child);
        }
        order
    }
}

fn parity_of_permutation(from: &[NodeIndex], to: &[NodeIndex]) -> bool {
    if from.len() != to.len() {
        return true;
    }
    let n = from.len();
    let perm: Vec<usize> = from
        .iter()
        .map(|f| to.iter().position(|t| t == f).unwrap_or(0))
        .collect();

    let mut visited = vec![false; n];
    let mut swaps = 0;
    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut cycle_len = 0;
        let mut j = i;
        while !visited[j] {
            visited[j] = true;
            j = perm[j];
            cycle_len += 1;
        }
        swaps += cycle_len - 1;
    }
    swaps % 2 == 0
}

fn resolve_chirality_for_smiles(
    mol: &Mol<Atom, Bond>,
    node: NodeIndex,
    ctx: &DfsContext,
) -> Chirality {
    let atom = mol.atom(node);
    if atom.chirality == Chirality::None {
        return Chirality::None;
    }

    let graph_neighbors: Vec<NodeIndex> = mol.neighbors(node).collect();
    let smiles_neighbors = ctx.smiles_neighbor_order(node);

    let even = parity_of_permutation(&graph_neighbors, &smiles_neighbors);

    if even {
        atom.chirality
    } else {
        match atom.chirality {
            Chirality::Cw => Chirality::Ccw,
            Chirality::Ccw => Chirality::Cw,
            Chirality::None => Chirality::None,
        }
    }
}

fn write_node(
    mol: &Mol<Atom, Bond>,
    node: NodeIndex,
    ctx: &DfsContext,
    out: &mut String,
) {
    let chirality = resolve_chirality_for_smiles(mol, node, ctx);
    write_atom_symbol(mol, node, chirality, out);

    for rc in &ctx.ring_opens[node.index()] {
        write_bond_between(mol, rc.order, node, rc.other, &ctx.bond_dirs, out);
        write_ring_digit(rc.ring_id, out);
    }

    for rc in &ctx.ring_closes[node.index()] {
        write_bond_between(mol, rc.order, node, rc.other, &ctx.bond_dirs, out);
        write_ring_digit(rc.ring_id, out);
    }

    let kids = &ctx.children[node.index()];
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
        write_bond_between(mol, mol.bond(edge).order, node, child, &ctx.bond_dirs, out);
        write_node(mol, child, ctx, out);
        if is_branch {
            out.push(')');
        }
    }
}

fn write_bond_between(
    mol: &Mol<Atom, Bond>,
    order: BondOrder,
    from: NodeIndex,
    to: NodeIndex,
    bond_dirs: &HashMap<(NodeIndex, NodeIndex), Direction>,
    out: &mut String,
) {
    if let Some(&dir) = bond_dirs.get(&(from, to)) {
        out.push(dir.as_char());
        return;
    }
    let from_atom = mol.atom(from);
    let to_atom = mol.atom(to);
    if from_atom.is_aromatic && to_atom.is_aromatic {
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

fn write_atom_symbol(
    mol: &Mol<Atom, Bond>,
    node: NodeIndex,
    chirality: Chirality,
    out: &mut String,
) {
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
        write_bracket_atom(atom, elem, chirality, out);
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

fn write_bracket_atom(
    atom: &Atom,
    elem: Option<Element>,
    chirality: Chirality,
    out: &mut String,
) {
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

    match chirality {
        Chirality::Ccw => out.push('@'),
        Chirality::Cw => out.push_str("@@"),
        Chirality::None => {}
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
    use crate::bond::BondStereo;
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

    fn has_chiral_atom(mol: &Mol<Atom, Bond>) -> bool {
        mol.atoms().any(|n| mol.atom(n).chirality != Chirality::None)
    }

    fn find_double_bond_stereo(mol: &Mol<Atom, Bond>) -> Option<bool> {
        mol.bonds().find_map(|e| {
            match mol.bond(e).stereo {
                BondStereo::Trans(_, _) => Some(true),
                BondStereo::Cis(_, _) => Some(false),
                BondStereo::None => None,
            }
        })
    }

    fn canonical_chirality(mol: &Mol<Atom, Bond>, node: NodeIndex) -> Chirality {
        let atom = mol.atom(node);
        if atom.chirality == Chirality::None {
            return Chirality::None;
        }
        let graph_neighbors: Vec<NodeIndex> = mol.neighbors(node).collect();
        let mut sorted_neighbors = graph_neighbors.clone();
        sorted_neighbors.sort_by_key(|n| mol.atom(*n).atomic_num);

        let n = graph_neighbors.len();
        let perm: Vec<usize> = graph_neighbors
            .iter()
            .map(|g| sorted_neighbors.iter().position(|s| s == g).unwrap())
            .collect();
        let mut visited = vec![false; n];
        let mut swaps = 0;
        for i in 0..n {
            if visited[i] {
                continue;
            }
            let mut cycle_len = 0;
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
                cycle_len += 1;
            }
            swaps += cycle_len - 1;
        }
        if swaps % 2 == 0 {
            atom.chirality
        } else {
            match atom.chirality {
                Chirality::Cw => Chirality::Ccw,
                Chirality::Ccw => Chirality::Cw,
                Chirality::None => Chirality::None,
            }
        }
    }

    fn assert_same_chirality(mol1: &Mol<Atom, Bond>, mol2: &Mol<Atom, Bond>, ctx: &str) {
        for n1 in mol1.atoms() {
            if mol1.atom(n1).chirality == Chirality::None {
                continue;
            }
            let n2 = mol2
                .atoms()
                .find(|&n2| mol2.atom(n2).atomic_num == mol1.atom(n1).atomic_num
                    && mol2.atom(n2).chirality != Chirality::None)
                .unwrap_or_else(|| panic!("{ctx}: no matching chiral atom in mol2"));
            let c1 = canonical_chirality(mol1, n1);
            let c2 = canonical_chirality(mol2, n2);
            assert_eq!(c1, c2, "{ctx}: chirality mismatch for atom {}", n1.index());
        }
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

    #[test]
    fn tetrahedral_with_h() {
        let (m1, m2, s) = round_trip("[C@@H](F)(Cl)Br");
        assert_same_structure(&m1, &m2, "tetrahedral_with_h");
        assert!(s.contains('@'), "expected chirality in '{s}'");
        assert!(has_chiral_atom(&m1));
        assert!(has_chiral_atom(&m2));
        assert_same_chirality(&m1, &m2, "tetrahedral_with_h");
    }

    #[test]
    fn tetrahedral_no_h() {
        let (m1, m2, s) = round_trip("[C@](F)(Cl)(Br)I");
        assert_same_structure(&m1, &m2, "tetrahedral_no_h");
        assert!(s.contains('@'), "expected chirality in '{s}'");
        assert!(has_chiral_atom(&m1));
        assert!(has_chiral_atom(&m2));
        assert_same_chirality(&m1, &m2, "tetrahedral_no_h");
    }

    #[test]
    fn ez_trans_round_trip() {
        let (m1, m2, s) = round_trip("F/C=C/F");
        assert_same_structure(&m1, &m2, "ez_trans");
        assert!(
            s.contains('/') || s.contains('\\'),
            "expected directional bonds in '{s}'"
        );
        let is_trans1 = find_double_bond_stereo(&m1).expect("m1 should have E/Z stereo");
        let is_trans2 = find_double_bond_stereo(&m2).expect("m2 should have E/Z stereo");
        assert_eq!(is_trans1, is_trans2);
    }

    #[test]
    fn ez_cis_round_trip() {
        let (m1, m2, s) = round_trip(r"F/C=C\F");
        assert_same_structure(&m1, &m2, "ez_cis");
        assert!(
            s.contains('/') || s.contains('\\'),
            "expected directional bonds in '{s}'"
        );
        let is_trans1 = find_double_bond_stereo(&m1).expect("m1 should have E/Z stereo");
        let is_trans2 = find_double_bond_stereo(&m2).expect("m2 should have E/Z stereo");
        assert_eq!(is_trans1, is_trans2);
    }

    #[test]
    fn ez_trans_chlorine_round_trip() {
        let (m1, m2, _s) = round_trip("Cl/C=C/Cl");
        assert_same_structure(&m1, &m2, "ez_trans_cl");
        let is_trans1 = find_double_bond_stereo(&m1).expect("m1 should have E/Z stereo");
        let is_trans2 = find_double_bond_stereo(&m2).expect("m2 should have E/Z stereo");
        assert_eq!(is_trans1, is_trans2);
    }

    #[test]
    fn combined_tetrahedral_and_ez() {
        let (m1, m2, s) = round_trip(r"F/C=C/[C@@H](Cl)Br");
        assert_same_structure(&m1, &m2, "combined_stereo");
        assert!(
            s.contains('/') || s.contains('\\'),
            "expected directional bonds in '{s}'"
        );
        assert!(has_chiral_atom(&m1));
        assert!(has_chiral_atom(&m2));
        assert_same_chirality(&m1, &m2, "combined_stereo");
        let is_trans1 = find_double_bond_stereo(&m1).expect("m1 should have E/Z stereo");
        let is_trans2 = find_double_bond_stereo(&m2).expect("m2 should have E/Z stereo");
        assert_eq!(is_trans1, is_trans2);
    }

    fn canonical(smiles: &str) -> String {
        let mol = from_smiles(smiles).unwrap();
        to_canonical_smiles(&mol)
    }

    fn canonical_round_trip(smiles: &str) -> (Mol<Atom, Bond>, Mol<Atom, Bond>, String) {
        let mol1 = from_smiles(smiles).unwrap();
        let written = to_canonical_smiles(&mol1);
        let mol2 = from_smiles(&written).unwrap_or_else(|e| {
            panic!("Failed to re-parse canonical '{written}' (from '{smiles}'): {e}");
        });
        (mol1, mol2, written)
    }

    // -- Determinism: same molecule from different SMILES -> same canonical SMILES --

    #[test]
    fn canonical_ethanol_determinism() {
        assert_eq!(canonical("OCC"), canonical("CCO"));
    }

    #[test]
    fn canonical_benzene_determinism() {
        let c1 = canonical("c1ccccc1");
        let c2 = canonical("c1ccccc1");
        assert_eq!(c1, c2);
    }

    #[test]
    fn canonical_acetic_acid_determinism() {
        assert_eq!(canonical("C(=O)(O)CC"), canonical("CCC(O)=O"));
    }

    #[test]
    fn canonical_propanol_orderings() {
        let a = canonical("CCCO");
        let b = canonical("OCCC");
        let c = canonical("C(CC)O");
        assert_eq!(a, b);
        assert_eq!(b, c);
    }

    #[test]
    fn canonical_methanol_orderings() {
        assert_eq!(canonical("CO"), canonical("OC"));
    }

    #[test]
    fn canonical_phenol_orderings() {
        let a = canonical("Oc1ccccc1");
        let b = canonical("c1ccccc1O");
        let c = canonical("c1ccc(O)cc1");
        assert_eq!(a, b);
        assert_eq!(b, c);
    }

    // -- Uniqueness: different molecules -> different canonical SMILES --

    #[test]
    fn canonical_ethanol_vs_methanol() {
        assert_ne!(canonical("CCO"), canonical("CO"));
    }

    #[test]
    fn canonical_benzene_vs_cyclohexane() {
        assert_ne!(canonical("c1ccccc1"), canonical("C1CCCCC1"));
    }

    #[test]
    fn canonical_methane_vs_ethane() {
        assert_ne!(canonical("C"), canonical("CC"));
    }

    // -- Stereochemistry: canonical SMILES preserves stereo --

    #[test]
    fn canonical_tetrahedral_round_trip() {
        let (m1, m2, s) = canonical_round_trip("[C@@H](F)(Cl)Br");
        assert_same_structure(&m1, &m2, "canonical_tetrahedral");
        assert!(s.contains('@'), "expected chirality in canonical '{s}'");
        assert!(has_chiral_atom(&m2));
        assert_same_chirality(&m1, &m2, "canonical_tetrahedral");
    }

    #[test]
    fn canonical_ez_trans_round_trip() {
        let (m1, m2, s) = canonical_round_trip("F/C=C/F");
        assert_same_structure(&m1, &m2, "canonical_ez_trans");
        assert!(
            s.contains('/') || s.contains('\\'),
            "expected stereo bonds in canonical '{s}'"
        );
        let t1 = find_double_bond_stereo(&m1).expect("m1 E/Z");
        let t2 = find_double_bond_stereo(&m2).expect("m2 E/Z");
        assert_eq!(t1, t2);
    }

    #[test]
    fn canonical_ez_cis_round_trip() {
        let (m1, m2, s) = canonical_round_trip(r"F/C=C\F");
        assert_same_structure(&m1, &m2, "canonical_ez_cis");
        assert!(
            s.contains('/') || s.contains('\\'),
            "expected stereo bonds in canonical '{s}'"
        );
        let t1 = find_double_bond_stereo(&m1).expect("m1 E/Z");
        let t2 = find_double_bond_stereo(&m2).expect("m2 E/Z");
        assert_eq!(t1, t2);
    }

    // -- Edge cases --

    #[test]
    fn canonical_single_atom_c() {
        let s = canonical("C");
        assert_eq!(s, "C");
    }

    #[test]
    fn canonical_single_atom_na() {
        let s = canonical("[Na+]");
        assert!(s.contains("Na"));
        assert!(s.contains('+'));
    }

    #[test]
    fn canonical_disconnected() {
        let s = canonical("[Na+].[Cl-]");
        assert!(s.contains('.'));
        let mol = from_smiles(&s).unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn canonical_symmetric_benzene() {
        let s = canonical("c1ccccc1");
        let mol = from_smiles(&s).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn canonical_charged_atoms() {
        let s = canonical("[NH4+]");
        let mol = from_smiles(&s).unwrap();
        assert_eq!(mol.atom(NodeIndex::new(0)).formal_charge, 1);
    }

    #[test]
    fn canonical_isotope() {
        let s = canonical("[13C]");
        assert!(s.contains("13"));
        let mol = from_smiles(&s).unwrap();
        assert_eq!(mol.atom(NodeIndex::new(0)).isotope, 13);
    }

    #[test]
    fn canonical_empty_mol() {
        let mol = Mol::<Atom, Bond>::new();
        assert_eq!(to_canonical_smiles(&mol), "");
    }

    #[test]
    fn canonical_round_trip_structure_preservation() {
        let cases = [
            "CC", "C=C", "C#C", "C1CCCCC1", "c1ccccc1", "CC(=O)O",
            "c1ccncc1", "c1ccc2ccccc2c1", "[Fe]", "[Na+].[Cl-]",
        ];
        for smiles in &cases {
            let (m1, m2, _) = canonical_round_trip(smiles);
            assert_same_structure(&m1, &m2, smiles);
        }
    }

    #[test]
    fn canonical_idempotent() {
        let cases = [
            "CCO", "c1ccccc1", "CC(=O)O", "[Na+].[Cl-]", "c1ccncc1",
        ];
        for smiles in &cases {
            let first = canonical(smiles);
            let second = canonical(&first);
            assert_eq!(first, second, "canonical not idempotent for {smiles}");
        }
    }

    #[test]
    fn canonical_three_fragments_determinism() {
        let s = canonical("[Na+].[Cl-].O");
        assert_eq!(s.matches('.').count(), 2);
        let mol = from_smiles(&s).unwrap();
        assert_eq!(mol.atom_count(), 3);
    }

    #[test]
    fn canonical_naphthalene_determinism() {
        let a = canonical("c1ccc2ccccc2c1");
        let b = canonical("c1cccc2ccccc12");
        assert_eq!(a, b);
    }

    #[test]
    fn canonical_toluene_orderings() {
        let a = canonical("Cc1ccccc1");
        let b = canonical("c1ccccc1C");
        let c = canonical("c1ccc(C)cc1");
        assert_eq!(a, b);
        assert_eq!(b, c);
    }
}
