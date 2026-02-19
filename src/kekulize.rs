use std::collections::VecDeque;
use std::fmt;

use petgraph::graph::{EdgeIndex, NodeIndex};

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, SmilesBond, SmilesBondOrder};
use crate::element::Element;
use crate::mol::Mol;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KekulizeError {
    Unkekulizable(Vec<NodeIndex>),
}

impl fmt::Display for KekulizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unkekulizable(atoms) => {
                write!(f, "cannot kekulize aromatic system: unmatched atoms [")?;
                for (i, idx) in atoms.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", idx.index())?;
                }
                write!(f, "]")
            }
        }
    }
}

impl std::error::Error for KekulizeError {}

pub fn kekulize(mol: Mol<Atom, SmilesBond>) -> Result<Mol<Atom, Bond>, KekulizeError> {
    let aromatic_edges: Vec<EdgeIndex> = mol
        .bonds()
        .filter(|&e| mol.bond(e).order == SmilesBondOrder::Aromatic)
        .collect();

    let n = mol.atom_count();

    let mut aromatic_adj: Vec<Vec<(NodeIndex, EdgeIndex)>> = vec![vec![]; n];
    for &e in &aromatic_edges {
        if let Some((a, b)) = mol.bond_endpoints(e) {
            aromatic_adj[a.index()].push((b, e));
            aromatic_adj[b.index()].push((a, e));
        }
    }

    let mut component_id: Vec<Option<usize>> = vec![None; n];
    let mut components: Vec<Vec<NodeIndex>> = Vec::new();
    for node in mol.atoms() {
        if aromatic_adj[node.index()].is_empty() || component_id[node.index()].is_some() {
            continue;
        }
        let cid = components.len();
        let mut stack = vec![node];
        let mut comp = Vec::new();
        while let Some(v) = stack.pop() {
            if component_id[v.index()].is_some() {
                continue;
            }
            component_id[v.index()] = Some(cid);
            comp.push(v);
            for &(w, _) in &aromatic_adj[v.index()] {
                if component_id[w.index()].is_none() {
                    stack.push(w);
                }
            }
        }
        components.push(comp);
    }

    let mut needs_double = vec![false; n];
    for comp in &components {
        for &node in comp {
            let atom = mol.atom(node);
            let elem = match Element::from_atomic_num(atom.atomic_num) {
                Some(e) => e,
                None => continue,
            };

            let bond_order_sum: u8 = mol
                .bonds_of(node)
                .map(|e| match mol.bond(e).order {
                    SmilesBondOrder::Single => 1,
                    SmilesBondOrder::Double => 2,
                    SmilesBondOrder::Triple => 3,
                    SmilesBondOrder::Aromatic => 1,
                    SmilesBondOrder::Implicit => 1,
                })
                .sum();

            let total_used = bond_order_sum + atom.hydrogen_count;

            let target = target_valence(elem, total_used);
            if let Some(tv) = target {
                if tv > total_used && tv - total_used == 1 {
                    needs_double[node.index()] = true;
                }
            }
        }
    }

    let mut matched_edge: Vec<Option<EdgeIndex>> = vec![None; n];

    for comp in &components {
        let candidates: Vec<NodeIndex> = comp
            .iter()
            .copied()
            .filter(|&v| needs_double[v.index()])
            .collect();

        for &start in &candidates {
            if matched_edge[start.index()].is_some() {
                continue;
            }
            augment(&mol, &aromatic_adj, &needs_double, &mut matched_edge, start);
        }

        let unmatched: Vec<NodeIndex> = candidates
            .iter()
            .copied()
            .filter(|&v| matched_edge[v.index()].is_none())
            .collect();

        if !unmatched.is_empty() {
            return Err(KekulizeError::Unkekulizable(unmatched));
        }
    }

    let mut result = Mol::new();
    let mut node_map: Vec<NodeIndex> = Vec::with_capacity(n);
    for node in mol.atoms() {
        let atom = mol.atom(node).clone();
        node_map.push(result.add_atom(atom));
    }

    let matched_edges: std::collections::HashSet<EdgeIndex> = matched_edge
        .iter()
        .filter_map(|e| *e)
        .collect();

    for edge in mol.bonds() {
        let (a, b) = match mol.bond_endpoints(edge) {
            Some(pair) => pair,
            None => continue,
        };
        let smiles_bond = mol.bond(edge);
        let order = match smiles_bond.order {
            SmilesBondOrder::Aromatic => {
                if matched_edges.contains(&edge) {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                }
            }
            SmilesBondOrder::Implicit | SmilesBondOrder::Single => BondOrder::Single,
            SmilesBondOrder::Double => BondOrder::Double,
            SmilesBondOrder::Triple => BondOrder::Triple,
        };
        result.add_bond(
            node_map[a.index()],
            node_map[b.index()],
            Bond {
                order,
                stereo: smiles_bond.stereo,
            },
        );
    }

    result.set_tetrahedral_stereo(mol.tetrahedral_stereo().to_vec());

    Ok(result)
}

fn target_valence(elem: Element, current_used: u8) -> Option<u8> {
    let valences = elem.default_valences();
    if valences.is_empty() {
        return None;
    }
    valences.iter().copied().find(|&v| v >= current_used)
}

fn augment(
    mol: &Mol<Atom, SmilesBond>,
    aromatic_adj: &[Vec<(NodeIndex, EdgeIndex)>],
    needs_double: &[bool],
    matched_edge: &mut [Option<EdgeIndex>],
    start: NodeIndex,
) -> bool {
    let n = mol.atom_count();
    let mut prev: Vec<Option<(NodeIndex, EdgeIndex)>> = vec![None; n];
    let mut visited = vec![false; n];
    let mut queue = VecDeque::new();

    visited[start.index()] = true;
    queue.push_back(start);

    while let Some(u) = queue.pop_front() {
        for &(v, e) in &aromatic_adj[u.index()] {
            if !needs_double[v.index()] || visited[v.index()] {
                continue;
            }
            if Some(e) == matched_edge[u.index()] {
                continue;
            }
            visited[v.index()] = true;
            prev[v.index()] = Some((u, e));

            if matched_edge[v.index()].is_none() {
                flip_path(mol, matched_edge, &prev, start, v);
                return true;
            }

            let matched_e = matched_edge[v.index()].expect("checked above");
            let (ea, eb) = mol.bond_endpoints(matched_e).expect("valid edge");
            let w = if ea == v { eb } else { ea };

            if !visited[w.index()] {
                visited[w.index()] = true;
                prev[w.index()] = Some((v, matched_e));
                queue.push_back(w);
            }
        }
    }
    false
}

fn flip_path(
    _mol: &Mol<Atom, SmilesBond>,
    matched_edge: &mut [Option<EdgeIndex>],
    prev: &[Option<(NodeIndex, EdgeIndex)>],
    start: NodeIndex,
    end: NodeIndex,
) {
    let mut cur = end;
    let mut is_new_match = true;
    while cur != start {
        let (p, e) = prev[cur.index()].expect("path exists");
        if is_new_match {
            matched_edge[cur.index()] = Some(e);
            matched_edge[p.index()] = Some(e);
        }
        is_new_match = !is_new_match;
        cur = p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bond::BondStereo;
    use crate::smiles::parse_smiles;

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    fn count_double_bonds(mol: &Mol<Atom, Bond>) -> usize {
        mol.bonds()
            .filter(|&e| mol.bond(e).order == BondOrder::Double)
            .count()
    }

    fn is_valid_kekulization(mol: &Mol<Atom, Bond>) -> bool {
        for node in mol.atoms() {
            let double_count = mol
                .bonds_of(node)
                .filter(|&e| mol.bond(e).order == BondOrder::Double)
                .count();
            if double_count > 1 && mol.atom(node).is_aromatic {
                return false;
            }
        }
        true
    }

    #[test]
    fn benzene() {
        let smiles_mol = parse_smiles("c1ccccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        assert_eq!(count_double_bonds(&mol), 3);
        assert!(is_valid_kekulization(&mol));
        for node in mol.atoms() {
            assert_eq!(mol.atom(node).hydrogen_count, 1);
            assert!(mol.atom(node).is_aromatic);
        }
    }

    #[test]
    fn naphthalene() {
        let smiles_mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 10);
        assert_eq!(mol.bond_count(), 11);
        assert_eq!(count_double_bonds(&mol), 5);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn pyridine() {
        let smiles_mol = parse_smiles("c1ccncc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        assert_eq!(count_double_bonds(&mol), 3);
        assert!(is_valid_kekulization(&mol));
        assert_eq!(mol.atom(n(3)).atomic_num, 7);
        assert_eq!(mol.atom(n(3)).hydrogen_count, 0);
    }

    #[test]
    fn pyrrole() {
        let smiles_mol = parse_smiles("[nH]1cccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(count_double_bonds(&mol), 2);
        assert!(is_valid_kekulization(&mol));
        assert_eq!(mol.atom(n(0)).hydrogen_count, 1);
    }

    #[test]
    fn furan() {
        let smiles_mol = parse_smiles("o1cccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(count_double_bonds(&mol), 2);
        assert!(is_valid_kekulization(&mol));
        assert_eq!(mol.atom(n(0)).hydrogen_count, 0);
    }

    #[test]
    fn thiophene() {
        let smiles_mol = parse_smiles("s1cccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(count_double_bonds(&mol), 2);
        assert!(is_valid_kekulization(&mol));
        assert_eq!(mol.atom(n(0)).hydrogen_count, 0);
    }

    #[test]
    fn cyclopentadienyl_anion() {
        let smiles_mol = parse_smiles("[c-]1cccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(count_double_bonds(&mol), 2);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn cyclopentadienyl_no_charge_fails() {
        let smiles_mol = parse_smiles("c1cccc1").unwrap();
        let result = kekulize(smiles_mol);
        assert!(result.is_err());
    }

    #[test]
    fn cyclobutadiene_kekulizes() {
        let smiles_mol = parse_smiles("c1ccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(count_double_bonds(&mol), 2);
    }

    #[test]
    fn phenol() {
        let smiles_mol = parse_smiles("Oc1ccccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 7);
        assert_eq!(mol.bond_count(), 7);
        let bond_o_c = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(bond_o_c).order, BondOrder::Single);
        assert_eq!(mol.atom(n(0)).hydrogen_count, 1);
        assert_eq!(count_double_bonds(&mol), 3);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn toluene() {
        let smiles_mol = parse_smiles("Cc1ccccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 7);
        assert_eq!(mol.bond_count(), 7);
        let bond_c_c = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(bond_c_c).order, BondOrder::Single);
        assert_eq!(count_double_bonds(&mol), 3);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn h_counts_preserved() {
        let smiles_mol = parse_smiles("c1ccccc1").unwrap();
        let expected: Vec<u8> = smiles_mol
            .atoms()
            .map(|n| smiles_mol.atom(n).hydrogen_count)
            .collect();
        let mol = kekulize(smiles_mol).unwrap();
        let actual: Vec<u8> = mol
            .atoms()
            .map(|n| mol.atom(n).hydrogen_count)
            .collect();
        assert_eq!(expected, actual);
    }

    #[test]
    fn non_aromatic_passthrough() {
        let smiles_mol = parse_smiles("C=CC").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        let e01 = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(e01).order, BondOrder::Double);
        let e12 = mol.bond_between(n(1), n(2)).unwrap();
        assert_eq!(mol.bond(e12).order, BondOrder::Single);
    }

    #[test]
    fn stereo_preserved() {
        let smiles_mol = parse_smiles("F/C=C/F").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        let e = mol.bond_between(n(1), n(2)).unwrap();
        assert_eq!(mol.bond(e).order, BondOrder::Double);
        assert!(matches!(mol.bond(e).stereo, BondStereo::Trans(_, _)));
    }

    #[test]
    fn error_display() {
        let err = KekulizeError::Unkekulizable(vec![NodeIndex::new(0), NodeIndex::new(2)]);
        let msg = format!("{}", err);
        assert!(msg.contains("0"));
        assert!(msg.contains("2"));
    }

    #[test]
    fn imidazole() {
        let smiles_mol = parse_smiles("c1c[nH]cn1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(count_double_bonds(&mol), 2);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn anthracene() {
        let smiles_mol = parse_smiles("c1ccc2cc3ccccc3cc2c1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 14);
        assert_eq!(count_double_bonds(&mol), 7);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn from_smiles_benzene() {
        let mol = crate::smiles::from_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        assert_eq!(count_double_bonds(&mol), 3);
    }

    #[test]
    fn from_smiles_non_aromatic() {
        let mol = crate::smiles::from_smiles("CC").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bond(mol.bonds().next().unwrap()).order, BondOrder::Single);
    }

    #[test]
    fn from_smiles_error_propagation() {
        let result = crate::smiles::from_smiles("");
        assert!(result.is_err());
    }

    #[test]
    fn pyridinium() {
        let smiles_mol = parse_smiles("[n+]1ccccc1").unwrap();
        let mol = kekulize(smiles_mol).unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(count_double_bonds(&mol), 3);
        assert!(is_valid_kekulization(&mol));
    }

    #[test]
    fn odd_ring_unkekulizable() {
        let smiles_mol = parse_smiles("c1cccc1").unwrap();
        let result = kekulize(smiles_mol);
        assert!(result.is_err());
        if let Err(KekulizeError::Unkekulizable(atoms)) = result {
            assert!(!atoms.is_empty());
        }
    }
}
