use std::collections::VecDeque;

use petgraph::graph::NodeIndex;

use crate::mol::Mol;

pub fn adjacency_matrix<A, B>(mol: &Mol<A, B>) -> Vec<Vec<bool>> {
    let n = mol.atom_count();
    let mut matrix = vec![vec![false; n]; n];
    for edge in mol.bonds() {
        if let Some((a, b)) = mol.bond_endpoints(edge) {
            matrix[a.index()][b.index()] = true;
            matrix[b.index()][a.index()] = true;
        }
    }
    matrix
}

pub fn distance_matrix<A, B>(mol: &Mol<A, B>) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    let mut dist = vec![vec![usize::MAX; n]; n];
    for start in mol.atoms() {
        let si = start.index();
        dist[si][si] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(current) = queue.pop_front() {
            let d = dist[si][current.index()];
            for neighbor in mol.neighbors(current) {
                if dist[si][neighbor.index()] == usize::MAX {
                    dist[si][neighbor.index()] = d + 1;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    dist
}

pub fn shortest_path<A, B>(
    mol: &Mol<A, B>,
    from: NodeIndex,
    to: NodeIndex,
) -> Option<Vec<NodeIndex>> {
    if from == to {
        return Some(vec![from]);
    }
    let n = mol.atom_count();
    let mut pred = vec![None; n];
    let mut visited = vec![false; n];
    visited[from.index()] = true;
    let mut queue = VecDeque::new();
    queue.push_back(from);
    while let Some(current) = queue.pop_front() {
        for neighbor in mol.neighbors(current) {
            if !visited[neighbor.index()] {
                visited[neighbor.index()] = true;
                pred[neighbor.index()] = Some(current);
                if neighbor == to {
                    let mut path = vec![to];
                    let mut node = to;
                    while let Some(p) = pred[node.index()] {
                        path.push(p);
                        node = p;
                    }
                    path.reverse();
                    return Some(path);
                }
                queue.push_back(neighbor);
            }
        }
    }
    None
}

pub fn connected_components<A, B>(mol: &Mol<A, B>) -> Vec<Vec<NodeIndex>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components = Vec::new();
    for node in mol.atoms() {
        if visited[node.index()] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![node];
        while let Some(current) = stack.pop() {
            if visited[current.index()] {
                continue;
            }
            visited[current.index()] = true;
            component.push(current);
            for neighbor in mol.neighbors(current) {
                if !visited[neighbor.index()] {
                    stack.push(neighbor);
                }
            }
        }
        component.sort();
        components.push(component);
    }
    components
}

pub fn num_components<A, B>(mol: &Mol<A, B>) -> usize {
    connected_components(mol).len()
}

pub fn get_fragments<A: Clone, B: Clone>(mol: &Mol<A, B>) -> Vec<Mol<A, B>> {
    let components = connected_components(mol);
    let mut fragments = Vec::with_capacity(components.len());
    for component in &components {
        let mut frag = Mol::new();
        let mut index_map = vec![NodeIndex::new(0); mol.atom_count()];
        for &old_idx in component {
            let new_idx = frag.add_atom(mol.atom(old_idx).clone());
            index_map[old_idx.index()] = new_idx;
        }
        for &old_idx in component {
            for edge in mol.bonds_of(old_idx) {
                if let Some((a, b)) = mol.bond_endpoints(edge) {
                    if a == old_idx && a.index() < b.index() {
                        frag.add_bond(
                            index_map[a.index()],
                            index_map[b.index()],
                            mol.bond(edge).clone(),
                        );
                    }
                }
            }
        }
        fragments.push(frag);
    }
    fragments
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    #[test]
    fn adjacency_ethane() {
        let mol = from_smiles("CC").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 2);
        assert!(!adj[0][0]);
        assert!(adj[0][1]);
        assert!(adj[1][0]);
        assert!(!adj[1][1]);
    }

    #[test]
    fn adjacency_methane() {
        let mol = from_smiles("C").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 1);
        assert!(!adj[0][0]);
    }

    #[test]
    fn adjacency_triangle() {
        let mol = from_smiles("C1CC1").unwrap();
        let adj = adjacency_matrix(&mol);
        assert_eq!(adj.len(), 3);
        assert!(adj[0][1]);
        assert!(adj[1][2]);
        assert!(adj[0][2]);
        assert!(adj[1][0]);
        assert!(adj[2][1]);
        assert!(adj[2][0]);
        for i in 0..3 {
            assert!(!adj[i][i]);
        }
    }

    #[test]
    fn distance_linear_chain() {
        let mol = from_smiles("CCCC").unwrap();
        let dist = distance_matrix(&mol);
        assert_eq!(dist[0][3], 3);
        assert_eq!(dist[0][1], 1);
        assert_eq!(dist[1][3], 2);
        assert_eq!(dist[0][0], 0);
    }

    #[test]
    fn distance_cyclohexane() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        let dist = distance_matrix(&mol);
        assert_eq!(dist.len(), 6);
        assert_eq!(dist[0][3], 3);
        assert_eq!(dist[0][1], 1);
        assert_eq!(dist[0][2], 2);
    }

    #[test]
    fn shortest_path_linear() {
        let mol = from_smiles("CCCCC").unwrap();
        let path = shortest_path(&mol, NodeIndex::new(0), NodeIndex::new(4));
        assert_eq!(
            path,
            Some(vec![
                NodeIndex::new(0),
                NodeIndex::new(1),
                NodeIndex::new(2),
                NodeIndex::new(3),
                NodeIndex::new(4),
            ])
        );
    }

    #[test]
    fn shortest_path_no_path() {
        let mol = from_smiles("[Na+].[Cl-]").unwrap();
        let path = shortest_path(&mol, NodeIndex::new(0), NodeIndex::new(1));
        assert_eq!(path, None);
    }

    #[test]
    fn components_nacl() {
        let mol = from_smiles("[Na+].[Cl-]").unwrap();
        let comps = connected_components(&mol);
        assert_eq!(comps.len(), 2);
    }

    #[test]
    fn components_single() {
        let mol = from_smiles("CCO").unwrap();
        assert_eq!(num_components(&mol), 1);
    }

    #[test]
    fn components_empty() {
        let mol: Mol<(), ()> = Mol::new();
        assert_eq!(num_components(&mol), 0);
    }

    #[test]
    fn fragments_three() {
        let mol = from_smiles("[Na+].[Cl-].O").unwrap();
        let frags = get_fragments(&mol);
        assert_eq!(frags.len(), 3);
        let mut counts: Vec<usize> = frags.iter().map(|f| f.atom_count()).collect();
        counts.sort();
        assert_eq!(counts, vec![1, 1, 1]);
    }

    #[test]
    fn fragments_single() {
        let mol = from_smiles("CCO").unwrap();
        let frags = get_fragments(&mol);
        assert_eq!(frags.len(), 1);
        assert_eq!(frags[0].atom_count(), mol.atom_count());
        assert_eq!(frags[0].bond_count(), mol.bond_count());
    }
}
