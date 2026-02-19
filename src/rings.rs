use std::collections::VecDeque;

use petgraph::algo::connected_components;
use petgraph::graph::NodeIndex;

use crate::mol::Mol;

#[derive(Debug, Clone)]
pub struct RingInfo {
    rings: Vec<Vec<NodeIndex>>,
}

impl RingInfo {
    pub fn sssr<A, B>(mol: &Mol<A, B>) -> Self {
        let num_expected = Self::expected_ring_count(mol);
        if num_expected == 0 {
            return Self { rings: vec![] };
        }

        let n = mol.atom_count();
        let num_edges = mol.bond_count();

        let dist = all_pairs_bfs(mol, n);
        let pred = all_pairs_predecessors(mol, n, &dist);

        let mut candidates: Vec<Vec<NodeIndex>> = Vec::new();

        for edge in mol.bonds() {
            let (u, v) = match mol.bond_endpoints(edge) {
                Some(pair) => pair,
                None => continue,
            };
            for w_idx in 0..n {
                let w = NodeIndex::new(w_idx);
                let du = dist[w.index()][u.index()];
                let dv = dist[w.index()][v.index()];
                if du == u32::MAX || dv == u32::MAX {
                    continue;
                }
                let ring_size = du as usize + dv as usize + 1;
                if ring_size < 3 {
                    continue;
                }
                let path_u = reconstruct_path(&pred, w, u);
                let path_v = reconstruct_path(&pred, w, v);
                if paths_share_internal_node(&path_u, &path_v) {
                    continue;
                }
                let mut ring = path_u;
                for &node in path_v[1..].iter().rev() {
                    ring.push(node);
                }
                candidates.push(ring);
            }
        }

        candidates.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        candidates.dedup();

        let rings =
            select_independent_rings(&candidates, num_expected, num_edges, mol);

        Self { rings }
    }

    pub fn num_rings(&self) -> usize {
        self.rings.len()
    }

    pub fn rings(&self) -> &[Vec<NodeIndex>] {
        &self.rings
    }

    pub fn is_ring_atom(&self, atom: NodeIndex) -> bool {
        self.rings.iter().any(|ring| ring.contains(&atom))
    }

    pub fn is_ring_bond(&self, a: NodeIndex, b: NodeIndex) -> bool {
        self.rings.iter().any(|ring| {
            let len = ring.len();
            (0..len).any(|i| {
                let j = (i + 1) % len;
                (ring[i] == a && ring[j] == b) || (ring[i] == b && ring[j] == a)
            })
        })
    }

    pub fn smallest_ring_size(&self, atom: NodeIndex) -> Option<usize> {
        self.rings
            .iter()
            .filter(|ring| ring.contains(&atom))
            .map(|ring| ring.len())
            .min()
    }

    pub fn atom_rings(&self, atom: NodeIndex) -> Vec<&Vec<NodeIndex>> {
        self.rings
            .iter()
            .filter(|ring| ring.contains(&atom))
            .collect()
    }

    pub fn expected_ring_count<A, B>(mol: &Mol<A, B>) -> usize {
        let v = mol.atom_count();
        let e = mol.bond_count();
        let c = connected_components(mol.graph());
        (e + c).saturating_sub(v)
    }
}

fn all_pairs_bfs<A, B>(mol: &Mol<A, B>, n: usize) -> Vec<Vec<u32>> {
    let mut dist = vec![vec![u32::MAX; n]; n];
    for (src_idx, row) in dist.iter_mut().enumerate() {
        let src = NodeIndex::new(src_idx);
        row[src_idx] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(src);
        while let Some(cur) = queue.pop_front() {
            let d = row[cur.index()];
            for nb in mol.neighbors(cur) {
                if row[nb.index()] == u32::MAX {
                    row[nb.index()] = d + 1;
                    queue.push_back(nb);
                }
            }
        }
    }
    dist
}

fn all_pairs_predecessors<A, B>(
    mol: &Mol<A, B>,
    n: usize,
    dist: &[Vec<u32>],
) -> Vec<Vec<Option<NodeIndex>>> {
    let mut pred = vec![vec![None; n]; n];
    for src_idx in 0..n {
        let src = NodeIndex::new(src_idx);
        let mut queue = VecDeque::new();
        queue.push_back(src);
        let mut visited = vec![false; n];
        visited[src_idx] = true;
        while let Some(cur) = queue.pop_front() {
            for nb in mol.neighbors(cur) {
                if !visited[nb.index()]
                    && dist[src_idx][nb.index()] == dist[src_idx][cur.index()] + 1
                {
                    visited[nb.index()] = true;
                    pred[src_idx][nb.index()] = Some(cur);
                    queue.push_back(nb);
                }
            }
        }
    }
    pred
}

fn reconstruct_path(
    pred: &[Vec<Option<NodeIndex>>],
    src: NodeIndex,
    dst: NodeIndex,
) -> Vec<NodeIndex> {
    let mut path = vec![dst];
    let mut cur = dst;
    while cur != src {
        match pred[src.index()][cur.index()] {
            Some(p) => {
                path.push(p);
                cur = p;
            }
            None => return vec![],
        }
    }
    path.reverse();
    path
}

fn paths_share_internal_node(path_u: &[NodeIndex], path_v: &[NodeIndex]) -> bool {
    if path_u.len() < 2 || path_v.len() < 2 {
        return false;
    }
    let internal_u = &path_u[1..];
    let internal_v = &path_v[1..];
    for node in internal_u {
        if internal_v.contains(node) {
            return true;
        }
    }
    false
}

fn ring_to_edge_bitvector<A, B>(
    ring: &[NodeIndex],
    num_edges: usize,
    mol: &Mol<A, B>,
) -> Vec<u64> {
    let num_words = num_edges.div_ceil(64);
    let mut bv = vec![0u64; num_words];
    let len = ring.len();
    for i in 0..len {
        let a = ring[i];
        let b = ring[(i + 1) % len];
        if let Some(edge) = mol.bond_between(a, b) {
            let idx = edge.index();
            bv[idx / 64] |= 1u64 << (idx % 64);
        }
    }
    bv
}

fn select_independent_rings<A, B>(
    candidates: &[Vec<NodeIndex>],
    num_needed: usize,
    num_edges: usize,
    mol: &Mol<A, B>,
) -> Vec<Vec<NodeIndex>> {
    let mut result = Vec::with_capacity(num_needed);
    let mut basis: Vec<Vec<u64>> = Vec::with_capacity(num_needed);

    for ring in candidates {
        if result.len() >= num_needed {
            break;
        }
        let bv = ring_to_edge_bitvector(ring, num_edges, mol);
        if bv.iter().all(|&w| w == 0) {
            continue;
        }
        if try_add_to_basis(&mut basis, bv) {
            result.push(normalize_ring(ring));
        }
    }

    result.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    result
}

fn try_add_to_basis(basis: &mut Vec<Vec<u64>>, candidate: Vec<u64>) -> bool {
    let mut v = candidate;
    for row in basis.iter() {
        let pivot = leading_bit(row);
        if let Some(p) = pivot {
            if v[p / 64] & (1u64 << (p % 64)) != 0 {
                xor_into(&mut v, row);
            }
        }
    }
    if v.iter().all(|&w| w == 0) {
        return false;
    }
    basis.push(v);
    true
}

fn leading_bit(bv: &[u64]) -> Option<usize> {
    for (i, &word) in bv.iter().enumerate() {
        if word != 0 {
            return Some(i * 64 + word.trailing_zeros() as usize);
        }
    }
    None
}

fn xor_into(a: &mut [u64], b: &[u64]) {
    for (aw, bw) in a.iter_mut().zip(b.iter()) {
        *aw ^= *bw;
    }
}

fn normalize_ring(ring: &[NodeIndex]) -> Vec<NodeIndex> {
    if ring.is_empty() {
        return vec![];
    }
    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, idx)| idx)
        .map(|(i, _)| i)
        .unwrap();

    let len = ring.len();
    let mut normalized = Vec::with_capacity(len);
    for i in 0..len {
        normalized.push(ring[(min_pos + i) % len]);
    }

    if len > 2 && normalized[1] > normalized[len - 1] {
        normalized[1..].reverse();
    }

    normalized
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    #[test]
    fn cyclohexane() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 1);
        assert_eq!(ri.rings()[0].len(), 6);
    }

    #[test]
    fn cyclopropane() {
        let mol = from_smiles("C1CC1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 1);
        assert_eq!(ri.rings()[0].len(), 3);
    }

    #[test]
    fn benzene() {
        let mol = from_smiles("c1ccccc1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 1);
        assert_eq!(ri.rings()[0].len(), 6);
    }

    #[test]
    fn acyclic() {
        let mol = from_smiles("CCCC").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 0);
    }

    #[test]
    fn naphthalene() {
        let mol = from_smiles("c1ccc2ccccc2c1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 2);
        for ring in ri.rings() {
            assert_eq!(ring.len(), 6);
        }
    }

    #[test]
    fn anthracene() {
        let mol = from_smiles("c1ccc2cc3ccccc3cc2c1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 3);
    }

    #[test]
    fn spiro_nonane() {
        let mol = from_smiles("C1CCC2(CC1)CCC2").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 2);
    }

    #[test]
    fn norbornane() {
        let mol = from_smiles("C1CC2CC1CC2").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 2);
    }

    #[test]
    fn benzene_all_atoms_in_ring() {
        let mol = from_smiles("c1ccccc1").unwrap();
        let ri = RingInfo::sssr(&mol);
        for i in 0..6 {
            assert!(ri.is_ring_atom(n(i)), "atom {} should be in ring", i);
        }
    }

    #[test]
    fn benzene_all_bonds_in_ring() {
        let mol = from_smiles("c1ccccc1").unwrap();
        let ri = RingInfo::sssr(&mol);
        for i in 0..6 {
            let j = (i + 1) % 6;
            assert!(
                ri.is_ring_bond(n(i), n(j)),
                "bond {}-{} should be in ring",
                i,
                j
            );
        }
    }

    #[test]
    fn phenol_oxygen_not_in_ring() {
        let mol = from_smiles("Oc1ccccc1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert!(!ri.is_ring_atom(n(0)));
        for i in 1..7 {
            assert!(ri.is_ring_atom(n(i)), "atom {} should be in ring", i);
        }
    }

    #[test]
    fn cubane_cyclomatic_number() {
        let mol = from_smiles("C12C3C4C1C5C3C4C25").unwrap();
        assert_eq!(RingInfo::expected_ring_count(&mol), 5);
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 5);
    }

    #[test]
    fn decalin() {
        let mol = from_smiles("C1CCC2CCCCC2C1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.num_rings(), 2);
        for ring in ri.rings() {
            assert_eq!(ring.len(), 6);
        }
    }

    #[test]
    fn smallest_ring_size_cyclohexane() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.smallest_ring_size(n(0)), Some(6));
    }

    #[test]
    fn smallest_ring_size_acyclic() {
        let mol = from_smiles("CCCC").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert_eq!(ri.smallest_ring_size(n(0)), None);
    }

    #[test]
    fn atom_rings_naphthalene_shared() {
        let mol = from_smiles("c1ccc2ccccc2c1").unwrap();
        let ri = RingInfo::sssr(&mol);
        let shared: Vec<NodeIndex> = mol
            .atoms()
            .filter(|&a| ri.atom_rings(a).len() == 2)
            .collect();
        assert_eq!(shared.len(), 2);
    }

    #[test]
    fn toluene_methyl_not_in_ring() {
        let mol = from_smiles("Cc1ccccc1").unwrap();
        let ri = RingInfo::sssr(&mol);
        assert!(!ri.is_ring_atom(n(0)));
    }

    #[test]
    fn cyclomatic_cyclohexane() {
        let mol = from_smiles("C1CCCCC1").unwrap();
        assert_eq!(RingInfo::expected_ring_count(&mol), 1);
    }

    #[test]
    fn cyclomatic_naphthalene() {
        let mol = from_smiles("c1ccc2ccccc2c1").unwrap();
        assert_eq!(RingInfo::expected_ring_count(&mol), 2);
    }
}

