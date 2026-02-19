use petgraph::graph::{EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomId {
    Node(NodeIndex),
    VirtualH(NodeIndex, u8),
}

pub struct Mol<A, B> {
    graph: UnGraph<A, B>,
    tetrahedral_stereo: Vec<[AtomId; 4]>,
}

impl<A, B> Mol<A, B> {
    pub fn new() -> Self {
        Self {
            graph: UnGraph::default(),
            tetrahedral_stereo: Vec::new(),
        }
    }

    pub fn graph(&self) -> &UnGraph<A, B> {
        &self.graph
    }

    pub fn atom(&self, idx: NodeIndex) -> &A {
        &self.graph[idx]
    }

    pub fn atom_mut(&mut self, idx: NodeIndex) -> &mut A {
        &mut self.graph[idx]
    }

    pub fn bond(&self, idx: EdgeIndex) -> &B {
        &self.graph[idx]
    }

    pub fn bond_mut(&mut self, idx: EdgeIndex) -> &mut B {
        &mut self.graph[idx]
    }

    pub fn add_atom(&mut self, atom: A) -> NodeIndex {
        self.graph.add_node(atom)
    }

    pub fn add_bond(&mut self, a: NodeIndex, b: NodeIndex, bond: B) -> EdgeIndex {
        self.graph.add_edge(a, b, bond)
    }

    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn neighbors(&self, idx: NodeIndex) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph.neighbors(idx)
    }

    pub fn bonds_of(&self, idx: NodeIndex) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph.edges(idx).map(|e| e.id())
    }

    pub fn atoms(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn bonds(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bond_between(&self, a: NodeIndex, b: NodeIndex) -> Option<EdgeIndex> {
        self.graph.find_edge(a, b)
    }

    pub fn bond_endpoints(&self, idx: EdgeIndex) -> Option<(NodeIndex, NodeIndex)> {
        self.graph.edge_endpoints(idx)
    }

    pub fn swap_remove_atom(&mut self, idx: NodeIndex) -> A {
        self.graph.remove_node(idx).expect("node index out of bounds")
    }

    pub fn swap_remove_bond(&mut self, idx: EdgeIndex) -> B {
        self.graph.remove_edge(idx).expect("edge index out of bounds")
    }

    pub fn tetrahedral_stereo(&self) -> &[[AtomId; 4]] {
        &self.tetrahedral_stereo
    }

    pub fn set_tetrahedral_stereo(&mut self, stereo: Vec<[AtomId; 4]>) {
        self.tetrahedral_stereo = stereo;
    }

    pub fn tetrahedral_stereo_for(&self, center: NodeIndex) -> Option<&[AtomId; 4]> {
        self.tetrahedral_stereo
            .iter()
            .find(|s| s[0] == AtomId::Node(center))
    }

    pub fn add_tetrahedral_stereo(&mut self, stereo: [AtomId; 4]) {
        self.tetrahedral_stereo.push(stereo);
    }

    pub fn remove_tetrahedral_stereo(&mut self, center: NodeIndex) {
        self.tetrahedral_stereo
            .retain(|s| s[0] != AtomId::Node(center));
    }
}

impl<A: Clone, B: Clone> Clone for Mol<A, B> {
    fn clone(&self) -> Self {
        Self {
            graph: self.graph.clone(),
            tetrahedral_stereo: self.tetrahedral_stereo.clone(),
        }
    }
}

impl<A, B> Default for Mol<A, B> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: PartialEq, B: PartialEq> PartialEq for Mol<A, B> {
    fn eq(&self, other: &Self) -> bool {
        if self.atom_count() != other.atom_count() || self.bond_count() != other.bond_count() {
            return false;
        }
        for idx in self.atoms() {
            if idx.index() >= other.atom_count() {
                return false;
            }
            if self.atom(idx) != other.atom(idx) {
                return false;
            }
        }
        for idx in self.bonds() {
            if idx.index() >= other.bond_count() {
                return false;
            }
            if self.bond(idx) != other.bond(idx) {
                return false;
            }
            if self.bond_endpoints(idx) != other.bond_endpoints(idx) {
                return false;
            }
        }
        self.tetrahedral_stereo == other.tetrahedral_stereo
    }
}

impl<A: std::fmt::Debug, B: std::fmt::Debug> std::fmt::Debug for Mol<A, B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mol")
            .field("atom_count", &self.atom_count())
            .field("bond_count", &self.bond_count())
            .field("tetrahedral_stereo", &self.tetrahedral_stereo)
            .finish()
    }
}

pub(crate) fn permutation_parity<T: Eq>(from: &[T], to: &[T]) -> bool {
    let n = from.len();
    if n != to.len() {
        return true;
    }
    let perm: Vec<usize> = from
        .iter()
        .map(|f| to.iter().position(|t| t == f).unwrap_or(0))
        .collect();
    let mut visited = vec![false; n];
    let mut swaps = 0usize;
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
    swaps.is_multiple_of(2)
}
