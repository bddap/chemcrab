use petgraph::graph::{EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;

pub struct Mol<A, B> {
    graph: UnGraph<A, B>,
}

impl<A, B> Mol<A, B> {
    pub fn new() -> Self {
        Self {
            graph: UnGraph::default(),
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
}

impl<A, B> Default for Mol<A, B> {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: std::fmt::Debug, B: std::fmt::Debug> std::fmt::Debug for Mol<A, B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mol")
            .field("atom_count", &self.atom_count())
            .field("bond_count", &self.bond_count())
            .finish()
    }
}
