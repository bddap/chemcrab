use petgraph::algo::connected_components;
use petgraph::graph::{EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;

/// Identifies an atom in a stereo context.
///
/// In cheminformatics, most hydrogen atoms are not stored as graph nodes —
/// they are implied by the heavy atom's valence and stored as an integer
/// count on the parent ([`Atom::hydrogen_count`](crate::Atom::hydrogen_count)). Stereo descriptors,
/// however, need to reference *all* neighbors of a stereocenter, including
/// those virtual hydrogens.
///
/// `AtomId` bridges the gap: [`Node`](AtomId::Node) refers to a real graph
/// node, while [`VirtualH`](AtomId::VirtualH) refers to a suppressed
/// hydrogen by naming its parent atom and a disambiguation index (for
/// parents with multiple virtual Hs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomId {
    /// A real atom in the molecular graph.
    Node(NodeIndex),
    /// A virtual (suppressed) hydrogen.
    ///
    /// The first field is the parent atom's [`NodeIndex`]; the second is a
    /// zero-based index distinguishing multiple virtual Hs on the same parent.
    VirtualH(NodeIndex, u8),
}

/// Tetrahedral stereochemistry descriptor.
///
/// A tetrahedral stereocenter — most commonly a carbon bonded to four
/// different substituents — has two non-superimposable mirror-image
/// arrangements (enantiomers). This struct encodes which arrangement is
/// intended.
///
/// **Plane-rule convention:** `above[0]` is above the plane defined by
/// `above[1..4]`, with orientation given by the right-hand rule applied to
/// the `[1, 2, 3]` sequence. Swapping any two elements of `above` flips
/// the handedness (converts one enantiomer to the other).
///
/// `center` is the stereocenter's node index, used for lookup — it does not
/// participate in the parity calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TetrahedralStereo {
    /// The stereocenter atom.
    pub center: NodeIndex,
    /// Four neighbors in parity-significant order.
    ///
    /// `above[0]` is above the plane of `above[1..4]` by the right-hand
    /// rule. Swapping any pair flips the handedness.
    pub above: [AtomId; 4],
}

/// E/Z (cis/trans) stereochemistry descriptor for a double bond.
///
/// Restricted rotation around a double bond locks substituents on one side
/// or the other. When substituents are on the same side, the configuration
/// is Z (from German *zusammen*, "together") or *cis*; when on opposite
/// sides, E (*entgegen*, "opposite") or *trans*.
///
/// This struct records a double bond and the pair of reference substituents
/// that are on the **same side** (cis). There is no explicit `Cis`/`Trans`
/// flag — the geometry is fully determined by which pair is recorded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EZStereo {
    /// The double bond's endpoint atoms, stored with the lower node index
    /// first.
    pub bond: (NodeIndex, NodeIndex),
    /// Two substituents that are on the same side (cis) of the double bond.
    ///
    /// `refs[0]` is a neighbor of `bond.0`; `refs[1]` is a neighbor of
    /// `bond.1`.
    pub refs: [AtomId; 2],
}

/// A molecular graph parameterized over atom type `A` and bond type `B`.
///
/// `Mol` is the central type in chemcrab. It wraps a petgraph
/// [`UnGraph`](petgraph::graph::UnGraph) and adds stereochemistry
/// descriptors. Stereochemistry lives here — not on individual atoms or
/// bonds — because it describes spatial *relationships* between neighbors.
///
/// The standard concrete type after SMILES parsing is `Mol<Atom, Bond>`,
/// where every bond has a definite Kekulé order. Before kekulization, the
/// intermediate type `Mol<Atom, AromaticBond>` preserves aromatic bonds.
///
/// Generic algorithms can accept `Mol<A, B>` with trait bounds on `A` and
/// `B` (e.g. `A: HasAtomicNum`) so they work with custom atom types too.
///
/// # Examples
///
/// ```
/// # use chemcrab::smiles;
/// let mol = smiles::from_smiles("CCO")?;
/// assert_eq!(mol.atom_count(), 3);
/// assert_eq!(mol.bond_count(), 2);
/// # Ok::<(), chemcrab::smiles::SmilesError>(())
/// ```
pub struct Mol<A, B> {
    graph: UnGraph<A, B>,
    tetrahedral_stereo: Vec<TetrahedralStereo>,
    ez_stereo: Vec<EZStereo>,
}

impl<A, B> Mol<A, B> {
    /// Creates an empty molecule with no atoms, bonds, or stereochemistry.
    pub fn new() -> Self {
        Self {
            graph: UnGraph::default(),
            tetrahedral_stereo: Vec::new(),
            ez_stereo: Vec::new(),
        }
    }

    /// Returns the number of connected components in the molecular graph.
    pub fn connected_component_count(&self) -> usize {
        connected_components(&self.graph)
    }

    /// Returns a reference to the atom at `idx`.
    pub fn atom(&self, idx: NodeIndex) -> &A {
        &self.graph[idx]
    }

    /// Returns a mutable reference to the atom at `idx`.
    pub fn atom_mut(&mut self, idx: NodeIndex) -> &mut A {
        &mut self.graph[idx]
    }

    /// Returns a reference to the bond at `idx`.
    pub fn bond(&self, idx: EdgeIndex) -> &B {
        &self.graph[idx]
    }

    /// Returns a mutable reference to the bond at `idx`.
    pub fn bond_mut(&mut self, idx: EdgeIndex) -> &mut B {
        &mut self.graph[idx]
    }

    /// Adds an atom to the molecule, returning its [`NodeIndex`].
    pub fn add_atom(&mut self, atom: A) -> NodeIndex {
        self.graph.add_node(atom)
    }

    /// Adds a bond between atoms `a` and `b`, returning its [`EdgeIndex`].
    pub fn add_bond(&mut self, a: NodeIndex, b: NodeIndex, bond: B) -> EdgeIndex {
        self.graph.add_edge(a, b, bond)
    }

    /// Returns the number of atoms (graph nodes) in the molecule.
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Returns the number of bonds (graph edges) in the molecule.
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Iterates over the graph neighbors of atom `idx`.
    pub fn neighbors(&self, idx: NodeIndex) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph.neighbors(idx)
    }

    /// Iterates over the [`EdgeIndex`]es of bonds incident to atom `idx`.
    pub fn bonds_of(&self, idx: NodeIndex) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph.edges(idx).map(|e| e.id())
    }

    /// Iterates over all atom indices in the molecule.
    pub fn atoms(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph.node_indices()
    }

    /// Iterates over all bond indices in the molecule.
    pub fn bonds(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph.edge_indices()
    }

    /// Returns the bond between atoms `a` and `b`, if one exists.
    pub fn bond_between(&self, a: NodeIndex, b: NodeIndex) -> Option<EdgeIndex> {
        self.graph.find_edge(a, b)
    }

    /// Returns the two endpoint atoms of a bond.
    pub fn bond_endpoints(&self, idx: EdgeIndex) -> Option<(NodeIndex, NodeIndex)> {
        self.graph.edge_endpoints(idx)
    }

    /// Returns a slice of all tetrahedral stereo descriptors.
    pub fn tetrahedral_stereo(&self) -> &[TetrahedralStereo] {
        &self.tetrahedral_stereo
    }

    /// Replaces all tetrahedral stereo descriptors.
    pub fn set_tetrahedral_stereo(&mut self, stereo: Vec<TetrahedralStereo>) {
        self.tetrahedral_stereo = stereo;
    }

    /// Looks up the tetrahedral stereo descriptor for a given center atom.
    pub fn tetrahedral_stereo_for(&self, center: NodeIndex) -> Option<&TetrahedralStereo> {
        self.tetrahedral_stereo.iter().find(|s| s.center == center)
    }

    /// Appends a tetrahedral stereo descriptor.
    pub fn add_tetrahedral_stereo(&mut self, stereo: TetrahedralStereo) {
        self.tetrahedral_stereo.push(stereo);
    }

    /// Removes the tetrahedral stereo descriptor for `center`, if present.
    pub fn remove_tetrahedral_stereo(&mut self, center: NodeIndex) {
        self.tetrahedral_stereo.retain(|s| s.center != center);
    }

    /// Returns a slice of all E/Z stereo descriptors.
    pub fn ez_stereo(&self) -> &[EZStereo] {
        &self.ez_stereo
    }

    /// Replaces all E/Z stereo descriptors.
    pub fn set_ez_stereo(&mut self, stereo: Vec<EZStereo>) {
        self.ez_stereo = stereo;
    }

    /// Looks up the E/Z stereo descriptor for the double bond between `a` and `b`.
    ///
    /// The argument order does not matter — endpoints are normalized internally.
    pub fn ez_stereo_for(&self, a: NodeIndex, b: NodeIndex) -> Option<&EZStereo> {
        let (lo, hi) = if a.index() < b.index() {
            (a, b)
        } else {
            (b, a)
        };
        self.ez_stereo.iter().find(|s| s.bond == (lo, hi))
    }

    /// Appends an E/Z stereo descriptor.
    pub fn add_ez_stereo(&mut self, stereo: EZStereo) {
        self.ez_stereo.push(stereo);
    }

    /// Removes the E/Z stereo descriptor for the bond between `a` and `b`, if present.
    pub fn remove_ez_stereo(&mut self, a: NodeIndex, b: NodeIndex) {
        let (lo, hi) = if a.index() < b.index() {
            (a, b)
        } else {
            (b, a)
        };
        self.ez_stereo.retain(|s| s.bond != (lo, hi));
    }
}

impl<A: Clone, B: Clone> Clone for Mol<A, B> {
    fn clone(&self) -> Self {
        Self {
            graph: self.graph.clone(),
            tetrahedral_stereo: self.tetrahedral_stereo.clone(),
            ez_stereo: self.ez_stereo.clone(),
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
        self.tetrahedral_stereo == other.tetrahedral_stereo && self.ez_stereo == other.ez_stereo
    }
}

impl<A: std::fmt::Debug, B: std::fmt::Debug> std::fmt::Debug for Mol<A, B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mol")
            .field("atom_count", &self.atom_count())
            .field("bond_count", &self.bond_count())
            .field("tetrahedral_stereo", &self.tetrahedral_stereo)
            .field("ez_stereo", &self.ez_stereo)
            .finish()
    }
}

/// Returns `true` if the permutation mapping `from` to `to` is even (an even
/// number of transpositions).
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
