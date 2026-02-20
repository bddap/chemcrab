/// Concrete bond order after kekulization.
///
/// Every bond in a [`Mol<Atom, Bond>`](crate::Mol) has one of these three
/// orders. There is no `Aromatic` variant — aromatic bonds only exist in the
/// intermediate [`SmilesBond`] type before kekulization resolves them to
/// alternating single and double bonds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondOrder {
    /// A single bond (bond order 1).
    #[default]
    Single,
    /// A double bond (bond order 2).
    Double,
    /// A triple bond (bond order 3).
    Triple,
}

/// Default bond type after kekulization.
///
/// Contains a single [`BondOrder`] field. This is the bond type used in the
/// standard `Mol<Atom, Bond>` that results from SMILES parsing.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    /// The bond order (single, double, or triple).
    pub order: BondOrder,
}

impl Default for Bond {
    fn default() -> Self {
        Self {
            order: BondOrder::Single,
        }
    }
}

impl crate::traits::HasBondOrder for Bond {
    fn bond_order(&self) -> BondOrder {
        self.order
    }
}

/// Bond order set used during SMILES parsing, before kekulization.
///
/// SMILES encodes bonds between lowercase (aromatic) atoms as [`Aromatic`](SmilesBondOrder::Aromatic),
/// and bonds written without an explicit symbol as [`Implicit`](SmilesBondOrder::Implicit)
/// (single between aliphatic atoms, aromatic between aromatic atoms).
/// Kekulization resolves these to concrete [`BondOrder`] values (single or
/// double).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SmilesBondOrder {
    /// Explicit single bond (`-` in SMILES).
    Single,
    /// Explicit double bond (`=` in SMILES).
    Double,
    /// Explicit triple bond (`#` in SMILES).
    Triple,
    /// Bond between two aromatic atoms (lowercase letters in SMILES).
    Aromatic,
    /// No explicit bond symbol between adjacent atoms.
    ///
    /// Resolves to single for aliphatic pairs, aromatic for aromatic pairs.
    #[default]
    Implicit,
}

/// Pre-kekulization bond type used during SMILES parsing.
///
/// This is the bond type in `Mol<Atom, SmilesBond>`, the intermediate
/// representation before kekulization converts it to `Mol<Atom, Bond>`.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesBond {
    /// The pre-kekulization bond order.
    pub order: SmilesBondOrder,
}

impl Default for SmilesBond {
    fn default() -> Self {
        Self {
            order: SmilesBondOrder::Implicit,
        }
    }
}
