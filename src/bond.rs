/// Concrete bond order after kekulization.
///
/// Every bond in a [`Mol<Atom, Bond>`](crate::Mol) has one of these three
/// orders. There is no `Aromatic` variant — aromatic bonds only exist in the
/// intermediate [`AromaticBond`] type before kekulization resolves them to
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

/// Bond order that may contain unresolved aromatic bonds.
///
/// Before kekulization, bonds between aromatic atoms are represented as
/// [`Aromatic`](Self::Aromatic). After kekulization, every bond has a
/// concrete [`BondOrder`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AromaticBondOrder {
    /// A resolved bond order (single, double, or triple).
    Known(BondOrder),
    /// An aromatic bond awaiting kekulization.
    Aromatic,
}

/// Pre-kekulization bond type.
///
/// This is the bond type in molecules that may still contain aromatic bonds.
/// Call [`kekulize`](crate::kekulize::kekulize) to resolve aromatic bonds
/// and obtain a `Mol<Atom, Bond>`.
#[derive(Debug, Clone, PartialEq)]
pub struct AromaticBond {
    /// The bond order, which may be aromatic.
    pub order: AromaticBondOrder,
}

impl Default for AromaticBond {
    fn default() -> Self {
        Self {
            order: AromaticBondOrder::Known(BondOrder::Single),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub(crate) enum SmilesBondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    #[default]
    Implicit,
}
