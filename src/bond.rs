#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondOrder {
    #[default]
    Single,
    Double,
    Triple,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SmilesBondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    #[default]
    Implicit,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SmilesBond {
    pub order: SmilesBondOrder,
}

impl Default for SmilesBond {
    fn default() -> Self {
        Self {
            order: SmilesBondOrder::Implicit,
        }
    }
}
