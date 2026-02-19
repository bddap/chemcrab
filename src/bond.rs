use petgraph::graph::NodeIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondOrder {
    #[default]
    Single,
    Double,
    Triple,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum BondStereo {
    #[default]
    None,
    Cis(NodeIndex, NodeIndex),
    Trans(NodeIndex, NodeIndex),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub order: BondOrder,
    pub stereo: BondStereo,
}

impl Default for Bond {
    fn default() -> Self {
        Self {
            order: BondOrder::Single,
            stereo: BondStereo::None,
        }
    }
}

impl crate::traits::HasBondOrder for Bond {
    fn bond_order(&self) -> BondOrder {
        self.order
    }
}

impl crate::traits::HasBondStereo for Bond {
    fn bond_stereo(&self) -> BondStereo {
        self.stereo
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
    pub stereo: BondStereo,
}

impl Default for SmilesBond {
    fn default() -> Self {
        Self {
            order: SmilesBondOrder::Implicit,
            stereo: BondStereo::None,
        }
    }
}

impl crate::traits::HasBondStereo for SmilesBond {
    fn bond_stereo(&self) -> BondStereo {
        self.stereo
    }
}
