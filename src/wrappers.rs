use crate::traits::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Hybridization {
    S,
    SP,
    SP2,
    #[default]
    SP3,
    SP3D,
    SP3D2,
    Other,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WithValence<T> {
    pub inner: T,
    pub valence: u8,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WithHybridization<T> {
    pub inner: T,
    pub hybridization: Hybridization,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WithPosition2D<T> {
    pub inner: T,
    pub position_2d: Option<[f64; 2]>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WithPosition3D<T> {
    pub inner: T,
    pub position_3d: Option<[f64; 3]>,
}

impl<T> HasValence for WithValence<T> {
    fn valence(&self) -> u8 {
        self.valence
    }
}

impl<T> HasHybridization for WithHybridization<T> {
    fn hybridization(&self) -> Hybridization {
        self.hybridization
    }
}

impl<T> HasPosition2D for WithPosition2D<T> {
    fn position_2d(&self) -> Option<[f64; 2]> {
        self.position_2d
    }
    fn set_position_2d(&mut self, pos: Option<[f64; 2]>) {
        self.position_2d = pos;
    }
}

impl<T> HasPosition3D for WithPosition3D<T> {
    fn position_3d(&self) -> Option<[f64; 3]> {
        self.position_3d
    }
    fn set_position_3d(&mut self, pos: Option<[f64; 3]>) {
        self.position_3d = pos;
    }
}

macro_rules! delegate_trait {
    ($wrapper:ident, $trait:ident, $method:ident, $ret:ty) => {
        impl<T: $trait> $trait for $wrapper<T> {
            fn $method(&self) -> $ret {
                self.inner.$method()
            }
        }
    };
}

macro_rules! delegate_position_2d {
    ($wrapper:ident) => {
        impl<T: HasPosition2D> HasPosition2D for $wrapper<T> {
            fn position_2d(&self) -> Option<[f64; 2]> {
                self.inner.position_2d()
            }
            fn set_position_2d(&mut self, pos: Option<[f64; 2]>) {
                self.inner.set_position_2d(pos);
            }
        }
    };
}

macro_rules! delegate_position_3d {
    ($wrapper:ident) => {
        impl<T: HasPosition3D> HasPosition3D for $wrapper<T> {
            fn position_3d(&self) -> Option<[f64; 3]> {
                self.inner.position_3d()
            }
            fn set_position_3d(&mut self, pos: Option<[f64; 3]>) {
                self.inner.set_position_3d(pos);
            }
        }
    };
}

macro_rules! delegate_common {
    ($wrapper:ident) => {
        delegate_trait!($wrapper, HasAtomicNum, atomic_num, u8);
        delegate_trait!($wrapper, HasFormalCharge, formal_charge, i8);
        delegate_trait!($wrapper, HasIsotope, isotope, u16);
        delegate_trait!($wrapper, HasHydrogenCount, hydrogen_count, u8);
        delegate_trait!($wrapper, HasAromaticity, is_aromatic, bool);
        delegate_trait!($wrapper, HasBondOrder, bond_order, crate::bond::BondOrder);
    };
}

delegate_common!(WithValence);
delegate_trait!(WithValence, HasHybridization, hybridization, Hybridization);
delegate_position_2d!(WithValence);
delegate_position_3d!(WithValence);

delegate_common!(WithHybridization);
delegate_trait!(WithHybridization, HasValence, valence, u8);
delegate_position_2d!(WithHybridization);
delegate_position_3d!(WithHybridization);

delegate_common!(WithPosition2D);
delegate_trait!(WithPosition2D, HasValence, valence, u8);
delegate_trait!(
    WithPosition2D,
    HasHybridization,
    hybridization,
    Hybridization
);
delegate_position_3d!(WithPosition2D);

delegate_common!(WithPosition3D);
delegate_trait!(WithPosition3D, HasValence, valence, u8);
delegate_trait!(
    WithPosition3D,
    HasHybridization,
    hybridization,
    Hybridization
);
delegate_position_2d!(WithPosition3D);
