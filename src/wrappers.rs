//! Newtype wrappers that enrich atom types with additional computed properties.
//!
//! Rather than storing every possible property on [`Atom`](crate::Atom),
//! chemcrab uses wrapper types: `WithValence<Atom>` adds a precomputed
//! valence, `WithHybridization<Atom>` adds hybridization, etc. Wrappers
//! compose freely and delegate all inner traits automatically via blanket
//! impls, so `WithHybridization<WithValence<Atom>>` implements every
//! trait that `Atom` does, plus `HasValence` and `HasHybridization`.

use crate::traits::*;

/// Orbital hybridization of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Hybridization {
    /// Pure s orbital (bare ion, no bonds).
    S,
    /// sp — linear geometry, 180° bond angles.
    SP,
    /// sp2 — trigonal planar, 120° bond angles.
    SP2,
    /// sp3 — tetrahedral, 109.5° bond angles.
    #[default]
    SP3,
    /// sp3d — trigonal bipyramidal, 5 electron domains.
    SP3D,
    /// sp3d2 — octahedral, 6 electron domains.
    SP3D2,
    /// Hybridization could not be determined.
    Other,
}

/// Wraps an atom type with a precomputed total valence.
#[derive(Debug, Clone, PartialEq)]
pub struct WithValence<T> {
    pub inner: T,
    pub valence: u8,
}

/// Wraps an atom type with a computed hybridization.
#[derive(Debug, Clone, PartialEq)]
pub struct WithHybridization<T> {
    pub inner: T,
    pub hybridization: Hybridization,
}

/// Wraps an atom type with optional 2D coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct WithPosition2D<T> {
    pub inner: T,
    pub position_2d: Option<[f64; 2]>,
}

/// Wraps an atom type with optional 3D coordinates.
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
