//! Atom and bond property traits for generic algorithms.
//!
//! These traits abstract over atom and bond properties so that algorithms
//! can be written generically. A function that needs to know atomic numbers
//! declares `A: HasAtomicNum`; one that needs bond orders declares
//! `B: HasBondOrder`. This lets custom atom and bond types work with every
//! algorithm in the library without modification.

use crate::bond::BondOrder;
use crate::wrappers::Hybridization;

/// Provides the atomic number (element identity).
pub trait HasAtomicNum {
    fn atomic_num(&self) -> u8;
}

/// Provides the formal charge in elementary charge units.
pub trait HasFormalCharge {
    fn formal_charge(&self) -> i8;
}

/// Provides the isotope mass number (0 means natural abundance).
pub trait HasIsotope {
    fn isotope(&self) -> u16;
}

/// Mutable access to the isotope label.
pub trait HasIsotopeMut: HasIsotope {
    fn isotope_mut(&mut self) -> &mut u16;
}

/// Provides the virtual (suppressed) hydrogen count.
pub trait HasHydrogenCount {
    fn hydrogen_count(&self) -> u8;
}

/// Whether the atom is aromatic.
pub trait HasAromaticity {
    fn is_aromatic(&self) -> bool;
}

/// 2D coordinates for depiction.
pub trait HasPosition2D {
    fn position_2d(&self) -> Option<[f64; 2]>;
    fn set_position_2d(&mut self, pos: Option<[f64; 2]>);
}

/// 3D coordinates for conformers.
pub trait HasPosition3D {
    fn position_3d(&self) -> Option<[f64; 3]>;
    fn set_position_3d(&mut self, pos: Option<[f64; 3]>);
}

/// Provides the bond order (single, double, or triple).
pub trait HasBondOrder {
    fn bond_order(&self) -> BondOrder;
}

/// Provides the precomputed total valence.
pub trait HasValence {
    fn valence(&self) -> u8;
}

/// Provides the orbital hybridization state.
pub trait HasHybridization {
    fn hybridization(&self) -> Hybridization;
}
