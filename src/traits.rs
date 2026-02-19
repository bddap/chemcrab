use crate::atom::Chirality;
use crate::bond::{BondOrder, BondStereo};
use crate::wrappers::Hybridization;

pub trait HasAtomicNum {
    fn atomic_num(&self) -> u8;
}

pub trait HasFormalCharge {
    fn formal_charge(&self) -> i8;
}

pub trait HasIsotope {
    fn isotope(&self) -> u16;
}

pub trait HasChirality {
    fn chirality(&self) -> Chirality;
}

pub trait HasHydrogenCount {
    fn hydrogen_count(&self) -> u8;
}

pub trait HasAromaticity {
    fn is_aromatic(&self) -> bool;
}

pub trait HasPosition2D {
    fn position_2d(&self) -> Option<[f64; 2]>;
    fn set_position_2d(&mut self, pos: Option<[f64; 2]>);
}

pub trait HasPosition3D {
    fn position_3d(&self) -> Option<[f64; 3]>;
    fn set_position_3d(&mut self, pos: Option<[f64; 3]>);
}

pub trait HasBondOrder {
    fn bond_order(&self) -> BondOrder;
}

pub trait HasBondStereo {
    fn bond_stereo(&self) -> BondStereo;
}

pub trait HasValence {
    fn valence(&self) -> u8;
}

pub trait HasHybridization {
    fn hybridization(&self) -> Hybridization;
}
