use crate::bond::BondOrder;
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

pub trait HasIsotopeMut: HasIsotope {
    fn isotope_mut(&mut self) -> &mut u16;
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

pub trait HasValence {
    fn valence(&self) -> u8;
}

pub trait HasHybridization {
    fn hybridization(&self) -> Hybridization;
}
