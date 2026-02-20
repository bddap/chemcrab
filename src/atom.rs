#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    #[default]
    None,
    Cw,
    Ccw,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Atom {
    pub atomic_num: u8,
    pub formal_charge: i8,
    /// 0 means natural abundance.
    pub isotope: u16,
    pub hydrogen_count: u8,
    pub is_aromatic: bool,
}

impl crate::traits::HasAtomicNum for Atom {
    fn atomic_num(&self) -> u8 {
        self.atomic_num
    }
}

impl crate::traits::HasFormalCharge for Atom {
    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

impl crate::traits::HasIsotope for Atom {
    fn isotope(&self) -> u16 {
        self.isotope
    }
}

impl crate::traits::HasIsotopeMut for Atom {
    fn isotope_mut(&mut self) -> &mut u16 {
        &mut self.isotope
    }
}

impl crate::traits::HasHydrogenCount for Atom {
    fn hydrogen_count(&self) -> u8 {
        self.hydrogen_count
    }
}

impl crate::traits::HasAromaticity for Atom {
    fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }
}
