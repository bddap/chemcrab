#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    #[default]
    None,
    Cw,
    Ccw,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub atomic_num: u8,
    pub formal_charge: i8,
    /// 0 means natural abundance.
    pub isotope: u16,
    pub chirality: Chirality,
    /// Virtual hydrogen count: the total number of suppressed (non-graph-node)
    /// hydrogens on this atom. Single source of truth — no implicit hydrogen
    /// rules are consulted at read time; the parser resolves hydrogen counts at
    /// parse time and stores the result.
    ///
    /// Credit to Richard L. Apodaca for articulating this model:
    /// <https://depth-first.com/articles/2019/11/06/virtual-hydrogens/>
    pub hydrogen_count: u8,
    pub is_aromatic: bool,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            atomic_num: 0,
            formal_charge: 0,
            isotope: 0,
            chirality: Chirality::None,
            hydrogen_count: 0,
            is_aromatic: false,
        }
    }
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

impl crate::traits::HasChirality for Atom {
    fn chirality(&self) -> Chirality {
        self.chirality
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
