/// Chirality tag used in SMARTS query patterns.
///
/// This enum appears in SMARTS atom expressions to match stereocenters by
/// handedness. It does **not** carry stereochemistry for actual molecules —
/// that is stored in [`TetrahedralStereo`](crate::TetrahedralStereo) on the
/// [`Mol`](crate::Mol).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Chirality {
    /// No chirality constraint.
    #[default]
    None,
    /// Clockwise (@@) arrangement.
    Cw,
    /// Counterclockwise (@) arrangement.
    Ccw,
}

/// Default atom type for a molecular graph node.
///
/// `Atom` stores intrinsic atomic properties — the things you would read off
/// a structural formula. It deliberately omits computed properties like
/// valence, hybridization, or coordinates; those are provided by wrapper
/// types in the [`wrappers`](crate::wrappers) module when needed.
///
/// # Examples
///
/// ```
/// use chemcrab::Atom;
///
/// let carbon = Atom {
///     atomic_num: 6,
///     formal_charge: 0,
///     isotope: 0,
///     hydrogen_count: 3,
///     is_aromatic: false,
/// };
/// assert_eq!(carbon.atomic_num, 6);
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Atom {
    /// Atomic number (1 = H, 6 = C, 7 = N, …). Identifies the element.
    pub atomic_num: u8,
    /// Formal charge in elementary charge units (e.g. −1 for a carboxylate oxygen).
    pub formal_charge: i8,
    /// Mass number. `0` means natural isotopic abundance (the common case).
    pub isotope: u16,
    /// Number of virtual (suppressed) hydrogens on this atom.
    ///
    /// These are not graph nodes — they are implied by the atom's valence.
    /// After SMILES parsing, this count is the single source of truth for
    /// how many Hs the atom carries.
    pub hydrogen_count: u8,
    /// Whether this atom is in an aromatic ring.
    ///
    /// Set during aromaticity perception. Does not affect bond orders — all
    /// bonds in a [`Mol<Atom, Bond>`](crate::Mol) have concrete Kekulé
    /// assignments regardless of this flag.
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
