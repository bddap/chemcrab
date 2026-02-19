use crate::atom::Chirality;
use crate::bond::BondStereo;
use crate::mol::Mol;
use crate::traits::{HasBondStereoMut, HasChiralityMut, HasIsotopeMut};

pub fn strip_chirality<A: HasChiralityMut, B>(mol: &mut Mol<A, B>) {
    for idx in mol.atoms().collect::<Vec<_>>() {
        *mol.atom_mut(idx).chirality_mut() = Chirality::None;
    }
}

pub fn strip_bond_stereo<A, B: HasBondStereoMut>(mol: &mut Mol<A, B>) {
    for idx in mol.bonds().collect::<Vec<_>>() {
        *mol.bond_mut(idx).bond_stereo_mut() = BondStereo::None;
    }
}

pub fn strip_isotope<A: HasIsotopeMut, B>(mol: &mut Mol<A, B>) {
    for idx in mol.atoms().collect::<Vec<_>>() {
        *mol.atom_mut(idx).isotope_mut() = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::{from_smiles, to_smiles};

    #[test]
    fn strip_chirality_removes_tetrahedral() {
        let mut mol = from_smiles("F[C@H](Cl)Br").unwrap();
        strip_chirality(&mut mol);
        for idx in mol.atoms() {
            assert_eq!(mol.atom(idx).chirality, Chirality::None);
        }
    }

    #[test]
    fn strip_bond_stereo_removes_ez() {
        let mut mol = from_smiles("F/C=C/F").unwrap();
        strip_bond_stereo(&mut mol);
        for idx in mol.bonds() {
            assert_eq!(mol.bond(idx).stereo, BondStereo::None);
        }
    }

    #[test]
    fn strip_isotope_clears_labels() {
        let mut mol = from_smiles("[13CH4]").unwrap();
        strip_isotope(&mut mol);
        for idx in mol.atoms() {
            assert_eq!(mol.atom(idx).isotope, 0);
        }
    }

    #[test]
    fn strip_all_produces_plain_smiles() {
        let mut mol = from_smiles("[13C@H](F)(Cl)Br").unwrap();
        strip_chirality(&mut mol);
        strip_bond_stereo(&mut mol);
        strip_isotope(&mut mol);
        let smiles = to_smiles(&mol);
        assert!(!smiles.contains('@'));
        assert!(!smiles.contains('/'));
        assert!(!smiles.contains('\\'));
        assert!(!smiles.contains("13"));
    }

    #[test]
    fn strip_chirality_on_achiral_is_noop() {
        let mut mol = from_smiles("CCO").unwrap();
        let before = to_smiles(&mol);
        strip_chirality(&mut mol);
        assert_eq!(to_smiles(&mol), before);
    }

    #[test]
    fn strip_bond_stereo_on_no_stereo_is_noop() {
        let mut mol = from_smiles("C=C").unwrap();
        let before = to_smiles(&mol);
        strip_bond_stereo(&mut mol);
        assert_eq!(to_smiles(&mol), before);
    }

    #[test]
    fn strip_isotope_on_natural_is_noop() {
        let mut mol = from_smiles("CCO").unwrap();
        let before = to_smiles(&mol);
        strip_isotope(&mut mol);
        assert_eq!(to_smiles(&mol), before);
    }

    #[test]
    fn strip_ez_from_complex() {
        let mut mol = from_smiles("F/C=C/C=C\\Cl").unwrap();
        strip_bond_stereo(&mut mol);
        let smiles = to_smiles(&mol);
        assert!(!smiles.contains('/'));
        assert!(!smiles.contains('\\'));
    }
}
