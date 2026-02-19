mod builder;
pub mod error;
mod parse_tree;
mod tokenizer;
mod writer;

use crate::atom::Atom;
use crate::bond::{Bond, SmilesBond};
use crate::kekulize;
use crate::mol::Mol;
pub use error::SmilesError;
pub use writer::to_smiles;

pub fn parse_smiles(s: &str) -> Result<Mol<Atom, SmilesBond>, SmilesError> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Err(SmilesError::EmptyInput);
    }
    let tokens = tokenizer::tokenize(trimmed)?;
    if tokens.is_empty() {
        return Err(SmilesError::EmptyInput);
    }
    let tree = parse_tree::build_parse_tree(&tokens)?;
    Ok(builder::build_mol(&tree))
}

pub fn from_smiles(s: &str) -> Result<Mol<Atom, Bond>, SmilesError> {
    let mol = parse_smiles(s)?;
    Ok(kekulize::kekulize(mol)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Chirality;
    use crate::bond::{BondStereo, SmilesBondOrder};
    use petgraph::graph::NodeIndex;

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    fn atom(mol: &Mol<Atom, SmilesBond>, i: usize) -> &Atom {
        mol.atom(n(i))
    }

    // ---- Simple molecules ----

    #[test]
    fn methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(atom(&mol, 0).atomic_num, 6);
        assert_eq!(atom(&mol, 0).hydrogen_count, 4);
    }

    #[test]
    fn ethane() {
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
        assert_eq!(atom(&mol, 1).hydrogen_count, 3);
    }

    #[test]
    fn ethene() {
        let mol = parse_smiles("C=C").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 2);
        assert_eq!(atom(&mol, 1).hydrogen_count, 2);
        let edge = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(edge).order, SmilesBondOrder::Double);
    }

    #[test]
    fn ethyne() {
        let mol = parse_smiles("C#C").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
        assert_eq!(atom(&mol, 1).hydrogen_count, 1);
        let edge = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(edge).order, SmilesBondOrder::Triple);
    }

    #[test]
    fn water_bare() {
        let mol = parse_smiles("O").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).atomic_num, 8);
        assert_eq!(atom(&mol, 0).hydrogen_count, 2);
    }

    #[test]
    fn ammonia_bare() {
        let mol = parse_smiles("N").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
    }

    #[test]
    fn hydrogen_fluoride() {
        let mol = parse_smiles("F").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    #[test]
    fn hydrogen_chloride() {
        let mol = parse_smiles("Cl").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).atomic_num, 17);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    #[test]
    fn hydrogen_bromide() {
        let mol = parse_smiles("Br").unwrap();
        assert_eq!(atom(&mol, 0).atomic_num, 35);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    #[test]
    fn methanol() {
        let mol = parse_smiles("CO").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
        assert_eq!(atom(&mol, 1).hydrogen_count, 1);
    }

    #[test]
    fn acetic_acid() {
        let mol = parse_smiles("CC(=O)O").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3); // CH3
        assert_eq!(atom(&mol, 1).hydrogen_count, 0); // C(=O)O
        assert_eq!(atom(&mol, 2).hydrogen_count, 0); // =O
        assert_eq!(atom(&mol, 3).hydrogen_count, 1); // OH
    }

    // ---- Branches ----

    #[test]
    fn isobutane() {
        let mol = parse_smiles("CC(C)C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
        assert_eq!(atom(&mol, 1).hydrogen_count, 1);
        assert_eq!(atom(&mol, 2).hydrogen_count, 3);
        assert_eq!(atom(&mol, 3).hydrogen_count, 3);
    }

    #[test]
    fn neopentane() {
        let mol = parse_smiles("CC(C)(C)C").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 4);
        assert_eq!(atom(&mol, 1).hydrogen_count, 0);
    }

    // ---- Ring closures ----

    #[test]
    fn cyclopropane() {
        let mol = parse_smiles("C1CC1").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 3);
        for i in 0..3 {
            assert_eq!(atom(&mol, i).hydrogen_count, 2);
        }
    }

    #[test]
    fn cyclohexane() {
        let mol = parse_smiles("C1CCCCC1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        for i in 0..6 {
            assert_eq!(atom(&mol, i).hydrogen_count, 2);
        }
    }

    #[test]
    fn multi_digit_ring() {
        let mol = parse_smiles("C%10CC%10").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 3);
    }

    // ---- Charges ----

    #[test]
    fn ammonium() {
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).atomic_num, 7);
        assert_eq!(atom(&mol, 0).formal_charge, 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 4);
    }

    #[test]
    fn oxide_anion() {
        let mol = parse_smiles("[O-]").unwrap();
        assert_eq!(atom(&mol, 0).formal_charge, -1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 0);
    }

    // ---- Isotopes ----

    #[test]
    fn carbon_13() {
        let mol = parse_smiles("[13C]").unwrap();
        assert_eq!(atom(&mol, 0).isotope, 13);
        assert_eq!(atom(&mol, 0).atomic_num, 6);
    }

    #[test]
    fn deuterium() {
        let mol = parse_smiles("[2H]").unwrap();
        assert_eq!(atom(&mol, 0).isotope, 2);
        assert_eq!(atom(&mol, 0).atomic_num, 1);
    }

    // ---- Bracket H counts ----

    #[test]
    fn bracket_ch4() {
        let mol = parse_smiles("[CH4]").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 4);
    }

    #[test]
    fn bracket_nh3() {
        let mol = parse_smiles("[NH3]").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
    }

    #[test]
    fn bracket_oh2() {
        let mol = parse_smiles("[OH2]").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 2);
    }

    // ---- Aromatic atoms ----

    #[test]
    fn benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        for i in 0..6 {
            assert!(atom(&mol, i).is_aromatic);
            assert_eq!(atom(&mol, i).hydrogen_count, 1);
        }
        for edge in mol.bonds() {
            assert_eq!(mol.bond(edge).order, SmilesBondOrder::Aromatic);
        }
    }

    #[test]
    fn pyridine() {
        let mol = parse_smiles("c1ccncc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        // N at index 3
        assert_eq!(atom(&mol, 3).atomic_num, 7);
        assert_eq!(atom(&mol, 3).hydrogen_count, 0);
        // C atoms each have 1 H
        for i in [0, 1, 2, 4, 5] {
            assert_eq!(atom(&mol, i).hydrogen_count, 1);
        }
    }

    #[test]
    fn furan() {
        let mol = parse_smiles("o1cccc1").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 5);
        assert_eq!(atom(&mol, 0).atomic_num, 8);
        assert_eq!(atom(&mol, 0).hydrogen_count, 0);
        for i in 1..5 {
            assert_eq!(atom(&mol, i).hydrogen_count, 1);
        }
    }

    #[test]
    fn pyrrole() {
        let mol = parse_smiles("[nH]1cccc1").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(atom(&mol, 0).atomic_num, 7);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
        for i in 1..5 {
            assert_eq!(atom(&mol, i).hydrogen_count, 1);
        }
    }

    // ---- Stereochemistry ----

    #[test]
    fn tetrahedral_ccw() {
        let mol = parse_smiles("[C@](F)(Cl)(Br)I").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(atom(&mol, 0).chirality, Chirality::Ccw);
    }

    #[test]
    fn tetrahedral_cw() {
        let mol = parse_smiles("[C@@](F)(Cl)(Br)I").unwrap();
        assert_eq!(atom(&mol, 0).chirality, Chirality::Cw);
    }

    #[test]
    fn tetrahedral_with_h() {
        let mol = parse_smiles("[C@@H](F)(Cl)Br").unwrap();
        assert_ne!(atom(&mol, 0).chirality, Chirality::None);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    #[test]
    fn ez_trans() {
        // F/C=C/F -> trans
        let mol = parse_smiles("F/C=C/F").unwrap();
        assert_eq!(mol.atom_count(), 4);
        let double_edge = mol.bond_between(n(1), n(2)).unwrap();
        assert_eq!(mol.bond(double_edge).order, SmilesBondOrder::Double);
        assert!(matches!(
            mol.bond(double_edge).stereo,
            BondStereo::Trans(_, _)
        ));
    }

    #[test]
    fn ez_cis() {
        // F/C=C\F -> cis
        let mol = parse_smiles(r"F/C=C\F").unwrap();
        let double_edge = mol.bond_between(n(1), n(2)).unwrap();
        assert!(matches!(
            mol.bond(double_edge).stereo,
            BondStereo::Cis(_, _)
        ));
    }

    #[test]
    fn ez_trans_chlorine() {
        let mol = parse_smiles(r"Cl/C=C/Cl").unwrap();
        let double_edge = mol.bond_between(n(1), n(2)).unwrap();
        assert!(matches!(
            mol.bond(double_edge).stereo,
            BondStereo::Trans(_, _)
        ));
    }

    #[test]
    fn ez_cis_chlorine() {
        let mol = parse_smiles(r"Cl/C=C\Cl").unwrap();
        let double_edge = mol.bond_between(n(1), n(2)).unwrap();
        assert!(matches!(
            mol.bond(double_edge).stereo,
            BondStereo::Cis(_, _)
        ));
    }

    // ---- Disconnected ----

    #[test]
    fn sodium_chloride() {
        let mol = parse_smiles("[Na+].[Cl-]").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(atom(&mol, 0).atomic_num, 11);
        assert_eq!(atom(&mol, 0).formal_charge, 1);
        assert_eq!(atom(&mol, 1).atomic_num, 17);
        assert_eq!(atom(&mol, 1).formal_charge, -1);
    }

    // ---- Error cases ----

    #[test]
    fn empty_string() {
        assert!(parse_smiles("").is_err());
    }

    #[test]
    fn whitespace_only() {
        assert!(parse_smiles("   ").is_err());
    }

    #[test]
    fn mismatched_paren_open() {
        assert!(parse_smiles("C(C").is_err());
    }

    #[test]
    fn mismatched_paren_close() {
        assert!(parse_smiles("C)C").is_err());
    }

    #[test]
    fn unclosed_ring() {
        assert!(parse_smiles("C1CC").is_err());
    }

    #[test]
    fn invalid_atom() {
        assert!(parse_smiles("X").is_err());
    }

    #[test]
    fn unclosed_bracket() {
        assert!(parse_smiles("[C").is_err());
    }

    // ---- H count correctness: more thorough tests ----

    #[test]
    fn phosphine() {
        let mol = parse_smiles("P").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
    }

    #[test]
    fn hydrogen_sulfide() {
        let mol = parse_smiles("S").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 2);
    }

    #[test]
    fn borane() {
        let mol = parse_smiles("B").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
    }

    #[test]
    fn molecular_hydrogen() {
        let mol = parse_smiles("[HH]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    #[test]
    fn explicit_single_bond() {
        let mol = parse_smiles("C-C").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        let edge = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(edge).order, SmilesBondOrder::Single);
        assert_eq!(atom(&mol, 0).hydrogen_count, 3);
    }

    // ---- Nitrogen with valence 5 ----

    #[test]
    fn nitro_group() {
        // C[N+](=O)[O-]
        let mol = parse_smiles("C[N+](=O)[O-]").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(atom(&mol, 1).atomic_num, 7);
        assert_eq!(atom(&mol, 1).formal_charge, 1);
    }

    // ---- Sulfur valences ----

    #[test]
    fn dmso() {
        // CS(=O)C
        let mol = parse_smiles("CS(=O)C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(atom(&mol, 1).hydrogen_count, 0); // S has 4 bonds used
    }

    // ---- Multiple ring closures on same atom ----

    #[test]
    fn bicyclo() {
        let mol = parse_smiles("C1CC2C1CC2").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 7);
    }

    // ---- Iodine ----

    #[test]
    fn iodine_hydrogen() {
        let mol = parse_smiles("I").unwrap();
        assert_eq!(atom(&mol, 0).hydrogen_count, 1);
    }

    // ---- Ring closure with bond type ----

    #[test]
    fn ring_with_double_bond() {
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    // ---- Complex molecules ----

    #[test]
    fn caffeine_atom_count() {
        let mol = parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
        assert_eq!(mol.atom_count(), 14);
    }

    #[test]
    fn naphthalene() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        assert_eq!(mol.atom_count(), 10);
        assert_eq!(mol.bond_count(), 11);
    }

    // ---- Bracket atom iron ----

    #[test]
    fn iron_bracket() {
        let mol = parse_smiles("[Fe]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(atom(&mol, 0).atomic_num, 26);
        assert_eq!(atom(&mol, 0).hydrogen_count, 0);
    }

    // ---- Atom class ----

    #[test]
    fn atom_class_preserved() {
        // We don't currently store atom_class on Atom, but parsing shouldn't fail
        let mol = parse_smiles("[C:1]").unwrap();
        assert_eq!(mol.atom_count(), 1);
    }

    // ---- Aromatic bond between aromatic and non-aromatic ----

    #[test]
    fn phenol() {
        let mol = parse_smiles("Oc1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 7);
        assert_eq!(atom(&mol, 0).hydrogen_count, 1); // OH
        let bond_o_c = mol.bond_between(n(0), n(1)).unwrap();
        assert_eq!(mol.bond(bond_o_c).order, SmilesBondOrder::Implicit);
    }

    // ---- Phosphorus with 5 bonds ----

    #[test]
    fn phosphorus_pentavalent() {
        let mol = parse_smiles("P(=O)(O)(O)O").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(atom(&mol, 0).hydrogen_count, 0);
    }

    // ---- Thiophene ----

    #[test]
    fn thiophene() {
        let mol = parse_smiles("s1cccc1").unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(atom(&mol, 0).atomic_num, 16);
        assert_eq!(atom(&mol, 0).hydrogen_count, 0);
        for i in 1..5 {
            assert_eq!(atom(&mol, i).hydrogen_count, 1);
        }
    }
}
