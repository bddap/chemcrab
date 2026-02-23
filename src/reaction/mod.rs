//! Reaction SMARTS — parse, apply, and write chemical reactions.
//!
//! Reaction SMARTS use `>>` to separate reactant and product templates,
//! with optional agents between single `>` separators. Atom map numbers
//! (`:1`, `:2`, …) link atoms across the transformation, so the engine
//! knows which reactant atom becomes which product atom.
//!
//! ```text
//! [C:1][Br:2].[OH-:3]>>[C:1][O:3].[Br-:2]
//! ```

pub mod error;
mod parser;
mod runner;
mod writer;

pub use error::{ReactionError, ReactionSmartsError};
pub use parser::parse_reaction_smarts;
pub use runner::extract_atom_map_num;
pub use writer::to_reaction_smarts;

use crate::mol::Mol;
use crate::smarts::{AtomExpr, BondExpr};

/// A parsed chemical reaction with reactant, product, and agent templates.
///
/// Construct via [`from_reaction_smarts`] and apply with [`Reaction::run`].
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    pub(crate) reactant_templates: Vec<Mol<AtomExpr, BondExpr>>,
    pub(crate) product_templates: Vec<Mol<AtomExpr, BondExpr>>,
    pub(crate) agent_templates: Vec<Mol<AtomExpr, BondExpr>>,
}

impl Reaction {
    /// Returns the reactant SMARTS templates.
    pub fn reactant_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.reactant_templates
    }

    /// Returns the product SMARTS templates.
    pub fn product_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.product_templates
    }

    /// Returns the agent SMARTS templates (catalysts, solvents, etc.).
    pub fn agent_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.agent_templates
    }
}

/// Parse a reaction SMARTS string into a [`Reaction`].
///
/// The string must contain `>>` (or `>agent>`) to separate reactants
/// from products. Components within each section are separated by `.`.
pub fn from_reaction_smarts(s: &str) -> Result<Reaction, ReactionSmartsError> {
    parse_reaction_smarts(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{AromaticBond, Bond};
    use crate::kekulize::kekulize;
    use crate::mol::Mol;
    use crate::smiles::from_smiles;

    fn mol(smiles: &str) -> Mol<Atom, Bond> {
        from_smiles(smiles).unwrap_or_else(|e| panic!("bad SMILES {smiles:?}: {e}"))
    }

    // --- Parsing tests ---

    #[test]
    fn parse_simple_reaction() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn parse_with_agents() {
        let rxn = from_reaction_smarts("[C:1]=[C:2]>[Pd]>[C:1][C:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 1);
    }

    #[test]
    fn parse_multi_reactant() {
        let rxn = from_reaction_smarts("[C:1]Br.[N:2]>>[C:1][N:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 2);
        assert_eq!(rxn.product_templates().len(), 1);
    }

    #[test]
    fn error_no_separator() {
        assert!(from_reaction_smarts("[C][Br]").is_err());
    }

    #[test]
    fn error_empty_products() {
        assert!(from_reaction_smarts("[C:1]>>").is_err());
    }

    #[test]
    fn error_empty_reactants() {
        assert!(from_reaction_smarts(">>[C:1]").is_err());
    }

    #[test]
    fn atom_map_class_parsed() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let tmpl = &rxn.reactant_templates()[0];
        let expr0 = tmpl.atom(petgraph::graph::NodeIndex::new(0));
        assert!(extract_atom_map_num(expr0).is_some());
        assert_eq!(extract_atom_map_num(expr0).unwrap(), 1);
    }

    #[test]
    fn multi_product_templates() {
        let rxn = from_reaction_smarts("[C:1][Br:2].[OH-:3]>>[C:1][O:3].[Br-:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 2);
        assert_eq!(rxn.product_templates().len(), 2);
    }

    // --- Application tests ---

    #[test]
    fn sn2_simple() {
        let rxn = from_reaction_smarts("[C:1][Br:2].[OH-:3]>>[C:1][O:3].[Br-:2]").unwrap();
        let r1 = mol("CBr");
        let r2 = mol("[OH-]");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert!(!products.is_empty());
        assert_eq!(products[0].len(), 2);
    }

    #[test]
    fn substituent_preservation() {
        let rxn = from_reaction_smarts("[C:1][Br:2].[OH-:3]>>[C:1][O:3].[Br-:2]").unwrap();
        let r1 = mol("CCBr");
        let r2 = mol("[OH-]");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert!(!products.is_empty());

        let product_mol = &products[0][0];
        // CCBr -> CC + OH- -> CCO (ethanol) - 3 heavy atoms
        assert_eq!(
            product_mol.atom_count(),
            3,
            "ethanol should have 3 heavy atoms"
        );
    }

    #[test]
    fn no_match_returns_empty() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert!(products.is_empty());
    }

    #[test]
    fn wrong_reactant_count_error() {
        let rxn = from_reaction_smarts("[C:1][Br:2].[OH-:3]>>[C:1][O:3]").unwrap();
        let r1 = mol("CBr");
        let result = rxn.run(&[&r1]);
        assert!(result.is_err());
    }

    #[test]
    fn simple_single_reactant_reaction() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let r = mol("CBr");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        assert_eq!(product.atom_count(), 2);
    }

    #[test]
    fn bond_formation() {
        let rxn = from_reaction_smarts("[C:1]Br.[N:2]>>[C:1][N:2]").unwrap();
        let r1 = mol("CBr");
        let r2 = mol("N");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        assert!(product.atom_count() >= 2);
    }

    #[test]
    fn double_bond_in_product() {
        let rxn = from_reaction_smarts("[C:1][C:2]>>[C:1]=[C:2]").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        assert_eq!(product.bond_count(), 1);
        let edge = product.bonds().next().unwrap();
        assert_eq!(
            product.bond(edge).order,
            crate::bond::AromaticBondOrder::Known(crate::bond::BondOrder::Double)
        );
    }

    // --- Writer tests ---

    #[test]
    fn writer_round_trip() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let s = to_reaction_smarts(&rxn);
        assert!(s.contains(">>"));
        let rxn2 = from_reaction_smarts(&s).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
    }

    #[test]
    fn writer_with_agents() {
        let rxn = from_reaction_smarts("[C:1]=[C:2]>[Pd]>[C:1][C:2]").unwrap();
        let s = to_reaction_smarts(&rxn);
        assert!(!s.contains(">>"));
        let parts: Vec<&str> = s.split('>').collect();
        assert_eq!(parts.len(), 3);
    }

    #[test]
    fn writer_multi_reactant_round_trip() {
        let rxn = from_reaction_smarts("[C:1]Br.[N:2]>>[C:1][N:2]").unwrap();
        let s = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&s).unwrap();
        assert_eq!(rxn2.reactant_templates().len(), 2);
        assert_eq!(rxn2.product_templates().len(), 1);
    }

    #[test]
    fn substituent_preservation_larger() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let r = mol("CCCBr");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        // CCCBr -> CCCO: C-C-C-O = 4 atoms
        assert_eq!(product.atom_count(), 4);
    }

    #[test]
    fn hydrogen_count_on_unmapped_product_atom() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let r = mol("CBr");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        let o_idx = product.atoms().find(|&i| product.atom(i).atomic_num == 8);
        assert!(o_idx.is_some());
    }

    #[test]
    fn charge_change_deprotonation() {
        let rxn = from_reaction_smarts("[OH:1]>>[O-:1]").unwrap();
        let r = mol("CO");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        let o_idx = product
            .atoms()
            .find(|&i| product.atom(i).atomic_num == 8)
            .unwrap();
        assert_eq!(product.atom(o_idx).formal_charge, -1);
    }

    #[test]
    fn unmapped_product_atom_created() {
        let rxn = from_reaction_smarts("[N:1]>>[N:1]C").unwrap();
        let r = mol("N");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        assert_eq!(product.atom_count(), 2);
        let c_idx = product.atoms().find(|&i| product.atom(i).atomic_num == 6);
        assert!(
            c_idx.is_some(),
            "unmapped carbon should be created in product"
        );
    }

    #[test]
    fn triple_bond_in_product() {
        let rxn = from_reaction_smarts("[C:1][C:2]>>[C:1]#[C:2]").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert!(!products.is_empty());
        let product = &products[0][0];
        let edge = product.bonds().next().unwrap();
        assert_eq!(
            product.bond(edge).order,
            crate::bond::AromaticBondOrder::Known(crate::bond::BondOrder::Triple)
        );
    }

    #[test]
    fn error_too_many_separators() {
        assert!(from_reaction_smarts("[C:1]>[A]>[B]>[C:1]").is_err());
    }

    #[test]
    fn multi_match_produces_multiple_results() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1][OH]").unwrap();
        let r = mol("BrCCBr");
        let products = rxn.run(&[&r]).unwrap();
        assert!(
            products.len() >= 2,
            "two C-Br bonds should give at least 2 match results"
        );
    }

    #[test]
    fn aromatic_substitution() {
        let rxn = from_reaction_smarts("[c:1][Br:2].[OH-:3]>>[c:1][O:3].[Br-:2]").unwrap();
        let reactant1 = mol("c1ccc(Br)cc1");
        let reactant2 = mol("[OH-]");
        let products = rxn.run(&[&reactant1, &reactant2]).unwrap();
        assert!(!products.is_empty());
    }

    #[test]
    fn aromatic_substitution_kekulized_bonds() {
        let rxn = from_reaction_smarts("[c:1][Br:2].[OH-:3]>>[c:1][O:3].[Br-:2]").unwrap();
        let reactant = mol("c1ccc(Br)cc1");
        let hydroxide = mol("[OH-]");
        let products = rxn.run(&[&reactant, &hydroxide]).unwrap();
        assert!(!products.is_empty());
        let phenol = kek(&products[0][0]);
        let double_count = phenol
            .bonds()
            .filter(|&e| phenol.bond(e).order == crate::bond::BondOrder::Double)
            .count();
        assert_eq!(
            double_count, 3,
            "phenol product ring should have 3 double bonds, not all single"
        );
    }

    #[test]
    fn aromatic_product_all_bonds_kekulized() {
        use crate::bond::BondOrder;
        let rxn = from_reaction_smarts("[c:1][Br:2].[OH-:3]>>[c:1][O:3].[Br-:2]").unwrap();
        let reactant = mol("c1ccc(Br)cc1");
        let hydroxide = mol("[OH-]");
        let products = rxn.run(&[&reactant, &hydroxide]).unwrap();
        assert!(!products.is_empty());
        let phenol = kek(&products[0][0]);
        for edge in phenol.bonds() {
            let order = phenol.bond(edge).order;
            assert!(
                order == BondOrder::Single
                    || order == BondOrder::Double
                    || order == BondOrder::Triple,
                "no bond should remain un-kekulized"
            );
        }
        let single_count = phenol
            .bonds()
            .filter(|&e| phenol.bond(e).order == BondOrder::Single)
            .count();
        let double_count = phenol
            .bonds()
            .filter(|&e| phenol.bond(e).order == BondOrder::Double)
            .count();
        assert!(
            double_count >= 3,
            "ring should have at least 3 double bonds, got {double_count}"
        );
        assert!(
            single_count >= 1,
            "should have at least 1 single bond (C-O), got {single_count}"
        );
    }

    // ================================================================
    // RDKit-generated test vectors
    // ================================================================

    use crate::smiles::to_canonical_smiles;

    fn canon(smiles: &str) -> String {
        to_canonical_smiles(&mol(smiles))
    }

    fn kek(m: &Mol<Atom, AromaticBond>) -> Mol<Atom, Bond> {
        let k = kekulize(m.clone()).expect("kekulization failed in test");
        crate::hydrogen::remove_hs(&k)
    }

    fn assert_product_sets(
        products: &[Vec<Mol<Atom, AromaticBond>>],
        expected_count: usize,
        expected_sets: &[&[&str]],
    ) {
        assert_eq!(
            products.len(),
            expected_count,
            "expected {expected_count} product sets, got {}",
            products.len()
        );
        for (i, (actual, expected)) in products.iter().zip(expected_sets.iter()).enumerate() {
            let mut actual_smiles: Vec<String> = actual
                .iter()
                .map(|m| to_canonical_smiles(&kek(m)))
                .collect();
            actual_smiles.sort();
            let mut expected_smiles: Vec<String> = expected.iter().map(|s| canon(s)).collect();
            expected_smiles.sort();
            assert_eq!(actual_smiles, expected_smiles, "product set {i} mismatch");
        }
    }

    // --- Parsing tests (RDKit vectors) ---

    #[test]
    fn test_rxn_parse_amide_formation() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 2);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_cyclization() {
        let rxn = from_reaction_smarts(
            "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]",
        )
        .unwrap();
        assert_eq!(rxn.reactant_templates().len(), 2);
        assert_eq!(rxn.product_templates().len(), 2);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_single_reactant_product() {
        let rxn = from_reaction_smarts("[C:1]-O>>[C:1]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_with_agents() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])O.[N:3]>CC>[C:1](=[O:2])[N:3]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 2);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 1);
    }

    #[test]
    fn test_rxn_parse_three_reactants() {
        let rxn = from_reaction_smarts("[A:1].[B:2].[C:3]>>[A:1][B:2][C:3]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 3);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_two_products() {
        let rxn = from_reaction_smarts("[C:1][C:2]>>[C:1].[C:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 2);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_aromatic_reactant() {
        let rxn = from_reaction_smarts("[c:1][c:2]>>[c:1][c:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_complex_smarts_atoms() {
        let rxn = from_reaction_smarts("[C;!$(C=O):1][OH:2]>>[C:1][O:2]C").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    #[test]
    fn test_rxn_parse_bond_not_in_reaction() {
        let rxn = from_reaction_smarts("[C:1]!@[C:2]>>[C:1][C:2]").unwrap();
        assert_eq!(rxn.reactant_templates().len(), 1);
        assert_eq!(rxn.product_templates().len(), 1);
        assert_eq!(rxn.agent_templates().len(), 0);
    }

    // --- Amide tests ---

    #[test]

    fn test_rxn_amide_simple() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        let r1 = mol("C(=O)O");
        let r2 = mol("CN");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["CNC=O"]]);
    }

    #[test]

    fn test_rxn_amide_acetic_acid() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        let r1 = mol("CC(=O)O");
        let r2 = mol("CN");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["CNC(C)=O"]]);
    }

    #[test]

    fn test_rxn_amide_two_acid_groups() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        let r1 = mol("CC(C(=O)O)C(=O)O");
        let r2 = mol("CN");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let mut smiles: Vec<String> =
                set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            smiles.sort();
            let mut expected: Vec<String> = vec![canon("CNC(=O)C(C)C(=O)O")];
            expected.sort();
            assert_eq!(smiles, expected);
        }
    }

    #[test]

    fn test_rxn_amide_two_amine_groups() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        let r1 = mol("CC(=O)O");
        let r2 = mol("NCN");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let mut smiles: Vec<String> =
                set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            smiles.sort();
            let mut expected: Vec<String> = vec![canon("CC(=O)NCN")];
            expected.sort();
            assert_eq!(smiles, expected);
        }
    }

    // --- Simple transform tests ---

    #[test]

    fn test_rxn_dehydroxylation() {
        let rxn = from_reaction_smarts("[C:1]-O>>[C:1]").unwrap();
        let r = mol("CCO");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC"]]);
    }

    #[test]

    fn test_rxn_oh_to_sh() {
        let rxn = from_reaction_smarts("[C:1]-O>>[C:1]-S").unwrap();
        let r = mol("CCO");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CCS"]]);
    }

    #[test]
    fn test_rxn_aldehyde_to_thioaldehyde() {
        let rxn = from_reaction_smarts("[C:2][C:1]=O>>[C:2][C:1]=S").unwrap();
        let r = mol("CC=O");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC=S"]]);
    }

    #[test]

    fn test_rxn_carboxylic_acid_reduction() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])[OH]>>[C:1][OH:2]").unwrap();
        let r = mol("CC(=O)O");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CCO"]]);
    }

    // --- Substituent tests ---

    #[test]

    fn test_rxn_carry_methyl() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[N:3]C").unwrap();
        let r = mol("CC(=O)NC");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC(=O)N(C)C"]]);
    }

    #[test]

    fn test_rxn_carry_ring() {
        let rxn = from_reaction_smarts("[c:1][NH2:2]>>[c:1][N:2](C)C").unwrap();
        let r = mol("c1ccc(N)cc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CN(C)c1ccccc1"]]);
    }

    #[test]

    fn test_rxn_carry_large_substituent() {
        let rxn = from_reaction_smarts("[C:1](=O)[OH]>>[C:1](=O)OCC").unwrap();
        let r = mol("c1ccc(C(=O)O)cc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CCOC(=O)c1ccccc1"]]);
    }

    // --- Bond tests ---

    #[test]

    fn test_rxn_bond_formation_c_n() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]").unwrap();
        let r1 = mol("CC(=O)O");
        let r2 = mol("NC");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["CNC(C)=O"]]);
    }

    #[test]

    fn test_rxn_bond_breaking_ester() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]").unwrap();
        let r = mol("CC(=O)OC");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC(=O)O", "CO"]]);
    }

    #[test]

    fn test_rxn_lactol_formation() {
        let rxn = from_reaction_smarts(
            "[OH:1][C:2][C:3][C:4][C:5]=[O:6]>>[O:1]1[C:2][C:3][C:4][C:5]1[OH:6]",
        )
        .unwrap();
        let r = mol("OCCCC=O");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["OC1CCCO1"]]);
    }

    // --- Charge tests ---

    #[test]
    fn test_rxn_deprotonation() {
        let rxn = from_reaction_smarts("[NH2:1]>>[NH-:1]").unwrap();
        let r = mol("CN");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["C[NH-]"]]);
    }

    #[test]
    fn test_rxn_protonation() {
        let rxn = from_reaction_smarts("[N:1]>>[NH+:1]").unwrap();
        let r = mol("NC");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 1);
        let product = &products[0][0];
        let n_idx = product
            .atoms()
            .find(|&i| product.atom(i).atomic_num == 7)
            .unwrap();
        assert_eq!(product.atom(n_idx).formal_charge, 1);
    }

    #[test]

    fn test_rxn_quaternization() {
        let rxn =
            from_reaction_smarts("[N:1]([C:2])([C:3])[C:4]>>[N+:1]([C:2])([C:3])([C:4])C").unwrap();
        let r = mol("CN(C)C");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 6);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("C[N+](C)(C)C"));
        }
    }

    // --- Aromatic tests ---

    #[test]
    fn test_rxn_aromatic_substitution_simple() {
        let rxn = from_reaction_smarts("[c:1][H]>>[c:1]Br").unwrap();
        let r = mol("c1ccccc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 6);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("Brc1ccccc1"));
        }
    }

    #[test]

    fn test_rxn_aromatic_nitration() {
        let rxn = from_reaction_smarts("[cH:1]>>[c:1][N+](=O)[O-]").unwrap();
        let r = mol("c1ccccc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 6);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("O=[N+]([O-])c1ccccc1"));
        }
    }

    #[test]

    fn test_rxn_phenol_methylation() {
        let rxn = from_reaction_smarts("[c:1][OH:2]>>[c:1][O:2]C").unwrap();
        let r = mol("c1ccc(O)cc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["COc1ccccc1"]]);
    }

    // --- Multi-product tests ---

    #[test]

    fn test_rxn_ester_hydrolysis_two_products() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[C:4][O:3]").unwrap();
        let r = mol("CC(=O)OCC");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC(=O)O", "CCO"]]);
    }

    // --- Multi-match tests ---

    #[test]

    fn test_rxn_disubstituted_benzene() {
        let rxn = from_reaction_smarts("[cH:1]>>[c:1]Cl").unwrap();
        let r = mol("c1ccccc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 6);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("Clc1ccccc1"));
        }
    }

    #[test]
    fn test_rxn_multiple_oh_groups() {
        let rxn = from_reaction_smarts("[C:1][OH:2]>>[C:1][O:2]C(=O)C").unwrap();
        let r = mol("OCC(O)CO");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 3);
    }

    // --- Isotope tests ---

    #[test]
    fn test_rxn_isotope_labeling() {
        let rxn = from_reaction_smarts("[C:1]>>[13C:1]").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("C[13CH3]"));
        }
    }

    #[test]
    fn test_rxn_deuteration() {
        let rxn = from_reaction_smarts("[C:1][H]>>[C:1]").unwrap();
        let r = mol("C");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 4);
    }

    // --- Unmapped atom tests ---

    #[test]

    fn test_rxn_add_methyl() {
        let rxn = from_reaction_smarts("[N:1]>>[N:1]C").unwrap();
        let r = mol("NC");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CNC"]]);
    }

    #[test]

    fn test_rxn_add_functional_group() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])O>>[C:1](=[O:2])NCCN").unwrap();
        let r = mol("CC(=O)O");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CC(=O)NCCN"]]);
    }

    // --- Hcount tests ---

    #[test]

    fn test_rxn_explicit_h_in_product() {
        let rxn = from_reaction_smarts("[NH2:1]>>[NH:1]C").unwrap();
        let r = mol("Nc1ccccc1");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CNc1ccccc1"]]);
    }

    #[test]

    fn test_rxn_remove_h() {
        let rxn = from_reaction_smarts("[CH3:1]>>[CH2:1]O").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("CCO"));
        }
    }

    // --- Roundtrip tests ---

    #[test]
    fn test_rxn_roundtrip_0() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])O.[N:3][C:4]>>[C:1](=[O:2])[N:3][C:4]").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    #[test]
    fn test_rxn_roundtrip_1() {
        let rxn = from_reaction_smarts("[C:1]-O>>[C:1]").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    #[test]
    fn test_rxn_roundtrip_2() {
        let rxn = from_reaction_smarts("[c:1][OH:2]>>[c:1][O:2]C").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    #[test]
    fn test_rxn_roundtrip_3() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    #[test]
    fn test_rxn_roundtrip_4() {
        let rxn = from_reaction_smarts("[N:1]>>[N:1]C").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    #[test]
    fn test_rxn_roundtrip_5() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])O.[N:3]>CC>[C:1](=[O:2])[N:3]").unwrap();
        let written = to_reaction_smarts(&rxn);
        let rxn2 = from_reaction_smarts(&written).unwrap();
        assert_eq!(
            rxn2.reactant_templates().len(),
            rxn.reactant_templates().len()
        );
        assert_eq!(
            rxn2.product_templates().len(),
            rxn.product_templates().len()
        );
        assert_eq!(rxn2.agent_templates().len(), rxn.agent_templates().len());
    }

    // --- Edge case tests ---

    #[test]
    fn test_rxn_identity_reaction() {
        let rxn = from_reaction_smarts("[C:1]>>[C:1]").unwrap();
        let r = mol("C");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["C"]]);
    }

    #[test]
    fn test_rxn_no_match() {
        let rxn = from_reaction_smarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]").unwrap();
        let r1 = mol("CCCCC");
        let r2 = mol("NC");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_eq!(products.len(), 0);
    }

    #[test]
    fn test_rxn_single_atom_reactant() {
        let rxn = from_reaction_smarts("[Br:1]>>[Cl:1]").unwrap();
        let r = mol("Br");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["Cl"]]);
    }

    #[test]

    fn test_rxn_empty_product_bonds() {
        let rxn = from_reaction_smarts("[C:1][C:2]>>[C:1].[C:2]").unwrap();
        let r = mol("CC");
        let products = rxn.run(&[&r]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let mut smiles: Vec<String> =
                set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            smiles.sort();
            let mut expected: Vec<String> = vec![canon("C"), canon("C")];
            expected.sort();
            assert_eq!(smiles, expected);
        }
    }

    // --- Named reaction tests ---

    #[test]

    fn test_rxn_fischer_esterification() {
        let rxn =
            from_reaction_smarts("[C:1](=[O:2])[OH].[OH:3][C:4]>>[C:1](=[O:2])[O:3][C:4]").unwrap();
        let r1 = mol("CC(=O)O");
        let r2 = mol("CO");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["COC(C)=O"]]);
    }

    #[test]

    fn test_rxn_williamson_ether() {
        let rxn = from_reaction_smarts("[C:1][OH:2].[Cl:3][C:4]>>[C:1][O:2][C:4]").unwrap();
        let r1 = mol("CO");
        let r2 = mol("CCl");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["COC"]]);
    }

    #[test]
    fn test_rxn_n_alkylation() {
        let rxn = from_reaction_smarts("[NH2:1].[C:2][Cl]>>[NH:1][C:2]").unwrap();
        let r1 = mol("NC");
        let r2 = mol("CCl");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_product_sets(&products, 1, &[&["CNC"]]);
    }

    #[test]

    fn test_rxn_suzuki_coupling_simplified() {
        let rxn = from_reaction_smarts("[c:1][Br].[c:2]B(O)O>>[c:1][c:2]").unwrap();
        let r1 = mol("c1ccc(Br)cc1");
        let r2 = mol("c1ccc(B(O)O)cc1");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert_eq!(products.len(), 2);
        for set in &products {
            let smiles: Vec<String> = set.iter().map(|m| to_canonical_smiles(&kek(m))).collect();
            assert_eq!(smiles.len(), 1);
            assert_eq!(smiles[0], canon("c1ccc(-c2ccccc2)cc1"));
        }
    }

    #[test]
    fn test_rxn_sn2_bromide_to_cyanide() {
        let rxn = from_reaction_smarts("[C:1][Br:2]>>[C:1]C#N").unwrap();
        let r = mol("CCBr");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CCC#N"]]);
    }

    #[test]

    fn test_rxn_boc_protection() {
        let rxn = from_reaction_smarts("[NH2:1]>>[NH:1]C(=O)OC(C)(C)C").unwrap();
        let r = mol("NCC");
        let products = rxn.run(&[&r]).unwrap();
        assert_product_sets(&products, 1, &[&["CCNC(=O)OC(C)(C)C"]]);
    }

    #[test]
    fn parenthesized_component_group_runs() {
        let rxn = from_reaction_smarts("([C:1](=O)O.[N:2])>>[C:1](=O)[N:2]").unwrap();
        let r1 = mol("CC(=O)O");
        let r2 = mol("N");
        let products = rxn.run(&[&r1, &r2]).unwrap();
        assert!(!products.is_empty());
    }
}
