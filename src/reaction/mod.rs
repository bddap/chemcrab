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

#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    pub(crate) reactant_templates: Vec<Mol<AtomExpr, BondExpr>>,
    pub(crate) product_templates: Vec<Mol<AtomExpr, BondExpr>>,
    pub(crate) agent_templates: Vec<Mol<AtomExpr, BondExpr>>,
}

impl Reaction {
    pub fn reactant_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.reactant_templates
    }

    pub fn product_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.product_templates
    }

    pub fn agent_templates(&self) -> &[Mol<AtomExpr, BondExpr>] {
        &self.agent_templates
    }
}

pub fn from_reaction_smarts(s: &str) -> Result<Reaction, ReactionSmartsError> {
    parse_reaction_smarts(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;
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
        assert_eq!(product.bond(edge).order, crate::bond::BondOrder::Double);
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
        assert_eq!(product.bond(edge).order, crate::bond::BondOrder::Triple);
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
        let phenol = &products[0][0];
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
        let phenol = &products[0][0];
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
}
