use crate::smarts::to_smarts;

use super::Reaction;

pub fn to_reaction_smarts(rxn: &Reaction) -> String {
    let reactants: Vec<String> = rxn
        .reactant_templates
        .iter()
        .map(to_smarts)
        .collect();
    let products: Vec<String> = rxn
        .product_templates
        .iter()
        .map(to_smarts)
        .collect();

    if rxn.agent_templates.is_empty() {
        format!("{}>>{}", reactants.join("."), products.join("."))
    } else {
        let agents: Vec<String> = rxn
            .agent_templates
            .iter()
            .map(to_smarts)
            .collect();
        format!(
            "{}>{}>{}", reactants.join("."), agents.join("."), products.join(".")
        )
    }
}
