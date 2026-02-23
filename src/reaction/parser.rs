use crate::mol::Mol;
use crate::smarts::{from_smarts, AtomExpr, BondExpr};

use super::error::ReactionSmartsError;
use super::Reaction;

pub fn parse_reaction_smarts(s: &str) -> Result<Reaction, ReactionSmartsError> {
    let (reactant_text, agent_text, product_text) = split_reaction(s)?;

    let reactant_templates = parse_section(reactant_text, "reactant")?;
    if reactant_templates.is_empty() {
        return Err(ReactionSmartsError::EmptyReactants);
    }

    let product_templates = parse_section(product_text, "product")?;
    if product_templates.is_empty() {
        return Err(ReactionSmartsError::EmptyProducts);
    }

    let agent_templates = if agent_text.is_empty() {
        Vec::new()
    } else {
        parse_section(agent_text, "agent")?
    };

    Ok(Reaction {
        reactant_templates,
        product_templates,
        agent_templates,
    })
}

fn split_reaction(s: &str) -> Result<(&str, &str, &str), ReactionSmartsError> {
    let gt_positions = find_gt_positions(s);

    if gt_positions.is_empty() {
        return Err(ReactionSmartsError::MissingSeparator);
    }

    // Look for >> (two consecutive >)
    for i in 0..gt_positions.len() - 1 {
        if gt_positions[i] + 1 == gt_positions[i + 1] {
            let reactants = &s[..gt_positions[i]];
            let products = &s[gt_positions[i + 1] + 1..];

            // Check if there are other > besides this >> pair
            let other_gt: Vec<_> = gt_positions
                .iter()
                .filter(|&&p| p != gt_positions[i] && p != gt_positions[i + 1])
                .collect();
            if !other_gt.is_empty() {
                return Err(ReactionSmartsError::TooManySeparators);
            }

            return Ok((reactants, "", products));
        }
    }

    // Look for exactly 2 > (the a>b>c form)
    if gt_positions.len() == 2 {
        let reactants = &s[..gt_positions[0]];
        let agents = &s[gt_positions[0] + 1..gt_positions[1]];
        let products = &s[gt_positions[1] + 1..];
        return Ok((reactants, agents, products));
    }

    if gt_positions.len() == 1 {
        return Err(ReactionSmartsError::MissingSeparator);
    }

    Err(ReactionSmartsError::TooManySeparators)
}

fn find_gt_positions(s: &str) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut bracket_depth = 0u32;
    let mut paren_depth = 0u32;

    for (i, ch) in s.chars().enumerate() {
        match ch {
            '[' => bracket_depth += 1,
            ']' => bracket_depth = bracket_depth.saturating_sub(1),
            '(' => paren_depth += 1,
            ')' => paren_depth = paren_depth.saturating_sub(1),
            '>' if bracket_depth == 0 && paren_depth == 0 => {
                positions.push(i);
            }
            _ => {}
        }
    }

    positions
}

fn split_on_dot(s: &str) -> Vec<&str> {
    if s.is_empty() {
        return Vec::new();
    }

    let mut parts = Vec::new();
    let mut start = 0;
    let mut bracket_depth = 0u32;
    let mut paren_depth = 0u32;

    for (i, ch) in s.char_indices() {
        match ch {
            '[' => bracket_depth += 1,
            ']' => bracket_depth = bracket_depth.saturating_sub(1),
            '(' => paren_depth += 1,
            ')' => paren_depth = paren_depth.saturating_sub(1),
            '.' if bracket_depth == 0 && paren_depth == 0 => {
                parts.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    parts.push(&s[start..]);
    parts.into_iter().filter(|p| !p.is_empty()).collect()
}

fn parse_section(
    text: &str,
    section: &'static str,
) -> Result<Vec<Mol<AtomExpr, BondExpr>>, ReactionSmartsError> {
    let components = split_on_dot(text);
    let mut mols = Vec::with_capacity(components.len());
    for comp in components {
        let stripped = strip_component_group(comp);
        for sub in split_on_dot(stripped) {
            let mol = from_smarts(sub)
                .map_err(|e| ReactionSmartsError::InvalidComponent { section, detail: e })?;
            mols.push(mol);
        }
    }
    Ok(mols)
}

fn strip_component_group(s: &str) -> &str {
    let bytes = s.as_bytes();
    if bytes.first() == Some(&b'(') && bytes.last() == Some(&b')') {
        let mut depth = 0i32;
        for (i, &b) in bytes.iter().enumerate() {
            match b {
                b'(' => depth += 1,
                b')' => {
                    depth -= 1;
                    if depth == 0 && i < bytes.len() - 1 {
                        return s;
                    }
                }
                _ => {}
            }
        }
        if depth == 0 {
            return &s[1..s.len() - 1];
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_simple_reaction() {
        let (r, a, p) = split_reaction("[C:1][Br:2]>>[C:1][OH]").unwrap();
        assert_eq!(r, "[C:1][Br:2]");
        assert_eq!(a, "");
        assert_eq!(p, "[C:1][OH]");
    }

    #[test]
    fn split_with_agents() {
        let (r, a, p) = split_reaction("[C:1]=[C:2]>[Pd]>[C:1][C:2]").unwrap();
        assert_eq!(r, "[C:1]=[C:2]");
        assert_eq!(a, "[Pd]");
        assert_eq!(p, "[C:1][C:2]");
    }

    #[test]
    fn no_separator_is_error() {
        assert!(split_reaction("[C][Br]").is_err());
    }

    #[test]
    fn single_gt_is_error() {
        assert!(split_reaction("[C]>[Br]").is_err());
    }

    #[test]
    fn split_on_dot_simple() {
        let parts = split_on_dot("[C:1]Br.[N:2]");
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0], "[C:1]Br");
        assert_eq!(parts[1], "[N:2]");
    }

    #[test]
    fn split_on_dot_with_brackets() {
        let parts = split_on_dot("[C.C]");
        assert_eq!(parts.len(), 1);
    }

    #[test]
    fn strip_outer_parens_simple() {
        assert_eq!(strip_component_group("(A.B)"), "A.B");
    }

    #[test]
    fn strip_no_parens() {
        assert_eq!(strip_component_group("A.B"), "A.B");
    }

    #[test]
    fn strip_inner_parens_untouched() {
        assert_eq!(strip_component_group("(A(=O).B)"), "A(=O).B");
    }

    #[test]
    fn strip_non_matching_outer_parens() {
        assert_eq!(strip_component_group("(A)(B)"), "(A)(B)");
    }

    #[test]
    fn parse_section_strips_component_group() {
        let result = parse_section("([C:1]=O.[N:2])", "reactant").unwrap();
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn parse_parenthesized_reaction() {
        let rxn = parse_reaction_smarts("([C:1]=O.[N:2])>>[C:1][N:2]").unwrap();
        assert_eq!(rxn.reactant_templates.len(), 2);
        assert_eq!(rxn.product_templates.len(), 1);
    }

    #[test]
    fn parse_complex_parenthesized_reaction() {
        let rxn =
            parse_reaction_smarts("([C:1](=O)[OH].[NH2:2][C:3])>>[C:1](=O)[N:2][C:3]").unwrap();
        assert_eq!(rxn.reactant_templates.len(), 2);
    }
}
