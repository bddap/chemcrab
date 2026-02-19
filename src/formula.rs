use std::collections::BTreeMap;
use std::fmt::Write;

use crate::element::Element;
use crate::mol::Mol;
use crate::traits::{HasAtomicNum, HasFormalCharge, HasHydrogenCount};

pub fn average_mol_weight<A: HasAtomicNum + HasHydrogenCount, B>(mol: &Mol<A, B>) -> f64 {
    let h_weight = Element::H.atomic_weight();
    mol.atoms().fold(0.0, |acc, idx| {
        let a = mol.atom(idx);
        let elem_weight = Element::from_atomic_num(a.atomic_num())
            .map_or(0.0, |e| e.atomic_weight());
        acc + elem_weight + a.hydrogen_count() as f64 * h_weight
    })
}

pub fn exact_mol_weight<A: HasAtomicNum + HasHydrogenCount, B>(mol: &Mol<A, B>) -> f64 {
    let h_mass = Element::H.exact_mass();
    mol.atoms().fold(0.0, |acc, idx| {
        let a = mol.atom(idx);
        let elem_mass = Element::from_atomic_num(a.atomic_num())
            .map_or(0.0, |e| e.exact_mass());
        acc + elem_mass + a.hydrogen_count() as f64 * h_mass
    })
}

pub fn mol_formula<A: HasAtomicNum + HasHydrogenCount + HasFormalCharge, B>(
    mol: &Mol<A, B>,
) -> String {
    let mut counts: BTreeMap<&'static str, u32> = BTreeMap::new();
    let mut net_charge: i32 = 0;

    for idx in mol.atoms() {
        let a = mol.atom(idx);
        if let Some(elem) = Element::from_atomic_num(a.atomic_num()) {
            *counts.entry(elem.symbol()).or_default() += 1;
        }
        let hc = a.hydrogen_count() as u32;
        if hc > 0 {
            *counts.entry("H").or_default() += hc;
        }
        net_charge += a.formal_charge() as i32;
    }

    let mut result = String::new();

    let has_carbon = counts.contains_key("C");
    if has_carbon {
        append_element(&mut result, "C", counts.remove("C").unwrap());
        if let Some(h) = counts.remove("H") {
            append_element(&mut result, "H", h);
        }
    }

    for (sym, count) in &counts {
        append_element(&mut result, sym, *count);
    }

    match net_charge.cmp(&0) {
        std::cmp::Ordering::Greater => {
            if net_charge > 1 {
                write!(result, "{net_charge}+").unwrap();
            } else {
                result.push('+');
            }
        }
        std::cmp::Ordering::Less => {
            if net_charge < -1 {
                write!(result, "{}−", net_charge.unsigned_abs()).unwrap();
            } else {
                result.push('-');
            }
        }
        std::cmp::Ordering::Equal => {}
    }

    result
}

fn append_element(buf: &mut String, symbol: &str, count: u32) {
    buf.push_str(symbol);
    if count > 1 {
        write!(buf, "{count}").unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    fn assert_approx(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "expected {expected} ± {tol}, got {actual}"
        );
    }

    #[test]
    fn methane_formula() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol_formula(&mol), "CH4");
    }

    #[test]
    fn methane_amw() {
        let mol = parse_smiles("C").unwrap();
        assert_approx(average_mol_weight(&mol), 16.043, 0.01);
    }

    #[test]
    fn methane_exact() {
        let mol = parse_smiles("C").unwrap();
        assert_approx(exact_mol_weight(&mol), 16.031, 0.01);
    }

    #[test]
    fn benzene_formula() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol_formula(&mol), "C6H6");
    }

    #[test]
    fn benzene_amw() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_approx(average_mol_weight(&mol), 78.112, 0.01);
    }

    #[test]
    fn benzene_exact() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_approx(exact_mol_weight(&mol), 78.047, 0.01);
    }

    #[test]
    fn water_formula() {
        let mol = parse_smiles("O").unwrap();
        assert_eq!(mol_formula(&mol), "H2O");
    }

    #[test]
    fn water_amw() {
        let mol = parse_smiles("O").unwrap();
        assert_approx(average_mol_weight(&mol), 18.015, 0.01);
    }

    #[test]
    fn ethanol_formula() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol_formula(&mol), "C2H6O");
    }

    #[test]
    fn nacl_formula() {
        let mol = parse_smiles("[Na+].[Cl-]").unwrap();
        assert_eq!(mol_formula(&mol), "ClNa");
    }

    #[test]
    fn ammonium_formula() {
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol_formula(&mol), "H4N+");
    }

    #[test]
    fn iron_formula() {
        let mol = parse_smiles("[Fe]").unwrap();
        assert_eq!(mol_formula(&mol), "Fe");
    }

    #[test]
    fn empty_mol_formula() {
        let mol: Mol<crate::atom::Atom, crate::bond::SmilesBond> = Mol::new();
        assert_eq!(mol_formula(&mol), "");
    }

    #[test]
    fn empty_mol_weight() {
        let mol: Mol<crate::atom::Atom, crate::bond::SmilesBond> = Mol::new();
        assert_eq!(average_mol_weight(&mol), 0.0);
        assert_eq!(exact_mol_weight(&mol), 0.0);
    }

    #[test]
    fn h2_formula() {
        let mol = parse_smiles("[HH]").unwrap();
        assert_eq!(mol_formula(&mol), "H2");
    }

    #[test]
    fn iron_amw() {
        let mol = parse_smiles("[Fe]").unwrap();
        assert_approx(average_mol_weight(&mol), 55.845, 0.01);
    }
}
