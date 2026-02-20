//! Molecular formula and molecular weight calculations.
//!
//! [`mol_formula`] produces a Hill system string, [`average_mol_weight`]
//! gives the average molecular weight in daltons, and [`exact_mol_weight`]
//! gives the monoisotopic exact mass.

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::element::isotope_exact_mass;
use crate::element::Element;
use crate::mol::Mol;
use crate::traits::{HasAtomicNum, HasFormalCharge, HasHydrogenCount, HasIsotope};

/// Compute the average molecular weight in daltons (Da).
///
/// Uses standard atomic weights averaged over natural isotopic abundance.
/// Atoms with an explicit isotope label use that isotope's exact mass
/// instead.
pub fn average_mol_weight<A: HasAtomicNum + HasHydrogenCount + HasIsotope, B>(
    mol: &Mol<A, B>,
) -> f64 {
    let h_weight = Element::H.atomic_weight();
    mol.atoms().fold(0.0, |acc, idx| {
        let a = mol.atom(idx);
        let elem = Element::from_atomic_num(a.atomic_num());
        let iso = a.isotope();
        let mass = if iso > 0 {
            isotope_exact_mass(a.atomic_num(), iso)
                .or_else(|| elem.map(|e| e.atomic_weight()))
                .unwrap_or(0.0)
        } else {
            elem.map_or(0.0, |e| e.atomic_weight())
        };
        acc + mass + a.hydrogen_count() as f64 * h_weight
    })
}

/// Compute the monoisotopic exact mass.
///
/// Uses the mass of the most abundant isotope of each element. Atoms
/// with an explicit isotope label use that isotope's exact mass.
pub fn exact_mol_weight<A: HasAtomicNum + HasHydrogenCount + HasIsotope, B>(
    mol: &Mol<A, B>,
) -> f64 {
    let h_mass = Element::H.exact_mass();
    mol.atoms().fold(0.0, |acc, idx| {
        let a = mol.atom(idx);
        let elem = Element::from_atomic_num(a.atomic_num());
        let iso = a.isotope();
        let mass = if iso > 0 {
            isotope_exact_mass(a.atomic_num(), iso)
                .or_else(|| elem.map(|e| e.exact_mass()))
                .unwrap_or(0.0)
        } else {
            elem.map_or(0.0, |e| e.exact_mass())
        };
        acc + mass + a.hydrogen_count() as f64 * h_mass
    })
}

/// Compute the molecular formula as a Hill system string.
///
/// The Hill system lists C first, then H, then remaining elements
/// alphabetically. This is the standard convention used in chemical
/// databases. Molecules without carbon list all elements alphabetically.
/// Net charge is appended as `+`, `2+`, `-`, `2-`, etc.
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
                write!(result, "{}-", net_charge.unsigned_abs()).unwrap();
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
    fn oxide_dianion_formula() {
        let mol = parse_smiles("[O-2]").unwrap();
        let f = mol_formula(&mol);
        assert_eq!(f, "O2-");
        assert!(
            !f.contains('\u{2212}'),
            "should use ASCII minus, not Unicode minus"
        );
    }

    #[test]
    fn iron_amw() {
        let mol = parse_smiles("[Fe]").unwrap();
        assert_approx(average_mol_weight(&mol), 55.845, 0.01);
    }

    #[test]
    fn deuterated_methane_exact() {
        let mol = parse_smiles("[2H]C([2H])([2H])[2H]").unwrap();
        let expected = 12.0 + 4.0 * 2.01410177812;
        assert_approx(exact_mol_weight(&mol), expected, 1e-6);
    }

    #[test]
    fn deuterated_methane_amw() {
        let mol = parse_smiles("[2H]C([2H])([2H])[2H]").unwrap();
        let expected = 12.011 + 4.0 * 2.01410177812;
        assert_approx(average_mol_weight(&mol), expected, 1e-6);
    }

    #[test]
    fn c13_toluene_exact() {
        // [13C]c1ccccc1 = C13-CH bonded to benzene ring: 7 carbons, 5 ring H + 0 on [13C]
        let mol = parse_smiles("[13C]c1ccccc1").unwrap();
        let expected = 13.00335483507 + 6.0 * 12.0 + 5.0 * 1.00782503207;
        assert_approx(exact_mol_weight(&mol), expected, 1e-6);
    }

    #[test]
    fn c13_toluene_amw() {
        let mol = parse_smiles("[13C]c1ccccc1").unwrap();
        let expected = 13.00335483507 + 6.0 * 12.011 + 5.0 * 1.008;
        assert_approx(average_mol_weight(&mol), expected, 1e-4);
    }

    #[test]
    fn non_isotopic_unchanged() {
        let mol = parse_smiles("CCO").unwrap();
        let amw = average_mol_weight(&mol);
        let exact = exact_mol_weight(&mol);
        assert_approx(amw, 2.0 * 12.011 + 6.0 * 1.008 + 15.999, 0.01);
        assert_approx(
            exact,
            2.0 * 12.0 + 6.0 * 1.00782503207 + 15.99491461957,
            1e-4,
        );
    }
}
