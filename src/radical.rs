//! Radical electron counting.
//!
//! A radical is an atom with unpaired electrons — for example, a methyl
//! radical `[CH3]` has one unpaired electron on carbon. The count is
//! inferred from the difference between the atom's valence shell capacity
//! and its actual bonding plus lone-pair electrons.

use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::element::{outer_shell_electrons, Element};
use crate::mol::Mol;
use crate::traits::{HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};

/// Count unpaired (radical) electrons on a specific atom.
///
/// Returns 0 for atoms in a normal bonding state. Returns 1 for
/// monoradicals (e.g., `[CH3]`), 2 for diradicals (e.g., `[CH2]`), etc.
pub fn num_radical_electrons<A, B>(mol: &Mol<A, B>, idx: NodeIndex) -> u8
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    let atomic_num = atom.atomic_num();
    let hcount = atom.hydrogen_count();
    let charge = atom.formal_charge() as i16;

    let elem = match Element::from_atomic_num(atomic_num) {
        Some(e) => e,
        None => return 0,
    };

    let default_valences = elem.default_valences();
    if default_valences.is_empty() {
        if mol.neighbors(idx).count() > 0 {
            return 0;
        }
        let n_outer = outer_shell_electrons(atomic_num) as i16;
        let n_valence = n_outer - charge;
        if n_valence < 0 {
            return 0;
        }
        return (n_valence % 2) as u8;
    }

    let bond_sum: i16 = mol
        .bonds_of(idx)
        .map(|ei| match mol.bond(ei).bond_order() {
            BondOrder::Single => 1i16,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        })
        .sum();

    let total_valence = bond_sum + hcount as i16;
    let n_outer = outer_shell_electrons(atomic_num) as i16;

    let base_count: i16 = if atomic_num <= 2 { 2 } else { 8 };

    let mut num_radicals = base_count - n_outer - total_valence + charge;
    if num_radicals < 0 {
        num_radicals = 0;
        if default_valences.len() > 1 {
            for &val in default_valences {
                let r = val as i16 - total_valence + charge;
                if r >= 0 {
                    num_radicals = r;
                    break;
                }
            }
        }
    }

    let num_radicals2 = n_outer - total_valence - charge;
    if num_radicals2 >= 0 && num_radicals2 < num_radicals {
        num_radicals = num_radicals2;
    }

    if num_radicals > 0 {
        num_radicals as u8
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;
    use petgraph::graph::NodeIndex;

    fn n(i: usize) -> NodeIndex {
        NodeIndex::new(i)
    }

    #[test]
    fn methane_no_radicals() {
        let mol = from_smiles("C").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }

    #[test]
    fn methyl_radical() {
        let mol = from_smiles("[CH3]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 1);
    }

    #[test]
    fn methylene_diradical() {
        let mol = from_smiles("[CH2]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 2);
    }

    #[test]
    fn ethane_no_radicals() {
        let mol = from_smiles("CC").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
        assert_eq!(num_radical_electrons(&mol, n(1)), 0);
    }

    #[test]
    fn amino_radical() {
        let mol = from_smiles("[NH2]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 1);
    }

    #[test]
    fn hydroxyl_radical() {
        let mol = from_smiles("[OH]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 1);
    }

    #[test]
    fn water_no_radicals() {
        let mol = from_smiles("O").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }

    #[test]
    fn ammonia_no_radicals() {
        let mol = from_smiles("N").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }

    #[test]
    fn charged_no_radicals() {
        let mol = from_smiles("[NH4+]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }

    #[test]
    fn benzene_no_radicals() {
        let mol = from_smiles("c1ccccc1").unwrap();
        for idx in mol.atoms() {
            assert_eq!(num_radical_electrons(&mol, idx), 0);
        }
    }

    #[test]
    fn metal_no_radicals() {
        let mol = from_smiles("[Fe]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }

    #[test]
    fn oxide_anion() {
        let mol = from_smiles("[O-]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 1);
    }

    #[test]
    fn chloride_anion() {
        let mol = from_smiles("[Cl-]").unwrap();
        assert_eq!(num_radical_electrons(&mol, n(0)), 0);
    }
}
