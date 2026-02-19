use petgraph::graph::NodeIndex;

use crate::bond::BondOrder;
use crate::element::Element;
use crate::mol::Mol;
use crate::traits::{HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceError {
    pub atom_idx: NodeIndex,
    pub atomic_num: u8,
    pub actual_valence: u8,
    pub allowed_valences: Vec<u8>,
}

impl std::fmt::Display for ValenceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let sym = Element::from_atomic_num(self.atomic_num)
            .map(|e| e.symbol())
            .unwrap_or("?");
        write!(
            f,
            "atom {} ({}): valence {} not in {:?}",
            self.atom_idx.index(),
            sym,
            self.actual_valence,
            self.allowed_valences,
        )
    }
}

impl std::error::Error for ValenceError {}

pub fn total_valence<A, B>(mol: &Mol<A, B>, atom: NodeIndex) -> u8
where
    A: HasHydrogenCount,
    B: HasBondOrder,
{
    let bond_sum: u8 = mol
        .bonds_of(atom)
        .map(|ei| match mol.bond(ei).bond_order() {
            BondOrder::Single => 1u8,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        })
        .sum();
    bond_sum + mol.atom(atom).hydrogen_count()
}

pub fn check_valence<A, B>(mol: &Mol<A, B>) -> Result<(), Vec<ValenceError>>
where
    A: HasAtomicNum + HasFormalCharge + HasHydrogenCount,
    B: HasBondOrder,
{
    let errors: Vec<ValenceError> = mol
        .atoms()
        .filter_map(|idx| {
            let atom = mol.atom(idx);
            if atom.formal_charge() != 0 {
                return None;
            }
            let elem = Element::from_atomic_num(atom.atomic_num())?;
            let allowed = elem.default_valences();
            if allowed.is_empty() {
                return None;
            }
            let v = total_valence(mol, idx);
            if allowed.contains(&v) {
                return None;
            }
            Some(ValenceError {
                atom_idx: idx,
                atomic_num: atom.atomic_num(),
                actual_valence: v,
                allowed_valences: allowed.to_vec(),
            })
        })
        .collect();

    if errors.is_empty() {
        Ok(())
    } else {
        Err(errors)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    #[test]
    fn methane_valid() {
        let mol = from_smiles("C").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn ethane_valid() {
        let mol = from_smiles("CC").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn benzene_valid() {
        let mol = from_smiles("c1ccccc1").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn water_valid() {
        let mol = from_smiles("O").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn ammonia_valid() {
        let mol = from_smiles("N").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn pentavalent_carbon_invalid() {
        use crate::atom::Atom;
        use crate::bond::Bond;
        use crate::mol::Mol;

        let mut mol = Mol::<Atom, Bond>::new();
        let c = mol.add_atom(Atom {
            atomic_num: 6,
            hydrogen_count: 5,
            ..Default::default()
        });
        let errs = check_valence(&mol).unwrap_err();
        assert_eq!(errs.len(), 1);
        assert_eq!(errs[0].atom_idx, c);
        assert_eq!(errs[0].actual_valence, 5);
        assert_eq!(errs[0].allowed_valences, vec![4]);
    }

    #[test]
    fn charged_ammonium_skipped() {
        let mol = from_smiles("[NH4+]").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn metal_skipped() {
        let mol = from_smiles("[Fe]").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn sulfur_hexafluoride_valid() {
        let mol = from_smiles("S(F)(F)(F)(F)(F)F").unwrap();
        assert!(check_valence(&mol).is_ok());
    }

    #[test]
    fn ethene_total_valence() {
        let mol = from_smiles("C=C").unwrap();
        for idx in mol.atoms() {
            assert_eq!(total_valence(&mol, idx), 4);
        }
    }
}
