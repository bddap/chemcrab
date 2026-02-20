use petgraph::graph::{EdgeIndex, NodeIndex};

use crate::bond::BondOrder;
use crate::element::{outer_shell_electrons, Element};
use crate::mol::Mol;
use crate::traits::{HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};
use crate::valence::total_valence;

fn count_atom_elec<A, B>(mol: &Mol<A, B>, idx: NodeIndex) -> i16
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    let anum = atom.atomic_num();

    let dv = Element::from_atomic_num(anum)
        .map(|e| e.default_valences())
        .unwrap_or(&[]);
    if dv.is_empty() || dv[0] <= 1 {
        return -1;
    }
    let default_val = dv[0];

    let degree = mol.neighbors(idx).count() as u8 + atom.hydrogen_count();
    if degree > 3 {
        return -1;
    }

    let nouter = outer_shell_electrons(anum);
    let nlp = (nouter as i16 - default_val as i16 - atom.formal_charge() as i16).max(0);
    let n_radicals = crate::radical::num_radical_electrons(mol, idx) as i16;

    (default_val as i16 - degree as i16) + nlp - n_radicals
}

fn is_conj_candidate<A, B>(mol: &Mol<A, B>, idx: NodeIndex) -> bool
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    let anum = atom.atomic_num();
    let charge = atom.formal_charge();

    let dv = Element::from_atomic_num(anum)
        .map(|e| e.default_valences())
        .unwrap_or(&[]);
    if dv.is_empty() || dv[0] <= 1 {
        return false;
    }

    let total_valence = total_valence(mol, idx) as i16;
    if charge == 0 && total_valence > dv[0] as i16 {
        return false;
    }

    let nouter = outer_shell_electrons(anum);
    let total_degree = mol.neighbors(idx).count() as u8 + atom.hydrogen_count();

    let row_check = anum <= 10
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && total_degree < 2);

    row_check && count_atom_elec(mol, idx) > 0
}

pub fn assign_conjugation<A, B>(mol: &Mol<A, B>) -> Vec<bool>
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + crate::traits::HasAromaticity,
    B: HasBondOrder,
{
    let mut conjugated = vec![false; mol.bond_count()];

    for edge in mol.bonds() {
        let (a, b) = mol.bond_endpoints(edge).unwrap();
        if mol.atom(a).is_aromatic() && mol.atom(b).is_aromatic() {
            conjugated[edge.index()] = true;
        }
    }

    for atom_idx in mol.atoms() {
        if !is_conj_candidate(mol, atom_idx) {
            continue;
        }

        let atom = mol.atom(atom_idx);
        let sbo = mol.neighbors(atom_idx).count() as u8 + atom.hydrogen_count();
        if !(2..=3).contains(&sbo) {
            continue;
        }

        let bonds: Vec<EdgeIndex> = mol.bonds_of(atom_idx).collect();

        for &bnd1 in &bonds {
            let bo1 = mol.bond(bnd1).bond_order();
            if bo1 != BondOrder::Double && bo1 != BondOrder::Triple {
                continue;
            }

            let (a1, b1) = mol.bond_endpoints(bnd1).unwrap();
            let other1 = if a1 == atom_idx { b1 } else { a1 };
            if !is_conj_candidate(mol, other1) {
                continue;
            }

            for &bnd2 in &bonds {
                if bnd1 == bnd2 {
                    continue;
                }
                let (a2, b2) = mol.bond_endpoints(bnd2).unwrap();
                let at2 = if a2 == atom_idx { b2 } else { a2 };
                let at2_atom = mol.atom(at2);
                let at2_sbo = mol.neighbors(at2).count() as u8 + at2_atom.hydrogen_count();
                if at2_sbo > 3 {
                    continue;
                }
                if is_conj_candidate(mol, at2) {
                    conjugated[bnd1.index()] = true;
                    conjugated[bnd2.index()] = true;
                }
            }
        }
    }

    conjugated
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    fn conj(smiles: &str) -> Vec<bool> {
        let mol = from_smiles(smiles).unwrap();
        assign_conjugation(&mol)
    }

    #[test]
    fn ethane_not_conjugated() {
        let c = conj("CC");
        assert_eq!(c, vec![false]);
    }

    #[test]
    fn ethene_not_conjugated() {
        // A single double bond with no adjacent unsaturation
        // A lone double bond is not conjugated; conjugation requires
        // adjacent unsaturation.
        let c = conj("C=C");
        assert_eq!(c, vec![false]);
    }

    #[test]
    fn butadiene_all_conjugated() {
        let c = conj("C=CC=C");
        assert_eq!(c, vec![true, true, true]);
    }

    #[test]
    fn propene_not_conjugated() {
        let c = conj("CC=C");
        assert_eq!(c, vec![false, false]);
    }

    #[test]
    fn acrolein_conjugated() {
        // C=CC=O
        let c = conj("C=CC=O");
        assert_eq!(c, vec![true, true, true]);
    }

    #[test]
    fn acetic_acid_oh_conjugated() {
        // CC(=O)O: bonds are C-C, C=O, C-OH
        // The C=O and C-OH should be conjugated (through C)
        let c = conj("CC(=O)O");
        assert!(!c[0]); // CH3-C not conjugated
        assert!(c[1]);  // C=O conjugated
        assert!(c[2]);  // C-OH conjugated
    }

    #[test]
    fn acetamide_conjugated() {
        // CC(N)=O: bonds are C-C, C-N, C=O
        let c = conj("CC(N)=O");
        assert!(!c[0]); // CH3-C
        assert!(c[1]);  // C-N conjugated
        assert!(c[2]);  // C=O conjugated
    }

    #[test]
    fn cyclohexane_not_conjugated() {
        let c = conj("C1CCCCC1");
        assert!(c.iter().all(|&x| !x));
    }

    #[test]
    fn aniline_n_c_conjugated() {
        // Nc1ccccc1: N-C(ring) bond should be conjugated
        let mol = from_smiles("Nc1ccccc1").unwrap();
        let c = assign_conjugation(&mol);
        assert!(c[0]); // N-C bond
    }

    #[test]
    fn phenol_o_c_conjugated() {
        let mol = from_smiles("Oc1ccccc1").unwrap();
        let c = assign_conjugation(&mol);
        assert!(c[0]); // O-C bond
    }
}
