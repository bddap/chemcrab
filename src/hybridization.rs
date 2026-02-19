use petgraph::graph::NodeIndex;

use crate::element::outer_shell_electrons;
use crate::mol::Mol;
use crate::traits::{HasAromaticity, HasAtomicNum, HasBondOrder, HasFormalCharge, HasHydrogenCount};
use crate::valence::total_valence;
use crate::wrappers::Hybridization;

fn num_bonds_plus_lone_pairs<A, B>(mol: &Mol<A, B>, idx: NodeIndex) -> i16
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    let atomic_num = atom.atomic_num();

    if atomic_num <= 1 {
        let degree = mol.neighbors(idx).count() as i16 + atom.hydrogen_count() as i16;
        return degree;
    }

    let degree = mol.neighbors(idx).count() as i16 + atom.hydrogen_count() as i16;
    let nouter = outer_shell_electrons(atomic_num) as i16;
    let total_valence = total_valence(mol, idx) as i16;
    let charge = atom.formal_charge() as i16;

    let free_electrons = nouter - (total_valence + charge);

    if total_valence + nouter - charge < 8 {
        let num_radicals = crate::radical::num_radical_electrons(mol, idx) as i16;
        let lone_pairs = (free_electrons - num_radicals) / 2;
        degree + lone_pairs + num_radicals
    } else {
        let lone_pairs = free_electrons / 2;
        degree + lone_pairs
    }
}

pub fn assign_hybridization_atom<A, B>(
    mol: &Mol<A, B>,
    idx: NodeIndex,
    has_conjugated_bond: bool,
) -> Hybridization
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge,
    B: HasBondOrder,
{
    let atom = mol.atom(idx);
    if atom.atomic_num() == 0 {
        return Hybridization::Other;
    }

    let norbs = if atom.atomic_num() < 89 {
        num_bonds_plus_lone_pairs(mol, idx)
    } else {
        mol.neighbors(idx).count() as i16 + atom.hydrogen_count() as i16
    };

    match norbs {
        i16::MIN..=0 => Hybridization::S,
        1 => Hybridization::S,
        2 => Hybridization::SP,
        3 => Hybridization::SP2,
        4 => {
            let total_degree = mol.neighbors(idx).count() as u8 + atom.hydrogen_count();
            if total_degree > 3 || !has_conjugated_bond {
                Hybridization::SP3
            } else {
                Hybridization::SP2
            }
        }
        5 => Hybridization::SP3D,
        6 => Hybridization::SP3D2,
        _ => Hybridization::Other,
    }
}

pub fn assign_hybridization<A, B>(mol: &Mol<A, B>) -> Vec<Hybridization>
where
    A: HasAtomicNum + HasHydrogenCount + HasFormalCharge + HasAromaticity,
    B: HasBondOrder,
{
    let conjugation = crate::conjugation::assign_conjugation(mol);

    mol.atoms()
        .map(|idx| {
            let has_conj = mol
                .bonds_of(idx)
                .any(|ei| ei.index() < conjugation.len() && conjugation[ei.index()]);
            assign_hybridization_atom(mol, idx, has_conj)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;
    use crate::wrappers::Hybridization::*;

    fn hyb(smiles: &str) -> Vec<Hybridization> {
        let mol = from_smiles(smiles).unwrap();
        assign_hybridization(&mol)
    }

    #[test]
    fn methane_sp3() {
        assert_eq!(hyb("C"), vec![SP3]);
    }

    #[test]
    fn ethane_sp3() {
        assert_eq!(hyb("CC"), vec![SP3, SP3]);
    }

    #[test]
    fn ethene_sp2() {
        assert_eq!(hyb("C=C"), vec![SP2, SP2]);
    }

    #[test]
    fn acetylene_sp() {
        assert_eq!(hyb("C#C"), vec![SP, SP]);
    }

    #[test]
    fn benzene_sp2() {
        let h = hyb("c1ccccc1");
        assert_eq!(h.len(), 6);
        assert!(h.iter().all(|&x| x == SP2));
    }

    #[test]
    fn water_sp3() {
        assert_eq!(hyb("O"), vec![SP3]);
    }

    #[test]
    fn ammonia_sp3() {
        assert_eq!(hyb("N"), vec![SP3]);
    }

    #[test]
    fn ethanol_sp3() {
        assert_eq!(hyb("CCO"), vec![SP3, SP3, SP3]);
    }

    #[test]
    fn acetaldehyde() {
        assert_eq!(hyb("CC=O"), vec![SP3, SP2, SP2]);
    }

    #[test]
    fn acetic_acid() {
        let h = hyb("CC(=O)O");
        assert_eq!(h[0], SP3);
        assert_eq!(h[1], SP2);
        assert_eq!(h[2], SP2);
        assert_eq!(h[3], SP2);
    }

    #[test]
    fn ammonium_sp3() {
        assert_eq!(hyb("[NH4+]"), vec![SP3]);
    }

    #[test]
    fn oxide_anion_sp3() {
        assert_eq!(hyb("[O-]"), vec![SP3]);
    }

    #[test]
    fn pyridine_sp2() {
        let h = hyb("c1ccncc1");
        assert_eq!(h.len(), 6);
        assert!(h.iter().all(|&x| x == SP2));
    }

    #[test]
    fn pyrrole_sp2() {
        let h = hyb("c1cc[nH]c1");
        assert_eq!(h.len(), 5);
        assert!(h.iter().all(|&x| x == SP2));
    }

    #[test]
    fn furan_sp2() {
        let h = hyb("c1ccoc1");
        assert_eq!(h.len(), 5);
        assert!(h.iter().all(|&x| x == SP2));
    }

    #[test]
    fn thiophene_sp2() {
        let h = hyb("c1ccsc1");
        assert_eq!(h.len(), 5);
        assert!(h.iter().all(|&x| x == SP2));
    }

    #[test]
    fn phenol_oxygen_sp2() {
        let h = hyb("Oc1ccccc1");
        assert_eq!(h[0], SP2);
        for hyb in &h[1..7] {
            assert_eq!(*hyb, SP2);
        }
    }

    #[test]
    fn aniline_nitrogen_sp2() {
        let h = hyb("Nc1ccccc1");
        assert_eq!(h[0], SP2);
    }

    #[test]
    fn thiophenol_sulfur_sp3() {
        let h = hyb("Sc1ccccc1");
        assert_eq!(h[0], SP3);
    }

    #[test]
    fn acetamide_nitrogen_sp2() {
        let h = hyb("CC(N)=O");
        assert_eq!(h[0], SP3);
        assert_eq!(h[1], SP2);
        assert_eq!(h[2], SP2);
        assert_eq!(h[3], SP2);
    }

    #[test]
    fn sodium_chloride() {
        let h = hyb("[Cl-].[Na+]");
        assert_eq!(h[0], SP3); // Cl-
        assert_eq!(h[1], S);   // Na+
    }

    #[test]
    fn cyanobenzene() {
        let h = hyb("N#Cc1ccccc1");
        assert_eq!(h[0], SP);
        assert_eq!(h[1], SP);
    }

    #[test]
    fn borane_sp2() {
        assert_eq!(hyb("B"), vec![SP2]);
    }

    #[test]
    fn phosphine_sp3() {
        assert_eq!(hyb("P"), vec![SP3]);
    }

    #[test]
    fn hydrogen_sulfide_sp3() {
        assert_eq!(hyb("S"), vec![SP3]);
    }

    #[test]
    fn dmso2_sulfur_sp3() {
        let h = hyb("CS(C)(=O)=O");
        assert_eq!(h[1], SP3);
    }

    #[test]
    fn phosphoric_acid_sp3() {
        let h = hyb("O=P(O)(O)O");
        assert_eq!(h[1], SP3);
    }

    #[test]
    fn nitrobenzene() {
        let h = hyb("O=[N+]([O-])c1ccccc1");
        assert_eq!(h[0], SP2);
        assert_eq!(h[1], SP2);
        assert_eq!(h[2], SP2);
    }

    #[test]
    fn deuterated_methane() {
        let h = hyb("[2H]C([2H])([2H])[2H]");
        assert_eq!(h[0], S);
        assert_eq!(h[1], SP3);
        assert_eq!(h[2], S);
        assert_eq!(h[3], S);
        assert_eq!(h[4], S);
    }

    #[test]
    fn methyl_radical() {
        assert_eq!(hyb("[CH3]"), vec![SP3]);
    }

    #[test]
    fn amino_radical() {
        assert_eq!(hyb("[NH2]"), vec![SP3]);
    }

    #[test]
    fn hydroxyl_radical() {
        assert_eq!(hyb("[OH]"), vec![SP3]);
    }
}
