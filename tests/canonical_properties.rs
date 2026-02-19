use chemcrab::{from_smiles, to_canonical_smiles, renumber_atoms};

fn canonical(smiles: &str) -> String {
    let mol = from_smiles(smiles).unwrap();
    to_canonical_smiles(&mol)
}

// Bug 1: Fragment ordering not canonical
#[test]
fn fragment_ordering_nacl() {
    let a = canonical("[Na+].[Cl-]");
    let b = canonical("[Cl-].[Na+]");
    assert_eq!(a, b, "fragment ordering: '{a}' vs '{b}'");
}

#[test]
fn fragment_ordering_three() {
    let a = canonical("[Na+].[Cl-].O");
    let b = canonical("O.[Na+].[Cl-]");
    assert_eq!(a, b, "fragment ordering: '{a}' vs '{b}'");
}

// Bug 2: Chirality not stable under different input orderings
#[test]
fn chirality_stability_1() {
    let a = canonical("[C@@H](F)(Cl)Br");
    let b = canonical("F[C@H](Cl)Br");
    assert_eq!(a, b, "chirality stability: '{a}' vs '{b}'");
}

#[test]
fn chirality_stability_2() {
    let a = canonical("[C@@H](F)(Cl)Br");
    let b = canonical("Cl[C@@H](F)Br");
    assert_eq!(a, b, "chirality stability: '{a}' vs '{b}'");
}

#[test]
fn chirality_stability_3() {
    let a = canonical("[C@@H](F)(Cl)Br");
    let b = canonical("Br[C@H](F)Cl");
    assert_eq!(a, b, "chirality stability: '{a}' vs '{b}'");
}

#[test]
fn chirality_stability_4() {
    let a = canonical("[C@](F)(Cl)(Br)I");
    let b = canonical("F[C@](Cl)(Br)I");
    assert_eq!(a, b, "chirality stability 4h: '{a}' vs '{b}'");
}

#[test]
fn chirality_stability_alanine() {
    let a = canonical("N[C@@H](C)C(=O)O");
    let b = canonical("[C@H](N)(C)C(=O)O");
    assert_eq!(a, b, "chirality stability alanine: '{a}' vs '{b}'");
}

// Bug 3: E/Z stereo lost on atom renumber
#[test]
fn ez_renumber_trans() {
    let mol = from_smiles("F/C=C/F").unwrap();
    let n = mol.atom_count();
    let reversed: Vec<usize> = (0..n).rev().collect();
    let renum = renumber_atoms(&mol, &reversed).unwrap();
    let s1 = to_canonical_smiles(&mol);
    let s2 = to_canonical_smiles(&renum);
    assert_eq!(s1, s2, "E/Z renumber trans: '{s1}' vs '{s2}'");
}

#[test]
fn ez_renumber_cis() {
    let mol = from_smiles(r"F/C=C\F").unwrap();
    let n = mol.atom_count();
    let reversed: Vec<usize> = (0..n).rev().collect();
    let renum = renumber_atoms(&mol, &reversed).unwrap();
    let s1 = to_canonical_smiles(&mol);
    let s2 = to_canonical_smiles(&renum);
    assert_eq!(s1, s2, "E/Z renumber cis: '{s1}' vs '{s2}'");
}

#[test]
fn ez_renumber_chlorine() {
    let mol = from_smiles("Cl/C=C/Cl").unwrap();
    let n = mol.atom_count();
    let reversed: Vec<usize> = (0..n).rev().collect();
    let renum = renumber_atoms(&mol, &reversed).unwrap();
    let s1 = to_canonical_smiles(&mol);
    let s2 = to_canonical_smiles(&renum);
    assert_eq!(s1, s2, "E/Z renumber chlorine: '{s1}' vs '{s2}'");
}

#[test]
fn ez_renumber_mixed() {
    let mol = from_smiles(r"F/C=C/[C@@H](Cl)Br").unwrap();
    let n = mol.atom_count();
    let reversed: Vec<usize> = (0..n).rev().collect();
    let renum = renumber_atoms(&mol, &reversed).unwrap();
    let s1 = to_canonical_smiles(&mol);
    let s2 = to_canonical_smiles(&renum);
    assert_eq!(s1, s2, "E/Z renumber mixed: '{s1}' vs '{s2}'");
}

#[test]
fn ez_renumber_shifted() {
    let mol = from_smiles("F/C=C/F").unwrap();
    let renum = renumber_atoms(&mol, &[1, 2, 3, 0]).unwrap();
    let s1 = to_canonical_smiles(&mol);
    let s2 = to_canonical_smiles(&renum);
    assert_eq!(s1, s2, "E/Z renumber shifted: '{s1}' vs '{s2}'");
}

// Bug 4: Caffeine aromatic input invariance
#[test]
fn caffeine_kekule_vs_aromatic() {
    let a = canonical("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
    let b = canonical("Cn1c(=O)c2c(ncn2C)n(C)c1=O");
    assert_eq!(a, b, "caffeine kekule vs aromatic: '{a}' vs '{b}'");
}

// Bug 5: Glucose idempotence — chirality flips on reparse
#[test]
fn glucose_idempotence() {
    let first = canonical("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O");
    let second = canonical(&first);
    assert_eq!(first, second, "glucose idempotence: '{first}' vs '{second}'");
}
