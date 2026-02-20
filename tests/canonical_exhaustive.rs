use chemcrab::graph::renumber_atoms;
use chemcrab::smiles::{from_smiles, to_canonical_smiles};

const MOLECULES: &[&str] = &[
    // Simple
    "C",
    "CC",
    "C=C",
    "C#C",
    "C=O",
    "O",
    "N",
    "[H][H]",
    // Heteroatoms
    "CCO",
    "CCN",
    "CCS",
    "CCF",
    "CCCl",
    "CCBr",
    // Branching
    "CC(C)C",
    "CC(C)(C)C",
    "CC(=O)O",
    "CC(=O)N",
    // Rings
    "C1CC1",
    "C1CCC1",
    "C1CCCCC1",
    "c1ccccc1",
    "c1ccncc1",
    "c1ccoc1",
    "c1ccsc1",
    "c1ccc2ccccc2c1",
    "C1CC2CCCC(C1)C2",
    // Charged
    "[NH4+]",
    "[O-]",
    "[Na+].[Cl-]",
    // Chirality
    "[C@@H](F)(Cl)Br",
    "[C@](F)(Cl)(Br)I",
    "N[C@@H](C)C(=O)O",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    // E/Z
    "F/C=C/F",
    r"F/C=C\F",
    r"Cl/C=C\Br",
    r"CC/C=C/CC",
    r"F/C=C/C=C/F",
    // Mixed stereo
    r"F/C=C/[C@@H](Cl)Br",
    // Highly symmetric chiral
    "[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    // Functional groups
    "CC(=O)OC",
    "c1ccc(cc1)O",
    "c1ccc(cc1)N",
    "CC(=O)Nc1ccccc1",
    "c1ccc2c(c1)[nH]cc2",
    // Multi-component
    "C.C",
    "[Na+].[Cl-].O",
    // Isotopes
    "[2H]C([2H])([2H])[2H]",
    "[13C]c1ccccc1",
    // Larger / drug-like
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "O=C(O)c1ccccc1O",
    "Nc1ccc(cc1)S(=O)(=O)Nc1ccccn1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "OC(=O)CC(O)(CC(=O)O)C(=O)O",
    // Additional edge cases
    "C#N",
    "[Cu+2]",
    "c1cc[nH]c1",
    "C1=CC=CC=C1",
    "c1cnc2ccccc2n1",
    "[O-][N+](=O)c1ccccc1",
    "c1cc2ccc3cccc4ccc(c1)c2c34",
    // Bug A regressions: symmetric chirality tie-breaking
    "[C@@H]1(O)C[C@H](O)C[C@@H](O)C1",
    "[C@H]1(O)C[C@@H](O)C[C@H](O)C1",
    "[C@@H]1(O)C[C@@H](O)C[C@@H](O)C1",
    "[C@H]1(O)C[C@H](O)C[C@H](O)C1",
    "[C@@H]1(F)C[C@H](F)C[C@@H](F)C1",
    "[C@@H]1(Cl)C[C@H](Cl)C[C@@H](Cl)C1",
    "[C@@H]1(CC1)[C@H]2CC2",
    "[C@H]1(CC1)[C@H]2CC2",
    "[C@@H]1(CCC1)[C@H]2CCC2",
    "[C@@H]1(CC1)[C@H]2CCC2",
    "[C@@H]1(CCC1)[C@H]2CC2",
    "[C@@H]1(CC1)[C@@H](F)CC",
    "[C@@H]1(CCCC1)[C@H]2CCCC2",
    // Bug B regressions: E/Z conjugated diene direction signs
    r"F/C=C/C=C\F",
    r"F/C=C\C=C\F",
    r"Cl/C=C\C=C\Cl",
    r"Br/C=C\C=C\Br",
    r"F/C=C\C=C/C=C\F",
    r"F/C=C/C=C\C=C/F",
    r"C1/C=C\CCC1",
    r"C/1=C\CCC/C1",
    r"C/1=C/C=C\CC1",
    r"C/1=C\C=C/CC1",
];

fn canonical(smiles: &str) -> String {
    let mol = from_smiles(smiles).unwrap_or_else(|e| panic!("parse failed for '{smiles}': {e}"));
    to_canonical_smiles(&mol)
}

fn all_permutations(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut state: Vec<usize> = (0..n).collect();
    result.push(state.clone());
    if n <= 1 {
        return result;
    }
    let mut c = vec![0usize; n];
    let mut i = 1;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                state.swap(0, i);
            } else {
                state.swap(c[i], i);
            }
            result.push(state.clone());
            c[i] += 1;
            i = 1;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    result
}

struct Xorshift64(u64);

impl Xorshift64 {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn shuffle(&mut self, slice: &mut [usize]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            slice.swap(i, j);
        }
    }
}

fn random_permutations(n: usize, count: usize) -> Vec<Vec<usize>> {
    let mut rng = Xorshift64(0xDEAD_BEEF_CAFE_BABE);
    (0..count)
        .map(|_| {
            let mut perm: Vec<usize> = (0..n).collect();
            rng.shuffle(&mut perm);
            perm
        })
        .collect()
}

const EXHAUSTIVE_THRESHOLD: usize = 8;
const RANDOM_SAMPLE_COUNT: usize = 100;

fn permutations_for(n: usize) -> Vec<Vec<usize>> {
    if n <= EXHAUSTIVE_THRESHOLD {
        all_permutations(n)
    } else {
        random_permutations(n, RANDOM_SAMPLE_COUNT)
    }
}

#[test]
fn determinism() {
    for &smiles in MOLECULES {
        let mol =
            from_smiles(smiles).unwrap_or_else(|e| panic!("parse failed for '{smiles}': {e}"));
        let a = to_canonical_smiles(&mol);
        let b = to_canonical_smiles(&mol);
        assert_eq!(a, b, "determinism failed for '{smiles}': '{a}' vs '{b}'");
    }
}

#[test]
fn round_trip_idempotence() {
    for &smiles in MOLECULES {
        let first = canonical(smiles);
        let reparsed = from_smiles(&first)
            .unwrap_or_else(|e| panic!("reparse failed for '{first}' (from '{smiles}'): {e}"));
        let second = to_canonical_smiles(&reparsed);
        assert_eq!(
            first, second,
            "round-trip failed for '{smiles}': first='{first}', second='{second}'"
        );
    }
}

#[test]
fn permutation_invariance() {
    for &smiles in MOLECULES {
        let mol =
            from_smiles(smiles).unwrap_or_else(|e| panic!("parse failed for '{smiles}': {e}"));
        let expected = to_canonical_smiles(&mol);
        let n = mol.atom_count();
        for perm in permutations_for(n) {
            let renum = renumber_atoms(&mol, &perm).unwrap_or_else(|e| {
                panic!("renumber failed for '{smiles}' with perm {perm:?}: {e}")
            });
            let got = to_canonical_smiles(&renum);
            assert_eq!(
                expected, got,
                "permutation invariance failed for '{smiles}' with perm {perm:?}: \
                 expected='{expected}', got='{got}'"
            );
        }
    }
}
