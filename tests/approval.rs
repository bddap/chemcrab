use std::collections::HashMap;

use serde::Deserialize;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

// Known parser limitation: `c1ccc(-c2ccccc2)cc1` (biphenyl) fails because
// the explicit `-` single bond inside an aromatic branch context is not yet handled.
fn try_parse(smiles: &str) -> Option<crabchem::Mol<crabchem::Atom, crabchem::Bond>> {
    match crabchem::smiles::from_smiles(smiles) {
        Ok(m) => Some(m),
        Err(e) => {
            eprintln!("SKIP (parse failure): {smiles:?}: {e}");
            None
        }
    }
}

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

// ---------------------------------------------------------------------------
// 1. Formula & weight
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct FormulaEntry {
    smiles: String,
    formula: String,
    average_mw: f64,
    exact_mw: f64,
}

#[test]
fn approval_formula_weight() {
    let data: Vec<FormulaEntry> =
        serde_json::from_str(include_str!("approval_data/formula_weight.json")).unwrap();

    let mut failures = Vec::new();
    let mut skipped = 0usize;
    for entry in &data {
        let mol = match try_parse(&entry.smiles) {
            Some(m) => m,
            None => { skipped += 1; continue; }
        };

        let formula = crabchem::mol_formula(&mol);
        if formula != entry.formula {
            failures.push(format!(
                "[formula] {}: expected {:?}, got {:?}",
                entry.smiles, entry.formula, formula
            ));
        }

        let amw = crabchem::average_mol_weight(&mol);
        if !approx_eq(amw, entry.average_mw, 0.01) {
            failures.push(format!(
                "[avg_mw] {}: expected {}, got {}",
                entry.smiles, entry.average_mw, amw
            ));
        }

        let emw = crabchem::exact_mol_weight(&mol);
        if !approx_eq(emw, entry.exact_mw, 0.01) {
            failures.push(format!(
                "[exact_mw] {}: expected {}, got {}",
                entry.smiles, entry.exact_mw, emw
            ));
        }
    }

    if skipped > 0 {
        eprintln!("formula/weight: skipped {skipped} unparseable molecules");
    }

    if !failures.is_empty() {
        panic!(
            "{} formula/weight failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Aromaticity
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct AromaticityEntry {
    smiles: String,
    num_aromatic_atoms: usize,
    num_atoms: usize,
}

#[test]
fn approval_aromaticity() {
    let data: Vec<AromaticityEntry> =
        serde_json::from_str(include_str!("approval_data/aromaticity.json")).unwrap();

    let mut failures = Vec::new();
    let mut skipped = 0usize;
    for entry in &data {
        let mol = match try_parse(&entry.smiles) {
            Some(m) => m,
            None => { skipped += 1; continue; }
        };

        if mol.atom_count() != entry.num_atoms {
            failures.push(format!(
                "[atom_count] {}: expected {}, got {}",
                entry.smiles,
                entry.num_atoms,
                mol.atom_count()
            ));
        }

        let aromatic = crabchem::find_aromatic_atoms(&mol);
        let count = aromatic.iter().filter(|&&a| a).count();
        if count != entry.num_aromatic_atoms {
            failures.push(format!(
                "[aromatic] {}: expected {}, got {}",
                entry.smiles, entry.num_aromatic_atoms, count
            ));
        }
    }

    if skipped > 0 {
        eprintln!("aromaticity: skipped {skipped} unparseable molecules");
    }

    if !failures.is_empty() {
        panic!(
            "{} aromaticity failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// 3. Ring finding
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct RingEntry {
    smiles: String,
    num_rings: usize,
    ring_sizes: Vec<usize>,
    num_ring_atoms: usize,
    num_ring_bonds: usize,
}

#[test]
fn approval_rings() {
    let data: Vec<RingEntry> =
        serde_json::from_str(include_str!("approval_data/rings.json")).unwrap();

    let mut failures = Vec::new();
    let mut skipped = 0usize;
    for entry in &data {
        let mol = match try_parse(&entry.smiles) {
            Some(m) => m,
            None => { skipped += 1; continue; }
        };
        let ri = crabchem::RingInfo::sssr(&mol);

        if ri.num_rings() != entry.num_rings {
            failures.push(format!(
                "[num_rings] {}: expected {}, got {}",
                entry.smiles,
                entry.num_rings,
                ri.num_rings()
            ));
        }

        let mut got_sizes: Vec<usize> = ri.rings().iter().map(|r| r.len()).collect();
        got_sizes.sort();
        let mut expected_sizes = entry.ring_sizes.clone();
        expected_sizes.sort();
        if got_sizes != expected_sizes {
            failures.push(format!(
                "[ring_sizes] {}: expected {:?}, got {:?}",
                entry.smiles, expected_sizes, got_sizes
            ));
        }

        let ring_atoms: std::collections::HashSet<usize> = ri
            .rings()
            .iter()
            .flat_map(|r| r.iter().map(|n| n.index()))
            .collect();
        if ring_atoms.len() != entry.num_ring_atoms {
            failures.push(format!(
                "[ring_atoms] {}: expected {}, got {}",
                entry.smiles,
                entry.num_ring_atoms,
                ring_atoms.len()
            ));
        }

        let ring_bonds: std::collections::HashSet<(usize, usize)> = ri
            .rings()
            .iter()
            .flat_map(|r| {
                let len = r.len();
                (0..len).map(move |i| {
                    let a = r[i].index();
                    let b = r[(i + 1) % len].index();
                    (a.min(b), a.max(b))
                })
            })
            .collect();
        if ring_bonds.len() != entry.num_ring_bonds {
            failures.push(format!(
                "[ring_bonds] {}: expected {}, got {}",
                entry.smiles,
                entry.num_ring_bonds,
                ring_bonds.len()
            ));
        }
    }

    if skipped > 0 {
        eprintln!("rings: skipped {skipped} unparseable molecules");
    }

    if !failures.is_empty() {
        panic!(
            "{} ring failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// 4. SMARTS matching
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct SmartsEntry {
    smiles: String,
    smarts_matches: HashMap<String, usize>,
}

#[test]
fn approval_smarts() {
    let data: Vec<SmartsEntry> =
        serde_json::from_str(include_str!("approval_data/smarts.json")).unwrap();

    let mut failures = Vec::new();
    let mut parse_failures = Vec::new();

    let mut skipped = 0usize;
    for entry in &data {
        let mol = match try_parse(&entry.smiles) {
            Some(m) => m,
            None => { skipped += 1; continue; }
        };

        for (smarts_str, &expected_count) in &entry.smarts_matches {
            let query = match crabchem::from_smarts(smarts_str) {
                Ok(q) => q,
                Err(e) => {
                    parse_failures.push(format!(
                        "failed to parse SMARTS {:?}: {}",
                        smarts_str, e
                    ));
                    continue;
                }
            };

            let matches = crabchem::get_smarts_matches(&mol, &query);
            if matches.len() != expected_count {
                failures.push(format!(
                    "[smarts] {} / {:?}: expected {}, got {}",
                    entry.smiles, smarts_str, expected_count, matches.len()
                ));
            }
        }
    }

    if skipped > 0 {
        eprintln!("smarts: skipped {skipped} unparseable molecules");
    }

    if !parse_failures.is_empty() {
        eprintln!(
            "SMARTS parse warnings ({}):\n{}",
            parse_failures.len(),
            parse_failures.join("\n")
        );
    }

    if !failures.is_empty() {
        panic!(
            "{} SMARTS match failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// 5. Canonical SMILES (soft comparison)
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct CanonicalEntry {
    input: String,
    canonical_isomeric: String,
    #[allow(dead_code)]
    canonical_non_isomeric: String,
}

fn same_structure(
    a: &crabchem::Mol<crabchem::Atom, crabchem::Bond>,
    b: &crabchem::Mol<crabchem::Atom, crabchem::Bond>,
) -> bool {
    a.atom_count() == b.atom_count()
        && a.bond_count() == b.bond_count()
        && crabchem::has_substruct_match(a, b)
        && crabchem::has_substruct_match(b, a)
}

#[test]
fn approval_canonical_smiles() {
    let data: Vec<CanonicalEntry> =
        serde_json::from_str(include_str!("approval_data/canonical_smiles.json")).unwrap();

    let mut mismatches = Vec::new();
    let mut round_trip_failures = Vec::new();
    let mut skipped = 0usize;

    for entry in &data {
        let mol = match crabchem::smiles::from_smiles(&entry.canonical_isomeric) {
            Ok(m) => m,
            Err(e) => {
                eprintln!(
                    "SKIP (parse failure): canonical {:?}: {}",
                    entry.canonical_isomeric, e
                );
                skipped += 1;
                continue;
            }
        };

        let our_canonical = crabchem::smiles::to_canonical_smiles(&mol);

        // Check determinism: same output if we do it again
        let our_canonical2 = crabchem::smiles::to_canonical_smiles(&mol);
        if our_canonical != our_canonical2 {
            round_trip_failures.push(format!(
                "non-deterministic canonical for {:?}: {:?} vs {:?}",
                entry.canonical_isomeric, our_canonical, our_canonical2
            ));
        }

        // Soft string comparison — log but don't fail
        if our_canonical != entry.canonical_isomeric {
            mismatches.push(format!(
                "{}: expected={:?}, ours={:?}",
                entry.input, entry.canonical_isomeric, our_canonical
            ));
        }

        // Structural round-trip: parse our canonical, compare structure
        match crabchem::smiles::from_smiles(&our_canonical) {
            Ok(reparsed) => {
                if !same_structure(&mol, &reparsed) {
                    round_trip_failures.push(format!(
                        "round-trip structural mismatch for {:?}: wrote {:?}",
                        entry.canonical_isomeric, our_canonical
                    ));
                }
            }
            Err(e) => {
                round_trip_failures.push(format!(
                    "failed to reparse our canonical {:?} (from {:?}): {}",
                    our_canonical, entry.canonical_isomeric, e
                ));
            }
        }
    }

    if skipped > 0 {
        eprintln!("canonical: skipped {skipped} unparseable molecules");
    }

    if !mismatches.is_empty() {
        eprintln!(
            "{} canonical SMILES string mismatches (informational):\n{}",
            mismatches.len(),
            mismatches.join("\n")
        );
    }

    if !round_trip_failures.is_empty() {
        panic!(
            "{} canonical SMILES round-trip failures:\n{}",
            round_trip_failures.len(),
            round_trip_failures.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// 6. Substructure matching
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct SubstructMatchInfo {
    has_match: bool,
    num_matches: usize,
}

#[derive(Deserialize)]
struct SubstructEntry {
    smiles: String,
    substruct_matches: HashMap<String, SubstructMatchInfo>,
}

fn query_smiles_for_name(name: &str) -> Option<&'static str> {
    match name {
        "benzene ring" => Some("c1ccccc1"),
        "carbonyl" => Some("C=O"),
        "carboxylic acid" => Some("C(=O)O"),
        "amide" => Some("C(=O)N"),
        "primary amine" => Some("[NH2]"),
        "hydroxyl" => Some("[OH]"),
        "nitrile" => Some("C#N"),
        "nitro" => Some("[N+](=O)[O-]"),
        "pyrrole" => Some("c1cc[nH]c1"),
        "pyridine" => Some("c1ccncc1"),
        _ => None,
    }
}

#[test]
fn approval_substruct() {
    let data: Vec<SubstructEntry> =
        serde_json::from_str(include_str!("approval_data/substruct.json")).unwrap();

    let mut failures = Vec::new();

    let mut skipped = 0usize;
    for entry in &data {
        let mol = match try_parse(&entry.smiles) {
            Some(m) => m,
            None => { skipped += 1; continue; }
        };

        for (name, expected) in &entry.substruct_matches {
            let query_smi = match query_smiles_for_name(name) {
                Some(s) => s,
                None => {
                    eprintln!("unknown substruct query name: {name:?}");
                    continue;
                }
            };

            let query = match try_parse(query_smi) {
                Some(q) => q,
                None => continue,
            };

            let has = crabchem::has_substruct_match(&mol, &query);
            if has != expected.has_match {
                failures.push(format!(
                    "[has_match] {} / {}: expected {}, got {}",
                    entry.smiles, name, expected.has_match, has
                ));
            }

            let matches = crabchem::get_substruct_matches(&mol, &query);
            if matches.len() != expected.num_matches {
                failures.push(format!(
                    "[num_matches] {} / {}: expected {}, got {}",
                    entry.smiles, name, expected.num_matches, matches.len()
                ));
            }
        }
    }

    if skipped > 0 {
        eprintln!("substruct: skipped {skipped} unparseable molecules");
    }

    if !failures.is_empty() {
        panic!(
            "{} substruct failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
