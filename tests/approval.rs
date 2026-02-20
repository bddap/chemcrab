use std::collections::HashMap;

use serde::Deserialize;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn parse(smiles: &str) -> chemcrab::Mol<chemcrab::Atom, chemcrab::Bond> {
    chemcrab::smiles::from_smiles(smiles)
        .unwrap_or_else(|e| panic!("failed to parse SMILES {smiles:?}: {e}"))
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
    for entry in &data {
        let mol = parse(&entry.smiles);

        let formula = chemcrab::formula::mol_formula(&mol);
        if formula != entry.formula {
            failures.push(format!(
                "[formula] {}: expected {:?}, got {:?}",
                entry.smiles, entry.formula, formula
            ));
        }

        let amw = chemcrab::formula::average_mol_weight(&mol);
        if !approx_eq(amw, entry.average_mw, 0.01) {
            failures.push(format!(
                "[avg_mw] {}: expected {}, got {}",
                entry.smiles, entry.average_mw, amw
            ));
        }

        let emw = chemcrab::formula::exact_mol_weight(&mol);
        if !approx_eq(emw, entry.exact_mw, 0.01) {
            failures.push(format!(
                "[exact_mw] {}: expected {}, got {}",
                entry.smiles, entry.exact_mw, emw
            ));
        }
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
    for entry in &data {
        let mol = parse(&entry.smiles);

        if mol.atom_count() != entry.num_atoms {
            failures.push(format!(
                "[atom_count] {}: expected {}, got {}",
                entry.smiles,
                entry.num_atoms,
                mol.atom_count()
            ));
        }

        let aromatic = chemcrab::aromaticity::find_aromatic_atoms(&mol);
        let count = aromatic.iter().filter(|&&a| a).count();
        if count != entry.num_aromatic_atoms {
            failures.push(format!(
                "[aromatic] {}: expected {}, got {}",
                entry.smiles, entry.num_aromatic_atoms, count
            ));
        }
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
    for entry in &data {
        let mol = parse(&entry.smiles);
        let ri = chemcrab::rings::RingInfo::symmetrized_sssr(&mol);

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

    if !failures.is_empty() {
        panic!("{} ring failures:\n{}", failures.len(), failures.join("\n"));
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
    for entry in &data {
        let mol = parse(&entry.smiles);

        for (smarts_str, &expected_count) in &entry.smarts_matches {
            let query = chemcrab::smarts::from_smarts(smarts_str)
                .unwrap_or_else(|e| panic!("failed to parse SMARTS {smarts_str:?}: {e}"));

            let matches = chemcrab::smarts::get_smarts_matches(&mol, &query);
            if matches.len() != expected_count {
                failures.push(format!(
                    "[smarts] {} / {:?}: expected {}, got {}",
                    entry.smiles,
                    smarts_str,
                    expected_count,
                    matches.len()
                ));
            }
        }
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
// 5. Substructure matching
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

fn query_smiles_for_name(name: &str) -> &'static str {
    match name {
        "benzene ring" => "c1ccccc1",
        "carbonyl" => "C=O",
        "carboxylic acid" => "C(=O)O",
        "amide" => "C(=O)N",
        "primary amine" => "[NH2]",
        "hydroxyl" => "[OH]",
        "nitrile" => "C#N",
        "nitro" => "[N+](=O)[O-]",
        "pyrrole" => "c1cc[nH]c1",
        "pyridine" => "c1ccncc1",
        _ => panic!("unknown substruct query name: {name:?}"),
    }
}

#[test]
fn approval_substruct() {
    let data: Vec<SubstructEntry> =
        serde_json::from_str(include_str!("approval_data/substruct.json")).unwrap();

    let mut failures = Vec::new();
    for entry in &data {
        let mol = parse(&entry.smiles);

        for (name, expected) in &entry.substruct_matches {
            let query = parse(query_smiles_for_name(name));

            let has = chemcrab::substruct::has_substruct_match(&mol, &query);
            if has != expected.has_match {
                failures.push(format!(
                    "[has_match] {} / {}: expected {}, got {}",
                    entry.smiles, name, expected.has_match, has
                ));
            }

            let matches = chemcrab::substruct::get_substruct_matches_unique(&mol, &query);
            if matches.len() != expected.num_matches {
                failures.push(format!(
                    "[num_matches] {} / {}: expected {}, got {}",
                    entry.smiles,
                    name,
                    expected.num_matches,
                    matches.len()
                ));
            }
        }
    }

    if !failures.is_empty() {
        panic!(
            "{} substruct failures:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
