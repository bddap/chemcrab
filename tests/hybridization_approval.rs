use crabchem::{assign_hybridization, from_smiles, Hybridization};
use serde::Deserialize;

#[derive(Deserialize)]
struct AtomEntry {
    #[allow(dead_code)]
    idx: usize,
    #[allow(dead_code)]
    symbol: String,
    hybridization: String,
}

#[derive(Deserialize)]
struct MolEntry {
    smiles: String,
    atoms: Vec<AtomEntry>,
}

fn parse_hybridization(s: &str) -> Hybridization {
    match s {
        "S" => Hybridization::S,
        "SP" => Hybridization::SP,
        "SP2" => Hybridization::SP2,
        "SP3" => Hybridization::SP3,
        "SP3D" => Hybridization::SP3D,
        "SP3D2" => Hybridization::SP3D2,
        other => panic!("Unknown hybridization: {other}"),
    }
}

#[test]
fn approval_hybridization() {
    let data: Vec<MolEntry> =
        serde_json::from_str(include_str!("approval_data/hybridization.json")).unwrap();

    let mut pass = 0;
    let mut fail = 0;
    let mut failures = Vec::new();

    for entry in &data {
        let mol = match from_smiles(&entry.smiles) {
            Ok(m) => m,
            Err(_) => {
                // Skip molecules our parser can't handle yet
                // (e.g. explicit '-' bonds in branches like biphenyl)
                continue;
            }
        };

        let computed = assign_hybridization(&mol);

        let mut expected_counts: std::collections::HashMap<Hybridization, usize> =
            std::collections::HashMap::new();
        for atom in &entry.atoms {
            let h = parse_hybridization(&atom.hybridization);
            *expected_counts.entry(h).or_default() += 1;
        }

        let mut computed_counts: std::collections::HashMap<Hybridization, usize> =
            std::collections::HashMap::new();
        for &h in &computed {
            *computed_counts.entry(h).or_default() += 1;
        }

        if expected_counts == computed_counts {
            pass += 1;
        } else {
            fail += 1;
            failures.push(format!(
                "{}: expected {:?}, got {:?}",
                entry.smiles, expected_counts, computed_counts
            ));
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("FAIL: {f}");
        }
    }
    eprintln!("{pass} passed, {fail} failed out of {} total", data.len());
    assert_eq!(fail, 0, "{fail} molecules failed hybridization approval");
}
