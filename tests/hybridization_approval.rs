use chemcrab::{assign_hybridization, from_smiles, Hybridization};
use serde::Deserialize;

#[derive(Deserialize)]
struct AtomEntry {
    #[allow(dead_code)]
    idx: usize,
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

fn hyb_label(h: Hybridization) -> &'static str {
    match h {
        Hybridization::S => "S",
        Hybridization::SP => "SP",
        Hybridization::SP2 => "SP2",
        Hybridization::SP3 => "SP3",
        Hybridization::SP3D => "SP3D",
        Hybridization::SP3D2 => "SP3D2",
        Hybridization::Other => "Other",
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
            Err(_) => continue,
        };

        let computed = assign_hybridization(&mol);

        let mut expected: Vec<(&str, Hybridization)> = entry
            .atoms
            .iter()
            .map(|a| (a.symbol.as_str(), parse_hybridization(&a.hybridization)))
            .collect();
        expected.sort_by(|a, b| a.0.cmp(&b.0).then(hyb_label(a.1).cmp(hyb_label(b.1))));

        let atom_symbols: Vec<&str> = (0..mol.atom_count())
            .map(|i| {
                use chemcrab::traits::HasAtomicNum;
                let anum = mol.atom(petgraph::graph::NodeIndex::new(i)).atomic_num();
                chemcrab::Element::from_atomic_num(anum)
                    .map(|e| e.symbol())
                    .unwrap_or("?")
            })
            .collect();

        let mut got: Vec<(&str, Hybridization)> = atom_symbols
            .iter()
            .zip(computed.iter())
            .map(|(&sym, &hyb)| (sym, hyb))
            .collect();
        got.sort_by(|a, b| a.0.cmp(&b.0).then(hyb_label(a.1).cmp(hyb_label(b.1))));

        if expected == got {
            pass += 1;
        } else {
            fail += 1;
            failures.push(format!(
                "{}: expected {:?}, got {:?}",
                entry.smiles, expected, got
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
