use criterion::{black_box, criterion_group, criterion_main, Criterion};

use chemcrab::reaction::from_reaction_smarts;
use chemcrab::smiles::from_smiles;

fn bench_sn2(c: &mut Criterion) {
    let rxn = from_reaction_smarts("[C:1][Br:2].[N:3]>>[C:1][N:3]").unwrap();
    let r1 = from_smiles("CBr").unwrap();
    let r2 = from_smiles("N").unwrap();

    c.bench_function("sn2_simple", |b| {
        b.iter(|| black_box(rxn.run(&[&r1, &r2]).unwrap()))
    });
}

fn bench_multi_match(c: &mut Criterion) {
    let rxn = from_reaction_smarts("[cH:1]>>[c:1]F").unwrap();
    let naphthalene = from_smiles("c1ccc2ccccc2c1").unwrap();

    c.bench_function("multi_match_naphthalene_8pos", |b| {
        b.iter(|| black_box(rxn.run(&[&naphthalene]).unwrap()))
    });
}

fn bench_recursive_smarts(c: &mut Criterion) {
    let rxn =
        from_reaction_smarts("[C:1][c;!$(c(C)(:c(-C):c):c(:c)-[!#1]):2]>>[C:1].[*:2]").unwrap();
    let toluene = from_smiles("c1ccccc1C").unwrap();

    c.bench_function("recursive_smarts_toluene", |b| {
        b.iter(|| black_box(rxn.run(&[&toluene]).unwrap()))
    });
}

fn bench_large_molecule(c: &mut Criterion) {
    let rxn = from_reaction_smarts("[NH2:1][C:2](=O)>>[NH:1][C:2](=O)C").unwrap();
    // atorvastatin (37 heavy atoms)
    let atorvastatin = from_smiles(
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    )
    .unwrap();

    c.bench_function("large_mol_atorvastatin", |b| {
        b.iter(|| black_box(rxn.run(&[&atorvastatin]).unwrap()))
    });
}

fn bench_ring_formation(c: &mut Criterion) {
    let rxn =
        from_reaction_smarts("[CH:1]=[CH:2].[CH:3]=[CH:4]>>[CH:1]1[CH:2][CH:3][CH:4]1").unwrap();
    let r1 = from_smiles("C=C").unwrap();
    let r2 = from_smiles("C=C").unwrap();

    c.bench_function("ring_formation_cyclobutane", |b| {
        b.iter(|| black_box(rxn.run(&[&r1, &r2]).unwrap()))
    });
}

criterion_group!(
    benches,
    bench_sn2,
    bench_multi_match,
    bench_recursive_smarts,
    bench_large_molecule,
    bench_ring_formation,
);
criterion_main!(benches);
