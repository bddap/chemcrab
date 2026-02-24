use criterion::{black_box, criterion_group, criterion_main, Criterion};

use chemcrab::smiles::{from_smiles, to_canonical_smiles, to_smiles};

const METHANE: &str = "C";
const CAFFEINE: &str = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
const ATORVASTATIN: &str =
    "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O";
const TAXOL: &str = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C";

fn bench_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");

    group.bench_function("methane", |b| {
        b.iter(|| black_box(from_smiles(black_box(METHANE)).unwrap()))
    });
    group.bench_function("caffeine", |b| {
        b.iter(|| black_box(from_smiles(black_box(CAFFEINE)).unwrap()))
    });
    group.bench_function("atorvastatin", |b| {
        b.iter(|| black_box(from_smiles(black_box(ATORVASTATIN)).unwrap()))
    });
    group.bench_function("taxol", |b| {
        b.iter(|| black_box(from_smiles(black_box(TAXOL)).unwrap()))
    });

    group.finish();
}

fn bench_write(c: &mut Criterion) {
    let methane = from_smiles(METHANE).unwrap();
    let caffeine = from_smiles(CAFFEINE).unwrap();
    let atorvastatin = from_smiles(ATORVASTATIN).unwrap();
    let taxol = from_smiles(TAXOL).unwrap();

    let mut group = c.benchmark_group("write");

    group.bench_function("methane", |b| {
        b.iter(|| black_box(to_smiles(black_box(&methane))))
    });
    group.bench_function("caffeine", |b| {
        b.iter(|| black_box(to_smiles(black_box(&caffeine))))
    });
    group.bench_function("atorvastatin", |b| {
        b.iter(|| black_box(to_smiles(black_box(&atorvastatin))))
    });
    group.bench_function("taxol", |b| {
        b.iter(|| black_box(to_smiles(black_box(&taxol))))
    });

    group.finish();
}

fn bench_canonical(c: &mut Criterion) {
    let methane = from_smiles(METHANE).unwrap();
    let caffeine = from_smiles(CAFFEINE).unwrap();
    let atorvastatin = from_smiles(ATORVASTATIN).unwrap();
    let taxol = from_smiles(TAXOL).unwrap();

    let mut group = c.benchmark_group("canonical");

    group.bench_function("methane", |b| {
        b.iter(|| black_box(to_canonical_smiles(black_box(&methane))))
    });
    group.bench_function("caffeine", |b| {
        b.iter(|| black_box(to_canonical_smiles(black_box(&caffeine))))
    });
    group.bench_function("atorvastatin", |b| {
        b.iter(|| black_box(to_canonical_smiles(black_box(&atorvastatin))))
    });
    group.bench_function("taxol", |b| {
        b.iter(|| black_box(to_canonical_smiles(black_box(&taxol))))
    });

    group.finish();
}

criterion_group!(benches, bench_parse, bench_write, bench_canonical);
criterion_main!(benches);
