# Sally Lightfoot Cheminformatics Library

## Substructure matching vs SMARTS matching

`get_substruct_matches` operates on `Mol<Atom, Bond>` pairs and compares
atoms by element, aromaticity, and formal charge. It does **not** check
hydrogen count or other SMARTS-level properties. This means a SMILES-derived
query like `C(=O)O` will match carboxylic acids, carboxylates, and esters
alike — it only asks "is there a carbon with two oxygen neighbors and a
double bond?"

For property-rich queries — hydrogen counts, ring membership, degree
constraints — use SMARTS matching (`get_smarts_matches`), which evaluates
the full `AtomExpr`/`BondExpr` query language.

`get_substruct_matches` returns every VF2 mapping including automorphisms
(benzene-on-benzene yields 12). `get_substruct_matches_unique` collapses
mappings that cover the same set of target atoms. SMARTS matching uniquifies
by default.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
