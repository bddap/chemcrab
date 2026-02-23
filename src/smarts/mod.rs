//! SMARTS (SMiles ARbitrary Target Specification) parsing and substructure search.
//!
//! SMARTS is a query language for substructure matching that extends SMILES
//! with logical operators and atom/bond primitives. Where SMILES describes a
//! single molecule, a SMARTS pattern describes a *set* of molecules (or
//! substructures) that match.
//!
//! Key extensions over SMILES:
//!
//! - **Atomic number**: `[#6]` matches any carbon regardless of aromaticity.
//! - **Ring membership**: `[R]` matches ring atoms, `[R0]` matches chain atoms.
//! - **Boolean logic**: `;` is AND (low precedence), `&` is AND (high
//!   precedence), `,` is OR, `!` is NOT. Example: `[C,N;H1]` matches carbon
//!   or nitrogen with exactly one hydrogen.
//! - **Recursive SMARTS**: `[$([CX3](=O))]` matches atoms that are part of a
//!   substructure matching the inner pattern.
//! - **Bond primitives**: `-` single, `=` double, `#` triple, `:` aromatic,
//!   `~` any, `@` ring bond.
//!
//! # Entry points
//!
//! | Function | Purpose |
//! |---|---|
//! | [`from_smarts`] | Parse a SMARTS string |
//! | [`has_smarts_match`] | Does target contain the pattern? |
//! | [`get_smarts_match`] | First match (atom mapping) |
//! | [`get_smarts_matches`] | All unique matches |
//! | [`has_smarts_match_chiral`] | Chirality-aware containment check |
//! | [`get_smarts_match_chiral`] | First chirality-aware match |
//! | [`get_smarts_matches_chiral`] | All unique chirality-aware matches |

mod error;
mod parser;
pub mod query;
mod writer;

pub use error::SmartsError;
pub use query::compute_hybridization;
pub use query::{AtomExpr, BondExpr, Hybridization, RangeKind};
pub use writer::to_smarts;

use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::atom::Chirality;
use crate::bond::Bond;
use crate::mol::{permutation_parity, AtomId, Mol};
use crate::rings::RingInfo;
use crate::substruct::{
    get_substruct_match_with, get_substruct_match_with_filter, get_substruct_matches_with,
    get_substruct_matches_with_filter, uniquify_atom_mappings, AtomMapping,
};

use query::MatchContext;

fn expr_references_hydrogen(expr: &AtomExpr) -> bool {
    match expr {
        AtomExpr::Element { atomic_num: 1, .. } => true,
        AtomExpr::And(parts) | AtomExpr::Or(parts) => parts.iter().any(expr_references_hydrogen),
        AtomExpr::Not(inner) => expr_references_hydrogen(inner),
        AtomExpr::Recursive(mol) => mol
            .atoms()
            .any(|idx| expr_references_hydrogen(mol.atom(idx))),
        _ => false,
    }
}

fn query_references_hydrogen(query: &Mol<AtomExpr, BondExpr>) -> bool {
    query
        .atoms()
        .any(|idx| expr_references_hydrogen(query.atom(idx)))
}

type RecursiveRef<'a> = (*const Mol<AtomExpr, BondExpr>, &'a Mol<AtomExpr, BondExpr>);

/// Parses a SMARTS pattern string into a query molecule.
///
/// The returned `Mol<AtomExpr, BondExpr>` uses expression trees for atoms
/// and bonds rather than concrete chemical types, enabling flexible
/// substructure matching.
///
/// # Errors
///
/// Returns [`SmartsError`] on invalid syntax, unclosed brackets, or
/// malformed recursive SMARTS.
pub fn from_smarts(s: &str) -> Result<Mol<AtomExpr, BondExpr>, SmartsError> {
    parser::parse(s)
}

/// Returns `true` if `target` contains a substructure matching `query`.
pub fn has_smarts_match(target: &Mol<Atom, Bond>, query: &Mol<AtomExpr, BondExpr>) -> bool {
    get_smarts_match(target, query).is_some()
}

/// Returns the first substructure match of `query` in `target`, if any.
///
/// The returned [`AtomMapping`] maps each query atom index to its
/// corresponding target atom index.
pub fn get_smarts_match(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Option<AtomMapping> {
    if query_references_hydrogen(query) {
        let explicit = crate::hydrogen::add_hs(target);
        return get_smarts_match_impl(&explicit, query);
    }
    get_smarts_match_impl(target, query)
}

fn get_smarts_match_impl(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Option<AtomMapping> {
    let ring_info = RingInfo::sssr(target);
    let recursive_matches = pre_evaluate_recursive(target, query, &ring_info);

    let ctx = MatchContext {
        mol: target,
        ring_info: &ring_info,
        recursive_matches,
    };

    get_substruct_match_with(
        target,
        query,
        |t_idx: NodeIndex, t_atom: &Atom, _q_idx: NodeIndex, q_expr: &AtomExpr| {
            match_atom_expr(q_expr, t_atom, &ctx, t_idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints, _q_endpoints| {
            let (t_a, t_atom_a) = t_endpoints.0;
            let (t_b, t_atom_b) = t_endpoints.1;
            q_bond.matches_with_ring_info(t_bond, (t_atom_a, t_atom_b), &ring_info, (t_a, t_b))
        },
    )
}

/// Returns all unique substructure matches of `query` in `target`.
///
/// Matches are deduplicated so that each set of matched target atoms appears
/// at most once, regardless of automorphic permutations.
pub fn get_smarts_matches(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Vec<AtomMapping> {
    if query_references_hydrogen(query) {
        let explicit = crate::hydrogen::add_hs(target);
        return uniquify_atom_mappings(&get_smarts_matches_all_impl(&explicit, query));
    }
    uniquify_atom_mappings(&get_smarts_matches_all_impl(target, query))
}

/// Returns all matches of `query` in `target`, including symmetric permutations.
///
/// Unlike [`get_smarts_matches`], this does not deduplicate matches that
/// map to the same set of target atoms in different order — useful for
/// reactions where each permutation produces a distinct product.
pub fn get_smarts_matches_all(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Vec<AtomMapping> {
    if query_references_hydrogen(query) {
        let explicit = crate::hydrogen::add_hs(target);
        return get_smarts_matches_all_impl(&explicit, query);
    }
    get_smarts_matches_all_impl(target, query)
}

fn get_smarts_matches_all_impl(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Vec<AtomMapping> {
    let ring_info = RingInfo::sssr(target);
    let recursive_matches = pre_evaluate_recursive(target, query, &ring_info);

    let ctx = MatchContext {
        mol: target,
        ring_info: &ring_info,
        recursive_matches,
    };

    get_substruct_matches_with(
        target,
        query,
        |t_idx: NodeIndex, t_atom: &Atom, _q_idx: NodeIndex, q_expr: &AtomExpr| {
            match_atom_expr(q_expr, t_atom, &ctx, t_idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints, _q_endpoints| {
            let (t_a, t_atom_a) = t_endpoints.0;
            let (t_b, t_atom_b) = t_endpoints.1;
            q_bond.matches_with_ring_info(t_bond, (t_atom_a, t_atom_b), &ring_info, (t_a, t_b))
        },
    )
}

/// Returns `true` if `target` contains a chirality-aware match for `query`.
///
/// In addition to the standard SMARTS atom and bond matching, this verifies
/// that tetrahedral stereocenters in the query have the same handedness in
/// the matched target atoms.
pub fn has_smarts_match_chiral(target: &Mol<Atom, Bond>, query: &Mol<AtomExpr, BondExpr>) -> bool {
    get_smarts_match_chiral(target, query).is_some()
}

/// Returns the first chirality-aware match of `query` in `target`, if any.
///
/// See [`has_smarts_match_chiral`] for details on chirality verification.
pub fn get_smarts_match_chiral(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Option<AtomMapping> {
    if query_references_hydrogen(query) {
        let explicit = crate::hydrogen::add_hs(target);
        return get_smarts_match_chiral_impl(&explicit, query);
    }
    get_smarts_match_chiral_impl(target, query)
}

fn get_smarts_match_chiral_impl(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Option<AtomMapping> {
    let ring_info = RingInfo::sssr(target);
    let recursive_matches = pre_evaluate_recursive(target, query, &ring_info);

    let ctx = MatchContext {
        mol: target,
        ring_info: &ring_info,
        recursive_matches,
    };

    let chiral_query_atoms = collect_chiral_query_atoms(query);

    get_substruct_match_with_filter(
        target,
        query,
        |t_idx: NodeIndex, t_atom: &Atom, _q_idx: NodeIndex, q_expr: &AtomExpr| {
            match_atom_expr(q_expr, t_atom, &ctx, t_idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints, _q_endpoints| {
            let (t_a, t_atom_a) = t_endpoints.0;
            let (t_b, t_atom_b) = t_endpoints.1;
            q_bond.matches_with_ring_info(t_bond, (t_atom_a, t_atom_b), &ring_info, (t_a, t_b))
        },
        |mapping: &AtomMapping| validate_chirality(mapping, target, query, &chiral_query_atoms),
    )
}

/// Returns all unique chirality-aware matches of `query` in `target`.
///
/// See [`has_smarts_match_chiral`] for details on chirality verification.
pub fn get_smarts_matches_chiral(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Vec<AtomMapping> {
    if query_references_hydrogen(query) {
        let explicit = crate::hydrogen::add_hs(target);
        return get_smarts_matches_chiral_impl(&explicit, query);
    }
    get_smarts_matches_chiral_impl(target, query)
}

fn get_smarts_matches_chiral_impl(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
) -> Vec<AtomMapping> {
    let ring_info = RingInfo::sssr(target);
    let recursive_matches = pre_evaluate_recursive(target, query, &ring_info);

    let ctx = MatchContext {
        mol: target,
        ring_info: &ring_info,
        recursive_matches,
    };

    let chiral_query_atoms = collect_chiral_query_atoms(query);

    let all = get_substruct_matches_with_filter(
        target,
        query,
        |t_idx: NodeIndex, t_atom: &Atom, _q_idx: NodeIndex, q_expr: &AtomExpr| {
            match_atom_expr(q_expr, t_atom, &ctx, t_idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints, _q_endpoints| {
            let (t_a, t_atom_a) = t_endpoints.0;
            let (t_b, t_atom_b) = t_endpoints.1;
            q_bond.matches_with_ring_info(t_bond, (t_atom_a, t_atom_b), &ring_info, (t_a, t_b))
        },
        |mapping: &AtomMapping| validate_chirality(mapping, target, query, &chiral_query_atoms),
    );
    uniquify_atom_mappings(&all)
}

struct ChiralQueryAtom {
    query_idx: NodeIndex,
    chirality: Chirality,
    has_implicit_h: bool,
}

fn collect_chiral_query_atoms(query: &Mol<AtomExpr, BondExpr>) -> Vec<ChiralQueryAtom> {
    let mut result = Vec::new();
    for q_idx in query.atoms() {
        if let Some((chiral, has_h)) = extract_chirality_and_h(query.atom(q_idx)) {
            if chiral != Chirality::None {
                result.push(ChiralQueryAtom {
                    query_idx: q_idx,
                    chirality: chiral,
                    has_implicit_h: has_h,
                });
            }
        }
    }
    result
}

fn extract_chirality_and_h(expr: &AtomExpr) -> Option<(Chirality, bool)> {
    match expr {
        AtomExpr::Chirality(c) => Some((*c, false)),
        AtomExpr::And(parts) => {
            let mut chiral = None;
            let mut has_h = false;
            for p in parts {
                match p {
                    AtomExpr::Chirality(c) => chiral = Some(*c),
                    AtomExpr::TotalHCount(n) if *n > 0 => has_h = true,
                    _ => {}
                }
            }
            chiral.map(|c| (c, has_h))
        }
        _ => None,
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum NeighborRef {
    Atom(NodeIndex),
    ImplicitH,
}

fn validate_chirality(
    mapping: &AtomMapping,
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
    chiral_query_atoms: &[ChiralQueryAtom],
) -> bool {
    for cqa in chiral_query_atoms {
        let q_idx = cqa.query_idx;
        let t_idx = match mapping.iter().find(|&&(q, _)| q == q_idx) {
            Some(&(_, t)) => t,
            None => continue,
        };

        let target_stereo = match target.tetrahedral_stereo_for(t_idx) {
            Some(s) => s,
            None => return false,
        };

        let q_neighbor_count =
            query.neighbors(q_idx).count() + if cqa.has_implicit_h { 1 } else { 0 };
        if q_neighbor_count < 3 {
            continue;
        }

        // Build query's neighbor ordering mapped to target NeighborRefs.
        // After SMARTS normalize_chirality, the stored chirality (Ccw/Cw) is
        // relative to petgraph neighbor iteration order (with implicit H first
        // if present).
        let mapped: Vec<NeighborRef> = {
            let mut v = Vec::new();
            if cqa.has_implicit_h {
                v.push(NeighborRef::ImplicitH);
            }
            for q_nb in query.neighbors(q_idx) {
                match mapping.iter().find(|&&(q, _)| q == q_nb) {
                    Some(&(_, t)) => v.push(NeighborRef::Atom(t)),
                    None => return false,
                }
            }
            v
        };

        if mapped.len() < 3 {
            continue;
        }

        let target_full: Vec<NeighborRef> = target_stereo
            .above
            .iter()
            .map(|aid| match aid {
                AtomId::Node(idx) => NeighborRef::Atom(*idx),
                AtomId::VirtualH(_, _) => NeighborRef::ImplicitH,
            })
            .collect();

        // Query's full tuple from its perspective.
        // Ccw (@): looking from first mapped neighbor, the rest wind CCW.
        //   → full tuple = [mapped[0], mapped[1], mapped[2], mapped[3]] for CCW.
        // Cw (@@): looking from first, the rest wind CW.
        //   → CW means the CCW order has last two swapped:
        //     full CCW tuple = [mapped[0], mapped[1], mapped[3], mapped[2]]
        //
        // For 3-neighbor case: no excluded, just [mapped[0], mapped[1], mapped[2]].
        let query_full: Vec<NeighborRef> = if mapped.len() == 4 {
            match cqa.chirality {
                Chirality::Ccw => vec![mapped[0], mapped[1], mapped[2], mapped[3]],
                Chirality::Cw => vec![mapped[0], mapped[1], mapped[3], mapped[2]],
                Chirality::None => continue,
            }
        } else {
            // 3 neighbors
            match cqa.chirality {
                Chirality::Ccw => vec![mapped[0], mapped[1], mapped[2]],
                Chirality::Cw => vec![mapped[0], mapped[2], mapped[1]],
                Chirality::None => continue,
            }
        };

        if target_full.len() != query_full.len() {
            if mapped.len() < target_full.len() {
                let target_overlap: Vec<NeighborRef> = target_full
                    .iter()
                    .filter(|n| query_full.contains(n))
                    .copied()
                    .collect();
                let query_overlap: Vec<NeighborRef> = query_full
                    .iter()
                    .filter(|n| target_full.contains(n))
                    .copied()
                    .collect();
                if target_overlap.len() != query_overlap.len() || target_overlap.len() < 3 {
                    continue;
                }
                let even = permutation_parity(&target_overlap, &query_overlap);
                if !even {
                    return false;
                }
                continue;
            }
            continue;
        }

        let even = permutation_parity(&target_full, &query_full);
        if !even {
            return false;
        }
    }
    true
}

fn match_atom_expr(expr: &AtomExpr, atom: &Atom, ctx: &MatchContext, idx: NodeIndex) -> bool {
    match expr {
        AtomExpr::Recursive(inner) => {
            let ptr = inner as *const Mol<AtomExpr, BondExpr> as usize;
            ctx.recursive_matches
                .get(&ptr)
                .is_some_and(|set| set.contains(&idx))
        }
        AtomExpr::And(exprs) => exprs.iter().all(|e| match_atom_expr(e, atom, ctx, idx)),
        AtomExpr::Or(exprs) => exprs.iter().any(|e| match_atom_expr(e, atom, ctx, idx)),
        AtomExpr::Not(inner) => !match_atom_expr(inner, atom, ctx, idx),
        _ => expr.matches(atom, ctx, idx),
    }
}

fn pre_evaluate_recursive(
    target: &Mol<Atom, Bond>,
    query: &Mol<AtomExpr, BondExpr>,
    ring_info: &RingInfo,
) -> HashMap<usize, HashSet<NodeIndex>> {
    let mut recursive_ptrs: Vec<RecursiveRef> = Vec::new();
    for atom_idx in query.atoms() {
        collect_recursive_refs(query.atom(atom_idx), &mut recursive_ptrs);
    }

    let mut results = HashMap::new();
    for &(ptr, inner_query) in &recursive_ptrs {
        let key = ptr as usize;
        if results.contains_key(&key) {
            continue;
        }

        let inner_recursive = pre_evaluate_recursive(target, inner_query, ring_info);
        let inner_ctx = MatchContext {
            mol: target,
            ring_info,
            recursive_matches: inner_recursive,
        };

        let mut matching_atoms = HashSet::new();
        let matches = get_substruct_matches_with(
            target,
            inner_query,
            |t_idx: NodeIndex, t_atom: &Atom, _q_idx: NodeIndex, q_expr: &AtomExpr| {
                match_atom_expr(q_expr, t_atom, &inner_ctx, t_idx)
            },
            |t_bond: &Bond, q_bond: &BondExpr, t_endpoints, _q_endpoints| {
                let (t_a, t_atom_a) = t_endpoints.0;
                let (t_b, t_atom_b) = t_endpoints.1;
                q_bond.matches_with_ring_info(t_bond, (t_atom_a, t_atom_b), ring_info, (t_a, t_b))
            },
        );

        let root = NodeIndex::new(0);
        for mapping in &matches {
            if let Some(&(_, t_root)) = mapping.iter().find(|&&(q, _)| q == root) {
                matching_atoms.insert(t_root);
            }
        }
        results.insert(key, matching_atoms);
    }

    results
}

fn collect_recursive_refs<'a>(expr: &'a AtomExpr, ptrs: &mut Vec<RecursiveRef<'a>>) {
    match expr {
        AtomExpr::Recursive(inner) => {
            ptrs.push((inner as *const Mol<AtomExpr, BondExpr>, inner));
        }
        AtomExpr::And(exprs) | AtomExpr::Or(exprs) => {
            for e in exprs {
                collect_recursive_refs(e, ptrs);
            }
        }
        AtomExpr::Not(inner) => {
            collect_recursive_refs(inner, ptrs);
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mol(smiles: &str) -> Mol<Atom, Bond> {
        crate::smiles::from_smiles(smiles).unwrap_or_else(|e| panic!("bad SMILES {smiles:?}: {e}"))
    }

    fn smarts(s: &str) -> Mol<AtomExpr, BondExpr> {
        from_smarts(s).unwrap_or_else(|e| panic!("bad SMARTS {s:?}: {e}"))
    }

    // ---- Parser tests ----

    #[test]
    fn parse_atomic_num() {
        let q = smarts("[#6]");
        assert_eq!(q.atom_count(), 1);
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 6,
                aromatic: None
            }
        ));
    }

    #[test]
    fn parse_aliphatic_carbon() {
        let q = smarts("[C]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 6,
                aromatic: Some(false)
            }
        ));
    }

    #[test]
    fn parse_aromatic_carbon() {
        let q = smarts("[c]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 6,
                aromatic: Some(true)
            }
        ));
    }

    #[test]
    fn parse_wildcard_bare() {
        let q = smarts("*");
        assert_eq!(q.atom_count(), 1);
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::True));
    }

    #[test]
    fn parse_wildcard_bracket() {
        let q = smarts("[*]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::True));
    }

    #[test]
    fn parse_or() {
        let q = smarts("[C,N]");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::Or(parts) => {
                assert_eq!(parts.len(), 2);
                assert!(matches!(
                    parts[0],
                    AtomExpr::Element {
                        atomic_num: 6,
                        aromatic: Some(false)
                    }
                ));
                assert!(matches!(
                    parts[1],
                    AtomExpr::Element {
                        atomic_num: 7,
                        aromatic: Some(false)
                    }
                ));
            }
            _ => panic!("expected Or, got {expr:?}"),
        }
    }

    #[test]
    fn parse_semicolon_and() {
        let q = smarts("[C;R]");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert_eq!(parts.len(), 2);
                assert!(matches!(
                    parts[0],
                    AtomExpr::Element {
                        atomic_num: 6,
                        aromatic: Some(false)
                    }
                ));
                assert!(matches!(parts[1], AtomExpr::InRing));
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    #[test]
    fn parse_not() {
        let q = smarts("[!C]");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::Not(inner) => {
                assert!(matches!(
                    **inner,
                    AtomExpr::Element {
                        atomic_num: 6,
                        aromatic: Some(false)
                    }
                ));
            }
            _ => panic!("expected Not, got {expr:?}"),
        }
    }

    #[test]
    fn parse_precedence_comma_semicolon() {
        let q = smarts("[C,N;H1]");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert_eq!(parts.len(), 2);
                match &parts[0] {
                    AtomExpr::Or(or_parts) => {
                        assert_eq!(or_parts.len(), 2);
                    }
                    _ => panic!("expected Or in first part"),
                }
                assert!(matches!(parts[1], AtomExpr::TotalHCount(1)));
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    #[test]
    fn parse_degree() {
        let q = smarts("[D2]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Degree(2)));
    }

    #[test]
    fn parse_valence() {
        let q = smarts("[v4]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Valence(4)));
    }

    #[test]
    fn parse_connectivity() {
        let q = smarts("[X4]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Connectivity(4)
        ));
    }

    #[test]
    fn parse_total_h_count() {
        let q = smarts("[H1]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::TotalHCount(1)
        ));
    }

    #[test]
    fn parse_implicit_h_count() {
        let q = smarts("[h1]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::ImplicitHCount(1)
        ));
    }

    #[test]
    fn parse_in_ring() {
        let q = smarts("[R]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::InRing));
    }

    #[test]
    fn parse_not_in_ring() {
        let q = smarts("[R0]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::NotInRing));
    }

    #[test]
    fn parse_ring_membership() {
        let q = smarts("[R2]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::RingMembership(2)
        ));
    }

    #[test]
    fn parse_smallest_ring_size() {
        let q = smarts("[r6]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::SmallestRingSize(6)
        ));
    }

    #[test]
    fn parse_ring_bond_count() {
        let q = smarts("[x2]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::RingBondCount(2)
        ));
    }

    #[test]
    fn parse_positive_charge() {
        let q = smarts("[+1]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Charge(1)));
    }

    #[test]
    fn parse_negative_charge() {
        let q = smarts("[-1]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Charge(-1)));
    }

    #[test]
    fn parse_default_bond() {
        let q = smarts("CC");
        assert_eq!(q.atom_count(), 2);
        assert_eq!(q.bond_count(), 1);
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert!(matches!(q.bond(edge), BondExpr::SingleOrAromatic));
    }

    #[test]
    fn parse_double_bond() {
        let q = smarts("C=C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert!(matches!(q.bond(edge), BondExpr::Double));
    }

    #[test]
    fn parse_triple_bond() {
        let q = smarts("C#C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert!(matches!(q.bond(edge), BondExpr::Triple));
    }

    #[test]
    fn parse_any_bond() {
        let q = smarts("C~C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert!(matches!(q.bond(edge), BondExpr::True));
    }

    #[test]
    fn parse_aromatic_bond() {
        let q = smarts("C:C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert!(matches!(q.bond(edge), BondExpr::Aromatic));
    }

    // ---- Matching tests ----

    #[test]
    fn match_carbon_in_ethane() {
        let target = mol("CC");
        let query = smarts("[#6]");
        assert!(has_smarts_match(&target, &query));
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_carbon_not_oxygen() {
        let target = mol("CCO");
        let query = smarts("[#6]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_nitrogen_in_pyridine() {
        let target = mol("c1ccncc1");
        let query = smarts("[#7]");
        assert!(has_smarts_match(&target, &query));
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_aromatic_in_benzene() {
        let target = mol("c1ccccc1");
        let query = smarts("[a]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_aromatic_not_in_cyclohexane() {
        let target = mol("C1CCCCC1");
        let query = smarts("[a]");
        assert!(!has_smarts_match(&target, &query));
    }

    #[test]
    fn match_aliphatic() {
        let target = mol("C1CCCCC1");
        let query = smarts("[A]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_wildcard() {
        let target = mol("CCO");
        let query = smarts("*");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 3);
    }

    #[test]
    fn match_degree_2() {
        let target = mol("CCC");
        let query = smarts("[D2]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_degree_1() {
        let target = mol("CCC");
        let query = smarts("[D1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_valence_4() {
        let target = mol("C");
        let query = smarts("[v4]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_connectivity_4() {
        let target = mol("C");
        let query = smarts("[X4]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_h_count_1() {
        let target = mol("CCO");
        let query = smarts("[H1]");
        let matches = get_smarts_matches(&target, &query);
        assert!(!matches.is_empty());
    }

    #[test]
    fn match_in_ring() {
        let target = mol("C1CCCCC1");
        let query = smarts("[R]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_not_in_ring() {
        let target = mol("CC1CCCCC1");
        let query = smarts("[R0]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_ring_size_6() {
        let target = mol("C1CCCCC1");
        let query = smarts("[r6]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_ring_membership_2() {
        let target = mol("c1ccc2ccccc2c1");
        let query = smarts("[R2]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_positive_charge() {
        let target = mol("[Na+]");
        let query = smarts("[+1]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_negative_charge() {
        let target = mol("[Cl-]");
        let query = smarts("[-1]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_or_c_n() {
        let target = mol("CN");
        let query = smarts("[C,N]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_not_carbon() {
        let target = mol("CCO");
        let query = smarts("[!C]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_carbon_in_ring() {
        let target = mol("CC1CCCCC1");
        let query = smarts("[C;R]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_aromatic_bond() {
        let target = mol("c1ccccc1");
        let query = smarts("c:c");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_any_bond() {
        let target = mol("CC");
        let query = smarts("C~C");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_single_or_aromatic_default() {
        let target = mol("CC");
        let query = smarts("CC");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_single_or_aromatic_aromatic() {
        let target = mol("c1ccccc1");
        let query = smarts("cc");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_explicit_single_only() {
        let target = mol("c1ccccc1");
        let query = smarts("c-c");
        assert!(!has_smarts_match(&target, &query));
    }

    #[test]
    fn match_explicit_single_on_single() {
        let target = mol("CC");
        let query = smarts("C-C");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_recursive_oh() {
        let target = mol("Oc1ccccc1");
        let query = smarts("[$([OH])]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_recursive_aromatic_cc() {
        let target = mol("c1ccccc1");
        let query = smarts("[$(cc)]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn match_multi_atom_cn() {
        let target = mol("CN");
        let query = smarts("[#6][#7]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_benzene_ring() {
        let target = mol("c1ccccc1");
        let query = smarts("c1ccccc1");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn match_cyclopropane() {
        let target = mol("C1CC1");
        let query = smarts("C1CC1");
        assert!(has_smarts_match(&target, &query));
    }

    // ---- Edge cases ----

    #[test]
    fn empty_string_error() {
        assert!(from_smarts("").is_err());
    }

    #[test]
    fn invalid_smarts_error() {
        assert!(from_smarts("[").is_err());
    }

    #[test]
    fn single_atom_query() {
        let target = mol("C");
        let query = smarts("C");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn disconnected_query() {
        let target = mol("[Na+].[Cl-]");
        let query = smarts("[Na+].[Cl-]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn parse_hydrogen_bracket() {
        let q = smarts("[H]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 1,
                aromatic: Some(false)
            }
        ));
    }

    #[test]
    fn parse_bare_aromatic_a() {
        let q = smarts("a");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Aromatic));
    }

    #[test]
    fn parse_bare_aliphatic_aa() {
        let q = smarts("A");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Aliphatic));
    }

    // ---- Writer tests ----

    #[test]
    fn writer_wildcard() {
        let q = smarts("*");
        assert_eq!(to_smarts(&q), "*");
    }

    #[test]
    fn writer_bare_element() {
        let q = smarts("C");
        assert_eq!(to_smarts(&q), "C");
    }

    #[test]
    fn writer_aromatic_element() {
        let q = smarts("c");
        assert_eq!(to_smarts(&q), "c");
    }

    #[test]
    fn writer_double_bond() {
        let q = smarts("C=C");
        assert_eq!(to_smarts(&q), "C=C");
    }

    #[test]
    fn writer_triple_bond() {
        let q = smarts("C#C");
        assert_eq!(to_smarts(&q), "C#C");
    }

    #[test]
    fn writer_any_bond() {
        let q = smarts("C~C");
        assert_eq!(to_smarts(&q), "C~C");
    }

    #[test]
    fn writer_default_bond_omitted() {
        let q = smarts("CC");
        assert_eq!(to_smarts(&q), "CC");
    }

    #[test]
    fn writer_aromatic_bond() {
        let q = smarts("c:c");
        assert_eq!(to_smarts(&q), "c:c");
    }

    #[test]
    fn writer_round_trip_bracket() {
        let cases = [
            "[#6]", "[D2]", "[v4]", "[X4]", "[H1]", "[h1]", "[R]", "[R0]", "[R2]", "[r6]", "[x2]",
            "[+1]", "[-1]",
        ];
        for s in &cases {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    #[test]
    fn writer_round_trip_logical() {
        let s = "[C,N]";
        let q = smarts(s);
        let written = to_smarts(&q);
        let reparsed = smarts(&written);
        assert_eq!(
            q.atom(NodeIndex::new(0)),
            reparsed.atom(NodeIndex::new(0)),
            "round-trip failed for {s}: wrote {written}"
        );
    }

    #[test]
    fn writer_disconnected() {
        let q = smarts("[Na].[Cl]");
        let written = to_smarts(&q);
        assert!(written.contains('.'));
    }

    #[test]
    fn writer_ring() {
        let q = smarts("C1CC1");
        let written = to_smarts(&q);
        assert!(written.contains('1'));
        let reparsed = smarts(&written);
        assert_eq!(reparsed.atom_count(), 3);
        assert_eq!(reparsed.bond_count(), 3);
    }

    #[test]
    fn match_bare_hydrogen() {
        let target = mol("C");
        let query = smarts("H");
        assert!(has_smarts_match(&target, &query));
        assert_eq!(get_smarts_matches(&target, &query).len(), 4);
    }

    #[test]
    fn match_bracket_na_charge() {
        let target = mol("[Na+]");
        let query = smarts("[Na+]");
        assert!(has_smarts_match(&target, &query));
        let expr = query.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts.iter().any(|p| matches!(
                    p,
                    AtomExpr::Element {
                        atomic_num: 11,
                        aromatic: Some(false)
                    }
                )));
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Charge(1))));
            }
            _ => panic!("expected And for [Na+], got {expr:?}"),
        }
    }

    #[test]
    fn match_bracket_cl_negative() {
        let target = mol("[Cl-]");
        let query = smarts("[Cl-]");
        assert!(has_smarts_match(&target, &query));
        let expr = query.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts.iter().any(|p| matches!(
                    p,
                    AtomExpr::Element {
                        atomic_num: 17,
                        aromatic: Some(false)
                    }
                )));
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Charge(-1))));
            }
            _ => panic!("expected And for [Cl-], got {expr:?}"),
        }
    }

    #[test]
    fn match_double_bond_ethene() {
        let target = mol("C=C");
        let query = smarts("C=C");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn no_match_double_bond_on_single() {
        let target = mol("CC");
        let query = smarts("C=C");
        assert!(!has_smarts_match(&target, &query));
    }

    #[test]
    fn parse_high_and_explicit() {
        let q = smarts("[c,n&H1]");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::Or(parts) => {
                assert_eq!(parts.len(), 2);
                assert!(matches!(
                    parts[0],
                    AtomExpr::Element {
                        atomic_num: 6,
                        aromatic: Some(true)
                    }
                ));
                match &parts[1] {
                    AtomExpr::And(inner) => {
                        assert_eq!(inner.len(), 2);
                        assert!(matches!(
                            inner[0],
                            AtomExpr::Element {
                                atomic_num: 7,
                                aromatic: Some(true)
                            }
                        ));
                        assert!(matches!(inner[1], AtomExpr::TotalHCount(1)));
                    }
                    _ => panic!("expected And in second Or part"),
                }
            }
            _ => panic!("expected Or, got {expr:?}"),
        }
    }

    #[test]
    fn parse_element_as_in_bracket() {
        let q = smarts("[As]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 33,
                aromatic: Some(false)
            }
        ));
    }

    #[test]
    fn parse_bare_bromine() {
        let q = smarts("Br");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 35,
                aromatic: Some(false)
            }
        ));
    }

    #[test]
    fn parse_bare_chlorine() {
        let q = smarts("Cl");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 17,
                aromatic: Some(false)
            }
        ));
    }

    #[test]
    fn match_ring_bond() {
        let target = mol("C1CCCCC1");
        let query = smarts("C@C");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn no_ring_bond_in_chain() {
        let target = mol("CCCC");
        let query = smarts("C@C");
        assert!(!has_smarts_match(&target, &query));
    }

    #[test]
    fn match_isotope() {
        let target = mol("[13C]");
        let query = smarts("[13]");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn unclosed_ring_error() {
        assert!(matches!(
            from_smarts("C1CC"),
            Err(SmartsError::UnclosedRing { .. })
        ));
    }

    #[test]
    fn unmatched_paren_error() {
        assert!(matches!(
            from_smarts("C(C"),
            Err(SmartsError::UnmatchedParen { .. })
        ));
    }

    #[test]
    fn bare_se_aromatic() {
        let q = smarts("[se]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(
            expr,
            AtomExpr::Element {
                atomic_num: 34,
                aromatic: Some(true)
            }
        ));
    }

    // ---- Non-H degree (d) ----

    #[test]
    fn parse_non_h_degree() {
        let q = smarts("[d2]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::NonHDegree(2)));
    }

    #[test]
    fn parse_non_h_degree_bare() {
        let q = smarts("[d]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::NonHDegree(1)));
    }

    #[test]
    fn match_non_h_degree_1_on_ccc() {
        let target = mol("CCC");
        let query = smarts("[d1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_non_h_degree_2_on_ccc() {
        let target = mol("CCC");
        let query = smarts("[d2]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_non_h_degree_3_branched() {
        let target = mol("CC(C)C");
        let query = smarts("[d3]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    // ---- Heteroatom neighbor count (z) ----

    #[test]
    fn parse_hetero_neighbor_count() {
        let q = smarts("[z1]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::HeteroNeighborCount(1)
        ));
    }

    #[test]
    fn parse_hetero_neighbor_bare() {
        let q = smarts("[z]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::HasHeteroNeighbor
        ));
    }

    #[test]
    fn match_z0_all_carbon() {
        let target = mol("CCC");
        let query = smarts("[z0]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 3);
    }

    #[test]
    fn match_z1_on_coc_c_n() {
        // COC(C)N: atom 0=C(bonded to O), 1=O(bonded to C,C), 2=C(bonded to O,C,N), 3=C(bonded to C), 4=N(bonded to C)
        let target = mol("COC(C)N");
        let query = smarts("[z1]");
        let matches = get_smarts_matches(&target, &query);
        // C(0): neighbor O → 1 heteroatom = yes
        // O(1): neighbors C,C → 0 heteroatoms = no
        // C(2): neighbors O,C,N → 2 heteroatoms = no
        // C(3): neighbor C → 0 heteroatoms = no
        // N(4): neighbor C → 0 heteroatoms = no
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_z2_on_coc_c_n() {
        let target = mol("COC(C)N");
        let query = smarts("[z2]");
        let matches = get_smarts_matches(&target, &query);
        // C(2) has O and N neighbors → 2 heteroatoms
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_bare_z_has_any_hetero() {
        let target = mol("COC(C)N");
        let query = smarts("[z]");
        let matches = get_smarts_matches(&target, &query);
        // C(0): 1 hetero neighbor (O) → yes
        // O(1): 0 hetero neighbors (C,C) → no
        // C(2): 2 hetero neighbors (O,N) → yes
        // C(3): 0 → no
        // N(4): 0 → no
        assert_eq!(matches.len(), 2);
    }

    // ---- Aliphatic heteroatom neighbor count (Z) ----

    #[test]
    fn parse_aliphatic_hetero_count() {
        let q = smarts("[Z1]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::AliphaticHeteroNeighborCount(1)
        ));
    }

    #[test]
    fn parse_aliphatic_hetero_bare() {
        let q = smarts("[Z]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::HasAliphaticHeteroNeighbor
        ));
    }

    #[test]
    fn match_big_z1_on_coc_c_n() {
        let target = mol("COC(C)N");
        let query = smarts("[Z1]");
        let matches = get_smarts_matches(&target, &query);
        // Same as z1 for all-aliphatic molecules
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn match_big_z0_aromatic_hetero() {
        // Pyridine: aromatic nitrogen — neighbors see aromatic N which doesn't count for Z
        let target = mol("c1ccncc1");
        let query = smarts("[Z0]");
        let matches = get_smarts_matches(&target, &query);
        // All atoms in pyridine: neighbors of aromatic carbons adjacent to N
        // see aromatic N, which is excluded from Z count.
        // c(0): neighbors c,c → Z=0
        // c(1): neighbors c,c → Z=0
        // c(2): neighbors c,n(aromatic) → Z=0 (n is aromatic, excluded)
        // n(3): neighbors c,c → Z=0
        // c(4): neighbors n(aromatic),c → Z=0
        // c(5): neighbors c,c → Z=0
        assert_eq!(matches.len(), 6);
    }

    // ---- Hybridization (^) ----

    #[test]
    fn parse_hybridization_sp3() {
        let q = smarts("[^3]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Hybridization(Hybridization::SP3)
        ));
    }

    #[test]
    fn parse_hybridization_sp2() {
        let q = smarts("[^2]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Hybridization(Hybridization::SP2)
        ));
    }

    #[test]
    fn parse_hybridization_sp() {
        let q = smarts("[^1]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Hybridization(Hybridization::SP)
        ));
    }

    #[test]
    fn match_sp3_ethane() {
        let target = mol("CC");
        let query = smarts("[^3]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_sp2_ethene() {
        let target = mol("C=C");
        let query = smarts("[^2]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_sp_ethyne() {
        let target = mol("C#C");
        let query = smarts("[^1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn match_sp2_benzene() {
        let target = mol("c1ccccc1");
        let query = smarts("[^2]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 6);
    }

    #[test]
    fn no_sp3_in_ethene() {
        let target = mol("C=C");
        let query = smarts("[^3]");
        assert!(!has_smarts_match(&target, &query));
    }

    // ---- Range syntax ({low-high}) ----

    #[test]
    fn parse_range_degree() {
        let q = smarts("[D{2-3}]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::Degree,
                low: Some(2),
                high: Some(3)
            }
        ));
    }

    #[test]
    fn parse_range_open_high() {
        let q = smarts("[D{2-}]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::Degree,
                low: Some(2),
                high: None
            }
        ));
    }

    #[test]
    fn parse_range_open_low() {
        let q = smarts("[D{-2}]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::Degree,
                low: None,
                high: Some(2)
            }
        ));
    }

    #[test]
    fn match_range_degree_2_3() {
        // COC(C)N: degrees: C(0)=1, O(1)=2, C(2)=3, C(3)=1, N(4)=1
        let target = mol("COC(C)N");
        let query = smarts("[D{2-3}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2); // O(D=2) and central C(D=3)
    }

    #[test]
    fn match_range_degree_2_up() {
        let target = mol("COC(C)N");
        let query = smarts("[D{2-}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2); // D>=2: O(2) and C(3)
    }

    #[test]
    fn match_range_degree_up_to_2() {
        let target = mol("COC(C)N");
        let query = smarts("[D{-2}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 4); // D<=2: C(1), O(2), C(1), N(1)
    }

    #[test]
    fn match_range_non_h_degree() {
        let target = mol("CCC");
        let query = smarts("[d{1-2}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 3); // all: d1, d2, d1
    }

    #[test]
    fn match_range_valence() {
        let target = mol("C=C");
        let query = smarts("[v{3-4}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2); // both have valence 4
    }

    #[test]
    fn match_range_hetero_neighbor() {
        let target = mol("COC(C)N");
        let query = smarts("[z{1-2}]");
        let matches = get_smarts_matches(&target, &query);
        // C(0): z=1, C(2): z=2
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn parse_range_total_h_count() {
        let q = smarts("[H{1-2}]");
        assert!(matches!(
            q.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::TotalHCount,
                low: Some(1),
                high: Some(2)
            }
        ));
    }

    #[test]
    fn match_range_total_h_count() {
        // methane CH4 has hydrogen_count=4, ethane C2H6 each C has hydrogen_count=3
        let target = mol("CC");
        let query = smarts("[H{2-3}]");
        let matches = get_smarts_matches(&target, &query);
        // Each carbon in ethane has 3 implicit H → both match H{2-3}
        assert_eq!(matches.len(), 2);

        let query_no_match = smarts("[H{0-0}]");
        let matches_none = get_smarts_matches(&target, &query_no_match);
        assert_eq!(matches_none.len(), 0);
    }

    #[test]
    fn range_h_vs_implicit_h_are_distinct() {
        let q_total = smarts("[H{1-2}]");
        let q_implicit = smarts("[h{1-2}]");
        assert!(matches!(
            q_total.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::TotalHCount,
                ..
            }
        ));
        assert!(matches!(
            q_implicit.atom(NodeIndex::new(0)),
            AtomExpr::Range {
                kind: RangeKind::ImplicitHCount,
                ..
            }
        ));
    }

    #[test]
    fn match_range_positive_charge() {
        let target = mol("[NH4+]");
        let query = smarts("[+{1-2}]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
    }

    // ---- Writer round-trips for new primitives ----

    #[test]
    fn writer_round_trip_non_h_degree() {
        let cases = ["[d1]", "[d2]", "[d3]"];
        for s in &cases {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    #[test]
    fn writer_round_trip_hetero_neighbor() {
        for s in &["[z0]", "[z1]", "[z2]", "[z]"] {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    #[test]
    fn writer_round_trip_aliphatic_hetero() {
        for s in &["[Z0]", "[Z1]", "[Z]"] {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    #[test]
    fn writer_round_trip_hybridization() {
        for s in &["[^0]", "[^1]", "[^2]", "[^3]", "[^4]", "[^5]"] {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    #[test]
    fn writer_round_trip_range() {
        for s in &[
            "[D{2-3}]", "[D{2-}]", "[D{-2}]", "[d{1-3}]", "[z{0-1}]", "[v{3-4}]", "[H{1-2}]",
            "[h{1-2}]",
        ] {
            let q = smarts(s);
            let written = to_smarts(&q);
            let reparsed = smarts(&written);
            assert_eq!(
                q.atom(NodeIndex::new(0)),
                reparsed.atom(NodeIndex::new(0)),
                "round-trip failed for {s}: wrote {written}"
            );
        }
    }

    // ---- Chirality parser tests ----

    #[test]
    fn parse_chirality_ccw() {
        let q = smarts("[C@](F)(Cl)Br");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(
                    parts
                        .iter()
                        .any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))),
                    "expected chirality variant, got {parts:?}"
                );
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    #[test]
    fn parse_chirality_cw() {
        let q = smarts("[C@@](F)(Cl)Br");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(
                    parts
                        .iter()
                        .any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))),
                    "expected chirality variant, got {parts:?}"
                );
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    #[test]
    fn parse_chirality_opposite_labels() {
        let q_at = smarts("[C@](F)(Cl)Br");
        let q_atat = smarts("[C@@](F)(Cl)Br");
        let get_chirality = |expr: &AtomExpr| -> Option<Chirality> {
            match expr {
                AtomExpr::And(parts) => parts.iter().find_map(|p| match p {
                    AtomExpr::Chirality(c) => Some(*c),
                    _ => None,
                }),
                _ => None,
            }
        };
        let c1 = get_chirality(q_at.atom(NodeIndex::new(0))).unwrap();
        let c2 = get_chirality(q_atat.atom(NodeIndex::new(0))).unwrap();
        assert_ne!(c1, c2, "@ and @@ must produce different chirality labels");
    }

    #[test]
    fn parse_chirality_with_hydrogen() {
        let q = smarts("[C@H](F)(Cl)Br");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts.iter().any(|p| matches!(
                    p,
                    AtomExpr::Element {
                        atomic_num: 6,
                        aromatic: Some(false)
                    }
                )));
                assert!(parts
                    .iter()
                    .any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))));
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::TotalHCount(1))));
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    #[test]
    fn parse_chirality_cw_with_hydrogen() {
        let q = smarts("[C@@H](F)(Cl)Br");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts
                    .iter()
                    .any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))));
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::TotalHCount(1))));
            }
            _ => panic!("expected And, got {expr:?}"),
        }
    }

    // ---- Chirality writer tests ----

    #[test]
    fn writer_chirality_ccw() {
        let q = smarts("[C@]");
        let written = to_smarts(&q);
        assert_eq!(written, "[C&@]");
        let reparsed = smarts(&written);
        assert_eq!(q.atom(NodeIndex::new(0)), reparsed.atom(NodeIndex::new(0)));
    }

    #[test]
    fn writer_chirality_cw() {
        let q = smarts("[C@@]");
        let written = to_smarts(&q);
        assert_eq!(written, "[C&@@]");
        let reparsed = smarts(&written);
        assert_eq!(q.atom(NodeIndex::new(0)), reparsed.atom(NodeIndex::new(0)));
    }

    #[test]
    fn writer_chirality_with_h_round_trip() {
        let q = smarts("[C@H]");
        let written = to_smarts(&q);
        let reparsed = smarts(&written);
        assert_eq!(q.atom(NodeIndex::new(0)), reparsed.atom(NodeIndex::new(0)));
    }

    // ---- Chirality matching tests ----

    #[test]
    fn chiral_query_matches_correct_enantiomer() {
        let target = mol("[C@](F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_query_rejects_wrong_enantiomer() {
        let target = mol("[C@@](F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(!has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn inverted_chirality_matches_opposite() {
        let target = mol("[C@@](F)(Cl)Br");
        let query = smarts("[C@@](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn achiral_query_matches_chiral_target() {
        let target = mol("[C@](F)(Cl)Br");
        let query = smarts("[C](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn achiral_query_matches_achiral_target() {
        let target = mol("C(F)(Cl)Br");
        let query = smarts("[C](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_query_rejects_achiral_target() {
        let target = mol("C(F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(!has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_query_few_neighbors_matches_any() {
        let target = mol("[C@@H](F)Cl");
        let query = smarts("[C@]F");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_with_implicit_h_matches() {
        let target = mol("[C@H](F)(Cl)Br");
        let query = smarts("[C@H](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_with_implicit_h_rejects_wrong() {
        let target = mol("[C@@H](F)(Cl)Br");
        let query = smarts("[C@H](F)(Cl)Br");
        assert!(!has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn non_chiral_matching_ignores_chirality() {
        let target = mol("[C@@](F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(has_smarts_match(&target, &query));
    }

    #[test]
    fn get_chiral_matches_returns_correct_count() {
        let target = mol("[C@](F)(Cl)Br");
        let query_match = smarts("[C@](F)(Cl)Br");
        let query_no_match = smarts("[C@@](F)(Cl)Br");
        let matches = get_smarts_matches_chiral(&target, &query_match);
        assert_eq!(matches.len(), 1);
        let no_matches = get_smarts_matches_chiral(&target, &query_no_match);
        assert_eq!(no_matches.len(), 0);
    }

    #[test]
    fn chiral_match_alanine_l() {
        let l_ala = mol("[C@@H](N)(C(=O)O)C");
        let query_l = smarts("[C@@H]([NH2])(C=O)C");
        assert!(has_smarts_match_chiral(&l_ala, &query_l));
    }

    #[test]
    fn chiral_match_alanine_d_vs_l() {
        let l_ala = mol("[C@@H](N)(C(=O)O)C");
        let query_d = smarts("[C@H]([NH2])(C=O)C");
        assert!(!has_smarts_match_chiral(&l_ala, &query_d));
    }

    #[test]
    fn chiral_ambiguous_neighbors_matches_either() {
        let l_ala = mol("[C@@H](N)(C(=O)O)C");
        let query = smarts("[C@H](N)(C)C");
        assert!(has_smarts_match_chiral(&l_ala, &query));
    }

    #[test]
    fn chiral_match_ring_atom() {
        let target = mol("[C@@H]1(C)CCCC[C@H]1C");
        let query = smarts("[C@@H](C)C");
        assert!(has_smarts_match_chiral(&target, &query));
    }

    #[test]
    fn chiral_match_two_centers() {
        let target = mol("[C@H](F)(Cl)[C@@](Br)(I)O");
        let query = smarts("[C@](F)(Cl)[C@@](Br)(I)O");
        assert!(has_smarts_match_chiral(&target, &query));

        let query_wrong = smarts("[C@@](F)(Cl)[C@@](Br)(I)O");
        assert!(!has_smarts_match_chiral(&target, &query_wrong));
    }

    #[test]
    fn debug_chirality_internals_simple() {
        let target = mol("[C@](F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));

        let query_opp = smarts("[C@@](F)(Cl)Br");
        assert!(!has_smarts_match_chiral(&target, &query_opp));
    }

    #[test]
    fn total_h_count_explicit_deuterium() {
        let target = mol("[2H]C([2H])([2H])[2H]");
        let query = smarts("[H0]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(
            matches.len(),
            4,
            "each deuterium has 0 total H, carbon has 4"
        );
        for m in &matches {
            let target_idx = m[0].1;
            assert_eq!(target.atom(target_idx).atomic_num, 1);
        }
    }

    #[test]
    fn total_h_count_virtual_only() {
        let target = mol("CCO");
        let query = smarts("[H1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1, "only oxygen has 1 total H");
        let target_idx = matches[0][0].1;
        assert_eq!(target.atom(target_idx).atomic_num, 8);
    }

    #[test]
    fn total_h_count_four_explicit() {
        let target = mol("[2H]C([2H])([2H])[2H]");
        let query = smarts("[H4]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1, "carbon has 4 total H (all explicit)");
        let target_idx = matches[0][0].1;
        assert_eq!(target.atom(target_idx).atomic_num, 6);
    }

    #[test]
    fn implicit_h_count_deuterium() {
        let target = mol("[2H]C([2H])([2H])[2H]");
        let query = smarts("[h0]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 5, "all atoms have 0 implicit H");
    }

    // ---- Compound bond expression parser tests ----

    #[test]
    fn parse_not_ring_bond() {
        let q = smarts("C!@C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(*q.bond(edge), BondExpr::Not(Box::new(BondExpr::Ring)));
    }

    #[test]
    fn parse_not_aromatic_bond() {
        let q = smarts("[#6]-!:[#6]");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(
            *q.bond(edge),
            BondExpr::And(vec![
                BondExpr::Single,
                BondExpr::Not(Box::new(BondExpr::Aromatic))
            ])
        );
    }

    #[test]
    fn parse_bond_or() {
        let q = smarts("C-,=C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(
            *q.bond(edge),
            BondExpr::Or(vec![BondExpr::Single, BondExpr::Double])
        );
    }

    #[test]
    fn parse_bond_and() {
        let q = smarts("C-&@C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(
            *q.bond(edge),
            BondExpr::And(vec![BondExpr::Single, BondExpr::Ring])
        );
    }

    #[test]
    fn parse_compound_bond_not_and() {
        let q = smarts("C-&!@C");
        let edge = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(
            *q.bond(edge),
            BondExpr::And(vec![
                BondExpr::Single,
                BondExpr::Not(Box::new(BondExpr::Ring))
            ])
        );
    }

    // ---- Compound bond expression match tests (ported from RDKit smatest.cpp) ----

    #[test]
    fn match_not_ring_bond() {
        let target = mol("C1CC1CC");
        let query = smarts("C!@C");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_not_ring_bond_chain() {
        let target = mol("CCCCC");
        let query = smarts("C!@C");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 4);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_not_ring_bond_hetero() {
        let target = mol("C1CC1CO");
        let query = smarts("[#6]!@[!#6]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_not_ring_bond_hetero_chain() {
        let target = mol("CCOCC");
        let query = smarts("[#6]!@[!#6]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 2);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_bond_or_double_single() {
        let target = mol("N=C");
        let query = smarts("N(=,-C)");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_bond_or_double_single_2() {
        let target = mol("NC");
        let query = smarts("N(=,-C)");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
        for m in &matches {
            assert_eq!(m.len(), 2);
        }
    }

    #[test]
    fn match_bond_or_no_triple() {
        let target = mol("N#C");
        let query = smarts("N(=,-C)");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn match_complex_not_ring() {
        let target = mol("N#CCC#N");
        let query = smarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn match_complex_not_ring_2() {
        let target = mol("OCCC#N");
        let query = smarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
        let matches = get_smarts_matches(&target, &query);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].len(), 2);
    }

    #[test]
    fn match_not_aromatic_bond() {
        let q = smarts("[#6]-!:[#6]");
        assert_eq!(q.atom_count(), 2);
        assert_eq!(q.bond_count(), 1);
        let target = mol("CC");
        assert!(has_smarts_match(&target, &q));
    }

    // ---- Compound bond expression writer round-trip tests ----

    #[test]
    fn writer_not_ring_bond() {
        let q = smarts("C!@C");
        let written = to_smarts(&q);
        assert!(written.contains("!@"), "expected !@ in {written}");
    }

    #[test]
    fn writer_bond_or() {
        let q = smarts("C-,=C");
        let written = to_smarts(&q);
        assert!(written.contains("-,="), "expected -,= in {written}");
    }

    #[test]
    fn writer_bond_and_not_round_trip() {
        let q = smarts("C-&!@C");
        let written = to_smarts(&q);
        assert!(written.contains("!@"), "expected !@ in {written}");
        let reparsed = smarts(&written);
        let edge_orig = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        let edge_re = reparsed
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(q.bond(edge_orig), reparsed.bond(edge_re));
    }

    #[test]
    fn writer_not_aromatic_round_trip() {
        let q = smarts("[#6]-!:[#6]");
        let written = to_smarts(&q);
        let reparsed = smarts(&written);
        let edge_orig = q
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        let edge_re = reparsed
            .bond_between(NodeIndex::new(0), NodeIndex::new(1))
            .unwrap();
        assert_eq!(q.bond(edge_orig), reparsed.bond(edge_re));
    }

    // ---- [#1] / hydrogen SMARTS matching ----

    #[test]
    fn hydrogen_by_atomic_num() {
        let methane = mol("C");
        let q = smarts("[#1]");
        assert_eq!(get_smarts_matches(&methane, &q).len(), 4);
    }

    #[test]
    fn oh_hydrogen_match() {
        let ethanol = mol("CCO");
        let q = smarts("[OX2:1][#1:2]");
        assert!(has_smarts_match(&ethanol, &q));
        let matches = get_smarts_matches(&ethanol, &q);
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn water_hydrogen_matches() {
        let water = mol("O");
        let q = smarts("[OX2:1][#1:2]");
        let matches = get_smarts_matches(&water, &q);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn ch_hydrogen_match() {
        let methane = mol("C");
        let q = smarts("[CX4:1][#1:2]");
        let matches = get_smarts_matches(&methane, &q);
        assert_eq!(matches.len(), 4);
    }

    #[test]
    fn no_hydrogen_query_unchanged() {
        let benzene = mol("c1ccccc1");
        let q = smarts("[cH]");
        assert_eq!(get_smarts_matches(&benzene, &q).len(), 6);
    }

    #[test]
    fn hydrogen_smarts_match_all() {
        let methane = mol("C");
        let q = smarts("[#1]");
        assert_eq!(get_smarts_matches_all(&methane, &q).len(), 4);
    }

    #[test]
    fn hydrogen_has_smarts_match() {
        let methane = mol("C");
        let q = smarts("[#1]");
        assert!(has_smarts_match(&methane, &q));
    }

    #[test]
    fn hydrogen_get_smarts_match() {
        let ethanol = mol("CCO");
        let q = smarts("[OX2][#1]");
        assert!(get_smarts_match(&ethanol, &q).is_some());
    }

    #[test]
    fn no_hydrogen_on_bare_metal() {
        let iron = mol("[Fe]");
        let q = smarts("[#1]");
        assert!(!has_smarts_match(&iron, &q));
    }

    #[test]
    fn bracket_h_matches_methane() {
        let methane = mol("C");
        let q = smarts("[H]");
        assert_eq!(get_smarts_matches(&methane, &q).len(), 4);
    }

    #[test]
    fn recursive_smarts_acetone_both_methyls() {
        let target = mol("CC(=O)C");
        let q = smarts("[$([CX4][CX3](=O))]");
        assert_eq!(get_smarts_matches(&target, &q).len(), 2);
    }

    #[test]
    fn recursive_smarts_mek_both_alphas() {
        let target = mol("CCC(=O)C");
        let q = smarts("[$([CX4][CX3](=O))]");
        assert_eq!(get_smarts_matches(&target, &q).len(), 2);
    }

    #[test]
    fn recursive_smarts_single_atom_unchanged() {
        let target = mol("c1ccccc1");
        let q = smarts("[$([cH1])]");
        assert_eq!(get_smarts_matches(&target, &q).len(), 6);
    }

    #[test]
    fn recursive_smarts_ethanol_oxygen() {
        let target = mol("CCO");
        let q = smarts("[$([OX2])]");
        assert_eq!(get_smarts_matches(&target, &q).len(), 1);
    }
}
