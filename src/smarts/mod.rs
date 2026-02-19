mod error;
mod parser;
pub mod query;
mod writer;

pub use error::SmartsError;
pub use query::{AtomExpr, BondExpr};
pub use writer::to_smarts;

use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::Bond;
use crate::mol::Mol;
use crate::rings::RingInfo;
use crate::atom::Chirality;
use crate::substruct::{
    get_substruct_match_with, get_substruct_match_with_filter, get_substruct_matches_with,
    get_substruct_matches_with_filter, AtomMapping,
};

use query::MatchContext;

type RecursiveRef<'a> = (*const Mol<AtomExpr, BondExpr>, &'a Mol<AtomExpr, BondExpr>);

pub fn from_smarts(s: &str) -> Result<Mol<AtomExpr, BondExpr>, SmartsError> {
    parser::parse(s)
}

pub fn has_smarts_match(target: &Mol<Atom, Bond>, query: &Mol<AtomExpr, BondExpr>) -> bool {
    get_smarts_match(target, query).is_some()
}

pub fn get_smarts_match(
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
        |t_atom: &Atom, q_expr: &AtomExpr| {
            let idx = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_atom))
                .unwrap();
            match_atom_expr(q_expr, t_atom, &ctx, idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints: (&Atom, &Atom), _q_endpoints| {
            let t_a = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.0))
                .unwrap();
            let t_b = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.1))
                .unwrap();
            q_bond.matches_with_ring_info(t_bond, t_endpoints, &ring_info, (t_a, t_b))
        },
    )
}

pub fn get_smarts_matches(
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
        |t_atom: &Atom, q_expr: &AtomExpr| {
            let idx = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_atom))
                .unwrap();
            match_atom_expr(q_expr, t_atom, &ctx, idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints: (&Atom, &Atom), _q_endpoints| {
            let t_a = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.0))
                .unwrap();
            let t_b = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.1))
                .unwrap();
            q_bond.matches_with_ring_info(t_bond, t_endpoints, &ring_info, (t_a, t_b))
        },
    )
}

pub fn has_smarts_match_chiral(target: &Mol<Atom, Bond>, query: &Mol<AtomExpr, BondExpr>) -> bool {
    get_smarts_match_chiral(target, query).is_some()
}

pub fn get_smarts_match_chiral(
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
        |t_atom: &Atom, q_expr: &AtomExpr| {
            let idx = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_atom))
                .unwrap();
            match_atom_expr(q_expr, t_atom, &ctx, idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints: (&Atom, &Atom), _q_endpoints| {
            let t_a = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.0))
                .unwrap();
            let t_b = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.1))
                .unwrap();
            q_bond.matches_with_ring_info(t_bond, t_endpoints, &ring_info, (t_a, t_b))
        },
        |mapping: &AtomMapping| validate_chirality(mapping, target, query, &chiral_query_atoms),
    )
}

pub fn get_smarts_matches_chiral(
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

    get_substruct_matches_with_filter(
        target,
        query,
        |t_atom: &Atom, q_expr: &AtomExpr| {
            let idx = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_atom))
                .unwrap();
            match_atom_expr(q_expr, t_atom, &ctx, idx)
        },
        |t_bond: &Bond, q_bond: &BondExpr, t_endpoints: (&Atom, &Atom), _q_endpoints| {
            let t_a = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.0))
                .unwrap();
            let t_b = target
                .atoms()
                .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.1))
                .unwrap();
            q_bond.matches_with_ring_info(t_bond, t_endpoints, &ring_info, (t_a, t_b))
        },
        |mapping: &AtomMapping| validate_chirality(mapping, target, query, &chiral_query_atoms),
    )
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
                    AtomExpr::And(inner) => {
                        if let Some((c, h)) = extract_chirality_and_h(&AtomExpr::And(inner.clone())) {
                            if c != Chirality::None {
                                chiral = Some(c);
                            }
                            if h {
                                has_h = true;
                            }
                        }
                    }
                    _ => {}
                }
            }
            chiral.map(|c| (c, has_h))
        }
        _ => None,
    }
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

        let t_atom = target.atom(t_idx);
        if t_atom.chirality == Chirality::None {
            return false;
        }

        let q_neighbor_count = query.neighbors(q_idx).count()
            + if cqa.has_implicit_h { 1 } else { 0 };
        if q_neighbor_count < 3 {
            continue;
        }

        let q_chiral = cqa.chirality;
        let t_chiral = t_atom.chirality;

        let mut mapped_target_neighbors: Vec<NodeIndex> = Vec::new();

        if cqa.has_implicit_h {
            let sentinel = NodeIndex::new(usize::MAX);
            mapped_target_neighbors.push(sentinel);
        }

        for q_nb in query.neighbors(q_idx) {
            let t_nb = match mapping.iter().find(|&&(q, _)| q == q_nb) {
                Some(&(_, t)) => t,
                None => return false,
            };
            mapped_target_neighbors.push(t_nb);
        }

        let target_stored_neighbors: Vec<NodeIndex> = {
            let mut nbs: Vec<NodeIndex> = target.neighbors(t_idx).collect();
            if t_atom.hydrogen_count > 0 && cqa.has_implicit_h {
                let sentinel = NodeIndex::new(usize::MAX);
                nbs.insert(0, sentinel);
            }
            nbs
        };

        let perm_count = count_inversions(&mapped_target_neighbors, &target_stored_neighbors);

        let parities_match = if perm_count.is_multiple_of(2) {
            q_chiral == t_chiral
        } else {
            q_chiral != t_chiral
        };

        if !parities_match {
            return false;
        }
    }
    true
}

fn count_inversions(mapped: &[NodeIndex], stored: &[NodeIndex]) -> usize {
    let position_of = |node: NodeIndex| -> Option<usize> {
        stored.iter().position(|&n| n == node)
    };

    let positions: Vec<usize> = mapped
        .iter()
        .filter_map(|&n| position_of(n))
        .collect();

    let mut inversions = 0;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            if positions[i] > positions[j] {
                inversions += 1;
            }
        }
    }
    inversions
}

fn match_atom_expr(expr: &AtomExpr, atom: &Atom, ctx: &MatchContext, idx: NodeIndex) -> bool {
    match expr {
        AtomExpr::Recursive(inner) => {
            let ptr = inner as *const Mol<AtomExpr, BondExpr> as usize;
            ctx.recursive_matches
                .get(&ptr)
                .is_some_and(|set| set.contains(&idx))
        }
        AtomExpr::And(exprs) => exprs
            .iter()
            .all(|e| match_atom_expr(e, atom, ctx, idx)),
        AtomExpr::Or(exprs) => exprs
            .iter()
            .any(|e| match_atom_expr(e, atom, ctx, idx)),
        AtomExpr::Chirality(_) => expr.matches(atom, ctx, idx),
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
            |t_atom: &Atom, q_expr: &AtomExpr| {
                let idx = target
                    .atoms()
                    .find(|&i| std::ptr::eq(target.atom(i), t_atom))
                    .unwrap();
                match_atom_expr(q_expr, t_atom, &inner_ctx, idx)
            },
            |t_bond: &Bond, q_bond: &BondExpr, t_endpoints: (&Atom, &Atom), _q_endpoints| {
                let t_a = target
                    .atoms()
                    .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.0))
                    .unwrap();
                let t_b = target
                    .atoms()
                    .find(|&i| std::ptr::eq(target.atom(i), t_endpoints.1))
                    .unwrap();
                q_bond.matches_with_ring_info(t_bond, t_endpoints, ring_info, (t_a, t_b))
            },
        );

        for mapping in &matches {
            if let Some(&(_q_first, t_first)) = mapping.first() {
                matching_atoms.insert(t_first);
            }
        }
        results.insert(key, matching_atoms);
    }

    results
}

fn collect_recursive_refs<'a>(
    expr: &'a AtomExpr,
    ptrs: &mut Vec<RecursiveRef<'a>>,
) {
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
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 6, aromatic: None }));
    }

    #[test]
    fn parse_aliphatic_carbon() {
        let q = smarts("[C]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 6, aromatic: Some(false) }));
    }

    #[test]
    fn parse_aromatic_carbon() {
        let q = smarts("[c]");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 6, aromatic: Some(true) }));
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
                assert!(matches!(parts[0], AtomExpr::Element { atomic_num: 6, aromatic: Some(false) }));
                assert!(matches!(parts[1], AtomExpr::Element { atomic_num: 7, aromatic: Some(false) }));
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
                assert!(matches!(parts[0], AtomExpr::Element { atomic_num: 6, aromatic: Some(false) }));
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
                assert!(matches!(**inner, AtomExpr::Element { atomic_num: 6, aromatic: Some(false) }));
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
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::Connectivity(4)));
    }

    #[test]
    fn parse_total_h_count() {
        let q = smarts("[H1]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::TotalHCount(1)));
    }

    #[test]
    fn parse_implicit_h_count() {
        let q = smarts("[h1]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::ImplicitHCount(1)));
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
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::RingMembership(2)));
    }

    #[test]
    fn parse_smallest_ring_size() {
        let q = smarts("[r6]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::SmallestRingSize(6)));
    }

    #[test]
    fn parse_ring_bond_count() {
        let q = smarts("[x2]");
        assert!(matches!(q.atom(NodeIndex::new(0)), AtomExpr::RingBondCount(2)));
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
        let edge = q.bond_between(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
        assert!(matches!(q.bond(edge), BondExpr::SingleOrAromatic));
    }

    #[test]
    fn parse_double_bond() {
        let q = smarts("C=C");
        let edge = q.bond_between(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
        assert!(matches!(q.bond(edge), BondExpr::Double));
    }

    #[test]
    fn parse_triple_bond() {
        let q = smarts("C#C");
        let edge = q.bond_between(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
        assert!(matches!(q.bond(edge), BondExpr::Triple));
    }

    #[test]
    fn parse_any_bond() {
        let q = smarts("C~C");
        let edge = q.bond_between(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
        assert!(matches!(q.bond(edge), BondExpr::True));
    }

    #[test]
    fn parse_aromatic_bond() {
        let q = smarts("C:C");
        let edge = q.bond_between(NodeIndex::new(0), NodeIndex::new(1)).unwrap();
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
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 1, aromatic: Some(false) }));
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
        let cases = ["[#6]", "[D2]", "[v4]", "[X4]", "[H1]", "[h1]", "[R]", "[R0]", "[R2]", "[r6]", "[x2]", "[+1]", "[-1]"];
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
        assert!(!has_smarts_match(&target, &query));
    }

    #[test]
    fn match_bracket_na_charge() {
        let target = mol("[Na+]");
        let query = smarts("[Na+]");
        assert!(has_smarts_match(&target, &query));
        let expr = query.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts
                    .iter()
                    .any(|p| matches!(p, AtomExpr::Element { atomic_num: 11, aromatic: Some(false) })));
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
                assert!(parts
                    .iter()
                    .any(|p| matches!(p, AtomExpr::Element { atomic_num: 17, aromatic: Some(false) })));
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
                assert!(matches!(parts[0], AtomExpr::Element { atomic_num: 6, aromatic: Some(true) }));
                match &parts[1] {
                    AtomExpr::And(inner) => {
                        assert_eq!(inner.len(), 2);
                        assert!(matches!(inner[0], AtomExpr::Element { atomic_num: 7, aromatic: Some(true) }));
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
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 33, aromatic: Some(false) }));
    }

    #[test]
    fn parse_bare_bromine() {
        let q = smarts("Br");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 35, aromatic: Some(false) }));
    }

    #[test]
    fn parse_bare_chlorine() {
        let q = smarts("Cl");
        let expr = q.atom(NodeIndex::new(0));
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 17, aromatic: Some(false) }));
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
        assert!(matches!(expr, AtomExpr::Element { atomic_num: 34, aromatic: Some(true) }));
    }

    // ---- Chirality parser tests ----

    #[test]
    fn parse_chirality_ccw() {
        let q = smarts("[C@](F)(Cl)Br");
        let expr = q.atom(NodeIndex::new(0));
        match expr {
            AtomExpr::And(parts) => {
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))),
                    "expected chirality variant, got {parts:?}");
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
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))),
                    "expected chirality variant, got {parts:?}");
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
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Element { atomic_num: 6, aromatic: Some(false) })));
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))));
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
                assert!(parts.iter().any(|p| matches!(p, AtomExpr::Chirality(Chirality::Cw | Chirality::Ccw))));
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
    fn debug_chirality_internals_simple() {
        let target = mol("[C@](F)(Cl)Br");
        let query = smarts("[C@](F)(Cl)Br");
        assert!(has_smarts_match_chiral(&target, &query));

        let query_opp = smarts("[C@@](F)(Cl)Br");
        assert!(!has_smarts_match_chiral(&target, &query_opp));
    }
}
