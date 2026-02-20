//! Substructure matching via the VF2 algorithm.
//!
//! Compares atoms by element, aromaticity, and formal charge. Compares bonds
//! by bond order (with aromatic bonds handled specially). For richer queries
//! involving hydrogen counts, ring membership, or degree constraints, use the
//! [`crate::smarts`] module instead.

use std::collections::HashSet;

use petgraph::graph::NodeIndex;

use crate::mol::Mol;
use crate::traits::{HasAromaticity, HasAtomicNum, HasBondOrder, HasFormalCharge};

/// A mapping from query atom indices to target atom indices.
pub type AtomMapping = Vec<(NodeIndex, NodeIndex)>;

/// Deduplicates atom mappings so each unique set of target atoms appears once.
///
/// Two mappings are considered duplicates if they map to the same set of
/// target atom indices, regardless of the query-side ordering.
pub fn uniquify_atom_mappings(mappings: &[AtomMapping]) -> Vec<AtomMapping> {
    let mut seen = HashSet::new();
    mappings
        .iter()
        .filter(|mapping| {
            let mut target_atoms: Vec<usize> = mapping.iter().map(|(_, t)| t.index()).collect();
            target_atoms.sort();
            seen.insert(target_atoms)
        })
        .cloned()
        .collect()
}

/// Returns `true` if `target` contains a substructure matching `query`.
pub fn has_substruct_match<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> bool
where
    A: HasAtomicNum + HasAromaticity + HasFormalCharge,
    B: HasBondOrder,
{
    get_substruct_match(target, query).is_some()
}

/// Returns the first substructure match of `query` in `target`, if any.
pub fn get_substruct_match<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> Option<AtomMapping>
where
    A: HasAtomicNum + HasAromaticity + HasFormalCharge,
    B: HasBondOrder,
{
    Vf2::new(target, query, default_atom_match, default_bond_match).find_first()
}

/// Returns all substructure matches of `query` in `target`.
///
/// This includes automorphic permutations. Use [`get_substruct_matches_unique`]
/// to deduplicate.
pub fn get_substruct_matches<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> Vec<AtomMapping>
where
    A: HasAtomicNum + HasAromaticity + HasFormalCharge,
    B: HasBondOrder,
{
    Vf2::new(target, query, default_atom_match, default_bond_match).find_all()
}

/// Returns all unique substructure matches, deduplicated by target atom set.
pub fn get_substruct_matches_unique<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> Vec<AtomMapping>
where
    A: HasAtomicNum + HasAromaticity + HasFormalCharge,
    B: HasBondOrder,
{
    uniquify_atom_mappings(&get_substruct_matches(target, query))
}

/// Like [`get_substruct_matches_unique`] but with custom match functions.
pub fn get_substruct_matches_with_unique<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
) -> Vec<AtomMapping> {
    uniquify_atom_mappings(&get_substruct_matches_with(
        target, query, atom_match, bond_match,
    ))
}

/// Like [`has_substruct_match`] but with custom atom and bond match functions.
pub fn has_substruct_match_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
) -> bool {
    get_substruct_match_with(target, query, atom_match, bond_match).is_some()
}

/// Like [`get_substruct_match`] but with custom atom and bond match functions.
pub fn get_substruct_match_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
) -> Option<AtomMapping> {
    Vf2::new(target, query, atom_match, bond_match).find_first()
}

/// Like [`get_substruct_matches`] but with custom atom and bond match functions.
pub fn get_substruct_matches_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
) -> Vec<AtomMapping> {
    Vf2::new(target, query, atom_match, bond_match).find_all()
}

/// Like [`get_substruct_match_with`] with an additional post-match filter.
pub fn get_substruct_match_with_filter<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
    filter: impl Fn(&AtomMapping) -> bool,
) -> Option<AtomMapping> {
    Vf2WithFilter::new(target, query, atom_match, bond_match, filter).find_first()
}

/// Like [`get_substruct_matches_with`] with an additional post-match filter.
pub fn get_substruct_matches_with_filter<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
    filter: impl Fn(&AtomMapping) -> bool,
) -> Vec<AtomMapping> {
    Vf2WithFilter::new(target, query, atom_match, bond_match, filter).find_all()
}

fn default_atom_match<A: HasAtomicNum + HasAromaticity + HasFormalCharge>(
    target: &A,
    query: &A,
) -> bool {
    if target.atomic_num() != query.atomic_num() {
        return false;
    }
    if query.is_aromatic() && !target.is_aromatic() {
        return false;
    }
    if target.formal_charge() != query.formal_charge() {
        return false;
    }
    true
}

fn default_bond_match<A: HasAromaticity, B: HasBondOrder>(
    target: &B,
    query: &B,
    target_endpoints: (&A, &A),
    query_endpoints: (&A, &A),
) -> bool {
    let both_target_aromatic = target_endpoints.0.is_aromatic() && target_endpoints.1.is_aromatic();
    let both_query_aromatic = query_endpoints.0.is_aromatic() && query_endpoints.1.is_aromatic();
    if both_query_aromatic && both_target_aromatic {
        return true;
    }
    target.bond_order() == query.bond_order()
}

struct Vf2<'a, A1, B1, A2, B2, FA, FB> {
    target: &'a Mol<A1, B1>,
    query: &'a Mol<A2, B2>,
    atom_match: FA,
    bond_match: FB,
    query_order: Vec<NodeIndex>,
    query_map: Vec<Option<NodeIndex>>,
    target_used: Vec<bool>,
}

impl<'a, A1, B1, A2, B2, FA, FB> Vf2<'a, A1, B1, A2, B2, FA, FB>
where
    FA: Fn(&A1, &A2) -> bool,
    FB: Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
{
    fn new(
        target: &'a Mol<A1, B1>,
        query: &'a Mol<A2, B2>,
        atom_match: FA,
        bond_match: FB,
    ) -> Self {
        let mut query_order: Vec<NodeIndex> = query.atoms().collect();
        query_order.sort_by(|&a, &b| query.neighbors(b).count().cmp(&query.neighbors(a).count()));
        Self {
            target,
            query,
            atom_match,
            bond_match,
            query_order,
            query_map: vec![None; query.atom_count()],
            target_used: vec![false; target.atom_count()],
        }
    }

    fn find_first(&mut self) -> Option<AtomMapping> {
        let mut results = Vec::new();
        self.recurse(0, &mut results, true);
        results.into_iter().next()
    }

    fn find_all(&mut self) -> Vec<AtomMapping> {
        let mut results = Vec::new();
        self.recurse(0, &mut results, false);
        results
    }

    fn recurse(&mut self, depth: usize, results: &mut Vec<AtomMapping>, first_only: bool) {
        if depth == self.query_order.len() {
            let mapping = self
                .query_order
                .iter()
                .map(|&qn| (qn, self.query_map[qn.index()].unwrap()))
                .collect();
            results.push(mapping);
            return;
        }

        if first_only && !results.is_empty() {
            return;
        }

        let query_node = self.query_order[depth];

        for t_idx in 0..self.target_used.len() {
            if self.target_used[t_idx] {
                continue;
            }

            let target_node = NodeIndex::new(t_idx);

            if !self.is_feasible(query_node, target_node) {
                continue;
            }

            self.query_map[query_node.index()] = Some(target_node);
            self.target_used[t_idx] = true;

            self.recurse(depth + 1, results, first_only);

            if first_only && !results.is_empty() {
                return;
            }

            self.query_map[query_node.index()] = None;
            self.target_used[t_idx] = false;
        }
    }

    fn is_feasible(&self, query_node: NodeIndex, target_node: NodeIndex) -> bool {
        if !(self.atom_match)(self.target.atom(target_node), self.query.atom(query_node)) {
            return false;
        }

        for q_neighbor in self.query.neighbors(query_node) {
            if let Some(t_mapped) = self.query_map[q_neighbor.index()] {
                let q_bond = self
                    .query
                    .bond_between(query_node, q_neighbor)
                    .expect("bond must exist between neighbors");
                match self.target.bond_between(target_node, t_mapped) {
                    Some(t_bond) => {
                        let target_endpoints =
                            (self.target.atom(target_node), self.target.atom(t_mapped));
                        let query_endpoints =
                            (self.query.atom(query_node), self.query.atom(q_neighbor));
                        if !(self.bond_match)(
                            self.target.bond(t_bond),
                            self.query.bond(q_bond),
                            target_endpoints,
                            query_endpoints,
                        ) {
                            return false;
                        }
                    }
                    None => return false,
                }
            }
        }

        true
    }
}

struct Vf2WithFilter<'a, A1, B1, A2, B2, FA, FB, FF> {
    target: &'a Mol<A1, B1>,
    query: &'a Mol<A2, B2>,
    atom_match: FA,
    bond_match: FB,
    filter: FF,
    query_order: Vec<NodeIndex>,
    query_map: Vec<Option<NodeIndex>>,
    target_used: Vec<bool>,
}

impl<'a, A1, B1, A2, B2, FA, FB, FF> Vf2WithFilter<'a, A1, B1, A2, B2, FA, FB, FF>
where
    FA: Fn(&A1, &A2) -> bool,
    FB: Fn(&B1, &B2, (&A1, &A1), (&A2, &A2)) -> bool,
    FF: Fn(&AtomMapping) -> bool,
{
    fn new(
        target: &'a Mol<A1, B1>,
        query: &'a Mol<A2, B2>,
        atom_match: FA,
        bond_match: FB,
        filter: FF,
    ) -> Self {
        let mut query_order: Vec<NodeIndex> = query.atoms().collect();
        query_order.sort_by(|&a, &b| query.neighbors(b).count().cmp(&query.neighbors(a).count()));
        Self {
            target,
            query,
            atom_match,
            bond_match,
            filter,
            query_order,
            query_map: vec![None; query.atom_count()],
            target_used: vec![false; target.atom_count()],
        }
    }

    fn find_first(&mut self) -> Option<AtomMapping> {
        let mut results = Vec::new();
        self.recurse(0, &mut results, true);
        results.into_iter().next()
    }

    fn find_all(&mut self) -> Vec<AtomMapping> {
        let mut results = Vec::new();
        self.recurse(0, &mut results, false);
        results
    }

    fn recurse(&mut self, depth: usize, results: &mut Vec<AtomMapping>, first_only: bool) {
        if depth == self.query_order.len() {
            let mapping: AtomMapping = self
                .query_order
                .iter()
                .map(|&qn| (qn, self.query_map[qn.index()].unwrap()))
                .collect();
            if (self.filter)(&mapping) {
                results.push(mapping);
            }
            return;
        }

        if first_only && !results.is_empty() {
            return;
        }

        let query_node = self.query_order[depth];

        for t_idx in 0..self.target_used.len() {
            if self.target_used[t_idx] {
                continue;
            }

            let target_node = NodeIndex::new(t_idx);

            if !self.is_feasible(query_node, target_node) {
                continue;
            }

            self.query_map[query_node.index()] = Some(target_node);
            self.target_used[t_idx] = true;

            self.recurse(depth + 1, results, first_only);

            if first_only && !results.is_empty() {
                return;
            }

            self.query_map[query_node.index()] = None;
            self.target_used[t_idx] = false;
        }
    }

    fn is_feasible(&self, query_node: NodeIndex, target_node: NodeIndex) -> bool {
        if !(self.atom_match)(self.target.atom(target_node), self.query.atom(query_node)) {
            return false;
        }

        for q_neighbor in self.query.neighbors(query_node) {
            if let Some(t_mapped) = self.query_map[q_neighbor.index()] {
                let q_bond = self
                    .query
                    .bond_between(query_node, q_neighbor)
                    .expect("bond must exist between neighbors");
                match self.target.bond_between(target_node, t_mapped) {
                    Some(t_bond) => {
                        let target_endpoints =
                            (self.target.atom(target_node), self.target.atom(t_mapped));
                        let query_endpoints =
                            (self.query.atom(query_node), self.query.atom(q_neighbor));
                        if !(self.bond_match)(
                            self.target.bond(t_bond),
                            self.query.bond(q_bond),
                            target_endpoints,
                            query_endpoints,
                        ) {
                            return false;
                        }
                    }
                    None => return false,
                }
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::from_smiles;

    fn mol(smiles: &str) -> Mol<crate::Atom, crate::Bond> {
        from_smiles(smiles).unwrap_or_else(|e| panic!("bad SMILES {smiles:?}: {e}"))
    }

    #[test]
    fn ethanol_contains_cc() {
        let target = mol("CCO");
        let query = mol("CC");
        assert!(has_substruct_match(&target, &query));
        let m = get_substruct_match(&target, &query).unwrap();
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn methane_does_not_contain_cc() {
        let target = mol("C");
        let query = mol("CC");
        assert!(!has_substruct_match(&target, &query));
        assert_eq!(get_substruct_match(&target, &query), None);
        assert!(get_substruct_matches(&target, &query).is_empty());
    }

    #[test]
    fn propane_cc_matches() {
        let target = mol("CCC");
        let query = mol("CC");
        let matches = get_substruct_matches(&target, &query);
        assert_eq!(matches.len(), 4);
    }

    #[test]
    fn cyclohexane_cc_single_bond_matches() {
        let target = mol("C1CCCCC1");
        let query = mol("CC");
        let matches = get_substruct_matches(&target, &query);
        assert_eq!(matches.len(), 12);
    }

    #[test]
    fn benzene_automorphisms() {
        let target = mol("c1ccccc1");
        let query = mol("c1ccccc1");
        let matches = get_substruct_matches(&target, &query);
        assert_eq!(matches.len(), 12);
    }

    #[test]
    fn empty_query_matches_anything() {
        let target = mol("CCO");
        let query = Mol::<crate::Atom, crate::Bond>::new();
        assert!(has_substruct_match(&target, &query));
        let m = get_substruct_match(&target, &query).unwrap();
        assert!(m.is_empty());
        let all = get_substruct_matches(&target, &query);
        assert_eq!(all.len(), 1);
        assert!(all[0].is_empty());
    }

    #[test]
    fn single_atom_query() {
        let target = mol("CCO");
        let query = mol("O");
        let matches = get_substruct_matches(&target, &query);
        assert_eq!(matches.len(), 1);
        let (q, t) = matches[0][0];
        assert_eq!(q, NodeIndex::new(0));
        assert_eq!(crate::traits::HasAtomicNum::atomic_num(target.atom(t)), 8);
    }

    #[test]
    fn query_larger_than_target_no_match() {
        let target = mol("C");
        let query = mol("CCCCCC");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn self_match() {
        let target = mol("CCO");
        let query = mol("CCO");
        assert!(has_substruct_match(&target, &query));
        let m = get_substruct_match(&target, &query).unwrap();
        assert_eq!(m.len(), 3);
    }

    #[test]
    fn bond_order_double_matches_double() {
        let target = mol("C=C");
        let query = mol("C=C");
        assert!(has_substruct_match(&target, &query));
    }

    #[test]
    fn bond_order_double_does_not_match_single() {
        let target = mol("CC");
        let query = mol("C=C");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn bond_order_single_does_not_match_double() {
        let target = mol("C=C");
        let query = mol("CC");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn aromatic_query_does_not_match_non_aromatic() {
        let target = mol("C1CCCCC1");
        let query = mol("c1ccccc1");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn aromatic_ring_in_naphthalene() {
        let target = mol("c1ccc2ccccc2c1");
        let query = mol("c1ccccc1");
        assert!(has_substruct_match(&target, &query));
        let matches = get_substruct_matches(&target, &query);
        assert!(!matches.is_empty());
    }

    #[test]
    fn phenol_oh_on_aromatic_ring() {
        let target = mol("c1ccc(O)cc1");
        let query = mol("OC");
        let matches = get_substruct_matches_with(
            &target,
            &query,
            |t: &crate::Atom, q: &crate::Atom| t.atomic_num == q.atomic_num,
            |t: &crate::Bond, _q, _, _| t.order == crate::BondOrder::Single,
        );
        assert!(!matches.is_empty());
    }

    #[test]
    fn custom_matchers_ignore_bond_order() {
        let target = mol("C=C");
        let query = mol("CC");
        let matches = get_substruct_matches_with(
            &target,
            &query,
            |t: &crate::Atom, q: &crate::Atom| t.atomic_num == q.atomic_num,
            |_t: &crate::Bond, _q, _, _| true,
        );
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn disconnected_target() {
        let target = mol("[Na+].[Cl-]");
        let query = mol("[Na+]");
        let matches = get_substruct_matches_with(
            &target,
            &query,
            |t: &crate::Atom, q: &crate::Atom| t.atomic_num == q.atomic_num,
            |_t: &crate::Bond, _q, _, _| true,
        );
        assert_eq!(matches.len(), 1);
    }

    #[test]
    fn triple_bond_match() {
        let target = mol("C#N");
        let query = mol("C#N");
        assert!(has_substruct_match(&target, &query));
    }

    #[test]
    fn triple_bond_does_not_match_single() {
        let target = mol("CN");
        let query = mol("C#N");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn empty_target_no_match_nonempty_query() {
        let target = Mol::<crate::Atom, crate::Bond>::new();
        let query = mol("C");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn both_empty() {
        let target = Mol::<crate::Atom, crate::Bond>::new();
        let query = Mol::<crate::Atom, crate::Bond>::new();
        assert!(has_substruct_match(&target, &query));
    }

    #[test]
    fn mapping_correctness() {
        let target = mol("CCO");
        let query = mol("CO");
        let m = get_substruct_match(&target, &query).unwrap();
        assert_eq!(m.len(), 2);
        for &(q, t) in &m {
            let q_num = crate::traits::HasAtomicNum::atomic_num(query.atom(q));
            let t_num = crate::traits::HasAtomicNum::atomic_num(target.atom(t));
            assert_eq!(q_num, t_num);
        }
    }

    #[test]
    fn all_mappings_are_valid() {
        let target = mol("c1ccccc1");
        let query = mol("c1ccccc1");
        let matches = get_substruct_matches(&target, &query);
        for mapping in &matches {
            assert_eq!(mapping.len(), query.atom_count());
            for &(q, t) in mapping {
                for q_neighbor in query.neighbors(q) {
                    let t_mapped = mapping
                        .iter()
                        .find(|&&(qn, _)| qn == q_neighbor)
                        .map(|&(_, tn)| tn)
                        .unwrap();
                    assert!(
                        target.bond_between(t, t_mapped).is_some(),
                        "mapped neighbors must be connected in target"
                    );
                }
            }
        }
    }

    #[test]
    fn no_duplicate_mappings() {
        let target = mol("c1ccccc1");
        let query = mol("c1ccccc1");
        let matches = get_substruct_matches(&target, &query);
        for (i, a) in matches.iter().enumerate() {
            for b in matches.iter().skip(i + 1) {
                assert_ne!(a, b, "duplicate mapping found");
            }
        }
    }

    #[test]
    fn cyclohexane_cc_matches_cover_all_edges() {
        let target = mol("C1CCCCC1");
        let query = mol("CC");
        let matches = get_substruct_matches(&target, &query);
        assert_eq!(matches.len(), 12);
        for mapping in &matches {
            assert_eq!(mapping.len(), 2);
            let t0 = mapping
                .iter()
                .find(|&&(q, _)| q == NodeIndex::new(0))
                .unwrap()
                .1;
            let t1 = mapping
                .iter()
                .find(|&&(q, _)| q == NodeIndex::new(1))
                .unwrap()
                .1;
            assert!(
                target.bond_between(t0, t1).is_some(),
                "matched atoms must be bonded"
            );
        }
    }

    #[test]
    fn benzene_on_benzene_raw_gives_12() {
        let target = mol("c1ccccc1");
        let query = mol("c1ccccc1");
        assert_eq!(get_substruct_matches(&target, &query).len(), 12);
    }

    #[test]
    fn benzene_on_benzene_unique_gives_1() {
        let target = mol("c1ccccc1");
        let query = mol("c1ccccc1");
        assert_eq!(get_substruct_matches_unique(&target, &query).len(), 1);
    }

    #[test]
    fn naphthalene_benzene_unique_gives_2() {
        let target = mol("c1ccc2ccccc2c1");
        let query = mol("c1ccccc1");
        assert_eq!(get_substruct_matches_unique(&target, &query).len(), 2);
    }

    #[test]
    fn charged_oxygen_does_not_match_neutral() {
        let target = mol("CCO");
        let query = mol("[O-]");
        assert!(!has_substruct_match(&target, &query));
    }

    #[test]
    fn nitro_group_matches_nitrobenzene() {
        let target = mol("[O-][N+](=O)c1ccccc1");
        let query = mol("[N+](=O)[O-]");
        assert!(has_substruct_match(&target, &query));
    }

    #[test]
    fn neutral_oxygen_does_not_match_charged_query() {
        let target = mol("CCO");
        let query = mol("[O-]");
        assert!(get_substruct_match(&target, &query).is_none());
    }
}
