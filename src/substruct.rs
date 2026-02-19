use petgraph::graph::NodeIndex;

use crate::mol::Mol;
use crate::traits::{HasAromaticity, HasAtomicNum, HasBondOrder};

pub type AtomMapping = Vec<(NodeIndex, NodeIndex)>;

pub fn has_substruct_match<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> bool
where
    A: HasAtomicNum + HasAromaticity,
    B: HasBondOrder,
{
    get_substruct_match(target, query).is_some()
}

pub fn get_substruct_match<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> Option<AtomMapping>
where
    A: HasAtomicNum + HasAromaticity,
    B: HasBondOrder,
{
    DefaultVf2::new(target, query).find_first()
}

pub fn get_substruct_matches<A, B>(target: &Mol<A, B>, query: &Mol<A, B>) -> Vec<AtomMapping>
where
    A: HasAtomicNum + HasAromaticity,
    B: HasBondOrder,
{
    DefaultVf2::new(target, query).find_all()
}

pub fn has_substruct_match_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2) -> bool,
) -> bool {
    get_substruct_match_with(target, query, atom_match, bond_match).is_some()
}

pub fn get_substruct_match_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2) -> bool,
) -> Option<AtomMapping> {
    Vf2::new(target, query, atom_match, bond_match).find_first()
}

pub fn get_substruct_matches_with<A1, B1, A2, B2>(
    target: &Mol<A1, B1>,
    query: &Mol<A2, B2>,
    atom_match: impl Fn(&A1, &A2) -> bool,
    bond_match: impl Fn(&B1, &B2) -> bool,
) -> Vec<AtomMapping> {
    Vf2::new(target, query, atom_match, bond_match).find_all()
}

fn default_atom_match<A: HasAtomicNum + HasAromaticity>(target: &A, query: &A) -> bool {
    if target.atomic_num() != query.atomic_num() {
        return false;
    }
    if query.is_aromatic() && !target.is_aromatic() {
        return false;
    }
    true
}

struct DefaultVf2<'a, A, B> {
    target: &'a Mol<A, B>,
    query: &'a Mol<A, B>,
    query_order: Vec<NodeIndex>,
    query_map: Vec<Option<NodeIndex>>,
    target_used: Vec<bool>,
}

impl<'a, A: HasAtomicNum + HasAromaticity, B: HasBondOrder> DefaultVf2<'a, A, B> {
    fn new(target: &'a Mol<A, B>, query: &'a Mol<A, B>) -> Self {
        let mut query_order: Vec<NodeIndex> = query.atoms().collect();
        query_order.sort_by(|&a, &b| {
            query.neighbors(b).count().cmp(&query.neighbors(a).count())
        });
        Self {
            target,
            query,
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

    fn recurse(
        &mut self,
        depth: usize,
        results: &mut Vec<AtomMapping>,
        first_only: bool,
    ) {
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
        if !default_atom_match(self.target.atom(target_node), self.query.atom(query_node)) {
            return false;
        }

        for q_neighbor in self.query.neighbors(query_node) {
            if let Some(t_mapped) = self.query_map[q_neighbor.index()] {
                let q_bond_idx = self
                    .query
                    .bond_between(query_node, q_neighbor)
                    .expect("bond must exist between neighbors");
                match self.target.bond_between(target_node, t_mapped) {
                    Some(t_bond_idx) => {
                        let both_query_aromatic = self.query.atom(query_node).is_aromatic()
                            && self.query.atom(q_neighbor).is_aromatic();
                        let both_target_aromatic = self.target.atom(target_node).is_aromatic()
                            && self.target.atom(t_mapped).is_aromatic();
                        if both_query_aromatic && both_target_aromatic {
                            continue;
                        }
                        if self.target.bond(t_bond_idx).bond_order()
                            != self.query.bond(q_bond_idx).bond_order()
                        {
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
    FB: Fn(&B1, &B2) -> bool,
{
    fn new(
        target: &'a Mol<A1, B1>,
        query: &'a Mol<A2, B2>,
        atom_match: FA,
        bond_match: FB,
    ) -> Self {
        let mut query_order: Vec<NodeIndex> = query.atoms().collect();
        query_order.sort_by(|&a, &b| {
            query.neighbors(b).count().cmp(&query.neighbors(a).count())
        });
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

    fn recurse(
        &mut self,
        depth: usize,
        results: &mut Vec<AtomMapping>,
        first_only: bool,
    ) {
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
                        if !(self.bond_match)(self.target.bond(t_bond), self.query.bond(q_bond)) {
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
            |t: &crate::Bond, _q: &crate::Bond| t.order == crate::BondOrder::Single,
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
            |_t: &crate::Bond, _q: &crate::Bond| true,
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
            |_t: &crate::Bond, _q: &crate::Bond| true,
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
}
