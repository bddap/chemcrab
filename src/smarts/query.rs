use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};
use crate::mol::Mol;
use crate::rings::RingInfo;

#[derive(Debug, Clone, PartialEq)]
pub enum AtomExpr {
    True,
    Element { atomic_num: u8, aromatic: Option<bool> },
    Aromatic,
    Aliphatic,
    Isotope(u16),
    Degree(u8),
    Valence(u8),
    Connectivity(u8),
    TotalHCount(u8),
    ImplicitHCount(u8),
    RingMembership(u8),
    SmallestRingSize(u8),
    RingBondCount(u8),
    Charge(i8),
    InRing,
    NotInRing,
    Recursive(Mol<AtomExpr, BondExpr>),
    And(Vec<AtomExpr>),
    Or(Vec<AtomExpr>),
    Chirality(crate::atom::Chirality),
    Not(Box<AtomExpr>),
}

#[derive(Debug, Clone, PartialEq)]
pub enum BondExpr {
    True,
    Single,
    Double,
    Triple,
    Aromatic,
    Ring,
    SingleOrAromatic,
    Up,
    Down,
    And(Vec<BondExpr>),
    Or(Vec<BondExpr>),
    Not(Box<BondExpr>),
}

pub struct MatchContext<'a> {
    pub mol: &'a Mol<Atom, Bond>,
    pub ring_info: &'a RingInfo,
    pub recursive_matches: HashMap<usize, HashSet<NodeIndex>>,
}

impl AtomExpr {
    pub fn matches(&self, atom: &Atom, ctx: &MatchContext, idx: NodeIndex) -> bool {
        match self {
            AtomExpr::True => true,
            AtomExpr::Element { atomic_num, aromatic } => {
                atom.atomic_num == *atomic_num
                    && aromatic.is_none_or(|a| atom.is_aromatic == a)
            }
            AtomExpr::Aromatic => atom.is_aromatic,
            AtomExpr::Aliphatic => !atom.is_aromatic,
            AtomExpr::Isotope(iso) => atom.isotope == *iso,
            AtomExpr::Degree(d) => {
                let degree = ctx.mol.neighbors(idx).count() as u8;
                degree == *d
            }
            AtomExpr::Valence(v) => {
                let bond_order_sum: u8 = ctx
                    .mol
                    .bonds_of(idx)
                    .map(|ei| match ctx.mol.bond(ei).order {
                        BondOrder::Single => 1u8,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                    })
                    .sum();
                let total = bond_order_sum + atom.hydrogen_count;
                total == *v
            }
            AtomExpr::Connectivity(x) => {
                let degree = ctx.mol.neighbors(idx).count() as u8;
                let total = degree + atom.hydrogen_count;
                total == *x
            }
            AtomExpr::TotalHCount(h) => {
                atom.hydrogen_count == *h
            }
            AtomExpr::ImplicitHCount(h) => {
                atom.hydrogen_count == *h
            }
            AtomExpr::RingMembership(n) => {
                let count = ctx.ring_info.atom_rings(idx).len() as u8;
                count == *n
            }
            AtomExpr::SmallestRingSize(r) => {
                match ctx.ring_info.smallest_ring_size(idx) {
                    Some(size) => size as u8 == *r,
                    None => *r == 0,
                }
            }
            AtomExpr::RingBondCount(x) => {
                let count = ctx
                    .mol
                    .neighbors(idx)
                    .filter(|&nb| ctx.ring_info.is_ring_bond(idx, nb))
                    .count() as u8;
                count == *x
            }
            AtomExpr::Charge(c) => atom.formal_charge == *c,
            AtomExpr::InRing => ctx.ring_info.is_ring_atom(idx),
            AtomExpr::NotInRing => !ctx.ring_info.is_ring_atom(idx),
            AtomExpr::Chirality(q_chiral) => {
                use crate::atom::Chirality;
                match q_chiral {
                    Chirality::None => true,
                    Chirality::Cw | Chirality::Ccw => atom.chirality != Chirality::None,
                }
            }
            AtomExpr::Recursive(ref _inner) => {
                unreachable!("Recursive SMARTS should be pre-evaluated")
            }
            AtomExpr::And(exprs) => exprs.iter().all(|e| e.matches(atom, ctx, idx)),
            AtomExpr::Or(exprs) => exprs.iter().any(|e| e.matches(atom, ctx, idx)),
            AtomExpr::Not(expr) => !expr.matches(atom, ctx, idx),
        }
    }


}

impl BondExpr {
    pub fn matches(&self, bond: &Bond, target_atoms: (&Atom, &Atom)) -> bool {
        match self {
            BondExpr::True => true,
            BondExpr::Single => {
                let both_aromatic =
                    target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
                bond.order == BondOrder::Single && !both_aromatic
            }
            BondExpr::Double => {
                let both_aromatic =
                    target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
                bond.order == BondOrder::Double && !both_aromatic
            }
            BondExpr::Triple => bond.order == BondOrder::Triple,
            BondExpr::Aromatic => {
                target_atoms.0.is_aromatic && target_atoms.1.is_aromatic
            }
            BondExpr::Ring => false,
            BondExpr::SingleOrAromatic => {
                let both_aromatic =
                    target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
                bond.order == BondOrder::Single || both_aromatic
            }
            BondExpr::Up => bond.order == BondOrder::Single,
            BondExpr::Down => bond.order == BondOrder::Single,
            BondExpr::And(exprs) => exprs.iter().all(|e| e.matches(bond, target_atoms)),
            BondExpr::Or(exprs) => exprs.iter().any(|e| e.matches(bond, target_atoms)),
            BondExpr::Not(expr) => !expr.matches(bond, target_atoms),
        }
    }

    pub fn matches_with_ring_info(
        &self,
        bond: &Bond,
        target_atoms: (&Atom, &Atom),
        ring_info: &RingInfo,
        endpoints: (NodeIndex, NodeIndex),
    ) -> bool {
        match self {
            BondExpr::Ring => ring_info.is_ring_bond(endpoints.0, endpoints.1),
            BondExpr::And(exprs) => exprs
                .iter()
                .all(|e| e.matches_with_ring_info(bond, target_atoms, ring_info, endpoints)),
            BondExpr::Or(exprs) => exprs
                .iter()
                .any(|e| e.matches_with_ring_info(bond, target_atoms, ring_info, endpoints)),
            BondExpr::Not(expr) => {
                !expr.matches_with_ring_info(bond, target_atoms, ring_info, endpoints)
            }
            _ => self.matches(bond, target_atoms),
        }
    }
}
