use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};
use crate::mol::Mol;
use crate::rings::RingInfo;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    S,
    SP,
    SP2,
    SP3,
    SP3D,
    SP3D2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RangeKind {
    Degree,
    NonHDegree,
    TotalHCount,
    ImplicitHCount,
    SmallestRingSize,
    RingMembership,
    Valence,
    RingBondCount,
    Connectivity,
    HeteroNeighborCount,
    AliphaticHeteroNeighborCount,
    PositiveCharge,
    NegativeCharge,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AtomExpr {
    True,
    Element { atomic_num: u8, aromatic: Option<bool> },
    Aromatic,
    Aliphatic,
    Isotope(u16),
    Degree(u8),
    NonHDegree(u8),
    Valence(u8),
    Connectivity(u8),
    TotalHCount(u8),
    ImplicitHCount(u8),
    RingMembership(u8),
    SmallestRingSize(u8),
    RingBondCount(u8),
    Charge(i8),
    HeteroNeighborCount(u8),
    AliphaticHeteroNeighborCount(u8),
    HasHeteroNeighbor,
    HasAliphaticHeteroNeighbor,
    Hybridization(Hybridization),
    Range { kind: RangeKind, low: Option<u8>, high: Option<u8> },
    InRing,
    NotInRing,
    Recursive(Mol<AtomExpr, BondExpr>),
    AtomMapClass(u16),
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

fn non_h_degree(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> u8 {
    mol.neighbors(idx)
        .filter(|&nb| mol.atom(nb).atomic_num != 1)
        .count() as u8
}

fn hetero_neighbor_count(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> u8 {
    mol.neighbors(idx)
        .filter(|&nb| {
            let a = mol.atom(nb).atomic_num;
            a != 6 && a != 1
        })
        .count() as u8
}

fn aliphatic_hetero_neighbor_count(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> u8 {
    mol.neighbors(idx)
        .filter(|&nb| {
            let nbr = mol.atom(nb);
            nbr.atomic_num != 6 && nbr.atomic_num != 1 && !nbr.is_aromatic
        })
        .count() as u8
}

fn explicit_h_count(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> u8 {
    mol.neighbors(idx)
        .filter(|&nb| mol.atom(nb).atomic_num == 1)
        .count() as u8
}

fn bond_order_sum(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> u8 {
    mol.bonds_of(idx)
        .map(|ei| match mol.bond(ei).order {
            BondOrder::Single => 1u8,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        })
        .sum()
}

pub fn compute_hybridization(mol: &Mol<Atom, Bond>, idx: NodeIndex) -> Hybridization {
    let atom = mol.atom(idx);
    let explicit_degree = mol.neighbors(idx).count() as u8;
    let sum_orders = bond_order_sum(mol, idx);
    let n_bonds = explicit_degree + atom.hydrogen_count;
    let total_orders = sum_orders + atom.hydrogen_count;
    let n_unsaturations = total_orders.saturating_sub(n_bonds);

    match n_unsaturations {
        0 => match n_bonds {
            0 | 1 => Hybridization::S,
            2..=4 => Hybridization::SP3,
            5 => Hybridization::SP3D,
            _ => Hybridization::SP3D2,
        },
        1 => match n_bonds {
            1 => Hybridization::SP,
            _ => Hybridization::SP2,
        },
        _ => Hybridization::SP,
    }
}

fn range_value(kind: RangeKind, atom: &Atom, ctx: &MatchContext, idx: NodeIndex) -> u8 {
    match kind {
        RangeKind::Degree => ctx.mol.neighbors(idx).count() as u8,
        RangeKind::NonHDegree => non_h_degree(ctx.mol, idx),
        RangeKind::TotalHCount => atom.hydrogen_count + explicit_h_count(ctx.mol, idx),
        RangeKind::ImplicitHCount => atom.hydrogen_count,
        RangeKind::SmallestRingSize => ctx
            .ring_info
            .smallest_ring_size(idx)
            .map(|s| s as u8)
            .unwrap_or(0),
        RangeKind::RingMembership => ctx.ring_info.atom_rings(idx).len() as u8,
        RangeKind::Valence => bond_order_sum(ctx.mol, idx) + atom.hydrogen_count,
        RangeKind::RingBondCount => ctx
            .mol
            .neighbors(idx)
            .filter(|&nb| ctx.ring_info.is_ring_bond(idx, nb))
            .count() as u8,
        RangeKind::Connectivity => ctx.mol.neighbors(idx).count() as u8 + atom.hydrogen_count,
        RangeKind::HeteroNeighborCount => hetero_neighbor_count(ctx.mol, idx),
        RangeKind::AliphaticHeteroNeighborCount => aliphatic_hetero_neighbor_count(ctx.mol, idx),
        RangeKind::PositiveCharge => {
            if atom.formal_charge > 0 {
                atom.formal_charge as u8
            } else {
                0
            }
        }
        RangeKind::NegativeCharge => {
            if atom.formal_charge < 0 {
                atom.formal_charge.unsigned_abs()
            } else {
                0
            }
        }
    }
}

fn in_range(val: u8, low: Option<u8>, high: Option<u8>) -> bool {
    if let Some(lo) = low {
        if val < lo {
            return false;
        }
    }
    if let Some(hi) = high {
        if val > hi {
            return false;
        }
    }
    true
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
            AtomExpr::NonHDegree(d) => non_h_degree(ctx.mol, idx) == *d,
            AtomExpr::Valence(v) => {
                let total = bond_order_sum(ctx.mol, idx) + atom.hydrogen_count;
                total == *v
            }
            AtomExpr::Connectivity(x) => {
                let degree = ctx.mol.neighbors(idx).count() as u8;
                let total = degree + atom.hydrogen_count;
                total == *x
            }
            AtomExpr::TotalHCount(h) => {
                (atom.hydrogen_count + explicit_h_count(ctx.mol, idx)) == *h
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
            AtomExpr::HeteroNeighborCount(n) => hetero_neighbor_count(ctx.mol, idx) == *n,
            AtomExpr::AliphaticHeteroNeighborCount(n) => {
                aliphatic_hetero_neighbor_count(ctx.mol, idx) == *n
            }
            AtomExpr::HasHeteroNeighbor => hetero_neighbor_count(ctx.mol, idx) > 0,
            AtomExpr::HasAliphaticHeteroNeighbor => {
                aliphatic_hetero_neighbor_count(ctx.mol, idx) > 0
            }
            AtomExpr::Hybridization(h) => compute_hybridization(ctx.mol, idx) == *h,
            AtomExpr::Range { kind, low, high } => {
                let val = range_value(*kind, atom, ctx, idx);
                in_range(val, *low, *high)
            }
            AtomExpr::AtomMapClass(_) => true,
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
