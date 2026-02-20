use std::collections::{HashMap, HashSet};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};
use crate::mol::Mol;
use crate::rings::RingInfo;

/// Hybridization state for SMARTS `^n` queries.
///
/// Computed from the bond order sum and neighbor count of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    /// s hybridization (no bonds or lone atoms).
    S,
    /// sp hybridization (e.g. alkynes, nitriles).
    SP,
    /// sp2 hybridization (e.g. alkenes, aromatic atoms).
    SP2,
    /// sp3 hybridization (e.g. saturated carbon).
    SP3,
    /// sp3d hybridization (e.g. pentacoordinate phosphorus).
    SP3D,
    /// sp3d2 hybridization (e.g. hexacoordinate sulfur).
    SP3D2,
}

/// Identifies the property being tested in a SMARTS range expression `{low-high}`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RangeKind {
    /// Heavy-atom degree (`D`).
    Degree,
    /// Non-hydrogen degree (`d`).
    NonHDegree,
    /// Total hydrogen count including explicit H neighbors (`H`).
    TotalHCount,
    /// Implicit (virtual) hydrogen count (`h`).
    ImplicitHCount,
    /// Size of the smallest SSSR ring containing the atom (`r`).
    SmallestRingSize,
    /// Number of SSSR rings containing the atom (`R`).
    RingMembership,
    /// Sum of bond orders plus implicit hydrogens (`v`).
    Valence,
    /// Number of ring bonds on the atom (`x`).
    RingBondCount,
    /// Total connectivity: degree plus implicit H count (`X`).
    Connectivity,
    /// Number of non-C, non-H neighbors (`z`).
    HeteroNeighborCount,
    /// Number of non-C, non-H, non-aromatic neighbors (`Z`).
    AliphaticHeteroNeighborCount,
    /// Positive formal charge value (`+{n-m}`).
    PositiveCharge,
    /// Negative formal charge magnitude (`-{n-m}`).
    NegativeCharge,
}

/// AST node for a SMARTS atom query expression.
///
/// Each variant represents a primitive test or a logical combination of tests.
/// During substructure search, [`AtomExpr::matches`] evaluates the expression
/// tree against a target atom.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomExpr {
    /// Matches any atom (wildcard `*`).
    True,
    /// Matches by element. `aromatic` is `None` for `#n` (either), `Some(true)`
    /// for lowercase (`c`), `Some(false)` for uppercase (`C`).
    Element {
        atomic_num: u8,
        aromatic: Option<bool>,
    },
    /// Matches any aromatic atom (`a`).
    Aromatic,
    /// Matches any aliphatic atom (`A`).
    Aliphatic,
    /// Matches a specific isotope number.
    Isotope(u16),
    /// Matches explicit degree — number of heavy-atom neighbors (`D`).
    Degree(u8),
    /// Matches non-hydrogen degree (`d`).
    NonHDegree(u8),
    /// Matches total valence: sum of bond orders plus implicit H count (`v`).
    Valence(u8),
    /// Matches total connectivity: degree plus implicit H count (`X`).
    Connectivity(u8),
    /// Matches total hydrogen count including explicit H neighbors (`H`).
    TotalHCount(u8),
    /// Matches implicit (virtual) hydrogen count (`h`).
    ImplicitHCount(u8),
    /// Matches the number of SSSR rings containing the atom (`R`).
    RingMembership(u8),
    /// Matches the size of the smallest SSSR ring containing the atom (`r`).
    SmallestRingSize(u8),
    /// Matches the number of ring bonds on the atom (`x`).
    RingBondCount(u8),
    /// Matches formal charge.
    Charge(i8),
    /// Matches the count of non-C, non-H neighbors (`z` with value).
    HeteroNeighborCount(u8),
    /// Matches the count of aliphatic non-C, non-H neighbors (`Z` with value).
    AliphaticHeteroNeighborCount(u8),
    /// Matches atoms with at least one heteroatom neighbor (`z` bare).
    HasHeteroNeighbor,
    /// Matches atoms with at least one aliphatic heteroatom neighbor (`Z` bare).
    HasAliphaticHeteroNeighbor,
    /// Matches computed hybridization state (`^n`).
    Hybridization(Hybridization),
    /// Matches a numeric property within a range (`{low-high}`).
    Range {
        kind: RangeKind,
        low: Option<u8>,
        high: Option<u8>,
    },
    /// Matches atoms in at least one ring (`R` bare).
    InRing,
    /// Matches atoms not in any ring (`R0`).
    NotInRing,
    /// A recursive SMARTS sub-query (`$(...)`).
    Recursive(Mol<AtomExpr, BondExpr>),
    /// Atom map class (`:n`) — always matches, used for reaction mapping.
    AtomMapClass(u16),
    /// Logical AND of sub-expressions.
    And(Vec<AtomExpr>),
    /// Logical OR of sub-expressions.
    Or(Vec<AtomExpr>),
    /// Matches tetrahedral chirality (`@` or `@@`).
    Chirality(crate::atom::Chirality),
    /// Logical NOT of a sub-expression.
    Not(Box<AtomExpr>),
}

/// AST node for a SMARTS bond query expression.
///
/// Implicit bonds in SMARTS default to [`BondExpr::SingleOrAromatic`], unlike
/// SMILES where implicit bonds are always single.
#[derive(Debug, Clone, PartialEq)]
pub enum BondExpr {
    /// Matches any bond (`~`).
    True,
    /// Matches a single bond (`-`), excluding aromatic bonds between aromatic atoms.
    Single,
    /// Matches a double bond (`=`), excluding aromatic bonds.
    Double,
    /// Matches a triple bond (`#`).
    Triple,
    /// Matches an aromatic bond (`:`) — both endpoints must be aromatic.
    Aromatic,
    /// Matches a ring bond (`@`).
    Ring,
    /// Default SMARTS bond: matches single or aromatic.
    SingleOrAromatic,
    /// Matches an up directional bond (`/`).
    Up,
    /// Matches a down directional bond (`\`).
    Down,
    /// Logical AND of sub-expressions.
    And(Vec<BondExpr>),
    /// Logical OR of sub-expressions.
    Or(Vec<BondExpr>),
    /// Logical NOT of a sub-expression.
    Not(Box<BondExpr>),
}

/// Context passed to [`AtomExpr::matches`] during SMARTS evaluation.
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

/// Computes the hybridization state of an atom from its bond orders and
/// neighbor count.
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
            AtomExpr::Element {
                atomic_num,
                aromatic,
            } => atom.atomic_num == *atomic_num && aromatic.is_none_or(|a| atom.is_aromatic == a),
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
            AtomExpr::ImplicitHCount(h) => atom.hydrogen_count == *h,
            AtomExpr::RingMembership(n) => {
                let count = ctx.ring_info.atom_rings(idx).len() as u8;
                count == *n
            }
            AtomExpr::SmallestRingSize(r) => match ctx.ring_info.smallest_ring_size(idx) {
                Some(size) => size as u8 == *r,
                None => *r == 0,
            },
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
                    Chirality::Cw | Chirality::Ccw => ctx.mol.tetrahedral_stereo_for(idx).is_some(),
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
                let both_aromatic = target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
                bond.order == BondOrder::Single && !both_aromatic
            }
            BondExpr::Double => {
                let both_aromatic = target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
                bond.order == BondOrder::Double && !both_aromatic
            }
            BondExpr::Triple => bond.order == BondOrder::Triple,
            BondExpr::Aromatic => target_atoms.0.is_aromatic && target_atoms.1.is_aromatic,
            BondExpr::Ring => false,
            BondExpr::SingleOrAromatic => {
                let both_aromatic = target_atoms.0.is_aromatic && target_atoms.1.is_aromatic;
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
