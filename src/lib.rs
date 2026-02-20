//! Cheminformatics library for Rust.
//!
//! Chemcrab provides a molecular representation, file format parsers, and
//! algorithms for working with chemical structures. If you are new to
//! cheminformatics, the module-level docs below introduce both the API and
//! the chemistry behind it.
//!
//! # Quick start
//!
//! ```
//! use chemcrab::Mol;
//! use chemcrab::smiles;
//!
//! // Parse a SMILES string into a molecule.
//! // SMILES encodes molecular graphs as text — atoms are nodes, bonds are
//! // edges, and ring closures are back-edges written with matching digits.
//! let caffeine = smiles::from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")?;
//!
//! assert_eq!(caffeine.atom_count(), 14);
//! assert_eq!(caffeine.bond_count(), 15);
//!
//! // Canonical SMILES: a unique string for each molecule, regardless of
//! // how the atoms were numbered in the input.
//! let canonical = smiles::to_canonical_smiles(&caffeine);
//! # Ok::<(), chemcrab::smiles::SmilesError>(())
//! ```
//!
//! # Molecular representation
//!
//! [`Mol<A, B>`] is a generic molecular graph parameterized over atom type `A`
//! and bond type `B`. The default concrete types are [`Atom`] and [`Bond`]:
//!
//! - **[`Atom`]** stores atomic number, formal charge, isotope, hydrogen count,
//!   and aromaticity. No cached properties — valence, hybridization, and
//!   coordinates are added via wrapper types when needed.
//!
//! - **[`Bond`]** stores a [`BondOrder`] (`Single`, `Double`, or `Triple`).
//!   There is no `Aromatic` bond order — every bond in a `Mol<Atom, Bond>` has
//!   a concrete Kekulé assignment.
//!
//! Algorithms declare their requirements via trait bounds
//! (`HasAtomicNum`, `HasBondOrder`, etc.) so they work with any atom/bond
//! types that implement the right traits. See the [`traits`] module.
//!
//! ## Hydrogen model
//!
//! Chemcrab uses a *virtual hydrogen* model: hydrogen atoms that are not
//! explicit nodes in the graph are represented as an integer count on their
//! parent atom ([`Atom::hydrogen_count`]). The SMILES parser resolves
//! implicit hydrogen counts at parse time using the element's
//! [default valences](Element::default_valences); after parsing, the stored
//! count is the single source of truth.
//!
//! Use [`hydrogen::add_hs`] and [`hydrogen::remove_hs`] to convert between
//! explicit and virtual representations.
//!
//! ## Stereochemistry
//!
//! Stereochemistry lives on [`Mol`], not on atoms or bonds:
//!
//! - **Tetrahedral**: [`TetrahedralStereo`] stores four neighbors in a
//!   specific order. The first neighbor is *above* the plane defined by the
//!   other three (right-hand rule). Swapping any pair flips the handedness.
//!
//! - **E/Z (cis/trans)**: [`EZStereo`] stores a double bond and two reference
//!   substituents that are on the *same side* (cis). There is no explicit
//!   `Cis`/`Trans` enum — the geometry is fully determined by which
//!   substituents are recorded as the cis pair.
//!
//! Both use [`AtomId`] to refer to neighbors, which can be either a graph
//! node or a virtual hydrogen.

// ── Internal modules ────────────────────────────────────────────────────
//
// These remain `pub` so that integration tests and internal code can reach
// them via `crate::` paths, but they are hidden from rustdoc.  The facade
// modules below re-export the curated public API.

#[doc(hidden)]
pub mod atom;
#[doc(hidden)]
pub mod bond;
#[doc(hidden)]
pub mod canonical;
#[doc(hidden)]
pub mod conjugation;
#[doc(hidden)]
pub mod graph_ops;
#[doc(hidden)]
pub mod mol;

// ── Public modules ──────────────────────────────────────────────────────

pub mod aromaticity;
pub mod chirality;
pub mod element;
pub mod formula;
pub mod hybridization;
pub mod hydrogen;
pub mod kekulize;
pub mod radical;
pub mod reaction;
pub mod rings;
pub mod smarts;
pub mod smiles;
pub mod strip;
pub mod substruct;
pub mod traits;
pub mod valence;
pub mod wrappers;

/// Graph algorithms on molecular graphs.
///
/// These operate on the topology of a [`Mol`] — connected components,
/// shortest paths, fragmentation, and atom renumbering. Higher-level
/// chemical algorithms (rings, aromaticity) live in their own modules.
pub mod graph {
    pub use crate::graph_ops::{
        connected_components, get_fragments, renumber_atoms, shortest_path, RenumberError,
    };
}

// ── Root re-exports ─────────────────────────────────────────────────────
//
// The most commonly used types are re-exported at the crate root so that
// `use chemcrab::{Mol, Atom, Bond}` works without diving into modules.

pub use atom::Atom;
pub use bond::{Bond, BondOrder};
pub use element::Element;
pub use mol::{AtomId, EZStereo, Mol, TetrahedralStereo};

#[cfg(test)]
mod tests;
