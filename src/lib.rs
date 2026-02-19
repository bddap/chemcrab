pub mod aromaticity;
pub mod atom;
pub mod bond;
pub mod canonical;
pub mod element;
pub mod formula;
pub mod graph_ops;
pub mod hydrogen;
pub mod kekulize;
pub mod mol;
pub mod reaction;
pub mod rings;
pub mod smarts;
pub mod smiles;
pub mod substruct;
pub mod traits;
pub mod valence;
pub mod wrappers;

pub use aromaticity::{find_aromatic_atoms, set_aromaticity, AromaticityModel};
pub use atom::{Atom, Chirality};
pub use canonical::canonical_ordering;
pub use bond::{Bond, BondOrder, BondStereo, SmilesBond, SmilesBondOrder};
pub use element::Element;
pub use formula::{average_mol_weight, exact_mol_weight, mol_formula};
pub use graph_ops::{
    adjacency_matrix, connected_components, distance_matrix, get_fragments, num_components,
    shortest_path,
};
pub use hydrogen::{add_hs, remove_hs};
pub use kekulize::{kekulize, KekulizeError};
pub use mol::Mol;
pub use reaction::{
    extract_atom_map_num, from_reaction_smarts, to_reaction_smarts, Reaction, ReactionError,
    ReactionSmartsError,
};
pub use rings::RingInfo;
pub use smarts::{
    from_smarts, get_smarts_match, get_smarts_matches, has_smarts_match, to_smarts, AtomExpr,
    BondExpr, SmartsError,
};
pub use smiles::{from_smiles, parse_smiles, to_canonical_smiles, to_smiles, SmilesError};
pub use substruct::{
    get_substruct_match, get_substruct_match_with, get_substruct_matches,
    get_substruct_matches_with, has_substruct_match, has_substruct_match_with, AtomMapping,
};
pub use traits::{
    HasAromaticity, HasAtomicNum, HasBondOrder, HasBondStereo, HasChirality, HasFormalCharge,
    HasHybridization, HasHydrogenCount, HasIsotope, HasPosition2D, HasPosition3D, HasValence,
};
pub use valence::{check_valence, total_valence, ValenceError};
pub use wrappers::{Hybridization, WithHybridization, WithPosition2D, WithPosition3D, WithValence};

#[cfg(test)]
mod tests;
