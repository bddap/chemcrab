use std::fmt;

use crate::kekulize::KekulizeError;
use crate::smarts::SmartsError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReactionSmartsError {
    MissingSeparator,
    TooManySeparators,
    EmptyReactants,
    EmptyProducts,
    InvalidComponent {
        section: &'static str,
        detail: SmartsError,
    },
}

impl fmt::Display for ReactionSmartsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingSeparator => write!(f, "no '>>' separator found in reaction SMARTS"),
            Self::TooManySeparators => write!(f, "too many '>' separators in reaction SMARTS"),
            Self::EmptyReactants => write!(f, "reaction has no reactant templates"),
            Self::EmptyProducts => write!(f, "reaction has no product templates"),
            Self::InvalidComponent { section, detail } => {
                write!(f, "invalid {section} component: {detail}")
            }
        }
    }
}

impl std::error::Error for ReactionSmartsError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReactionError {
    WrongReactantCount { expected: usize, got: usize },
    TooManyCombinations,
    DuplicateAtomMap { map_num: u16 },
    Kekulize(KekulizeError),
}

impl fmt::Display for ReactionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongReactantCount { expected, got } => {
                write!(f, "expected {expected} reactants, got {got}")
            }
            Self::TooManyCombinations => {
                write!(f, "match combination count exceeds limit")
            }
            Self::DuplicateAtomMap { map_num } => {
                write!(
                    f,
                    "duplicate atom map number {map_num} in reactant templates"
                )
            }
            Self::Kekulize(e) => write!(f, "product kekulization failed: {e}"),
        }
    }
}

impl std::error::Error for ReactionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Kekulize(e) => Some(e),
            _ => None,
        }
    }
}

impl From<KekulizeError> for ReactionError {
    fn from(e: KekulizeError) -> Self {
        Self::Kekulize(e)
    }
}
