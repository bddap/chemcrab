use std::fmt;

use crate::kekulize::KekulizeError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    UnexpectedEnd,
    UnexpectedChar { pos: usize, ch: char },
    InvalidElement { pos: usize, text: String },
    UnclosedBracket { pos: usize },
    UnclosedRing { digit: u16 },
    UnmatchedParen { pos: usize },
    InvalidCharge { pos: usize },
    InvalidRingBond { digit: u16, pos: usize },
    EmptyInput,
    RingBondConflict { digit: u16 },
    Kekulize(KekulizeError),
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEnd => write!(f, "unexpected end of SMILES"),
            Self::UnexpectedChar { pos, ch } => {
                write!(f, "unexpected character '{}' at position {}", ch, pos)
            }
            Self::InvalidElement { pos, text } => {
                write!(f, "invalid element '{}' at position {}", text, pos)
            }
            Self::UnclosedBracket { pos } => {
                write!(f, "unclosed bracket atom starting at position {}", pos)
            }
            Self::UnclosedRing { digit } => write!(f, "unclosed ring {}", digit),
            Self::UnmatchedParen { pos } => {
                write!(f, "unmatched parenthesis at position {}", pos)
            }
            Self::InvalidCharge { pos } => {
                write!(f, "invalid charge at position {}", pos)
            }
            Self::InvalidRingBond { digit, pos } => {
                write!(f, "invalid ring bond {} at position {}", digit, pos)
            }
            Self::EmptyInput => write!(f, "empty SMILES string"),
            Self::RingBondConflict { digit } => {
                write!(f, "conflicting bond types on ring closure {}", digit)
            }
            Self::Kekulize(e) => write!(f, "{}", e),
        }
    }
}

impl std::error::Error for SmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Kekulize(e) => Some(e),
            _ => None,
        }
    }
}

impl From<KekulizeError> for SmilesError {
    fn from(e: KekulizeError) -> Self {
        Self::Kekulize(e)
    }
}
