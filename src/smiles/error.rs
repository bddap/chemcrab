use std::fmt;

use crate::kekulize::KekulizeError;

/// Errors produced when parsing a SMILES string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Input ended before a complete token could be read.
    UnexpectedEnd,
    /// An unexpected character was encountered at the given position.
    UnexpectedChar { pos: usize, ch: char },
    /// An unrecognized element symbol was found.
    InvalidElement { pos: usize, text: String },
    /// A bracket atom `[` was opened but never closed with `]`.
    UnclosedBracket { pos: usize },
    /// A ring-opening digit was never matched by a ring-closing digit.
    UnclosedRing { digit: u16 },
    /// A parenthesis was opened without a matching close, or vice versa.
    UnmatchedParen { pos: usize },
    /// A charge specifier inside a bracket atom could not be parsed.
    InvalidCharge { pos: usize },
    /// An isotope number overflowed or was otherwise invalid.
    InvalidIsotope { pos: usize },
    /// An atom class (`:n`) could not be parsed.
    InvalidAtomClass { pos: usize },
    /// A bond was specified on a ring-closure digit that is inconsistent.
    InvalidRingBond { digit: u16, pos: usize },
    /// The input string was empty or contained only whitespace.
    EmptyInput,
    /// Two ring-closure bonds on the same digit specify conflicting bond types.
    RingBondConflict { digit: u16 },
    /// Kekulization of the aromatic system failed.
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
            Self::InvalidIsotope { pos } => {
                write!(f, "isotope overflow at position {}", pos)
            }
            Self::InvalidAtomClass { pos } => {
                write!(f, "atom class overflow at position {}", pos)
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
