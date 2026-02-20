use std::fmt;

/// Errors produced when parsing a SMARTS pattern string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmartsError {
    /// The input string was empty.
    EmptyInput,
    /// An unexpected character was encountered at the given position.
    UnexpectedChar { pos: usize, ch: char },
    /// A bracket atom `[` was opened but never closed with `]`.
    UnclosedBracket { pos: usize },
    /// A ring-opening digit was never matched by a ring-closing digit.
    UnclosedRing { digit: u16 },
    /// A parenthesis was opened without a matching close, or vice versa.
    UnmatchedParen { pos: usize },
    /// An `#n` atomic number specifier could not be parsed.
    InvalidAtomicNum { pos: usize },
    /// A recursive SMARTS `$( ... )` was opened but never closed.
    UnclosedRecursive { pos: usize },
    /// A catch-all for other SMARTS parse errors.
    InvalidSmarts { pos: usize, msg: String },
}

impl fmt::Display for SmartsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty SMARTS string"),
            Self::UnexpectedChar { pos, ch } => {
                write!(f, "unexpected character '{ch}' at position {pos}")
            }
            Self::UnclosedBracket { pos } => {
                write!(f, "unclosed bracket starting at position {pos}")
            }
            Self::UnclosedRing { digit } => write!(f, "unclosed ring {digit}"),
            Self::UnmatchedParen { pos } => {
                write!(f, "unmatched parenthesis at position {pos}")
            }
            Self::InvalidAtomicNum { pos } => {
                write!(f, "invalid atomic number at position {pos}")
            }
            Self::UnclosedRecursive { pos } => {
                write!(f, "unclosed recursive SMARTS at position {pos}")
            }
            Self::InvalidSmarts { pos, msg } => {
                write!(f, "invalid SMARTS at position {pos}: {msg}")
            }
        }
    }
}

impl std::error::Error for SmartsError {}
