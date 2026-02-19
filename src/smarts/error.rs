use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmartsError {
    EmptyInput,
    UnexpectedChar { pos: usize, ch: char },
    UnclosedBracket { pos: usize },
    UnclosedRing { digit: u16 },
    UnmatchedParen { pos: usize },
    InvalidAtomicNum { pos: usize },
    UnclosedRecursive { pos: usize },
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
