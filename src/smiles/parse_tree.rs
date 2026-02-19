use crate::element::Element;
use crate::smiles::error::SmilesError;
use crate::smiles::tokenizer::{AtomToken, BondToken, ChiralityToken, Token};

#[derive(Debug, Clone)]
pub struct ParseAtom {
    pub element: Element,
    pub is_aromatic: bool,
    pub isotope: u16,
    pub chirality: ChiralityToken,
    pub hcount: Option<u8>,
    pub charge: i8,
    #[allow(dead_code)]
    pub atom_class: u16,
    pub is_bracket: bool,
    pub neighbors: Vec<Neighbor>,
}

#[derive(Debug, Clone)]
pub struct Neighbor {
    pub bond: Option<BondToken>,
    pub atom_idx: usize,
}

#[derive(Debug, Clone)]
pub struct ParseTree {
    pub atoms: Vec<ParseAtom>,
}

pub fn build_parse_tree(tokens: &[Token]) -> Result<ParseTree, SmilesError> {
    let mut atoms: Vec<ParseAtom> = Vec::new();
    let mut stack: Vec<usize> = Vec::new(); // stack of current atom indices for branches
    let mut current: Option<usize> = None;
    let mut pending_bond: Option<BondToken> = None;
    let mut ring_opens: Vec<Option<(usize, Option<BondToken>, usize)>> = vec![None; 100];

    for token in tokens {
        match token {
            Token::Atom(atom_tok) => {
                let idx = atoms.len();
                atoms.push(parse_atom_from_token(atom_tok));

                if let Some(cur) = current {
                    let bond = pending_bond.take();
                    atoms[cur].neighbors.push(Neighbor {
                        bond,
                        atom_idx: idx,
                    });
                    atoms[idx].neighbors.push(Neighbor {
                        bond,
                        atom_idx: cur,
                    });
                } else {
                    pending_bond = None;
                }

                current = Some(idx);
            }
            Token::Bond(b) => {
                pending_bond = Some(*b);
            }
            Token::RingClosure { bond, digit, pos } => {
                let d = *digit as usize;
                let cur = current.ok_or(SmilesError::InvalidRingBond {
                    digit: *digit,
                    pos: *pos,
                })?;

                if let Some((open_idx, open_bond, _open_pos)) = ring_opens[d].take() {
                    let ring_bond = match (bond.or(pending_bond.take()), open_bond) {
                        (None, None) => None,
                        (Some(b), None) | (None, Some(b)) => Some(b),
                        (Some(b1), Some(b2)) => {
                            if b1 == b2 {
                                Some(b1)
                            } else {
                                return Err(SmilesError::RingBondConflict { digit: *digit });
                            }
                        }
                    };

                    atoms[open_idx].neighbors.push(Neighbor {
                        bond: ring_bond,
                        atom_idx: cur,
                    });
                    atoms[cur].neighbors.push(Neighbor {
                        bond: ring_bond,
                        atom_idx: open_idx,
                    });
                } else {
                    ring_opens[d] = Some((cur, bond.or(pending_bond.take()), *pos));
                }
            }
            Token::OpenParen(pos) => {
                let cur = current.ok_or(SmilesError::UnmatchedParen { pos: *pos })?;
                stack.push(cur);
            }
            Token::CloseParen(pos) => {
                current = Some(stack.pop().ok_or(SmilesError::UnmatchedParen { pos: *pos })?);
                pending_bond = None;
            }
            Token::Dot(_) => {
                current = None;
                pending_bond = None;
            }
        }
    }

    if !stack.is_empty() {
        return Err(SmilesError::UnmatchedParen {
            pos: 0,
        });
    }

    for (digit, entry) in ring_opens.iter().enumerate() {
        if entry.is_some() {
            return Err(SmilesError::UnclosedRing {
                digit: digit as u16,
            });
        }
    }

    Ok(ParseTree { atoms })
}

fn parse_atom_from_token(tok: &AtomToken) -> ParseAtom {
    ParseAtom {
        element: tok.element,
        is_aromatic: tok.is_aromatic,
        isotope: tok.isotope,
        chirality: tok.chirality,
        hcount: tok.hcount,
        charge: tok.charge,
        atom_class: tok.atom_class,
        is_bracket: tok.is_bracket,
        neighbors: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::tokenizer::tokenize;

    #[test]
    fn ethane_tree() {
        let tokens = tokenize("CC").unwrap();
        let tree = build_parse_tree(&tokens).unwrap();
        assert_eq!(tree.atoms.len(), 2);
        assert_eq!(tree.atoms[0].neighbors.len(), 1);
        assert_eq!(tree.atoms[0].neighbors[0].atom_idx, 1);
    }

    #[test]
    fn cyclohexane_tree() {
        let tokens = tokenize("C1CCCCC1").unwrap();
        let tree = build_parse_tree(&tokens).unwrap();
        assert_eq!(tree.atoms.len(), 6);
        for atom in &tree.atoms {
            assert_eq!(atom.neighbors.len(), 2);
        }
    }

    #[test]
    fn branch_tree() {
        let tokens = tokenize("CC(C)C").unwrap();
        let tree = build_parse_tree(&tokens).unwrap();
        assert_eq!(tree.atoms.len(), 4);
        assert_eq!(tree.atoms[1].neighbors.len(), 3);
    }

    #[test]
    fn unclosed_ring_error() {
        let tokens = tokenize("C1CC").unwrap();
        let result = build_parse_tree(&tokens);
        assert!(result.is_err());
    }

    #[test]
    fn unmatched_paren_error() {
        let tokens = tokenize("C(C").unwrap();
        let result = build_parse_tree(&tokens);
        assert!(result.is_err());
    }

    #[test]
    fn disconnected() {
        let tokens = tokenize("[Na+].[Cl-]").unwrap();
        let tree = build_parse_tree(&tokens).unwrap();
        assert_eq!(tree.atoms.len(), 2);
        assert_eq!(tree.atoms[0].neighbors.len(), 0);
        assert_eq!(tree.atoms[1].neighbors.len(), 0);
    }
}
