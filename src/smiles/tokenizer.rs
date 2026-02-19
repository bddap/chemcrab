use crate::element::Element;
use crate::smiles::error::SmilesError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token {
    Atom(AtomToken),
    Bond(BondToken),
    RingClosure {
        bond: Option<BondToken>,
        digit: u16,
        pos: usize,
    },
    OpenParen(usize),
    CloseParen(usize),
    Dot(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomToken {
    pub element: Element,
    pub is_aromatic: bool,
    pub isotope: u16,
    pub chirality: ChiralityToken,
    pub hcount: Option<u8>,
    pub charge: i8,
    pub atom_class: u16,
    pub is_bracket: bool,
    pub pos: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChiralityToken {
    None,
    CounterClockwise,
    Clockwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondToken {
    Single,
    Double,
    Triple,
    Aromatic,
    Up,
    Down,
}

pub fn tokenize(input: &str) -> Result<Vec<Token>, SmilesError> {
    let chars: Vec<char> = input.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < chars.len() {
        match chars[i] {
            ' ' | '\t' | '\r' | '\n' => {
                i += 1;
            }
            '[' => {
                let (tok, next) = parse_bracket_atom(&chars, i)?;
                tokens.push(Token::Atom(tok));
                i = next;
            }
            'B' => {
                if i + 1 < chars.len() && chars[i + 1] == 'r' {
                    tokens.push(Token::Atom(bare_atom(Element::Br, false, i)));
                    i += 2;
                } else {
                    tokens.push(Token::Atom(bare_atom(Element::B, false, i)));
                    i += 1;
                }
            }
            'C' => {
                if i + 1 < chars.len() && chars[i + 1] == 'l' {
                    tokens.push(Token::Atom(bare_atom(Element::Cl, false, i)));
                    i += 2;
                } else {
                    tokens.push(Token::Atom(bare_atom(Element::C, false, i)));
                    i += 1;
                }
            }
            'N' => {
                tokens.push(Token::Atom(bare_atom(Element::N, false, i)));
                i += 1;
            }
            'O' => {
                tokens.push(Token::Atom(bare_atom(Element::O, false, i)));
                i += 1;
            }
            'P' => {
                tokens.push(Token::Atom(bare_atom(Element::P, false, i)));
                i += 1;
            }
            'S' => {
                tokens.push(Token::Atom(bare_atom(Element::S, false, i)));
                i += 1;
            }
            'F' => {
                tokens.push(Token::Atom(bare_atom(Element::F, false, i)));
                i += 1;
            }
            'I' => {
                tokens.push(Token::Atom(bare_atom(Element::I, false, i)));
                i += 1;
            }
            'b' => {
                tokens.push(Token::Atom(bare_atom(Element::B, true, i)));
                i += 1;
            }
            'c' => {
                tokens.push(Token::Atom(bare_atom(Element::C, true, i)));
                i += 1;
            }
            'n' => {
                tokens.push(Token::Atom(bare_atom(Element::N, true, i)));
                i += 1;
            }
            'o' => {
                tokens.push(Token::Atom(bare_atom(Element::O, true, i)));
                i += 1;
            }
            'p' => {
                tokens.push(Token::Atom(bare_atom(Element::P, true, i)));
                i += 1;
            }
            's' => {
                tokens.push(Token::Atom(bare_atom(Element::S, true, i)));
                i += 1;
            }
            '-' => {
                if looks_like_bond(&tokens) {
                    tokens.push(Token::Bond(BondToken::Single));
                    i += 1;
                } else {
                    return Err(SmilesError::UnexpectedChar { pos: i, ch: '-' });
                }
            }
            '=' => {
                tokens.push(Token::Bond(BondToken::Double));
                i += 1;
            }
            '#' => {
                tokens.push(Token::Bond(BondToken::Triple));
                i += 1;
            }
            ':' => {
                tokens.push(Token::Bond(BondToken::Aromatic));
                i += 1;
            }
            '/' => {
                tokens.push(Token::Bond(BondToken::Up));
                i += 1;
            }
            '\\' => {
                tokens.push(Token::Bond(BondToken::Down));
                i += 1;
            }
            '(' => {
                tokens.push(Token::OpenParen(i));
                i += 1;
            }
            ')' => {
                tokens.push(Token::CloseParen(i));
                i += 1;
            }
            '.' => {
                tokens.push(Token::Dot(i));
                i += 1;
            }
            '%' => {
                let (bond_tok, digit, next) = parse_percent_ring(&chars, i, &tokens)?;
                tokens.push(Token::RingClosure {
                    bond: bond_tok,
                    digit,
                    pos: i,
                });
                i = next;
            }
            d @ '0'..='9' => {
                let pending_bond = try_consume_pending_bond(&mut tokens);
                tokens.push(Token::RingClosure {
                    bond: pending_bond,
                    digit: (d as u16) - b'0' as u16,
                    pos: i,
                });
                i += 1;
            }
            ch => return Err(SmilesError::UnexpectedChar { pos: i, ch }),
        }
    }

    Ok(tokens)
}

fn bare_atom(element: Element, aromatic: bool, pos: usize) -> AtomToken {
    AtomToken {
        element,
        is_aromatic: aromatic,
        isotope: 0,
        chirality: ChiralityToken::None,
        hcount: None,
        charge: 0,
        atom_class: 0,
        is_bracket: false,
        pos,
    }
}

fn looks_like_bond(tokens: &[Token]) -> bool {
    matches!(
        tokens.last(),
        Some(Token::Atom(_))
            | Some(Token::RingClosure { .. })
            | Some(Token::CloseParen(_))
            | None
    )
}

fn try_consume_pending_bond(tokens: &mut Vec<Token>) -> Option<BondToken> {
    if let Some(Token::Bond(_)) = tokens.last() {
        if let Some(Token::Bond(b)) = tokens.pop() {
            return Some(b);
        }
    }
    None
}

fn parse_percent_ring(
    chars: &[char],
    start: usize,
    _tokens: &[Token],
) -> Result<(Option<BondToken>, u16, usize), SmilesError> {
    let i = start + 1;
    if i + 1 >= chars.len() || !chars[i].is_ascii_digit() || !chars[i + 1].is_ascii_digit() {
        return Err(SmilesError::UnexpectedChar {
            pos: start,
            ch: '%',
        });
    }
    let d1 = (chars[i] as u16) - b'0' as u16;
    let d2 = (chars[i + 1] as u16) - b'0' as u16;
    let digit = d1 * 10 + d2;

    Ok((None, digit, i + 2))
}

fn parse_bracket_atom(chars: &[char], start: usize) -> Result<(AtomToken, usize), SmilesError> {
    let mut i = start + 1; // skip '['

    let isotope = parse_isotope(chars, &mut i);

    let (element, is_aromatic) = parse_bracket_element(chars, &mut i, start)?;

    let chirality = parse_chirality(chars, &mut i);

    let hcount = parse_hcount(chars, &mut i);

    let charge = parse_charge(chars, &mut i, start)?;

    let atom_class = parse_atom_class(chars, &mut i);

    if i >= chars.len() || chars[i] != ']' {
        return Err(SmilesError::UnclosedBracket { pos: start });
    }
    i += 1; // skip ']'

    Ok((
        AtomToken {
            element,
            is_aromatic,
            isotope,
            chirality,
            hcount: Some(hcount.unwrap_or(0)),
            charge,
            atom_class,
            is_bracket: true,
            pos: start,
        },
        i,
    ))
}

fn parse_isotope(chars: &[char], i: &mut usize) -> u16 {
    let mut val: u16 = 0;
    let mut found = false;
    while *i < chars.len() && chars[*i].is_ascii_digit() {
        found = true;
        val = val * 10 + (chars[*i] as u16 - b'0' as u16);
        *i += 1;
    }
    if found {
        val
    } else {
        0
    }
}

fn parse_bracket_element(
    chars: &[char],
    i: &mut usize,
    bracket_start: usize,
) -> Result<(Element, bool), SmilesError> {
    if *i >= chars.len() {
        return Err(SmilesError::UnclosedBracket {
            pos: bracket_start,
        });
    }

    let aromatic_map: &[(&str, Element)] = &[
        ("se", Element::Se),
        ("te", Element::Te),
        ("b", Element::B),
        ("c", Element::C),
        ("n", Element::N),
        ("o", Element::O),
        ("p", Element::P),
        ("s", Element::S),
    ];

    for &(pat, elem) in aromatic_map {
        if *i + pat.len() <= chars.len() {
            let slice: String = chars[*i..*i + pat.len()].iter().collect();
            if slice == pat {
                let after = *i + pat.len();
                let next_is_lower = after < chars.len() && chars[after].is_ascii_lowercase();
                if !next_is_lower || pat.len() == 2 {
                    *i += pat.len();
                    return Ok((elem, true));
                }
            }
        }
    }

    // Try two-char uppercase element first, then one-char
    if *i + 1 < chars.len() && chars[*i].is_ascii_uppercase() && chars[*i + 1].is_ascii_lowercase()
    {
        let sym: String = chars[*i..=*i + 1].iter().collect();
        if let Some(e) = Element::from_symbol(&sym) {
            *i += 2;
            return Ok((e, false));
        }
    }

    if chars[*i].is_ascii_uppercase() {
        let sym: String = chars[*i..=*i].iter().collect();
        if let Some(e) = Element::from_symbol(&sym) {
            *i += 1;
            return Ok((e, false));
        }
    }

    Err(SmilesError::InvalidElement {
        pos: *i,
        text: chars.get(*i).map(|c| c.to_string()).unwrap_or_default(),
    })
}

fn parse_chirality(chars: &[char], i: &mut usize) -> ChiralityToken {
    if *i < chars.len() && chars[*i] == '@' {
        *i += 1;
        if *i < chars.len() && chars[*i] == '@' {
            *i += 1;
            ChiralityToken::Clockwise
        } else {
            ChiralityToken::CounterClockwise
        }
    } else {
        ChiralityToken::None
    }
}

fn parse_hcount(chars: &[char], i: &mut usize) -> Option<u8> {
    if *i < chars.len() && chars[*i] == 'H' {
        *i += 1;
        let mut count: u8 = 1;
        if *i < chars.len() && chars[*i].is_ascii_digit() {
            count = chars[*i] as u8 - b'0';
            *i += 1;
        }
        Some(count)
    } else {
        None
    }
}

fn parse_charge(chars: &[char], i: &mut usize, bracket_start: usize) -> Result<i8, SmilesError> {
    if *i >= chars.len() {
        return Ok(0);
    }

    match chars[*i] {
        '+' => {
            *i += 1;
            if *i < chars.len() && chars[*i] == '+' {
                let mut count: i8 = 1;
                while *i < chars.len() && chars[*i] == '+' {
                    count = count
                        .checked_add(1)
                        .ok_or(SmilesError::InvalidCharge { pos: bracket_start })?;
                    *i += 1;
                }
                Ok(count)
            } else if *i < chars.len() && chars[*i].is_ascii_digit() {
                let mut val: i8 = 0;
                while *i < chars.len() && chars[*i].is_ascii_digit() {
                    val = val
                        .checked_mul(10)
                        .and_then(|v| v.checked_add((chars[*i] as i8) - b'0' as i8))
                        .ok_or(SmilesError::InvalidCharge { pos: bracket_start })?;
                    *i += 1;
                }
                Ok(val)
            } else {
                Ok(1)
            }
        }
        '-' => {
            *i += 1;
            if *i < chars.len() && chars[*i] == '-' {
                let mut count: i8 = -1;
                while *i < chars.len() && chars[*i] == '-' {
                    count = count
                        .checked_sub(1)
                        .ok_or(SmilesError::InvalidCharge { pos: bracket_start })?;
                    *i += 1;
                }
                Ok(count)
            } else if *i < chars.len() && chars[*i].is_ascii_digit() {
                let mut val: i8 = 0;
                while *i < chars.len() && chars[*i].is_ascii_digit() {
                    val = val
                        .checked_mul(10)
                        .and_then(|v| v.checked_add((chars[*i] as i8) - b'0' as i8))
                        .ok_or(SmilesError::InvalidCharge { pos: bracket_start })?;
                    *i += 1;
                }
                Ok(-val)
            } else {
                Ok(-1)
            }
        }
        _ => Ok(0),
    }
}

fn parse_atom_class(chars: &[char], i: &mut usize) -> u16 {
    if *i < chars.len() && chars[*i] == ':' {
        *i += 1;
        let mut val: u16 = 0;
        while *i < chars.len() && chars[*i].is_ascii_digit() {
            val = val * 10 + (chars[*i] as u16 - b'0' as u16);
            *i += 1;
        }
        val
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenize_methane() {
        let tokens = tokenize("C").unwrap();
        assert_eq!(tokens.len(), 1);
        match &tokens[0] {
            Token::Atom(a) => {
                assert_eq!(a.element, Element::C);
                assert!(!a.is_bracket);
                assert!(!a.is_aromatic);
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn tokenize_ethene() {
        let tokens = tokenize("C=C").unwrap();
        assert_eq!(tokens.len(), 3);
    }

    #[test]
    fn tokenize_bracket_atom() {
        let tokens = tokenize("[NH4+]").unwrap();
        assert_eq!(tokens.len(), 1);
        match &tokens[0] {
            Token::Atom(a) => {
                assert_eq!(a.element, Element::N);
                assert!(a.is_bracket);
                assert_eq!(a.hcount, Some(4));
                assert_eq!(a.charge, 1);
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn tokenize_isotope() {
        let tokens = tokenize("[13C]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => {
                assert_eq!(a.isotope, 13);
                assert_eq!(a.element, Element::C);
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn tokenize_ring_closure() {
        let tokens = tokenize("C1CC1").unwrap();
        assert_eq!(tokens.len(), 5);
        assert!(matches!(
            &tokens[1],
            Token::RingClosure { digit: 1, .. }
        ));
    }

    #[test]
    fn tokenize_percent_ring() {
        let tokens = tokenize("C%10CC%10").unwrap();
        assert!(matches!(
            &tokens[1],
            Token::RingClosure { digit: 10, .. }
        ));
    }

    #[test]
    fn tokenize_chirality() {
        let tokens = tokenize("[C@@H](F)(Cl)Br").unwrap();
        match &tokens[0] {
            Token::Atom(a) => {
                assert_eq!(a.chirality, ChiralityToken::Clockwise);
                assert_eq!(a.hcount, Some(1));
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn tokenize_aromatic() {
        let tokens = tokenize("c1ccccc1").unwrap();
        assert_eq!(tokens.len(), 8);
        match &tokens[0] {
            Token::Atom(a) => {
                assert!(a.is_aromatic);
                assert_eq!(a.element, Element::C);
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn bracket_aromatic_se() {
        let tokens = tokenize("[se]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => {
                assert!(a.is_aromatic);
                assert_eq!(a.element, Element::Se);
            }
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn negative_charge_variants() {
        let tokens = tokenize("[O-]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => assert_eq!(a.charge, -1),
            _ => panic!("expected atom"),
        }

        let tokens = tokenize("[O-2]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => assert_eq!(a.charge, -2),
            _ => panic!("expected atom"),
        }

        let tokens = tokenize("[O--]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => assert_eq!(a.charge, -2),
            _ => panic!("expected atom"),
        }
    }

    #[test]
    fn atom_class() {
        let tokens = tokenize("[C:1]").unwrap();
        match &tokens[0] {
            Token::Atom(a) => assert_eq!(a.atom_class, 1),
            _ => panic!("expected atom"),
        }
    }
}
