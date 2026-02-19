use crate::element::Element;
use crate::mol::Mol;

use super::error::SmartsError;
use super::query::{AtomExpr, BondExpr};

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    input: &'a str,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            chars: input.chars().collect(),
            pos: 0,
            input,
        }
    }

    fn advance(&mut self) -> Option<char> {
        let ch = self.chars.get(self.pos).copied()?;
        self.pos += 1;
        Some(ch)
    }

    fn expect(&mut self, ch: char) -> Result<(), SmartsError> {
        match self.advance() {
            Some(c) if c == ch => Ok(()),
            Some(c) => Err(SmartsError::UnexpectedChar { pos: self.pos - 1, ch: c }),
            None => Err(SmartsError::InvalidSmarts {
                pos: self.pos,
                msg: format!("expected '{ch}', got end of input"),
            }),
        }
    }

    fn parse_number(&mut self) -> Option<u32> {
        let start = self.pos;
        while self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        if self.pos > start {
            let s: String = self.chars[start..self.pos].iter().collect();
            s.parse().ok()
        } else {
            None
        }
    }

    fn parse_smarts(&mut self) -> Result<Mol<AtomExpr, BondExpr>, SmartsError> {
        let mut mol = Mol::new();
        let mut stack = Vec::new();
        let mut current = None;
        let mut pending_bond: Option<BondExpr> = None;
        let mut ring_map: std::collections::HashMap<u16, (petgraph::graph::NodeIndex, BondExpr)> =
            std::collections::HashMap::new();
        let mut insertion_order: std::collections::HashMap<usize, Vec<petgraph::graph::NodeIndex>> =
            std::collections::HashMap::new();

        while self.pos < self.chars.len() {
            let ch = self.chars[self.pos];

            match ch {
                '[' => {
                    let atom_expr = self.parse_bracket_atom()?;
                    let idx = mol.add_atom(atom_expr);
                    if let Some(prev) = current {
                        let bond = pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic);
                        mol.add_bond(prev, idx, bond);
                        insertion_order.entry(prev.index()).or_default().push(idx);
                        insertion_order.entry(idx.index()).or_default().push(prev);
                    }
                    current = Some(idx);
                }
                '(' => {
                    self.pos += 1;
                    if let Some(cur) = current {
                        stack.push((cur, pending_bond.take()));
                    } else {
                        return Err(SmartsError::UnmatchedParen { pos: self.pos - 1 });
                    }
                    continue;
                }
                ')' => {
                    self.pos += 1;
                    if let Some((prev, saved_bond)) = stack.pop() {
                        current = Some(prev);
                        pending_bond = saved_bond;
                    } else {
                        return Err(SmartsError::UnmatchedParen { pos: self.pos - 1 });
                    }
                    continue;
                }
                '.' => {
                    self.pos += 1;
                    current = None;
                    pending_bond = None;
                    continue;
                }
                '-' | '=' | '#' | '~' | ':' | '/' | '\\' | '@' => {
                    if pending_bond.is_some() {
                        return Err(SmartsError::InvalidSmarts {
                            pos: self.pos,
                            msg: "consecutive bond expressions".into(),
                        });
                    }
                    pending_bond = Some(self.parse_bond_expr()?);
                    continue;
                }
                '0'..='9' | '%' => {
                    let (digit, _bond_pos) = self.parse_ring_closure()?;
                    if let Some(cur) = current {
                        if let Some((other, saved_bond)) = ring_map.remove(&digit) {
                            let bond = pending_bond
                                .take()
                                .or(Some(saved_bond))
                                .unwrap_or(BondExpr::SingleOrAromatic);
                            mol.add_bond(cur, other, bond);
                            insertion_order.entry(cur.index()).or_default().push(other);
                            insertion_order.entry(other.index()).or_default().push(cur);
                        } else {
                            ring_map.insert(
                                digit,
                                (cur, pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic)),
                            );
                        }
                    } else {
                        return Err(SmartsError::InvalidSmarts {
                            pos: self.pos,
                            msg: "ring closure without preceding atom".into(),
                        });
                    }
                    continue;
                }
                _ => {
                    let atom_expr = self.parse_bare_atom()?;
                    let idx = mol.add_atom(atom_expr);
                    if let Some(prev) = current {
                        let bond = pending_bond.take().unwrap_or(BondExpr::SingleOrAromatic);
                        mol.add_bond(prev, idx, bond);
                        insertion_order.entry(prev.index()).or_default().push(idx);
                        insertion_order.entry(idx.index()).or_default().push(prev);
                    }
                    current = Some(idx);
                }
            }
        }

        if !stack.is_empty() {
            return Err(SmartsError::UnmatchedParen { pos: self.pos });
        }

        if let Some((&digit, _)) = ring_map.iter().next() {
            return Err(SmartsError::UnclosedRing { digit });
        }

        normalize_chirality(&mut mol, &insertion_order);

        Ok(mol)
    }

    fn parse_ring_closure(&mut self) -> Result<(u16, usize), SmartsError> {
        let start = self.pos;
        if self.chars[self.pos] == '%' {
            self.pos += 1;
            if self.pos + 1 < self.chars.len()
                && self.chars[self.pos].is_ascii_digit()
                && self.chars[self.pos + 1].is_ascii_digit()
            {
                let d1 = self.chars[self.pos].to_digit(10).unwrap() as u16;
                let d2 = self.chars[self.pos + 1].to_digit(10).unwrap() as u16;
                self.pos += 2;
                Ok((d1 * 10 + d2, start))
            } else {
                Err(SmartsError::InvalidSmarts {
                    pos: start,
                    msg: "expected two digits after %".into(),
                })
            }
        } else {
            let d = self.chars[self.pos].to_digit(10).unwrap() as u16;
            self.pos += 1;
            Ok((d, start))
        }
    }

    fn parse_bond_expr(&mut self) -> Result<BondExpr, SmartsError> {
        let ch = self.chars[self.pos];
        self.pos += 1;
        match ch {
            '-' => Ok(BondExpr::Single),
            '=' => Ok(BondExpr::Double),
            '#' => Ok(BondExpr::Triple),
            '~' => Ok(BondExpr::True),
            ':' => Ok(BondExpr::Aromatic),
            '/' => Ok(BondExpr::Up),
            '\\' => Ok(BondExpr::Down),
            '@' => Ok(BondExpr::Ring),
            _ => Err(SmartsError::UnexpectedChar { pos: self.pos - 1, ch }),
        }
    }

    fn parse_bare_atom(&mut self) -> Result<AtomExpr, SmartsError> {
        let start = self.pos;
        let ch = self.chars[self.pos];

        if ch == '*' {
            self.pos += 1;
            return Ok(AtomExpr::True);
        }
        if ch == 'A' {
            if self.pos + 1 < self.chars.len() {
                let next = self.chars[self.pos + 1];
                if next.is_ascii_lowercase() {
                    return self.parse_bare_element();
                }
            }
            self.pos += 1;
            return Ok(AtomExpr::Aliphatic);
        }
        if ch == 'a' {
            self.pos += 1;
            return Ok(AtomExpr::Aromatic);
        }

        self.parse_bare_element().map_err(|_| SmartsError::UnexpectedChar {
            pos: start,
            ch,
        })
    }

    fn parse_bare_element(&mut self) -> Result<AtomExpr, SmartsError> {
        let start = self.pos;
        let ch = self.chars[self.pos];

        let aromatic_atoms = [
            ('c', 6u8),
            ('n', 7),
            ('o', 8),
            ('s', 16),
            ('p', 15),
        ];

        for &(sym, num) in &aromatic_atoms {
            if ch == sym {
                self.pos += 1;
                return Ok(AtomExpr::Element {
                    atomic_num: num,
                    aromatic: Some(true),
                });
            }
        }

        if ch == 'H' {
            self.pos += 1;
            return Ok(AtomExpr::Element {
                atomic_num: 1,
                aromatic: Some(false),
            });
        }

        if ch.is_ascii_uppercase() {
            let mut symbol = String::new();
            symbol.push(ch);
            self.pos += 1;
            if self.pos < self.chars.len() && self.chars[self.pos].is_ascii_lowercase() {
                symbol.push(self.chars[self.pos]);
                if let Some(elem) = Element::from_symbol(&symbol) {
                    self.pos += 1;
                    return Ok(AtomExpr::Element {
                        atomic_num: elem.atomic_num(),
                        aromatic: Some(false),
                    });
                }
                symbol.pop();
            }
            if let Some(elem) = Element::from_symbol(&symbol) {
                return Ok(AtomExpr::Element {
                    atomic_num: elem.atomic_num(),
                    aromatic: Some(false),
                });
            }
            self.pos = start;
        }

        Err(SmartsError::UnexpectedChar { pos: start, ch })
    }

    fn parse_bracket_atom(&mut self) -> Result<AtomExpr, SmartsError> {
        let bracket_start = self.pos;
        self.expect('[')?;

        let expr = self.parse_semicolon_expr()?;

        if self.pos >= self.chars.len() || self.chars[self.pos] != ']' {
            return Err(SmartsError::UnclosedBracket { pos: bracket_start });
        }
        self.pos += 1;

        Ok(expr)
    }

    fn parse_semicolon_expr(&mut self) -> Result<AtomExpr, SmartsError> {
        let mut parts = vec![self.parse_comma_expr()?];
        while self.pos < self.chars.len() && self.chars[self.pos] == ';' {
            self.pos += 1;
            parts.push(self.parse_comma_expr()?);
        }
        Ok(flatten_and(parts))
    }

    fn parse_comma_expr(&mut self) -> Result<AtomExpr, SmartsError> {
        let mut parts = vec![self.parse_high_and_expr()?];
        while self.pos < self.chars.len() && self.chars[self.pos] == ',' {
            self.pos += 1;
            parts.push(self.parse_high_and_expr()?);
        }
        Ok(flatten_or(parts))
    }

    fn parse_high_and_expr(&mut self) -> Result<AtomExpr, SmartsError> {
        let mut parts = Vec::new();
        loop {
            if self.pos >= self.chars.len() {
                break;
            }
            let ch = self.chars[self.pos];
            if ch == ']' || ch == ',' || ch == ';' {
                break;
            }
            if ch == '&' {
                self.pos += 1;
                continue;
            }
            parts.push(self.parse_not_expr()?);
        }
        if parts.is_empty() {
            Ok(AtomExpr::True)
        } else {
            Ok(flatten_and(parts))
        }
    }

    fn parse_not_expr(&mut self) -> Result<AtomExpr, SmartsError> {
        if self.pos < self.chars.len() && self.chars[self.pos] == '!' {
            self.pos += 1;
            let inner = self.parse_primitive()?;
            Ok(AtomExpr::Not(Box::new(inner)))
        } else {
            self.parse_primitive()
        }
    }

    fn parse_primitive(&mut self) -> Result<AtomExpr, SmartsError> {
        if self.pos >= self.chars.len() {
            return Err(SmartsError::InvalidSmarts {
                pos: self.pos,
                msg: "expected atom primitive".into(),
            });
        }

        let ch = self.chars[self.pos];

        match ch {
            '*' => {
                self.pos += 1;
                Ok(AtomExpr::True)
            }
            'A' => {
                if self.pos + 1 < self.chars.len() && self.chars[self.pos + 1].is_ascii_lowercase()
                {
                    self.parse_bracket_element()
                } else {
                    self.pos += 1;
                    Ok(AtomExpr::Aliphatic)
                }
            }
            'a' => {
                if self.pos + 1 < self.chars.len() && self.chars[self.pos + 1] == 's' {
                    self.parse_bracket_element()
                } else {
                    self.pos += 1;
                    Ok(AtomExpr::Aromatic)
                }
            }
            '#' => {
                self.pos += 1;
                let num = self.parse_number().ok_or(SmartsError::InvalidAtomicNum { pos: self.pos })?;
                if num == 0 || num > 118 {
                    return Err(SmartsError::InvalidAtomicNum { pos: self.pos });
                }
                Ok(AtomExpr::Element {
                    atomic_num: num as u8,
                    aromatic: None,
                })
            }
            'D' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::Degree(n as u8))
            }
            'v' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::Valence(n as u8))
            }
            'X' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::Connectivity(n as u8))
            }
            'H' => {
                if self.is_hydrogen_element_context() {
                    self.parse_bracket_element()
                } else {
                    self.pos += 1;
                    let n = self.parse_number().unwrap_or(1);
                    Ok(AtomExpr::TotalHCount(n as u8))
                }
            }
            'h' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::ImplicitHCount(n as u8))
            }
            'R' => {
                self.pos += 1;
                if self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
                    let n = self.parse_number().unwrap();
                    if n == 0 {
                        Ok(AtomExpr::NotInRing)
                    } else {
                        Ok(AtomExpr::RingMembership(n as u8))
                    }
                } else {
                    Ok(AtomExpr::InRing)
                }
            }
            'r' => {
                self.pos += 1;
                if self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
                    let n = self.parse_number().unwrap();
                    Ok(AtomExpr::SmallestRingSize(n as u8))
                } else {
                    Ok(AtomExpr::InRing)
                }
            }
            'x' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::RingBondCount(n as u8))
            }
            '@' => {
                self.pos += 1;
                if self.pos < self.chars.len() && self.chars[self.pos] == '@' {
                    self.pos += 1;
                    Ok(AtomExpr::Chirality(crate::atom::Chirality::Cw))
                } else {
                    Ok(AtomExpr::Chirality(crate::atom::Chirality::Ccw))
                }
            }
            '+' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::Charge(n as i8))
            }
            '-' => {
                self.pos += 1;
                let n = self.parse_number().unwrap_or(1);
                Ok(AtomExpr::Charge(-(n as i8)))
            }
            '$' => {
                self.pos += 1;
                if self.pos >= self.chars.len() || self.chars[self.pos] != '(' {
                    return Err(SmartsError::UnclosedRecursive { pos: self.pos });
                }
                self.pos += 1;
                let inner_start = self.pos;
                let inner = self.extract_balanced_parens(inner_start)?;
                let inner_mol = parse(inner)?;
                Ok(AtomExpr::Recursive(inner_mol))
            }
            _ if ch.is_ascii_uppercase() || ch.is_ascii_lowercase() => {
                self.parse_bracket_element()
            }
            _ if ch.is_ascii_digit() => {
                let n = self.parse_number().unwrap();
                Ok(AtomExpr::Isotope(n as u16))
            }
            _ => Err(SmartsError::UnexpectedChar { pos: self.pos, ch }),
        }
    }

    fn is_hydrogen_element_context(&self) -> bool {
        if self.pos >= self.chars.len() || self.chars[self.pos] != 'H' {
            return false;
        }
        let next_pos = self.pos + 1;
        if next_pos >= self.chars.len() {
            return false;
        }
        let next = self.chars[next_pos];
        if next == ']' {
            let before_h = self.bracket_content_before_pos();
            return before_h.is_empty();
        }
        false
    }

    fn bracket_content_before_pos(&self) -> &[char] {
        let mut start = self.pos;
        while start > 0 && self.chars[start - 1] != '[' {
            start -= 1;
        }
        &self.chars[start..self.pos]
    }

    fn parse_bracket_element(&mut self) -> Result<AtomExpr, SmartsError> {
        let start = self.pos;
        let ch = self.chars[self.pos];

        let aromatic_bracket = [
            ("se", 34u8),
            ("as", 33u8),
            ("c", 6u8),
            ("n", 7),
            ("o", 8),
            ("s", 16),
            ("p", 15),
            ("te", 52),
        ];

        if ch.is_ascii_lowercase() {
            for &(sym, num) in &aromatic_bracket {
                if self.matches_str(sym) {
                    self.pos += sym.len();
                    return Ok(AtomExpr::Element {
                        atomic_num: num,
                        aromatic: Some(true),
                    });
                }
            }
            return Err(SmartsError::UnexpectedChar { pos: start, ch });
        }

        let mut symbol = String::new();
        symbol.push(ch);
        self.pos += 1;

        if self.pos < self.chars.len() && self.chars[self.pos].is_ascii_lowercase() {
            let next = self.chars[self.pos];
            symbol.push(next);
            if let Some(elem) = Element::from_symbol(&symbol) {
                self.pos += 1;
                return Ok(AtomExpr::Element {
                    atomic_num: elem.atomic_num(),
                    aromatic: Some(false),
                });
            }
            symbol.pop();
        }

        if let Some(elem) = Element::from_symbol(&symbol) {
            return Ok(AtomExpr::Element {
                atomic_num: elem.atomic_num(),
                aromatic: Some(false),
            });
        }

        self.pos = start;
        Err(SmartsError::UnexpectedChar { pos: start, ch })
    }

    fn matches_str(&self, s: &str) -> bool {
        let chars: Vec<char> = s.chars().collect();
        if self.pos + chars.len() > self.chars.len() {
            return false;
        }
        for (i, &c) in chars.iter().enumerate() {
            if self.chars[self.pos + i] != c {
                return false;
            }
        }
        if self.pos + chars.len() < self.chars.len() {
            let after = self.chars[self.pos + chars.len()];
            if after.is_ascii_lowercase() && chars.len() < 2 {
                if let Some(elem) = Element::from_symbol(&format!("{}{}", s, after)) {
                    let _ = elem;
                    return false;
                }
            }
        }
        true
    }

    fn extract_balanced_parens(&mut self, _start: usize) -> Result<&str, SmartsError> {
        let start_pos = self.pos;
        let mut depth = 1;
        let begin_byte = self.chars[..self.pos]
            .iter()
            .map(|c| c.len_utf8())
            .sum::<usize>();

        while self.pos < self.chars.len() {
            let ch = self.chars[self.pos];
            match ch {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        let end_byte = begin_byte
                            + self.chars[start_pos..self.pos]
                                .iter()
                                .map(|c| c.len_utf8())
                                .sum::<usize>();
                        let inner = &self.input[begin_byte..end_byte];
                        self.pos += 1;
                        return Ok(inner);
                    }
                }
                _ => {}
            }
            self.pos += 1;
        }

        Err(SmartsError::UnclosedRecursive { pos: start_pos })
    }
}

fn flatten_and(mut parts: Vec<AtomExpr>) -> AtomExpr {
    let mut flattened = Vec::new();
    for p in parts.drain(..) {
        match p {
            AtomExpr::And(inner) => flattened.extend(inner),
            other => flattened.push(other),
        }
    }
    if flattened.len() == 1 {
        flattened.pop().unwrap()
    } else {
        AtomExpr::And(flattened)
    }
}

fn flatten_or(mut parts: Vec<AtomExpr>) -> AtomExpr {
    let mut flattened = Vec::new();
    for p in parts.drain(..) {
        match p {
            AtomExpr::Or(inner) => flattened.extend(inner),
            other => flattened.push(other),
        }
    }
    if flattened.len() == 1 {
        flattened.pop().unwrap()
    } else {
        AtomExpr::Or(flattened)
    }
}

fn normalize_chirality(
    mol: &mut Mol<AtomExpr, BondExpr>,
    insertion_order: &std::collections::HashMap<usize, Vec<petgraph::graph::NodeIndex>>,
) {
    let atoms: Vec<petgraph::graph::NodeIndex> = mol.atoms().collect();
    for idx in atoms {
        let needs_flip = {
            let expr = mol.atom(idx);
            if !has_chirality(expr) {
                continue;
            }
            let Some(input_order) = insertion_order.get(&idx.index()) else {
                continue;
            };
            let graph_order: Vec<petgraph::graph::NodeIndex> = mol.neighbors(idx).collect();
            if input_order.len() != graph_order.len() || input_order.len() < 2 {
                continue;
            }
            !parity_even(input_order, &graph_order)
        };
        if needs_flip {
            flip_chirality(mol.atom_mut(idx));
        }
    }
}

fn has_chirality(expr: &AtomExpr) -> bool {
    match expr {
        AtomExpr::Chirality(c) => *c != crate::atom::Chirality::None,
        AtomExpr::And(parts) => parts.iter().any(has_chirality),
        _ => false,
    }
}

fn flip_chirality(expr: &mut AtomExpr) {
    use crate::atom::Chirality;
    match expr {
        AtomExpr::Chirality(c) => {
            *c = match *c {
                Chirality::Cw => Chirality::Ccw,
                Chirality::Ccw => Chirality::Cw,
                Chirality::None => Chirality::None,
            };
        }
        AtomExpr::And(parts) => {
            for p in parts {
                flip_chirality(p);
            }
        }
        _ => {}
    }
}

fn parity_even(from: &[petgraph::graph::NodeIndex], to: &[petgraph::graph::NodeIndex]) -> bool {
    let n = from.len();
    let perm: Vec<usize> = from
        .iter()
        .map(|f| to.iter().position(|t| t == f).unwrap_or(0))
        .collect();

    let mut visited = vec![false; n];
    let mut swaps = 0;
    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut cycle_len = 0;
        let mut j = i;
        while !visited[j] {
            visited[j] = true;
            j = perm[j];
            cycle_len += 1;
        }
        swaps += cycle_len - 1;
    }
    swaps % 2 == 0
}

pub fn parse(input: &str) -> Result<Mol<AtomExpr, BondExpr>, SmartsError> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(SmartsError::EmptyInput);
    }
    let mut parser = Parser::new(trimmed);
    parser.parse_smarts()
}
