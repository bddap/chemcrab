use std::collections::{HashMap, HashSet, VecDeque};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, SmilesBond, SmilesBondOrder};
use crate::kekulize::kekulize;
use crate::mol::Mol;
use crate::smarts::{get_smarts_matches, AtomExpr, BondExpr};

use super::error::ReactionError;
use super::Reaction;

const MAX_COMBINATIONS: usize = 1000;

impl Reaction {
    /// Apply this reaction to a set of reactant molecules.
    ///
    /// Returns one result per match combination. Each result is a `Vec`
    /// of product molecules (one per product template). Returns an empty
    /// `Vec` if no template matches the input.
    pub fn run(
        &self,
        reactants: &[&Mol<Atom, Bond>],
    ) -> Result<Vec<Vec<Mol<Atom, Bond>>>, ReactionError> {
        if reactants.len() != self.reactant_templates.len() {
            return Err(ReactionError::WrongReactantCount {
                expected: self.reactant_templates.len(),
                got: reactants.len(),
            });
        }

        let per_template_matches: Vec<Vec<HashMap<NodeIndex, NodeIndex>>> = self
            .reactant_templates
            .iter()
            .zip(reactants.iter())
            .map(|(tmpl, mol)| {
                get_smarts_matches(mol, tmpl)
                    .into_iter()
                    .map(|mapping| mapping.into_iter().collect())
                    .collect()
            })
            .collect();

        if per_template_matches.iter().any(|m| m.is_empty()) {
            return Ok(Vec::new());
        }

        let combinations = cartesian_product(&per_template_matches, MAX_COMBINATIONS)?;

        let mut results = Vec::new();
        for combo in &combinations {
            let products = self.generate_products(combo, reactants)?;
            results.push(products);
        }

        Ok(results)
    }

    fn generate_products(
        &self,
        match_combo: &[&HashMap<NodeIndex, NodeIndex>],
        reactants: &[&Mol<Atom, Bond>],
    ) -> Result<Vec<Mol<Atom, Bond>>, ReactionError> {
        let reactant_atom_map = build_reactant_atom_map(&self.reactant_templates, match_combo)?;

        let matched_target_atoms = collect_matched_atoms(match_combo);

        let reactant_template_bonds = collect_template_bonds(&self.reactant_templates);

        self.product_templates
            .iter()
            .map(|product_tmpl| {
                generate_single_product(
                    product_tmpl,
                    &reactant_atom_map,
                    &matched_target_atoms,
                    &reactant_template_bonds,
                    reactants,
                )
            })
            .collect()
    }
}

fn build_reactant_atom_map(
    reactant_templates: &[Mol<AtomExpr, BondExpr>],
    match_combo: &[&HashMap<NodeIndex, NodeIndex>],
) -> Result<HashMap<u16, (usize, NodeIndex)>, ReactionError> {
    let mut map = HashMap::new();
    for (ri, tmpl) in reactant_templates.iter().enumerate() {
        let mapping = match_combo[ri];
        for q_idx in tmpl.atoms() {
            if let Some(map_num) = extract_atom_map_num(tmpl.atom(q_idx)) {
                if let Some(&t_idx) = mapping.get(&q_idx) {
                    if map.contains_key(&map_num) {
                        return Err(ReactionError::DuplicateAtomMap { map_num });
                    }
                    map.insert(map_num, (ri, t_idx));
                }
            }
        }
    }
    Ok(map)
}

fn collect_matched_atoms(
    match_combo: &[&HashMap<NodeIndex, NodeIndex>],
) -> Vec<HashSet<NodeIndex>> {
    match_combo
        .iter()
        .map(|mapping| mapping.values().copied().collect())
        .collect()
}

fn collect_template_bonds(templates: &[Mol<AtomExpr, BondExpr>]) -> Vec<HashSet<(u16, u16)>> {
    templates
        .iter()
        .map(|tmpl| {
            let mut bonds = HashSet::new();
            for edge in tmpl.bonds() {
                if let Some((a, b)) = tmpl.bond_endpoints(edge) {
                    if let (Some(ma), Some(mb)) = (
                        extract_atom_map_num(tmpl.atom(a)),
                        extract_atom_map_num(tmpl.atom(b)),
                    ) {
                        bonds.insert((ma.min(mb), ma.max(mb)));
                    }
                }
            }
            bonds
        })
        .collect()
}

fn bond_to_smiles_bond(bond: &Bond) -> SmilesBond {
    SmilesBond {
        order: match bond.order {
            BondOrder::Single => SmilesBondOrder::Single,
            BondOrder::Double => SmilesBondOrder::Double,
            BondOrder::Triple => SmilesBondOrder::Triple,
        },
    }
}

fn generate_single_product(
    product_tmpl: &Mol<AtomExpr, BondExpr>,
    reactant_atom_map: &HashMap<u16, (usize, NodeIndex)>,
    matched_target_atoms: &[HashSet<NodeIndex>],
    reactant_template_bonds: &[HashSet<(u16, u16)>],
    reactants: &[&Mol<Atom, Bond>],
) -> Result<Mol<Atom, Bond>, ReactionError> {
    let mut product: Mol<Atom, SmilesBond> = Mol::new();

    let mut product_node_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    let mut map_num_to_product_node: HashMap<u16, NodeIndex> = HashMap::new();
    let mut carried_atom_map: HashMap<(usize, NodeIndex), NodeIndex> = HashMap::new();

    for p_idx in product_tmpl.atoms() {
        let p_expr = product_tmpl.atom(p_idx);
        let map_num = extract_atom_map_num(p_expr);

        let atom = if let Some(mn) = map_num {
            if let Some(&(ri, t_idx)) = reactant_atom_map.get(&mn) {
                build_mapped_atom(reactants[ri].atom(t_idx), p_expr)
            } else {
                build_unmapped_atom(p_expr)
            }
        } else {
            build_unmapped_atom(p_expr)
        };

        let new_idx = product.add_atom(atom);
        product_node_map.insert(p_idx, new_idx);
        if let Some(mn) = map_num {
            map_num_to_product_node.insert(mn, new_idx);
        }
    }

    for edge in product_tmpl.bonds() {
        if let Some((a, b)) = product_tmpl.bond_endpoints(edge) {
            let bond_expr = product_tmpl.bond(edge);
            let bond = smiles_bond_from_expr(bond_expr);
            let pa = product_node_map[&a];
            let pb = product_node_map[&b];
            product.add_bond(pa, pb, bond);
        }
    }

    let product_template_mapped_pairs = collect_mapped_bond_pairs(product_tmpl);

    for p_idx in product_tmpl.atoms() {
        let p_expr = product_tmpl.atom(p_idx);
        let map_num = match extract_atom_map_num(p_expr) {
            Some(mn) => mn,
            None => continue,
        };
        let &(ri, t_idx) = match reactant_atom_map.get(&map_num) {
            Some(v) => v,
            None => continue,
        };

        let product_node = product_node_map[&p_idx];
        let reactant_mol = reactants[ri];
        let matched = &matched_target_atoms[ri];

        for neighbor in reactant_mol.neighbors(t_idx) {
            if matched.contains(&neighbor) {
                let neighbor_map = find_map_num_for_target(
                    &self_reactant_template_ref(reactant_atom_map, ri),
                    neighbor,
                );

                if let Some(n_mn) = neighbor_map {
                    let pair = (map_num.min(n_mn), map_num.max(n_mn));
                    let in_reactant_template = reactant_template_bonds
                        .get(ri)
                        .is_some_and(|bonds| bonds.contains(&pair));
                    let in_product_template = product_template_mapped_pairs.contains(&pair);

                    if in_reactant_template && !in_product_template {
                        continue;
                    }

                    if !in_reactant_template && !in_product_template {
                        if let Some(&neighbor_product_node) = map_num_to_product_node.get(&n_mn) {
                            if product
                                .bond_between(product_node, neighbor_product_node)
                                .is_none()
                            {
                                if let Some(edge) = reactant_mol.bond_between(t_idx, neighbor) {
                                    let bond = bond_to_smiles_bond(reactant_mol.bond(edge));
                                    product.add_bond(product_node, neighbor_product_node, bond);
                                }
                            }
                        }
                    }
                }
                continue;
            }

            carry_substituents(
                &mut product,
                reactant_mol,
                matched,
                t_idx,
                neighbor,
                product_node,
                &mut carried_atom_map,
                ri,
            );
        }
    }

    Ok(kekulize(product)?)
}

#[allow(clippy::too_many_arguments)]
fn carry_substituents(
    product: &mut Mol<Atom, SmilesBond>,
    reactant_mol: &Mol<Atom, Bond>,
    matched: &HashSet<NodeIndex>,
    anchor: NodeIndex,
    start_neighbor: NodeIndex,
    product_anchor: NodeIndex,
    carried_atom_map: &mut HashMap<(usize, NodeIndex), NodeIndex>,
    reactant_idx: usize,
) {
    let mut queue = VecDeque::new();

    let key = (reactant_idx, start_neighbor);
    if carried_atom_map.contains_key(&key) {
        let carried_node = carried_atom_map[&key];
        if product.bond_between(product_anchor, carried_node).is_none() {
            if let Some(edge) = reactant_mol.bond_between(anchor, start_neighbor) {
                product.add_bond(
                    product_anchor,
                    carried_node,
                    bond_to_smiles_bond(reactant_mol.bond(edge)),
                );
            }
        }
        return;
    }

    let new_node = product.add_atom(reactant_mol.atom(start_neighbor).clone());
    carried_atom_map.insert(key, new_node);

    if let Some(edge) = reactant_mol.bond_between(anchor, start_neighbor) {
        product.add_bond(
            product_anchor,
            new_node,
            bond_to_smiles_bond(reactant_mol.bond(edge)),
        );
    }

    queue.push_back((start_neighbor, new_node));

    while let Some((r_node, p_node)) = queue.pop_front() {
        for nb in reactant_mol.neighbors(r_node) {
            if matched.contains(&nb) {
                continue;
            }
            let nb_key = (reactant_idx, nb);
            if carried_atom_map.contains_key(&nb_key) {
                continue;
            }

            let nb_product = product.add_atom(reactant_mol.atom(nb).clone());
            carried_atom_map.insert(nb_key, nb_product);

            if let Some(edge) = reactant_mol.bond_between(r_node, nb) {
                product.add_bond(
                    p_node,
                    nb_product,
                    bond_to_smiles_bond(reactant_mol.bond(edge)),
                );
            }

            queue.push_back((nb, nb_product));
        }
    }
}

fn collect_mapped_bond_pairs(tmpl: &Mol<AtomExpr, BondExpr>) -> HashSet<(u16, u16)> {
    let mut pairs = HashSet::new();
    for edge in tmpl.bonds() {
        if let Some((a, b)) = tmpl.bond_endpoints(edge) {
            if let (Some(ma), Some(mb)) = (
                extract_atom_map_num(tmpl.atom(a)),
                extract_atom_map_num(tmpl.atom(b)),
            ) {
                pairs.insert((ma.min(mb), ma.max(mb)));
            }
        }
    }
    pairs
}

fn self_reactant_template_ref(
    reactant_atom_map: &HashMap<u16, (usize, NodeIndex)>,
    reactant_idx: usize,
) -> HashMap<u16, NodeIndex> {
    reactant_atom_map
        .iter()
        .filter(|(_, &(ri, _))| ri == reactant_idx)
        .map(|(&mn, &(_, t_idx))| (mn, t_idx))
        .collect()
}

fn find_map_num_for_target(
    map_num_to_target: &HashMap<u16, NodeIndex>,
    target_node: NodeIndex,
) -> Option<u16> {
    map_num_to_target
        .iter()
        .find(|(_, &t)| t == target_node)
        .map(|(&mn, _)| mn)
}

fn build_mapped_atom(reactant_atom: &Atom, product_expr: &AtomExpr) -> Atom {
    let mut atom = reactant_atom.clone();
    apply_expr_to_atom(&mut atom, product_expr);
    atom
}

fn build_unmapped_atom(expr: &AtomExpr) -> Atom {
    let mut atom = Atom::default();
    apply_expr_to_atom(&mut atom, expr);
    atom
}

fn apply_expr_to_atom(atom: &mut Atom, expr: &AtomExpr) {
    match expr {
        AtomExpr::Element {
            atomic_num,
            aromatic,
        } => {
            atom.atomic_num = *atomic_num;
            if let Some(arom) = aromatic {
                atom.is_aromatic = *arom;
            }
        }
        AtomExpr::Charge(c) => {
            atom.formal_charge = *c;
        }
        AtomExpr::Isotope(i) => {
            atom.isotope = *i;
        }
        AtomExpr::TotalHCount(h) => {
            atom.hydrogen_count = *h;
        }
        AtomExpr::And(parts) => {
            for p in parts {
                apply_expr_to_atom(atom, p);
            }
        }
        _ => {}
    }
}

fn smiles_bond_from_expr(expr: &BondExpr) -> SmilesBond {
    let order = match expr {
        BondExpr::Single | BondExpr::SingleOrAromatic => SmilesBondOrder::Single,
        BondExpr::Double => SmilesBondOrder::Double,
        BondExpr::Triple => SmilesBondOrder::Triple,
        BondExpr::Aromatic => SmilesBondOrder::Aromatic,
        BondExpr::And(_)
        | BondExpr::Or(_)
        | BondExpr::Not(_)
        | BondExpr::Ring
        | BondExpr::True
        | BondExpr::Up
        | BondExpr::Down => SmilesBondOrder::Single,
    };
    SmilesBond { order }
}

/// Extract the atom map number from a SMARTS atom expression, if present.
pub fn extract_atom_map_num(expr: &AtomExpr) -> Option<u16> {
    match expr {
        AtomExpr::AtomMapClass(0) => None,
        AtomExpr::AtomMapClass(n) => Some(*n),
        AtomExpr::And(parts) => parts.iter().find_map(extract_atom_map_num),
        _ => None,
    }
}

fn cartesian_product<'a, T>(
    sets: &'a [Vec<T>],
    max: usize,
) -> Result<Vec<Vec<&'a T>>, ReactionError> {
    let mut result: Vec<Vec<&'a T>> = vec![vec![]];
    for set in sets {
        let mut new_result = Vec::new();
        for combo in &result {
            for item in set {
                let mut new_combo = combo.clone();
                new_combo.push(item);
                new_result.push(new_combo);
                if new_result.len() > max {
                    return Err(ReactionError::TooManyCombinations);
                }
            }
        }
        result = new_result;
    }
    Ok(result)
}
