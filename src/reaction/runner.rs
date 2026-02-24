use std::borrow::Cow;
use std::collections::{HashMap, HashSet, VecDeque};

use petgraph::graph::NodeIndex;

use crate::atom::Atom;
use crate::bond::{AromaticBond, AromaticBondOrder, Bond, BondOrder};
use crate::element::Element;
use crate::hydrogen;
use crate::mol::Mol;
use crate::smarts::{get_smarts_matches_all_impl, query_references_hydrogen, AtomExpr, BondExpr};

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
    ) -> Result<Vec<Vec<Mol<Atom, AromaticBond>>>, ReactionError> {
        if reactants.len() != self.reactant_templates.len() {
            return Err(ReactionError::WrongReactantCount {
                expected: self.reactant_templates.len(),
                got: reactants.len(),
            });
        }

        let explicit_mols: Vec<Cow<'_, Mol<Atom, Bond>>> = self
            .reactant_templates
            .iter()
            .zip(reactants.iter())
            .map(|(tmpl, mol)| {
                if query_references_hydrogen(tmpl) {
                    Cow::Owned(hydrogen::add_hs(mol))
                } else {
                    Cow::Borrowed(*mol)
                }
            })
            .collect();

        let per_template_matches: Vec<Vec<HashMap<NodeIndex, NodeIndex>>> = self
            .reactant_templates
            .iter()
            .zip(explicit_mols.iter())
            .map(|(tmpl, mol)| {
                get_smarts_matches_all_impl(mol, tmpl)
                    .into_iter()
                    .map(|mapping| mapping.into_iter().collect())
                    .collect()
            })
            .collect();

        if per_template_matches.iter().any(|m| m.is_empty()) {
            return Ok(Vec::new());
        }

        let effective_reactants: Vec<&Mol<Atom, Bond>> =
            explicit_mols.iter().map(|c| c.as_ref()).collect();

        let combinations = cartesian_product(&per_template_matches, MAX_COMBINATIONS)?;

        let mut results = Vec::new();
        for combo in &combinations {
            let products = self.generate_products(combo, &effective_reactants)?;
            results.push(products);
        }

        Ok(results)
    }

    fn generate_products(
        &self,
        match_combo: &[&HashMap<NodeIndex, NodeIndex>],
        reactants: &[&Mol<Atom, Bond>],
    ) -> Result<Vec<Mol<Atom, AromaticBond>>, ReactionError> {
        let reactant_atom_map = build_reactant_atom_map(&self.reactant_templates, match_combo)?;

        let target_to_map_num = build_target_to_map_num(&reactant_atom_map, reactants.len());

        let matched_target_atoms = collect_matched_atoms(match_combo);

        let reactant_template_bonds = collect_template_bonds(&self.reactant_templates);

        self.product_templates
            .iter()
            .map(|product_tmpl| {
                generate_single_product(
                    product_tmpl,
                    &reactant_atom_map,
                    &target_to_map_num,
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

fn build_target_to_map_num(
    reactant_atom_map: &HashMap<u16, (usize, NodeIndex)>,
    num_reactants: usize,
) -> Vec<HashMap<NodeIndex, u16>> {
    let mut maps = vec![HashMap::new(); num_reactants];
    for (&mn, &(ri, t_idx)) in reactant_atom_map {
        maps[ri].insert(t_idx, mn);
    }
    maps
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

fn bond_to_aromatic_bond(bond: &Bond) -> AromaticBond {
    AromaticBond {
        order: AromaticBondOrder::Known(bond.order),
    }
}

fn bond_to_aromatic_bond_aromatic_aware(
    bond: &Bond,
    src_aromatic: bool,
    dst_aromatic: bool,
) -> AromaticBond {
    if src_aromatic && dst_aromatic {
        AromaticBond {
            order: AromaticBondOrder::Aromatic,
        }
    } else {
        bond_to_aromatic_bond(bond)
    }
}

fn generate_single_product(
    product_tmpl: &Mol<AtomExpr, BondExpr>,
    reactant_atom_map: &HashMap<u16, (usize, NodeIndex)>,
    target_to_map_num: &[HashMap<NodeIndex, u16>],
    matched_target_atoms: &[HashSet<NodeIndex>],
    reactant_template_bonds: &[HashSet<(u16, u16)>],
    reactants: &[&Mol<Atom, Bond>],
) -> Result<Mol<Atom, AromaticBond>, ReactionError> {
    let mut product: Mol<Atom, AromaticBond> = Mol::new();

    let mut product_node_map: HashMap<NodeIndex, NodeIndex> = HashMap::new();
    let mut map_num_to_product_node: HashMap<u16, NodeIndex> = HashMap::new();
    let mut carried_atom_map: HashMap<(usize, NodeIndex), NodeIndex> = HashMap::new();
    let mut explicit_h_atoms: HashSet<NodeIndex> = HashSet::new();

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
        if expr_has_explicit_h(p_expr) {
            explicit_h_atoms.insert(new_idx);
        }
    }

    for edge in product_tmpl.bonds() {
        if let Some((a, b)) = product_tmpl.bond_endpoints(edge) {
            let bond_expr = product_tmpl.bond(edge);
            let bond = aromatic_bond_from_expr(bond_expr);
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
        let reverse_map = &target_to_map_num[ri];

        for neighbor in reactant_mol.neighbors(t_idx) {
            if matched.contains(&neighbor) {
                let neighbor_map = reverse_map.get(&neighbor).copied();

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
                                    let src_arom = reactant_mol.atom(t_idx).is_aromatic;
                                    let dst_arom = reactant_mol.atom(neighbor).is_aromatic;
                                    let bond = bond_to_aromatic_bond_aromatic_aware(
                                        reactant_mol.bond(edge),
                                        src_arom,
                                        dst_arom,
                                    );
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

    let carried_product_nodes: HashSet<NodeIndex> = carried_atom_map.values().copied().collect();
    recalculate_hydrogen_counts(&mut product, &explicit_h_atoms, &carried_product_nodes);

    Ok(product)
}

#[allow(clippy::too_many_arguments)]
fn carry_substituents(
    product: &mut Mol<Atom, AromaticBond>,
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
                let src_arom = reactant_mol.atom(anchor).is_aromatic;
                let dst_arom = reactant_mol.atom(start_neighbor).is_aromatic;
                product.add_bond(
                    product_anchor,
                    carried_node,
                    bond_to_aromatic_bond_aromatic_aware(
                        reactant_mol.bond(edge),
                        src_arom,
                        dst_arom,
                    ),
                );
            }
        }
        return;
    }

    let new_node = product.add_atom(reactant_mol.atom(start_neighbor).clone());
    carried_atom_map.insert(key, new_node);

    if let Some(edge) = reactant_mol.bond_between(anchor, start_neighbor) {
        let src_arom = reactant_mol.atom(anchor).is_aromatic;
        let dst_arom = reactant_mol.atom(start_neighbor).is_aromatic;
        product.add_bond(
            product_anchor,
            new_node,
            bond_to_aromatic_bond_aromatic_aware(reactant_mol.bond(edge), src_arom, dst_arom),
        );
    }

    queue.push_back((start_neighbor, new_node));

    while let Some((r_node, p_node)) = queue.pop_front() {
        for nb in reactant_mol.neighbors(r_node) {
            if matched.contains(&nb) {
                continue;
            }
            let nb_key = (reactant_idx, nb);
            if let Some(&existing) = carried_atom_map.get(&nb_key) {
                if product.bond_between(p_node, existing).is_none() {
                    if let Some(edge) = reactant_mol.bond_between(r_node, nb) {
                        let src_arom = reactant_mol.atom(r_node).is_aromatic;
                        let dst_arom = reactant_mol.atom(nb).is_aromatic;
                        product.add_bond(
                            p_node,
                            existing,
                            bond_to_aromatic_bond_aromatic_aware(
                                reactant_mol.bond(edge),
                                src_arom,
                                dst_arom,
                            ),
                        );
                    }
                }
                continue;
            }

            let nb_product = product.add_atom(reactant_mol.atom(nb).clone());
            carried_atom_map.insert(nb_key, nb_product);

            if let Some(edge) = reactant_mol.bond_between(r_node, nb) {
                let src_arom = reactant_mol.atom(r_node).is_aromatic;
                let dst_arom = reactant_mol.atom(nb).is_aromatic;
                product.add_bond(
                    p_node,
                    nb_product,
                    bond_to_aromatic_bond_aromatic_aware(
                        reactant_mol.bond(edge),
                        src_arom,
                        dst_arom,
                    ),
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

fn aromatic_bond_from_expr(expr: &BondExpr) -> AromaticBond {
    let order = match expr {
        BondExpr::Single | BondExpr::SingleOrAromatic => {
            AromaticBondOrder::Known(BondOrder::Single)
        }
        BondExpr::Double => AromaticBondOrder::Known(BondOrder::Double),
        BondExpr::Triple => AromaticBondOrder::Known(BondOrder::Triple),
        BondExpr::Aromatic => AromaticBondOrder::Aromatic,
        BondExpr::And(_)
        | BondExpr::Or(_)
        | BondExpr::Not(_)
        | BondExpr::Ring
        | BondExpr::True
        | BondExpr::Up
        | BondExpr::Down => AromaticBondOrder::Known(BondOrder::Single),
    };
    AromaticBond { order }
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

fn expr_has_explicit_h(expr: &AtomExpr) -> bool {
    match expr {
        AtomExpr::TotalHCount(_) => true,
        AtomExpr::And(parts) => parts.iter().any(expr_has_explicit_h),
        _ => false,
    }
}

fn recalculate_hydrogen_counts(
    product: &mut Mol<Atom, AromaticBond>,
    explicit_h_atoms: &HashSet<NodeIndex>,
    carried_atoms: &HashSet<NodeIndex>,
) {
    let atoms: Vec<NodeIndex> = product.atoms().collect();
    for idx in atoms {
        if explicit_h_atoms.contains(&idx) || carried_atoms.contains(&idx) {
            continue;
        }

        let atom = product.atom(idx);
        let elem = match Element::from_atomic_num(atom.atomic_num) {
            Some(e) => e,
            None => {
                product.atom_mut(idx).hydrogen_count = 0;
                continue;
            }
        };

        let default_valences = elem.default_valences();
        if default_valences.is_empty() {
            product.atom_mut(idx).hydrogen_count = 0;
            continue;
        }

        let bond_order_sum: u16 = product
            .bonds_of(idx)
            .map(|ei| match product.bond(ei).order {
                AromaticBondOrder::Known(BondOrder::Single) => 1u16,
                AromaticBondOrder::Known(BondOrder::Double) => 2,
                AromaticBondOrder::Known(BondOrder::Triple) => 3,
                AromaticBondOrder::Aromatic => 1,
            })
            .sum();

        let charge = product.atom(idx).formal_charge as i16;
        let is_aromatic = product.atom(idx).is_aromatic;

        let h = default_valences
            .iter()
            .filter_map(|&v| {
                let adjusted = v as i16 + charge;
                if adjusted > 0 {
                    Some(adjusted as u16)
                } else {
                    None
                }
            })
            .find(|&adj| adj >= bond_order_sum)
            .map(|adj| {
                let mut h = adj - bond_order_sum;
                if is_aromatic && h > 0 {
                    h -= 1;
                }
                h as u8
            })
            .unwrap_or(0);

        product.atom_mut(idx).hydrogen_count = h;
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
