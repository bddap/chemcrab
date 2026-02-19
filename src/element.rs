pub fn outer_shell_electrons(atomic_num: u8) -> u8 {
    OUTER_ELECTRONS
        .get(atomic_num as usize)
        .copied()
        .unwrap_or(0)
}

static OUTER_ELECTRONS: [u8; 119] = [
    0,  // dummy
    1, 2,                                                       // H  He
    1, 2, 3, 4, 5, 6, 7, 8,                                    // Li Be B  C  N  O  F  Ne
    1, 2, 3, 4, 5, 6, 7, 8,                                    // Na Mg Al Si P  S  Cl Ar
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8, // K  Ca Sc..Zn Ga Ge As Se Br Kr
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8, // Rb Sr Y ..Cd In Sn Sb Te I  Xe
    1, 2,                                                       // Cs Ba
    3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,            // La Ce..Yb
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,       // Lu Hf..Hg Tl Pb Bi Po At Rn
    1, 2,                                                       // Fr Ra
    3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,            // Ac Th..No
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8,       // Lr Rf..Cn Nh Fl Mc Lv Ts Og
];

/// Periodic table data for elements 1–118.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Element {
    H = 1,
    He = 2,
    Li = 3,
    Be = 4,
    B = 5,
    C = 6,
    N = 7,
    O = 8,
    F = 9,
    Ne = 10,
    Na = 11,
    Mg = 12,
    Al = 13,
    Si = 14,
    P = 15,
    S = 16,
    Cl = 17,
    Ar = 18,
    K = 19,
    Ca = 20,
    Sc = 21,
    Ti = 22,
    V = 23,
    Cr = 24,
    Mn = 25,
    Fe = 26,
    Co = 27,
    Ni = 28,
    Cu = 29,
    Zn = 30,
    Ga = 31,
    Ge = 32,
    As = 33,
    Se = 34,
    Br = 35,
    Kr = 36,
    Rb = 37,
    Sr = 38,
    Y = 39,
    Zr = 40,
    Nb = 41,
    Mo = 42,
    Tc = 43,
    Ru = 44,
    Rh = 45,
    Pd = 46,
    Ag = 47,
    Cd = 48,
    In = 49,
    Sn = 50,
    Sb = 51,
    Te = 52,
    I = 53,
    Xe = 54,
    Cs = 55,
    Ba = 56,
    La = 57,
    Ce = 58,
    Pr = 59,
    Nd = 60,
    Pm = 61,
    Sm = 62,
    Eu = 63,
    Gd = 64,
    Tb = 65,
    Dy = 66,
    Ho = 67,
    Er = 68,
    Tm = 69,
    Yb = 70,
    Lu = 71,
    Hf = 72,
    Ta = 73,
    W = 74,
    Re = 75,
    Os = 76,
    Ir = 77,
    Pt = 78,
    Au = 79,
    Hg = 80,
    Tl = 81,
    Pb = 82,
    Bi = 83,
    Po = 84,
    At = 85,
    Rn = 86,
    Fr = 87,
    Ra = 88,
    Ac = 89,
    Th = 90,
    Pa = 91,
    U = 92,
    Np = 93,
    Pu = 94,
    Am = 95,
    Cm = 96,
    Bk = 97,
    Cf = 98,
    Es = 99,
    Fm = 100,
    Md = 101,
    No = 102,
    Lr = 103,
    Rf = 104,
    Db = 105,
    Sg = 106,
    Bh = 107,
    Hs = 108,
    Mt = 109,
    Ds = 110,
    Rg = 111,
    Cn = 112,
    Nh = 113,
    Fl = 114,
    Mc = 115,
    Lv = 116,
    Ts = 117,
    Og = 118,
}

impl Element {
    pub fn from_atomic_num(n: u8) -> Option<Element> {
        if (1..=118).contains(&n) {
            // SAFETY: Element is repr(u8) with variants 1..=118, and we checked bounds.
            Some(unsafe { std::mem::transmute::<u8, Element>(n) })
        } else {
            None
        }
    }

    pub fn from_symbol(s: &str) -> Option<Element> {
        SYMBOL_TABLE.iter().find(|(sym, _)| *sym == s).map(|(_, e)| *e)
    }

    pub fn atomic_num(self) -> u8 {
        self as u8
    }

    pub fn symbol(self) -> &'static str {
        SYMBOLS[self as usize - 1]
    }

    pub fn name(self) -> &'static str {
        NAMES[self as usize - 1]
    }

    pub fn atomic_weight(self) -> f64 {
        ATOMIC_WEIGHTS[self as usize - 1]
    }

    pub fn exact_mass(self) -> f64 {
        EXACT_MASSES[self as usize - 1]
    }

    pub fn covalent_radius(self) -> Option<f64> {
        let v = COVALENT_RADII[self as usize - 1];
        if v < 0.0 { None } else { Some(v) }
    }

    pub fn vdw_radius(self) -> Option<f64> {
        let v = VDW_RADII[self as usize - 1];
        if v < 0.0 { None } else { Some(v) }
    }

    pub fn electronegativity(self) -> Option<f64> {
        let v = ELECTRONEGATIVITIES[self as usize - 1];
        if v < 0.0 { None } else { Some(v) }
    }

    pub fn default_valences(self) -> &'static [u8] {
        match self {
            Element::H => &[1],
            Element::B => &[3],
            Element::C => &[4],
            Element::N => &[3, 5],
            Element::O => &[2],
            Element::F | Element::Cl | Element::Br | Element::At => &[1],
            Element::Si | Element::Ge => &[4],
            Element::P | Element::As => &[3, 5],
            Element::S | Element::Se | Element::Te => &[2, 4, 6],
            Element::I => &[1, 3, 5, 7],
            _ => &[],
        }
    }

    pub fn is_organic_subset(self) -> bool {
        matches!(
            self,
            Element::B
                | Element::C
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::F
                | Element::Cl
                | Element::Br
                | Element::I
        )
    }
}

// symbol, Element pairs for from_symbol lookup
const SYMBOL_TABLE: [(&str, Element); 118] = [
    ("H", Element::H), ("He", Element::He), ("Li", Element::Li), ("Be", Element::Be),
    ("B", Element::B), ("C", Element::C), ("N", Element::N), ("O", Element::O),
    ("F", Element::F), ("Ne", Element::Ne), ("Na", Element::Na), ("Mg", Element::Mg),
    ("Al", Element::Al), ("Si", Element::Si), ("P", Element::P), ("S", Element::S),
    ("Cl", Element::Cl), ("Ar", Element::Ar), ("K", Element::K), ("Ca", Element::Ca),
    ("Sc", Element::Sc), ("Ti", Element::Ti), ("V", Element::V), ("Cr", Element::Cr),
    ("Mn", Element::Mn), ("Fe", Element::Fe), ("Co", Element::Co), ("Ni", Element::Ni),
    ("Cu", Element::Cu), ("Zn", Element::Zn), ("Ga", Element::Ga), ("Ge", Element::Ge),
    ("As", Element::As), ("Se", Element::Se), ("Br", Element::Br), ("Kr", Element::Kr),
    ("Rb", Element::Rb), ("Sr", Element::Sr), ("Y", Element::Y), ("Zr", Element::Zr),
    ("Nb", Element::Nb), ("Mo", Element::Mo), ("Tc", Element::Tc), ("Ru", Element::Ru),
    ("Rh", Element::Rh), ("Pd", Element::Pd), ("Ag", Element::Ag), ("Cd", Element::Cd),
    ("In", Element::In), ("Sn", Element::Sn), ("Sb", Element::Sb), ("Te", Element::Te),
    ("I", Element::I), ("Xe", Element::Xe), ("Cs", Element::Cs), ("Ba", Element::Ba),
    ("La", Element::La), ("Ce", Element::Ce), ("Pr", Element::Pr), ("Nd", Element::Nd),
    ("Pm", Element::Pm), ("Sm", Element::Sm), ("Eu", Element::Eu), ("Gd", Element::Gd),
    ("Tb", Element::Tb), ("Dy", Element::Dy), ("Ho", Element::Ho), ("Er", Element::Er),
    ("Tm", Element::Tm), ("Yb", Element::Yb), ("Lu", Element::Lu), ("Hf", Element::Hf),
    ("Ta", Element::Ta), ("W", Element::W), ("Re", Element::Re), ("Os", Element::Os),
    ("Ir", Element::Ir), ("Pt", Element::Pt), ("Au", Element::Au), ("Hg", Element::Hg),
    ("Tl", Element::Tl), ("Pb", Element::Pb), ("Bi", Element::Bi), ("Po", Element::Po),
    ("At", Element::At), ("Rn", Element::Rn), ("Fr", Element::Fr), ("Ra", Element::Ra),
    ("Ac", Element::Ac), ("Th", Element::Th), ("Pa", Element::Pa), ("U", Element::U),
    ("Np", Element::Np), ("Pu", Element::Pu), ("Am", Element::Am), ("Cm", Element::Cm),
    ("Bk", Element::Bk), ("Cf", Element::Cf), ("Es", Element::Es), ("Fm", Element::Fm),
    ("Md", Element::Md), ("No", Element::No), ("Lr", Element::Lr), ("Rf", Element::Rf),
    ("Db", Element::Db), ("Sg", Element::Sg), ("Bh", Element::Bh), ("Hs", Element::Hs),
    ("Mt", Element::Mt), ("Ds", Element::Ds), ("Rg", Element::Rg), ("Cn", Element::Cn),
    ("Nh", Element::Nh), ("Fl", Element::Fl), ("Mc", Element::Mc), ("Lv", Element::Lv),
    ("Ts", Element::Ts), ("Og", Element::Og),
];

static SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

static NAMES: [&str; 118] = [
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
    "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon",
    "Sodium", "Magnesium", "Aluminium", "Silicon", "Phosphorus",
    "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium",
    "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese",
    "Iron", "Cobalt", "Nickel", "Copper", "Zinc",
    "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine",
    "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium",
    "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium",
    "Palladium", "Silver", "Cadmium", "Indium", "Tin",
    "Antimony", "Tellurium", "Iodine", "Xenon", "Caesium",
    "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
    "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium",
    "Dysprosium", "Holmium", "Erbium", "Thulium", "Ytterbium",
    "Lutetium", "Hafnium", "Tantalum", "Tungsten", "Rhenium",
    "Osmium", "Iridium", "Platinum", "Gold", "Mercury",
    "Thallium", "Lead", "Bismuth", "Polonium", "Astatine",
    "Radon", "Francium", "Radium", "Actinium", "Thorium",
    "Protactinium", "Uranium", "Neptunium", "Plutonium", "Americium",
    "Curium", "Berkelium", "Californium", "Einsteinium", "Fermium",
    "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium", "Dubnium",
    "Seaborgium", "Bohrium", "Hassium", "Meitnerium", "Darmstadtium",
    "Roentgenium", "Copernicium", "Nihonium", "Flerovium", "Moscovium",
    "Livermorium", "Tennessine", "Oganesson",
];

// IUPAC CIAAW 2021 standard atomic weights.
// For radioactive elements without stable isotopes, mass number of longest-lived isotope.
static ATOMIC_WEIGHTS: [f64; 118] = [
    1.008,    // H
    4.002602, // He
    6.941,    // Li
    9.0121831,// Be
    10.81,    // B
    12.011,   // C
    14.007,   // N
    15.999,   // O
    18.998403163, // F
    20.1797,  // Ne
    22.98976928, // Na
    24.305,   // Mg
    26.9815384, // Al
    28.085,   // Si
    30.973761998, // P
    32.06,    // S
    35.45,    // Cl
    39.948,   // Ar
    39.0983,  // K
    40.078,   // Ca
    44.955908, // Sc
    47.867,   // Ti
    50.9415,  // V
    51.9961,  // Cr
    54.938043, // Mn
    55.845,   // Fe
    58.933194, // Co
    58.6934,  // Ni
    63.546,   // Cu
    65.38,    // Zn
    69.723,   // Ga
    72.630,   // Ge
    74.921595, // As
    78.971,   // Se
    79.904,   // Br
    83.798,   // Kr
    85.4678,  // Rb
    87.62,    // Sr
    88.90584, // Y
    91.224,   // Zr
    92.90637, // Nb
    95.95,    // Mo
    97.0,     // Tc (longest-lived isotope: 97)
    101.07,   // Ru
    102.90549, // Rh
    106.42,   // Pd
    107.8682, // Ag
    112.414,  // Cd
    114.818,  // In
    118.710,  // Sn
    121.760,  // Sb
    127.60,   // Te
    126.90447, // I
    131.293,  // Xe
    132.90545196, // Cs
    137.327,  // Ba
    138.90547, // La
    140.116,  // Ce
    140.90766, // Pr
    144.242,  // Nd
    145.0,    // Pm (longest-lived isotope: 145)
    150.36,   // Sm
    151.964,  // Eu
    157.25,   // Gd
    158.925354, // Tb
    162.500,  // Dy
    164.930328, // Ho
    167.259,  // Er
    168.934218, // Tm
    173.045,  // Yb
    174.9668, // Lu
    178.486,  // Hf
    180.94788, // Ta
    183.84,   // W
    186.207,  // Re
    190.23,   // Os
    192.217,  // Ir
    195.084,  // Pt
    196.966570, // Au
    200.592,  // Hg
    204.38,   // Tl
    207.2,    // Pb
    208.98040, // Bi
    209.0,    // Po (longest-lived: 209)
    210.0,    // At (longest-lived: 210)
    222.0,    // Rn (longest-lived: 222)
    223.0,    // Fr (longest-lived: 223)
    226.0,    // Ra (longest-lived: 226)
    227.0,    // Ac (longest-lived: 227)
    232.0377, // Th
    231.03588, // Pa
    238.02891, // U
    237.0,    // Np (longest-lived: 237)
    244.0,    // Pu (longest-lived: 244)
    243.0,    // Am (longest-lived: 243)
    247.0,    // Cm (longest-lived: 247)
    247.0,    // Bk (longest-lived: 247)
    251.0,    // Cf (longest-lived: 251)
    252.0,    // Es (longest-lived: 252)
    257.0,    // Fm (longest-lived: 257)
    258.0,    // Md (longest-lived: 258)
    259.0,    // No (longest-lived: 259)
    266.0,    // Lr (longest-lived: 266)
    267.0,    // Rf (longest-lived: 267)
    268.0,    // Db (longest-lived: 268)
    269.0,    // Sg (longest-lived: 269)
    270.0,    // Bh (longest-lived: 270)
    277.0,    // Hs (longest-lived: 277)
    278.0,    // Mt (longest-lived: 278)
    281.0,    // Ds (longest-lived: 281)
    282.0,    // Rg (longest-lived: 282)
    285.0,    // Cn (longest-lived: 285)
    286.0,    // Nh (longest-lived: 286)
    289.0,    // Fl (longest-lived: 289)
    290.0,    // Mc (longest-lived: 290)
    293.0,    // Lv (longest-lived: 293)
    294.0,    // Ts (longest-lived: 294)
    294.0,    // Og (longest-lived: 294)
];

// Monoisotopic exact masses (most abundant isotope), in daltons.
static EXACT_MASSES: [f64; 118] = [
    1.00782503207,   // H-1
    4.00260325413,   // He-4
    7.0160034366,    // Li-7
    9.012183065,     // Be-9
    11.00930536,     // B-11
    12.0,            // C-12
    14.00307400443,  // N-14
    15.99491461957,  // O-16
    18.99840316273,  // F-19
    19.9924401762,   // Ne-20
    22.9897692820,   // Na-23
    23.985041697,    // Mg-24
    26.98153853,     // Al-27
    27.97692653465,  // Si-28
    30.97376199842,  // P-31
    31.9720711744,   // S-32
    34.96885268,     // Cl-35
    39.9623831237,   // Ar-40
    38.9637064864,   // K-39
    39.962590863,    // Ca-40
    44.95590828,     // Sc-45
    47.94794198,     // Ti-48
    50.94395704,     // V-51
    51.94050623,     // Cr-52
    54.93804391,     // Mn-55
    55.93493633,     // Fe-56
    58.93319429,     // Co-59
    57.93534241,     // Ni-58
    62.92959772,     // Cu-63
    63.92914201,     // Zn-64
    68.9255735,      // Ga-69
    73.921177761,    // Ge-74
    74.92159457,     // As-75
    79.9165218,      // Se-80
    78.9183376,      // Br-79
    83.9114977282,   // Kr-84
    84.9117897379,   // Rb-85
    87.9056125,      // Sr-88
    88.9058403,      // Y-89
    89.9046977,      // Zr-90
    92.9063730,      // Nb-93
    97.90540482,     // Mo-98
    96.9063667,      // Tc-97
    101.9043441,     // Ru-102
    102.905498,      // Rh-103
    105.903483,      // Pd-106
    106.905092,      // Ag-107
    113.903365,      // Cd-114
    114.903878776,   // In-115
    119.902202,      // Sn-120
    120.903812,      // Sb-121
    129.906222748,   // Te-130
    126.904473,      // I-127
    131.904155086,   // Xe-132
    132.905451961,   // Cs-133
    137.905247,      // Ba-138
    138.906353,      // La-139
    139.905439,      // Ce-140
    140.907657,      // Pr-141
    141.907729,      // Nd-142
    144.912756,      // Pm-145
    151.919739,      // Sm-152
    152.921238,      // Eu-153
    157.924112,      // Gd-158
    158.925354,      // Tb-159
    163.929181,      // Dy-164
    164.930328,      // Ho-165
    165.930299,      // Er-166
    168.934218,      // Tm-169
    173.938867,      // Yb-174
    174.940777,      // Lu-175
    179.946557,      // Hf-180
    180.947999,      // Ta-181
    183.950933,      // W-184
    186.955752,      // Re-187
    191.961477,      // Os-192
    192.962942,      // Ir-193
    195.965836,      // Pt-195 (most abundant isotope of Pt)
    196.966570,      // Au-197
    201.970644,      // Hg-202
    204.974427,      // Tl-205
    207.976653,      // Pb-208
    208.980399,      // Bi-209
    208.982430,      // Po-209
    209.987148,      // At-210
    222.017578,      // Rn-222
    223.019736,      // Fr-223
    226.025410,      // Ra-226
    227.027752,      // Ac-227
    232.038055,      // Th-232
    231.035884,      // Pa-231
    238.050788,      // U-238
    237.048174,      // Np-237
    244.064205,      // Pu-244
    243.061381,      // Am-243
    247.070354,      // Cm-247
    247.070307,      // Bk-247
    251.079587,      // Cf-251
    252.082980,      // Es-252
    257.095106,      // Fm-257
    258.098431,      // Md-258
    259.101030,      // No-259
    266.120,         // Lr-266
    267.122,         // Rf-267
    268.126,         // Db-268
    269.129,         // Sg-269
    270.133,         // Bh-270
    277.150,         // Hs-277
    278.156,         // Mt-278
    281.165,         // Ds-281
    282.169,         // Rg-282
    285.177,         // Cn-285
    286.183,         // Nh-286
    289.190,         // Fl-289
    290.196,         // Mc-290
    293.205,         // Lv-293
    294.211,         // Ts-294
    294.214,         // Og-294
];

// Covalent radii in angstroms (Cordero 2008). -1.0 = unavailable.
static COVALENT_RADII: [f64; 118] = [
    0.31,  // H
    0.28,  // He
    1.28,  // Li
    0.96,  // Be
    0.84,  // B
    0.76,  // C(sp3)
    0.71,  // N
    0.66,  // O
    0.57,  // F
    0.58,  // Ne
    1.66,  // Na
    1.41,  // Mg
    1.21,  // Al
    1.11,  // Si
    1.07,  // P
    1.05,  // S
    1.02,  // Cl
    1.06,  // Ar
    2.03,  // K
    1.76,  // Ca
    1.70,  // Sc
    1.60,  // Ti
    1.53,  // V
    1.39,  // Cr
    1.39,  // Mn (low spin)
    1.32,  // Fe (low spin)
    1.26,  // Co (low spin)
    1.24,  // Ni
    1.32,  // Cu
    1.22,  // Zn
    1.22,  // Ga
    1.20,  // Ge
    1.19,  // As
    1.20,  // Se
    1.20,  // Br
    1.16,  // Kr
    2.20,  // Rb
    1.95,  // Sr
    1.90,  // Y
    1.75,  // Zr
    1.64,  // Nb
    1.54,  // Mo
    1.47,  // Tc
    1.46,  // Ru
    1.42,  // Rh
    1.39,  // Pd
    1.45,  // Ag
    1.44,  // Cd
    1.42,  // In
    1.39,  // Sn
    1.39,  // Sb
    1.38,  // Te
    1.39,  // I
    1.40,  // Xe
    2.44,  // Cs
    2.15,  // Ba
    2.07,  // La
    2.04,  // Ce
    2.03,  // Pr
    2.01,  // Nd
    1.99,  // Pm
    1.98,  // Sm
    1.98,  // Eu
    1.96,  // Gd
    1.94,  // Tb
    1.92,  // Dy
    1.92,  // Ho
    1.89,  // Er
    1.90,  // Tm
    1.87,  // Yb
    1.87,  // Lu
    1.75,  // Hf
    1.70,  // Ta
    1.62,  // W
    1.51,  // Re
    1.44,  // Os
    1.41,  // Ir
    1.36,  // Pt
    1.36,  // Au
    1.32,  // Hg
    1.45,  // Tl
    1.46,  // Pb
    1.48,  // Bi
    1.40,  // Po
    1.50,  // At
    1.50,  // Rn
    2.60,  // Fr
    2.21,  // Ra
    2.15,  // Ac
    2.06,  // Th
    2.00,  // Pa
    1.96,  // U
    1.90,  // Np
    1.87,  // Pu
    1.80,  // Am
    1.69,  // Cm
    -1.0,  // Bk
    -1.0,  // Cf
    -1.0,  // Es
    -1.0,  // Fm
    -1.0,  // Md
    -1.0,  // No
    -1.0,  // Lr
    -1.0,  // Rf
    -1.0,  // Db
    -1.0,  // Sg
    -1.0,  // Bh
    -1.0,  // Hs
    -1.0,  // Mt
    -1.0,  // Ds
    -1.0,  // Rg
    -1.0,  // Cn
    -1.0,  // Nh
    -1.0,  // Fl
    -1.0,  // Mc
    -1.0,  // Lv
    -1.0,  // Ts
    -1.0,  // Og
];

// Van der Waals radii in angstroms (Bondi 1964 + Rowland/Taylor 1996). -1.0 = unavailable.
static VDW_RADII: [f64; 118] = [
    1.20,  // H
    1.40,  // He
    1.82,  // Li
    1.53,  // Be
    1.92,  // B
    1.70,  // C
    1.55,  // N
    1.52,  // O
    1.47,  // F
    1.54,  // Ne
    2.27,  // Na
    1.73,  // Mg
    1.84,  // Al
    2.10,  // Si
    1.80,  // P
    1.80,  // S
    1.75,  // Cl
    1.88,  // Ar
    2.75,  // K
    2.31,  // Ca
    -1.0,  // Sc
    -1.0,  // Ti
    -1.0,  // V
    -1.0,  // Cr
    -1.0,  // Mn
    -1.0,  // Fe
    -1.0,  // Co
    1.63,  // Ni
    1.40,  // Cu
    1.39,  // Zn
    1.87,  // Ga
    2.11,  // Ge
    1.85,  // As
    1.90,  // Se
    1.85,  // Br
    2.02,  // Kr
    3.03,  // Rb
    2.49,  // Sr
    -1.0,  // Y
    -1.0,  // Zr
    -1.0,  // Nb
    -1.0,  // Mo
    -1.0,  // Tc
    -1.0,  // Ru
    -1.0,  // Rh
    1.63,  // Pd
    1.72,  // Ag
    1.58,  // Cd
    1.93,  // In
    2.17,  // Sn
    2.06,  // Sb
    2.06,  // Te
    1.98,  // I
    2.16,  // Xe
    3.43,  // Cs
    2.68,  // Ba
    -1.0,  // La
    -1.0,  // Ce
    -1.0,  // Pr
    -1.0,  // Nd
    -1.0,  // Pm
    -1.0,  // Sm
    -1.0,  // Eu
    -1.0,  // Gd
    -1.0,  // Tb
    -1.0,  // Dy
    -1.0,  // Ho
    -1.0,  // Er
    -1.0,  // Tm
    -1.0,  // Yb
    -1.0,  // Lu
    -1.0,  // Hf
    -1.0,  // Ta
    -1.0,  // W
    -1.0,  // Re
    -1.0,  // Os
    -1.0,  // Ir
    1.75,  // Pt
    1.66,  // Au
    1.55,  // Hg
    1.96,  // Tl
    2.02,  // Pb
    2.07,  // Bi
    1.97,  // Po
    2.02,  // At
    2.20,  // Rn
    3.48,  // Fr
    2.83,  // Ra
    -1.0,  // Ac
    -1.0,  // Th
    -1.0,  // Pa
    1.86,  // U
    -1.0,  // Np
    -1.0,  // Pu
    -1.0,  // Am
    -1.0,  // Cm
    -1.0,  // Bk
    -1.0,  // Cf
    -1.0,  // Es
    -1.0,  // Fm
    -1.0,  // Md
    -1.0,  // No
    -1.0,  // Lr
    -1.0,  // Rf
    -1.0,  // Db
    -1.0,  // Sg
    -1.0,  // Bh
    -1.0,  // Hs
    -1.0,  // Mt
    -1.0,  // Ds
    -1.0,  // Rg
    -1.0,  // Cn
    -1.0,  // Nh
    -1.0,  // Fl
    -1.0,  // Mc
    -1.0,  // Lv
    -1.0,  // Ts
    -1.0,  // Og
];

// Pauling electronegativity. -1.0 = no reliable value.
static ELECTRONEGATIVITIES: [f64; 118] = [
    2.20,  // H
    -1.0,  // He
    0.98,  // Li
    1.57,  // Be
    2.04,  // B
    2.55,  // C
    3.04,  // N
    3.44,  // O
    3.98,  // F
    -1.0,  // Ne
    0.93,  // Na
    1.31,  // Mg
    1.61,  // Al
    1.90,  // Si
    2.19,  // P
    2.58,  // S
    3.16,  // Cl
    -1.0,  // Ar
    0.82,  // K
    1.00,  // Ca
    1.36,  // Sc
    1.54,  // Ti
    1.63,  // V
    1.66,  // Cr
    1.55,  // Mn
    1.83,  // Fe
    1.88,  // Co
    1.91,  // Ni
    1.90,  // Cu
    1.65,  // Zn
    1.81,  // Ga
    2.01,  // Ge
    2.18,  // As
    2.55,  // Se
    2.96,  // Br
    3.00,  // Kr
    0.82,  // Rb
    0.95,  // Sr
    1.22,  // Y
    1.33,  // Zr
    1.6,   // Nb
    2.16,  // Mo
    1.9,   // Tc
    2.2,   // Ru
    2.28,  // Rh
    2.20,  // Pd
    1.93,  // Ag
    1.69,  // Cd
    1.78,  // In
    1.96,  // Sn
    2.05,  // Sb
    2.1,   // Te
    2.66,  // I
    2.60,  // Xe
    0.79,  // Cs
    0.89,  // Ba
    1.10,  // La
    1.12,  // Ce
    1.13,  // Pr
    1.14,  // Nd
    -1.0,  // Pm
    1.17,  // Sm
    -1.0,  // Eu
    1.20,  // Gd
    -1.0,  // Tb
    1.22,  // Dy
    1.23,  // Ho
    1.24,  // Er
    1.25,  // Tm
    -1.0,  // Yb
    1.27,  // Lu
    1.3,   // Hf
    1.5,   // Ta
    2.36,  // W
    1.9,   // Re
    2.2,   // Os
    2.20,  // Ir
    2.28,  // Pt
    2.54,  // Au
    2.00,  // Hg
    1.62,  // Tl
    2.33,  // Pb
    2.02,  // Bi
    2.0,   // Po
    2.2,   // At
    -1.0,  // Rn
    0.7,   // Fr
    0.9,   // Ra
    1.1,   // Ac
    1.3,   // Th
    1.5,   // Pa
    1.38,  // U
    1.36,  // Np
    1.28,  // Pu
    1.3,   // Am
    1.3,   // Cm
    1.3,   // Bk
    1.3,   // Cf
    1.3,   // Es
    1.3,   // Fm
    1.3,   // Md
    1.3,   // No
    1.3,   // Lr
    -1.0,  // Rf
    -1.0,  // Db
    -1.0,  // Sg
    -1.0,  // Bh
    -1.0,  // Hs
    -1.0,  // Mt
    -1.0,  // Ds
    -1.0,  // Rg
    -1.0,  // Cn
    -1.0,  // Nh
    -1.0,  // Fl
    -1.0,  // Mc
    -1.0,  // Lv
    -1.0,  // Ts
    -1.0,  // Og
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_atomic_num_round_trip() {
        for n in 1u8..=118 {
            let e = Element::from_atomic_num(n).unwrap();
            assert_eq!(e.atomic_num(), n);
        }
    }

    #[test]
    fn from_atomic_num_boundaries() {
        assert!(Element::from_atomic_num(0).is_none());
        assert!(Element::from_atomic_num(119).is_none());
        assert!(Element::from_atomic_num(255).is_none());
        assert_eq!(Element::from_atomic_num(1), Some(Element::H));
        assert_eq!(Element::from_atomic_num(118), Some(Element::Og));
    }

    #[test]
    fn from_symbol_exact_match() {
        assert_eq!(Element::from_symbol("He"), Some(Element::He));
        assert_eq!(Element::from_symbol("Fe"), Some(Element::Fe));
        assert_eq!(Element::from_symbol("Og"), Some(Element::Og));
    }

    #[test]
    fn from_symbol_case_sensitive() {
        assert!(Element::from_symbol("he").is_none());
        assert!(Element::from_symbol("HE").is_none());
        assert!(Element::from_symbol("hE").is_none());
        assert!(Element::from_symbol("").is_none());
        assert!(Element::from_symbol("Xx").is_none());
    }

    #[test]
    fn symbol_round_trip() {
        for n in 1u8..=118 {
            let e = Element::from_atomic_num(n).unwrap();
            assert_eq!(Element::from_symbol(e.symbol()), Some(e));
        }
    }

    #[test]
    fn names_spot_check() {
        assert_eq!(Element::H.name(), "Hydrogen");
        assert_eq!(Element::C.name(), "Carbon");
        assert_eq!(Element::Fe.name(), "Iron");
        assert_eq!(Element::Au.name(), "Gold");
        assert_eq!(Element::Og.name(), "Oganesson");
    }

    #[test]
    fn atomic_weight_spot_check() {
        assert!((Element::H.atomic_weight() - 1.008).abs() < 0.001);
        assert!((Element::C.atomic_weight() - 12.011).abs() < 0.001);
        assert!((Element::O.atomic_weight() - 15.999).abs() < 0.001);
        assert!((Element::Fe.atomic_weight() - 55.845).abs() < 0.001);
        assert!((Element::U.atomic_weight() - 238.02891).abs() < 0.001);
    }

    #[test]
    fn exact_mass_spot_check() {
        assert!((Element::C.exact_mass() - 12.0).abs() < 1e-6);
        assert!((Element::H.exact_mass() - 1.00782503207).abs() < 1e-6);
        assert!((Element::O.exact_mass() - 15.99491461957).abs() < 1e-6);
    }

    #[test]
    fn covalent_radius_available() {
        assert!((Element::C.covalent_radius().unwrap() - 0.76).abs() < 0.01);
        assert!((Element::H.covalent_radius().unwrap() - 0.31).abs() < 0.01);
        assert!(Element::Og.covalent_radius().is_none());
    }

    #[test]
    fn vdw_radius_available() {
        assert!((Element::C.vdw_radius().unwrap() - 1.70).abs() < 0.01);
        assert!((Element::H.vdw_radius().unwrap() - 1.20).abs() < 0.01);
        assert!(Element::Sc.vdw_radius().is_none());
    }

    #[test]
    fn electronegativity_available() {
        assert!((Element::F.electronegativity().unwrap() - 3.98).abs() < 0.01);
        assert!((Element::Cs.electronegativity().unwrap() - 0.79).abs() < 0.01);
        assert!(Element::He.electronegativity().is_none());
        assert!(Element::Ne.electronegativity().is_none());
    }

    #[test]
    fn default_valences_smiles() {
        assert_eq!(Element::B.default_valences(), &[3]);
        assert_eq!(Element::C.default_valences(), &[4]);
        assert_eq!(Element::N.default_valences(), &[3, 5]);
        assert_eq!(Element::O.default_valences(), &[2]);
        assert_eq!(Element::P.default_valences(), &[3, 5]);
        assert_eq!(Element::S.default_valences(), &[2, 4, 6]);
        assert_eq!(Element::F.default_valences(), &[1]);
        assert_eq!(Element::Cl.default_valences(), &[1]);
        assert_eq!(Element::Br.default_valences(), &[1]);
        assert_eq!(Element::I.default_valences(), &[1, 3, 5, 7]);
    }

    #[test]
    fn default_valences_hydrogen() {
        assert_eq!(Element::H.default_valences(), &[1]);
    }

    #[test]
    fn default_valences_extended() {
        assert_eq!(Element::Si.default_valences(), &[4]);
        assert_eq!(Element::Ge.default_valences(), &[4]);
        assert_eq!(Element::As.default_valences(), &[3, 5]);
        assert_eq!(Element::Se.default_valences(), &[2, 4, 6]);
        assert_eq!(Element::Te.default_valences(), &[2, 4, 6]);
        assert_eq!(Element::I.default_valences(), &[1, 3, 5, 7]);
        assert_eq!(Element::At.default_valences(), &[1]);
    }

    #[test]
    fn default_valences_non_organic_empty() {
        assert_eq!(Element::He.default_valences(), &[] as &[u8]);
        assert_eq!(Element::Fe.default_valences(), &[] as &[u8]);
        assert_eq!(Element::Og.default_valences(), &[] as &[u8]);
    }

    #[test]
    fn organic_subset() {
        assert!(Element::C.is_organic_subset());
        assert!(Element::Br.is_organic_subset());
        assert!(!Element::Fe.is_organic_subset());
        assert!(!Element::H.is_organic_subset());
    }

    #[test]
    fn enum_is_copy_and_hashable() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        let c = Element::C;
        let c2 = c;
        set.insert(c);
        set.insert(c2);
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn debug_format() {
        assert_eq!(format!("{:?}", Element::H), "H");
        assert_eq!(format!("{:?}", Element::Og), "Og");
    }

    #[test]
    fn all_symbols_unique() {
        use std::collections::HashSet;
        let symbols: HashSet<&str> = (1u8..=118).map(|n| Element::from_atomic_num(n).unwrap().symbol()).collect();
        assert_eq!(symbols.len(), 118);
    }

    #[test]
    fn all_names_unique() {
        use std::collections::HashSet;
        let names: HashSet<&str> = (1u8..=118).map(|n| Element::from_atomic_num(n).unwrap().name()).collect();
        assert_eq!(names.len(), 118);
    }

    #[test]
    fn weights_positive() {
        for n in 1u8..=118 {
            assert!(Element::from_atomic_num(n).unwrap().atomic_weight() > 0.0);
        }
    }

    #[test]
    fn exact_masses_positive() {
        for n in 1u8..=118 {
            assert!(Element::from_atomic_num(n).unwrap().exact_mass() > 0.0);
        }
    }
}
