import numpy as np
import os

restypes = [
    "A",  # 0
    "R",  # 1
    "N",  # 2
    "D",  # 3
    "C",  # 4
    "Q",  # 5
    "E",  # 6
    "G",  # 7
    "H",  # 8
    "I",  # 9
    "L",  # 10
    "K",  # 11
    "M",  # 12
    "F",  # 13
    "P",  # 14
    "S",  # 15
    "T",  # 16
    "W",  # 17
    "Y",  # 18
    "V",  # 19
]
restypes_with_x = restypes + ["X"]
restype_order = {restype: i for i, restype in enumerate(restypes)}
restype_num = len(restypes)  # := 20.

restype_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}
restype_3to1 = {v: k for k, v in restype_1to3.items()}

# Atoms positions relative to the 8 rigid groups, defined by the pre-omega, phi,
# psi and chi angles:
# 0: 'backbone group',
# 1: 'pre-omega-group', (empty)
# 2: 'phi-group', (currently empty, because it defines only hydrogens)
# 3: 'psi-group',
# 4,5,6,7: 'chi1,2,3,4-group'
# The atom positions are relative to the axis-end-atom of the corresponding
# rotation axis. The x-axis is in direction of the rotation axis, and the y-axis
# is defined such that the dihedral-angle-definiting atom (the last entry in
# chi_angles_atoms above) is in the xy-plane (with a positive y-coordinate).
# format: [atomname, group_idx, rel_position]
rigid_group_atom_positions = {
    "ALA": [
        ["N", 0, (-0.525, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.529, -0.774, -1.205)],
        ["O", 3, (0.627, 1.062, 0.000)],
    ],
    "ARG": [
        ["N", 0, (-0.524, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.524, -0.778, -1.209)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG", 4, (0.616, 1.390, -0.000)],
        ["CD", 5, (0.564, 1.414, 0.000)],
        ["NE", 6, (0.539, 1.357, -0.000)],
        ["NH1", 7, (0.206, 2.301, 0.000)],
        ["NH2", 7, (2.078, 0.978, -0.000)],
        ["CZ", 7, (0.758, 1.093, -0.000)],
    ],
    "ASN": [
        ["N", 0, (-0.536, 1.357, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.531, -0.787, -1.200)],
        ["O", 3, (0.625, 1.062, 0.000)],
        ["CG", 4, (0.584, 1.399, 0.000)],
        ["ND2", 5, (0.593, -1.188, 0.001)],
        ["OD1", 5, (0.633, 1.059, 0.000)],
    ],
    "ASP": [
        ["N", 0, (-0.525, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, 0.000, -0.000)],
        ["CB", 0, (-0.526, -0.778, -1.208)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.593, 1.398, -0.000)],
        ["OD1", 5, (0.610, 1.091, 0.000)],
        ["OD2", 5, (0.592, -1.101, -0.003)],
    ],
    "CYS": [
        ["N", 0, (-0.522, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, 0.000, 0.000)],
        ["CB", 0, (-0.519, -0.773, -1.212)],
        ["O", 3, (0.625, 1.062, -0.000)],
        ["SG", 4, (0.728, 1.653, 0.000)],
    ],
    "GLN": [
        ["N", 0, (-0.526, 1.361, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, 0.000)],
        ["CB", 0, (-0.525, -0.779, -1.207)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.615, 1.393, 0.000)],
        ["CD", 5, (0.587, 1.399, -0.000)],
        ["NE2", 6, (0.593, -1.189, -0.001)],
        ["OE1", 6, (0.634, 1.060, 0.000)],
    ],
    "GLU": [
        ["N", 0, (-0.528, 1.361, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, -0.000, -0.000)],
        ["CB", 0, (-0.526, -0.781, -1.207)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG", 4, (0.615, 1.392, 0.000)],
        ["CD", 5, (0.600, 1.397, 0.000)],
        ["OE1", 6, (0.607, 1.095, -0.000)],
        ["OE2", 6, (0.589, -1.104, -0.001)],
    ],
    "GLY": [
        ["N", 0, (-0.572, 1.337, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.517, -0.000, -0.000)],
        ["O", 3, (0.626, 1.062, -0.000)],
    ],
    "HIS": [
        ["N", 0, (-0.527, 1.360, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, 0.000, 0.000)],
        ["CB", 0, (-0.525, -0.778, -1.208)],
        ["O", 3, (0.625, 1.063, 0.000)],
        ["CG", 4, (0.600, 1.370, -0.000)],
        ["CD2", 5, (0.889, -1.021, 0.003)],
        ["ND1", 5, (0.744, 1.160, -0.000)],
        ["CE1", 5, (2.030, 0.851, 0.002)],
        ["NE2", 5, (2.145, -0.466, 0.004)],
    ],
    "ILE": [
        ["N", 0, (-0.493, 1.373, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, -0.000)],
        ["CB", 0, (-0.536, -0.793, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG1", 4, (0.534, 1.437, -0.000)],
        ["CG2", 4, (0.540, -0.785, -1.199)],
        ["CD1", 5, (0.619, 1.391, 0.000)],
    ],
    "LEU": [
        ["N", 0, (-0.520, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.522, -0.773, -1.214)],
        ["O", 3, (0.625, 1.063, -0.000)],
        ["CG", 4, (0.678, 1.371, 0.000)],
        ["CD1", 5, (0.530, 1.430, -0.000)],
        ["CD2", 5, (0.535, -0.774, 1.200)],
    ],
    "LYS": [
        ["N", 0, (-0.526, 1.362, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, 0.000)],
        ["CB", 0, (-0.524, -0.778, -1.208)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.619, 1.390, 0.000)],
        ["CD", 5, (0.559, 1.417, 0.000)],
        ["CE", 6, (0.560, 1.416, 0.000)],
        ["NZ", 7, (0.554, 1.387, 0.000)],
    ],
    "MET": [
        ["N", 0, (-0.521, 1.364, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, 0.000, 0.000)],
        ["CB", 0, (-0.523, -0.776, -1.210)],
        ["O", 3, (0.625, 1.062, -0.000)],
        ["CG", 4, (0.613, 1.391, -0.000)],
        ["SD", 5, (0.703, 1.695, 0.000)],
        ["CE", 6, (0.320, 1.786, -0.000)],
    ],
    "PHE": [
        ["N", 0, (-0.518, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, 0.000, -0.000)],
        ["CB", 0, (-0.525, -0.776, -1.212)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["CG", 4, (0.607, 1.377, 0.000)],
        ["CD1", 5, (0.709, 1.195, -0.000)],
        ["CD2", 5, (0.706, -1.196, 0.000)],
        ["CE1", 5, (2.102, 1.198, -0.000)],
        ["CE2", 5, (2.098, -1.201, -0.000)],
        ["CZ", 5, (2.794, -0.003, -0.001)],
    ],
    "PRO": [
        ["N", 0, (-0.566, 1.351, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, 0.000)],
        ["CB", 0, (-0.546, -0.611, -1.293)],
        ["O", 3, (0.621, 1.066, 0.000)],
        ["CG", 4, (0.382, 1.445, 0.0)],
        # ['CD', 5, (0.427, 1.440, 0.0)],
        ["CD", 5, (0.477, 1.424, 0.0)],  # manually made angle 2 degrees larger
    ],
    "SER": [
        ["N", 0, (-0.529, 1.360, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, -0.000)],
        ["CB", 0, (-0.518, -0.777, -1.211)],
        ["O", 3, (0.626, 1.062, -0.000)],
        ["OG", 4, (0.503, 1.325, 0.000)],
    ],
    "THR": [
        ["N", 0, (-0.517, 1.364, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.526, 0.000, -0.000)],
        ["CB", 0, (-0.516, -0.793, -1.215)],
        ["O", 3, (0.626, 1.062, 0.000)],
        ["CG2", 4, (0.550, -0.718, -1.228)],
        ["OG1", 4, (0.472, 1.353, 0.000)],
    ],
    "TRP": [
        ["N", 0, (-0.521, 1.363, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.525, -0.000, 0.000)],
        ["CB", 0, (-0.523, -0.776, -1.212)],
        ["O", 3, (0.627, 1.062, 0.000)],
        ["CG", 4, (0.609, 1.370, -0.000)],
        ["CD1", 5, (0.824, 1.091, 0.000)],
        ["CD2", 5, (0.854, -1.148, -0.005)],
        ["CE2", 5, (2.186, -0.678, -0.007)],
        ["CE3", 5, (0.622, -2.530, -0.007)],
        ["NE1", 5, (2.140, 0.690, -0.004)],
        ["CH2", 5, (3.028, -2.890, -0.013)],
        ["CZ2", 5, (3.283, -1.543, -0.011)],
        ["CZ3", 5, (1.715, -3.389, -0.011)],
    ],
    "TYR": [
        ["N", 0, (-0.522, 1.362, 0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.524, -0.000, -0.000)],
        ["CB", 0, (-0.522, -0.776, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG", 4, (0.607, 1.382, -0.000)],
        ["CD1", 5, (0.716, 1.195, -0.000)],
        ["CD2", 5, (0.713, -1.194, -0.001)],
        ["CE1", 5, (2.107, 1.200, -0.002)],
        ["CE2", 5, (2.104, -1.201, -0.003)],
        ["OH", 5, (4.168, -0.002, -0.005)],
        ["CZ", 5, (2.791, -0.001, -0.003)],
    ],
    "VAL": [
        ["N", 0, (-0.494, 1.373, -0.000)],
        ["CA", 0, (0.000, 0.000, 0.000)],
        ["C", 0, (1.527, -0.000, -0.000)],
        ["CB", 0, (-0.533, -0.795, -1.213)],
        ["O", 3, (0.627, 1.062, -0.000)],
        ["CG1", 4, (0.540, 1.429, -0.000)],
        ["CG2", 4, (0.533, -0.776, 1.203)],
    ],
    "UNK": [],
}

# A list of atoms (excluding hydrogen) for each AA type. PDB naming convention.
residue_atoms = {
    "ALA": ["C", "CA", "CB", "N", "O"],
    "ARG": ["C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "O", "NH1", "NH2"],
    "ASP": ["C", "CA", "CB", "CG", "N", "O", "OD1", "OD2"],
    "ASN": ["C", "CA", "CB", "CG", "N", "ND2", "O", "OD1"],
    "CYS": ["C", "CA", "CB", "N", "O", "SG"],
    "GLU": ["C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2"],
    "GLN": ["C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1"],
    "GLY": ["C", "CA", "N", "O"],
    "HIS": ["C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O"],
    "ILE": ["C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O"],
    "LEU": ["C", "CA", "CB", "CG", "CD1", "CD2", "N", "O"],
    "LYS": ["C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O"],
    "MET": ["C", "CA", "CB", "CG", "CE", "N", "O", "SD"],
    "PHE": ["C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O"],
    "PRO": ["C", "CA", "CB", "CG", "CD", "N", "O"],
    "SER": ["C", "CA", "CB", "N", "O", "OG"],
    "THR": ["C", "CA", "CB", "CG2", "N", "O", "OG1"],
    "TRP": [
        "C",
        "CA",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
        "N",
        "NE1",
        "O",
    ],
    "TYR": [
        "C",
        "CA",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "N",
        "O",
        "OH",
    ],
    "VAL": ["C", "CA", "CB", "CG1", "CG2", "N", "O"],
}


residue_atom_renaming_swaps = {
    "PHE": [["CD1", "CD2"], ["CE1", "CE2"]],
    "TYR": [["CD1", "CD2"], ["CE1", "CE2"]],
    "ARG": [["NH1", "NH2"]],
    "ASP": [["OD1", "OD2"]],
    "GLU": [["OE1", "OE2"]],
}


# Van der Waals radii [Angstroem] of the atoms (from Wikipedia)
van_der_waals_radius = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
}

# Sidechain bond lengths (from Rosetta database)
sc_bond_lengths = {
    "ARG": {
        ("CB", "CG"): 1.5204,
        ("CG", "CD"): 1.4854,
        ("CD", "NE"): 1.4541,
        ("NE", "CZ"): 1.3473,
    },
    "ASN": {("CB", "CG"): 1.5035, ("CG", "OD1"): 1.2364},
    "ASP": {("CB", "CG"): 1.5228, ("CG", "OD1"): 1.2082},
    "CYS": {("CB", "SG"): 1.8088},
    "GLU": {("CB", "CG"): 1.5221, ("CG", "CD"): 1.5034, ("CD", "OE1"): 1.2076},
    "GLN": {("CB", "CG"): 1.5191, ("CG", "CD"): 1.5169, ("CD", "OE1"): 1.2342},
    "HIS": {("CB", "CG"): 1.4972, ("CG", "ND1"): 1.3792},
    "ILE": {("CB", "CG1"): 1.5309, ("CG1", "CD1"): 1.5117},
    "LEU": {("CB", "CG"): 1.5340, ("CG", "CD1"): 1.5227},
    "LYS": {
        ("CB", "CG"): 1.5229,
        ("CG", "CD"): 1.5213,
        ("CD", "CE"): 1.5216,
        ("CE", "NZ"): 1.4881,
    },
    "MET": {("CB", "CG"): 1.5222, ("CG", "SD"): 1.8038, ("SD", "CE"): 1.7904},
    "PHE": {("CB", "CG"): 1.5022, ("CG", "CD1"): 1.3870},
    "PRO": {("CB", "CG"): 1.4906, ("CG", "CD"): 1.5055},
    "SER": {("CB", "OG"): 1.4012},
    "THR": {("CB", "OG1"): 1.4335},
    "TRP": {("CB", "CG"): 1.4987, ("CG", "CD1"): 1.3627},
    "TYR": {("CB", "CG"): 1.5127, ("CG", "CD1"): 1.3872},
    "VAL": {("CB", "CG1"): 1.5214},
}


# This mapping is used when we need to store atom data in a format that requires
# fixed atom data size for every residue (e.g. a numpy array).
atom_types = [
    "N",
    "CA",
    "C",
    "CB",
    "O",
    "CG",
    "CG1",
    "CG2",
    "OG",
    "OG1",
    "SG",
    "CD",
    "CD1",
    "CD2",
    "ND1",
    "ND2",
    "OD1",
    "OD2",
    "SD",
    "CE",
    "CE1",
    "CE2",
    "CE3",
    "NE",
    "NE1",
    "NE2",
    "OE1",
    "OE2",
    "CH2",
    "NH1",
    "NH2",
    "OH",
    "CZ",
    "CZ2",
    "CZ3",
    "NZ",
    "OXT",
]
atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
atom_type_num = len(atom_types)  # := 37.

# A compact atom encoding with 14 columns
# pylint: disable=line-too-long
# pylint: disable=bad-whitespace
restype_name_to_atom14_names = {
    "ALA": ["N", "CA", "C", "O", "CB", "", "", "", "", "", "", "", "", ""],
    "ARG": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "NE",
        "CZ",
        "NH1",
        "NH2",
        "",
        "",
        "",
    ],
    "ASN": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "OD1",
        "ND2",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "ASP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "OD1",
        "OD2",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "CYS": ["N", "CA", "C", "O", "CB", "SG", "", "", "", "", "", "", "", ""],
    "GLN": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "OE1",
        "NE2",
        "",
        "",
        "",
        "",
        "",
    ],
    "GLU": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "OE1",
        "OE2",
        "",
        "",
        "",
        "",
        "",
    ],
    "GLY": ["N", "CA", "C", "O", "", "", "", "", "", "", "", "", "", ""],
    "HIS": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "ND1",
        "CD2",
        "CE1",
        "NE2",
        "",
        "",
        "",
        "",
    ],
    "ILE": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "LEU": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "LYS": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "CE",
        "NZ",
        "",
        "",
        "",
        "",
        "",
    ],
    "MET": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "SD",
        "CE",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "PHE": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "",
        "",
        "",
    ],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD", "", "", "", "", "", "", ""],
    "SER": ["N", "CA", "C", "O", "CB", "OG", "", "", "", "", "", "", "", ""],
    "THR": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "OG1",
        "CG2",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "TRP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],
    "TYR": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "OH",
        "",
        "",
    ],
    "VAL": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG1",
        "CG2",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "UNK": ["", "", "", "", "", "", "", "", "", "", "", "", "", ""],
}

# An atom encoding with 27 columns (includes heavy atoms and hydrogens)
# Note for HIS, the V1 atom is a virtual atom representing the center of the imidazole ring
restype_name_to_atom27_names = {
    "ALA": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB1",
        "HB2",
        "HB3",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "ARG": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "NE",
        "CZ",
        "NH1",
        "NH2",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HD2",
        "HD3",
        "HE",
        "HH11",
        "HH12",
        "HH21",
        "HH22",
    ],
    "ASN": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "OD1",
        "ND2",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HD21",
        "HD22",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "ASP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "OD1",
        "OD2",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "CYS": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "SG",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "GLN": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "OE1",
        "NE2",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HE21",
        "HE22",
        "",
        "",
        "",
        "",
        "",
    ],
    "GLU": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "OE1",
        "OE2",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "GLY": [
        "N",
        "CA",
        "C",
        "O",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA2",
        "HA3",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "HIS": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "ND1",
        "CD2",
        "CE1",
        "NE2",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HD1",
        "HD2",
        "HE1",
        "HE2",
        "V1",
        "",
        "",
        "",
        "",
    ],
    "ILE": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG1",
        "CG2",
        "CD1",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB",
        "HG12",
        "HG13",
        "HG21",
        "HG22",
        "HG23",
        "HD11",
        "HD12",
        "HD13",
        "",
        "",
    ],
    "LEU": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG",
        "HD11",
        "HD12",
        "HD13",
        "HD21",
        "HD22",
        "HD23",
        "",
        "",
    ],
    "LYS": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "CE",
        "NZ",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HD2",
        "HD3",
        "HE2",
        "HE3",
        "HZ1",
        "HZ2",
        "HZ3",
    ],
    "MET": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "SD",
        "CE",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HE1",
        "HE2",
        "HE3",
        "",
        "",
        "",
        "",
    ],
    "PHE": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HD1",
        "HD2",
        "HE1",
        "HE2",
        "HZ",
        "",
        "",
        "",
        "",
    ],
    "PRO": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "HA",
        "HB2",
        "HB3",
        "HG2",
        "HG3",
        "HD2",
        "HD3",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "SER": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "OG",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HG",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "THR": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "OG1",
        "CG2",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB",
        "HG1",
        "HG21",
        "HG22",
        "HG23",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
    "TRP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HD1",
        "HE1",
        "HE3",
        "HZ2",
        "HZ3",
        "HH2",
        "",
        "",
        "",
    ],
    "TYR": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "OH",
        "",
        "",
        "H",
        "HA",
        "HB2",
        "HB3",
        "HD1",
        "HD2",
        "HE1",
        "HE2",
        "HH",
        "",
        "",
        "",
        "",
    ],
    "VAL": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG1",
        "CG2",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "H",
        "HA",
        "HB",
        "HG11",
        "HG12",
        "HG13",
        "HG21",
        "HG22",
        "HG23",
        "",
        "",
        "",
        "",
    ],
    "UNK": [
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ],
}


def get_atom27_index(res: str, atom: str) -> int:
    if atom in restype_name_to_atom27_names[res]:
        return restype_name_to_atom27_names[res].index(atom)
    else:
        return None


def _make_rigid_transformation_4x4(ex, ey, translation):
    """Create a rigid 4x4 transformation matrix from two axes and transl."""
    # Normalize ex.
    ex_normalized = ex / np.linalg.norm(ex)

    # make ey perpendicular to ex
    ey_normalized = ey - np.dot(ey, ex_normalized) * ex_normalized
    ey_normalized /= np.linalg.norm(ey_normalized)

    # compute ez as cross product
    eznorm = np.cross(ex_normalized, ey_normalized)
    m = np.stack([ex_normalized, ey_normalized, eznorm, translation]).transpose()
    m = np.concatenate([m, [[0.0, 0.0, 0.0, 1.0]]], axis=0)
    return m


hbond_donor_atoms = [
    "OG",
    "OG1",
    "NE2",
    "ND1",
    "ND2",
    "NZ",
    "NE",
    "NH1",
    "NH2",
    "NE1",
    "OH",
    "N",
]
hbond_acceptor_atoms = [
    "ND1",
    "NE2",
    "OE1",
    "OE2",
    "OD1",
    "OD2",
    "OH",
    "OG",
    "OG1",
    "O",
]

# Can't just use H atom names b/c some overlap between polar/nonpolar H atoms
# e.g., HG can mean CYS thiol (polar) or LEU CG (nonpolar)
restype_name_to_atom27_polar_h_atoms = {
    "ALA": ["H", "H1", "H2", "H3"],
    "ARG": ["H", "H1", "H2", "H3", "HE", "HH11", "HH12", "HH21", "HH22"],
    "ASN": ["H", "H1", "H2", "H3", "HD21", "HD22"],
    "ASP": ["H", "H1", "H2", "H3"],
    "CYS": ["H", "H1", "H2", "H3", "HG"],
    "GLN": ["H", "H1", "H2", "H3", "HE21", "HE22"],
    "GLU": ["H", "H1", "H2", "H3"],
    "GLY": ["H", "H1", "H2", "H3"],
    "HIS": ["H", "H1", "H2", "H3", "HD1", "HE2"],
    "ILE": ["H", "H1", "H2", "H3"],
    "LEU": ["H", "H1", "H2", "H3"],
    "LYS": ["H", "H1", "H2", "H3", "HZ1", "HZ2", "HZ3"],
    "MET": ["H", "H1", "H2", "H3"],
    "PHE": ["H", "H1", "H2", "H3"],
    "PRO": ["H", "H1", "H2", "H3"],
    "SER": ["H", "H1", "H2", "H3", "HG"],
    "THR": ["H", "H1", "H2", "H3", "HG1"],
    "TRP": ["H", "H1", "H2", "H3", "HE1"],
    "TYR": ["H", "H1", "H2", "H3", "HH"],
    "VAL": ["H", "H1", "H2", "H3"],
    "UNK": ["H", "H1", "H2", "H3"],
}


hbond_donors = np.zeros(atom_type_num, dtype=np.float32)
hbond_acceptors = np.zeros(atom_type_num, dtype=np.float32)
for atom in hbond_donor_atoms:
    hbond_donors[atom_order[atom]] = 1.0
for atom in hbond_acceptor_atoms:
    hbond_acceptors[atom_order[atom]] = 1.0


def _get_restype_atom14_hbond_donors_and_acceptors():
    restype_hbond_donors = []
    restype_hbond_acceptors = []
    for res_name in restypes:
        res_name = restype_1to3[res_name]

        res_hbond_donors = [
            1.0 if atom in hbond_donor_atoms else 0.0
            for atom in restype_name_to_atom14_names[res_name]
        ]
        res_hbond_acceptors = [
            1.0 if atom in hbond_acceptor_atoms else 0.0
            for atom in restype_name_to_atom14_names[res_name]
        ]

        restype_hbond_donors.append(res_hbond_donors)
        restype_hbond_acceptors.append(res_hbond_acceptors)

    # Update for unknown restype
    restype_hbond_donors.append([0.0] * 14)
    restype_hbond_acceptors.append([0.0] * 14)

    return restype_hbond_donors, restype_hbond_acceptors


restype_hbond_donors_atom14, restype_hbond_acceptors_atom14 = (
    _get_restype_atom14_hbond_donors_and_acceptors()
)

# Format: The list for each AA type contains chi1, chi2, chi3, chi4 in
# this order (or a relevant subset from chi1 onwards). ALA and GLY don't have
# chi angles so their chi angle lists are empty.
chi_angles_atoms = {
    "ALA": [],
    # Chi5 in arginine is always 0 +- 5 degrees, so ignore it.
    "ARG": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "NE"],
        ["CG", "CD", "NE", "CZ"],
    ],
    "ASN": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    "ASP": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    "CYS": [["N", "CA", "CB", "SG"]],
    "GLN": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "OE1"],
    ],
    "GLU": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "OE1"],
    ],
    "GLY": [],
    "HIS": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    "ILE": [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD1"]],
    "LEU": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "LYS": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "CD"],
        ["CB", "CG", "CD", "CE"],
        ["CG", "CD", "CE", "NZ"],
    ],
    "MET": [
        ["N", "CA", "CB", "CG"],
        ["CA", "CB", "CG", "SD"],
        ["CB", "CG", "SD", "CE"],
    ],
    "PHE": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "PRO": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"]],
    "SER": [["N", "CA", "CB", "OG"]],
    "THR": [["N", "CA", "CB", "OG1"]],
    "TRP": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "TYR": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    "VAL": [["N", "CA", "CB", "CG1"]],
}

# If chi angles given in fixed-length array, this matrix determines how to mask
# them for each AA type. The order is as per restype_order (see below).
chi_angles_mask = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [1.0, 1.0, 1.0, 1.0],  # ARG
    [1.0, 1.0, 0.0, 0.0],  # ASN
    [1.0, 1.0, 0.0, 0.0],  # ASP
    [1.0, 0.0, 0.0, 0.0],  # CYS
    [1.0, 1.0, 1.0, 0.0],  # GLN
    [1.0, 1.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [1.0, 1.0, 0.0, 0.0],  # HIS
    [1.0, 1.0, 0.0, 0.0],  # ILE
    [1.0, 1.0, 0.0, 0.0],  # LEU
    [1.0, 1.0, 1.0, 1.0],  # LYS
    [1.0, 1.0, 1.0, 0.0],  # MET
    [1.0, 1.0, 0.0, 0.0],  # PHE
    [1.0, 1.0, 0.0, 0.0],  # PRO
    [1.0, 0.0, 0.0, 0.0],  # SER
    [1.0, 0.0, 0.0, 0.0],  # THR
    [1.0, 1.0, 0.0, 0.0],  # TRP
    [1.0, 1.0, 0.0, 0.0],  # TYR
    [1.0, 0.0, 0.0, 0.0],  # VAL
]

# The following chi angles are pi periodic: they can be rotated by a multiple
# of pi without affecting the structure.
chi_pi_periodic = [
    [0.0, 0.0, 0.0, 0.0],  # ALA
    [0.0, 0.0, 0.0, 0.0],  # ARG
    [0.0, 0.0, 0.0, 0.0],  # ASN
    [0.0, 1.0, 0.0, 0.0],  # ASP
    [0.0, 0.0, 0.0, 0.0],  # CYS
    [0.0, 0.0, 0.0, 0.0],  # GLN
    [0.0, 0.0, 1.0, 0.0],  # GLU
    [0.0, 0.0, 0.0, 0.0],  # GLY
    [0.0, 0.0, 0.0, 0.0],  # HIS
    [0.0, 0.0, 0.0, 0.0],  # ILE
    [0.0, 0.0, 0.0, 0.0],  # LEU
    [0.0, 0.0, 0.0, 0.0],  # LYS
    [0.0, 0.0, 0.0, 0.0],  # MET
    [0.0, 1.0, 0.0, 0.0],  # PHE
    [0.0, 0.0, 0.0, 0.0],  # PRO
    [0.0, 0.0, 0.0, 0.0],  # SER
    [0.0, 0.0, 0.0, 0.0],  # THR
    [0.0, 0.0, 0.0, 0.0],  # TRP
    [0.0, 1.0, 0.0, 0.0],  # TYR
    [0.0, 0.0, 0.0, 0.0],  # VAL
    [0.0, 0.0, 0.0, 0.0],  # UNK
]


# create an array with (restype, atomtype) --> rigid_group_idx
# and an array with (restype, atomtype, coord) for the atom positions
# and compute affine transformation matrices (4,4) from one rigid group to the
# previous group
restype_atom37_to_rigid_group = np.zeros([21, 37], dtype=np.int64)
restype_atom37_mask = np.zeros([21, 37], dtype=np.float32)
restype_atom37_rigid_group_positions = np.zeros([21, 37, 3], dtype=np.float32)
restype_atom14_to_rigid_group = np.zeros([21, 14], dtype=np.int64)
restype_atom14_mask = np.zeros([21, 14], dtype=np.float32)
restype_atom14_rigid_group_positions = np.zeros([21, 14, 3], dtype=np.float32)
restype_rigid_group_default_frame = np.zeros([21, 8, 4, 4], dtype=np.float32)


def _make_rigid_group_constants():
    """Fill the arrays above."""
    for restype, restype_letter in enumerate(restypes):
        resname = restype_1to3[restype_letter]
        for atomname, group_idx, atom_position in rigid_group_atom_positions[resname]:
            atomtype = atom_order[atomname]
            restype_atom37_to_rigid_group[restype, atomtype] = group_idx
            restype_atom37_mask[restype, atomtype] = 1
            restype_atom37_rigid_group_positions[restype, atomtype, :] = atom_position

            atom14idx = restype_name_to_atom14_names[resname].index(atomname)
            restype_atom14_to_rigid_group[restype, atom14idx] = group_idx
            restype_atom14_mask[restype, atom14idx] = 1
            restype_atom14_rigid_group_positions[restype, atom14idx, :] = atom_position

    for restype, restype_letter in enumerate(restypes):
        resname = restype_1to3[restype_letter]
        atom_positions = {
            name: np.array(pos) for name, _, pos in rigid_group_atom_positions[resname]
        }

        # backbone to backbone is the identity transform
        restype_rigid_group_default_frame[restype, 0, :, :] = np.eye(4)

        # pre-omega-frame to backbone (currently dummy identity matrix)
        restype_rigid_group_default_frame[restype, 1, :, :] = np.eye(4)

        # phi-frame to backbone
        mat = _make_rigid_transformation_4x4(
            ex=atom_positions["N"] - atom_positions["CA"],
            ey=np.array([1.0, 0.0, 0.0]),
            translation=atom_positions["N"],
        )
        restype_rigid_group_default_frame[restype, 2, :, :] = mat

        # psi-frame to backbone
        mat = _make_rigid_transformation_4x4(
            ex=atom_positions["C"] - atom_positions["CA"],
            ey=atom_positions["CA"] - atom_positions["N"],
            translation=atom_positions["C"],
        )
        restype_rigid_group_default_frame[restype, 3, :, :] = mat

        # chi1-frame to backbone
        if chi_angles_mask[restype][0]:
            base_atom_names = chi_angles_atoms[resname][0]
            base_atom_positions = [atom_positions[name] for name in base_atom_names]
            mat = _make_rigid_transformation_4x4(
                ex=base_atom_positions[2] - base_atom_positions[1],
                ey=base_atom_positions[0] - base_atom_positions[1],
                translation=base_atom_positions[2],
            )
            restype_rigid_group_default_frame[restype, 4, :, :] = mat

        # chi2-frame to chi1-frame
        # chi3-frame to chi2-frame
        # chi4-frame to chi3-frame
        # luckily all rotation axes for the next frame start at (0,0,0) of the
        # previous frame
        for chi_idx in range(1, 4):
            if chi_angles_mask[restype][chi_idx]:
                axis_end_atom_name = chi_angles_atoms[resname][chi_idx][2]
                axis_end_atom_position = atom_positions[axis_end_atom_name]
                mat = _make_rigid_transformation_4x4(
                    ex=axis_end_atom_position,
                    ey=np.array([-1.0, 0.0, 0.0]),
                    translation=axis_end_atom_position,
                )
                restype_rigid_group_default_frame[restype, 4 + chi_idx, :, :] = mat


_make_rigid_group_constants()


cwd = os.path.dirname(os.path.realpath(__file__))
stereo_chemical_props_path = os.path.join(cwd, "stereo_chemical_props.txt")


def restype_bonded_atoms(self_bonds=False, atom14=True):
    with open(stereo_chemical_props_path, "r") as f:
        stereo_chemical_props = f.read()
    lines_iter = iter(stereo_chemical_props.splitlines())

    # Determine bonded residues
    if atom14:
        restype_bonded_atoms = np.zeros([21, 14, 14], dtype=np.float32)
    else:
        restype_bonded_atoms = np.zeros([21, 37, 37], dtype=np.float32)

    next(lines_iter)  # Skip header line.
    for line in lines_iter:
        if line.strip() == "-":
            break
        bond, resname, _, _ = line.split()
        atom1, atom2 = bond.split("-")

        # Get residue and atom indices
        res_idx = restype_order[restype_3to1[resname]]
        if atom14:
            atom1_idx = restype_name_to_atom14_names[resname].index(atom1)
            atom2_idx = restype_name_to_atom14_names[resname].index(atom2)
        else:
            atom1_idx = atom_order[atom1]
            atom2_idx = atom_order[atom2]

        # Symmetrically mark each bonded atom
        restype_bonded_atoms[res_idx, atom1_idx, atom2_idx] = 1.0
        restype_bonded_atoms[res_idx, atom2_idx, atom1_idx] = 1.0

    if self_bonds:
        for restype in restypes:
            res_idx = restype_order[restype]
            for atom in atom_types:
                if atom14:
                    if atom not in residue_atoms[restype_1to3[restype]]:
                        continue
                    atom_idx = restype_name_to_atom14_names[restype].index(atom)
                else:
                    atom_idx = atom_order[atom]
                restype_bonded_atoms[res_idx, atom_idx, atom_idx] = 1.0

    return restype_bonded_atoms


def _get_chi_atom_indices_and_mask(use_atom14=True):
    chi_atom_indices = []
    chi_mask = []
    for res_name in restypes:
        res_name = restype_1to3[res_name]
        res_chi_angles = chi_angles_atoms[res_name]

        # Chi mask where 1 for existing chi angle and 0 for nonexistent chi angle
        chi_mask.append([1] * len(res_chi_angles) + [0] * (4 - len(res_chi_angles)))

        # All unique atoms for chi angles
        atoms = [atom for chi in res_chi_angles for atom in chi]
        atoms = sorted(set(atoms), key=lambda x: atoms.index(x))

        # Indices of unique atoms
        if use_atom14:
            atom_indices = [
                restype_name_to_atom14_names[res_name].index(atom) for atom in atoms
            ]
        else:
            atom_indices = [atom_order[atom] for atom in atoms]

        for _ in range(7 - len(atom_indices)):
            atom_indices.append(0)

        chi_atom_indices.append(atom_indices)

    # Update for unknown restype
    chi_atom_indices.append([0] * 7)
    chi_mask.append([0] * 4)

    return chi_atom_indices, chi_mask


chi_atom_indices_atom14, chi_mask_atom14 = _get_chi_atom_indices_and_mask(
    use_atom14=True
)
chi_atom_indices_atom37, chi_mask_atom37 = _get_chi_atom_indices_and_mask(
    use_atom14=False
)


def _get_restype_atom_radius_atom14():
    restype_atom_radius = []
    for res_name in restypes:
        res_name = restype_1to3[res_name]
        atom_radius = [
            van_der_waals_radius[name[0]]
            for name in restype_name_to_atom14_names[res_name]
            if name != ""
        ]

        for _ in range(14 - len(atom_radius)):
            atom_radius.append(0)

        restype_atom_radius.append(atom_radius)

    # Update for unknown restype
    restype_atom_radius.append([0] * 14)

    return restype_atom_radius


restype_atom_radius_atom14 = _get_restype_atom_radius_atom14()