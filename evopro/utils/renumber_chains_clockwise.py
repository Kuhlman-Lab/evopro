from Bio import PDB
import os, sys
import numpy as np
import dataclasses
import io
from typing import Any, Mapping, Optional, Dict
from scipy.spatial.transform import Rotation as R
from Bio.PDB import PDBParser


sys.path.append("/proj/kuhl_lab/evopro/")
#from evopro.utils.pdb_parser import get_coordinates_pdb

def get_coordinates_pdb(pdb, fil=False):
    lines = []
    chains = []
    residues = {}
    residueindices = {}
    if fil:
        with open(pdb,"r") as f:
            i=0
            for lin in f:
                x = lin[30:38].strip(' ')
                y = lin[38:46].strip(' ')
                z = lin[46:54].strip(' ')
                l = lin.strip().split()
                if l:
                    if 'ATOM' in l[0] or 'HETATM' in l[0]:
                        resid = lin[21] + lin[22:26].strip()
                        #resid = l[4]+l[5]
                        atominfo = (l[1], l[2], (x, y, z))
                        if lin[21] not in chains:
                            print(lin, lin[21])
                            chains.append(lin[21])
                        if resid not in residues:
                            residues[resid] = [atominfo]
                        else:
                            residues[resid].append(atominfo)
                        if resid not in residueindices:
                            residueindices[resid] = i
                            i = i+1

    else:
        pdb_split = pdb.split("\n")
        i=0
        pdb_split = [x for x in pdb_split if x]
        for lin in pdb_split:
            x = lin[30:38].strip(' ')
            y = lin[38:46].strip(' ')
            z = lin[46:54].strip(' ')
            l = lin.strip().split()
            if 'ATOM' in l[0] or 'HETATM' in l[0]:
                resid = l[4]+l[5]
                atominfo = (l[1], l[2], (x, y, z))
                if l[4] not in chains:
                    chains.append(l[4])
                if resid not in residues:
                    residues[resid] = [atominfo]
                else:
                    residues[resid].append(atominfo)
                if resid not in residueindices:
                    residueindices[resid] = i
                    i = i+1 

    return chains, residues, residueindices

FeatureDict = Mapping[str, np.ndarray]
ModelOutput = Mapping[str, Any]  # Is a nested dict.

# Complete sequence of chain IDs supported by the PDB format.
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.

restypes = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
restypes_with_x = restypes + ["X"]
restype_order = {restype: i for i, restype in enumerate(restypes)}
restype_num = len(restypes)  # := 20.
restype_1to3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
                "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
                "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
restype_3to1 = {v: k for k, v in restype_1to3.items()}
restype_name_to_atom14_names = {
    "ALA": ["N", "CA", "C", "O", "CB", "", "", "", "", "", "", "", "", ""],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "", "", ""],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "", "", "", "", "", ""],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "", "", "", "", "", ""],
    "CYS": ["N", "CA", "C", "O", "CB", "SG", "", "", "", "", "", "", "", ""],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "", "", "", "", ""],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "", "", "", "", ""],
    "GLY": ["N", "CA", "C", "O", "", "", "", "", "", "", "", "", "", ""],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "", "", "", ""],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "", "", "", "", "", ""],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "", "", "", "", "", ""],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "", "", "", "", ""],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE", "", "", "", "", "", ""],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "", "", ""],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD", "", "", "", "", "", "", ""],
    "SER": ["N", "CA", "C", "O", "CB", "OG", "", "", "", "", "", "", "", ""],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2", "", "", "", "", "", "", ""],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "", ""],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "", "", "", "", "", "", ""],
    "UNK": ["", "", "", "", "", "", "", "", "", "", "", "", "", ""],
}


@dataclasses.dataclass(frozen=False)
class Protein:
    """Protein structure representation."""

    # Cartesian coordinates of atoms in angstroms. The atom types correspond to
    # residue_constants.atom_types, i.e. the first three are N, CA, CB.
    atom_positions: np.ndarray  # [num_res, num_atom_type, 3]

    # Amino-acid type for each residue represented as an integer between 0 and
    # 20, where 20 is 'X'.
    aatype: np.ndarray  # [num_res]

    # Binary float mask to indicate presence of a particular atom. 1.0 if an atom
    # is present and 0.0 if not. This should be used for loss masking.
    atom_mask: np.ndarray  # [num_res, num_atom_type]

    # Residue index as used in PDB. It is not necessarily continuous or 0-indexed.
    residue_index: np.ndarray  # [num_res]

    # 0-indexed number corresponding to the chain in the protein that this residue
    # belongs to.
    chain_index: np.ndarray  # [num_res]

    # B-factors, or temperature factors, of each residue (in sq. angstroms units),
    # representing the displacement of the residue from its ground truth mean
    # value.
    b_factors: np.ndarray  # [num_res, num_atom_type]
    
    # Mapping connecting chain ID to integer based chain index.
    chain_id_mapping: Dict[str, int]

    def __post_init__(self):
        if len(np.unique(self.chain_index)) > PDB_MAX_CHAINS:
            raise ValueError(
                f'Cannot build an instance with more than {PDB_MAX_CHAINS} chains '
                'because these cannot be written to PDB format.'
            )
def from_pdb_str(pdb_str: str, **kwargs) -> Protein:
    return from_pdb_string(pdb_str, **kwargs)
            
def from_pdb_file(pdb_file: str, **kwargs) -> Protein:
    # Obtain PDB string from PDB file.
    with open(pdb_file, 'r') as f:
        pdb_str = f.read()
        
    # Parse the string and get Protein.
    return from_pdb_string(pdb_str, **kwargs)

def to_pdb(prot: Protein) -> str:
    """Converts a `Protein` instance to a PDB string.

    Args:
        prot: The protein to convert to PDB.

    Returns:
        PDB string.
    """
    res_1to3 = lambda r: restype_1to3.get(restypes_with_x[r], 'UNK')

    pdb_lines = []

    atom_mask = prot.atom_mask
    aatype = prot.aatype
    atom_positions = prot.atom_positions
    residue_index = prot.residue_index.astype(np.int32)
    chain_index = prot.chain_index.astype(np.int32)
    b_factors = prot.b_factors

    if np.any(aatype > restype_num):
        raise ValueError('Invalid aatypes.')

    # Construct a mapping from chain integer indices to chain ID strings.
    chain_ids = {}
    for i in np.unique(chain_index):  # np.unique gives sorted output.
        if i >= PDB_MAX_CHAINS:
            raise ValueError(
                f'The PDB format supports at most {PDB_MAX_CHAINS} chains.')
        chain_ids[i] = PDB_CHAIN_IDS[i]

    pdb_lines.append('MODEL     1')
    atom_index = 1
    last_chain_index = chain_index[0]
    # Add all atom sites.
    for i in range(aatype.shape[0]):
        # Close the previous chain if in a multichain PDB.
        if last_chain_index != chain_index[i]:
            pdb_lines.append(_chain_end(
                atom_index, res_1to3(aatype[i - 1]), chain_ids[chain_index[i - 1]],
                residue_index[i - 1]))
            last_chain_index = chain_index[i]
            atom_index += 1  # Atom index increases at the TER symbol.

        res_name_3 = res_1to3(aatype[i])
        for atom_name, pos, mask, b_factor in zip(
                restype_name_to_atom14_names[res_name_3], atom_positions[i], atom_mask[i], b_factors[i]):
            if mask < 0.5:
                continue

            record_type = 'ATOM'
            name = atom_name if len(atom_name) == 4 else f' {atom_name}'
            alt_loc = ''
            insertion_code = ''
            occupancy = 1.00
            element = atom_name[0]  # Protein supports only C, N, O, S, this works.
            charge = ''

            # PDB is a columnar format, every space matters here!
            atom_line = (f'{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}'
                         f'{res_name_3:>3} {chain_ids[chain_index[i]]:>1}'
                         f'{residue_index[i]:>4}{insertion_code:>1}   '
                         f'{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}'
                         f'{occupancy:>6.2f}{b_factor:>6.2f}          '
                         f'{element:>2}{charge:>2}')
            pdb_lines.append(atom_line)
            atom_index += 1

    # Close the final chain.
    pdb_lines.append(_chain_end(atom_index, res_1to3(aatype[-1]),
                                chain_ids[chain_index[-1]], residue_index[-1]))
    pdb_lines.append('ENDMDL')
    pdb_lines.append('END')

    # Pad all lines to 80 characters.
    pdb_lines = [line.ljust(80) for line in pdb_lines]
    return '\n'.join(pdb_lines) + '\n'  # Add terminating newline.

def calculate_chain_angle(chain, centroid):
    angles = np.arctan2([atom.coord[1] - centroid[1] for atom in chain.get_atoms()],
                       [atom.coord[0] - centroid[0] for atom in chain.get_atoms()])
    return np.mean(angles)

def orient_chains(pdb_file, output_filename, clockwise=True, start_char='A'):
    # Load the PDB structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Get the coordinates of the atoms to calculate the centroid
    coords = []
    for model in structure:
        for chain in model:
            for atom in chain.get_atoms():
                coords.append(atom.coord)

    coords = np.array(coords)
    centroid = np.mean(coords, axis=0)

    # Sort chains based on the average angle of all atoms in each chain
    chains = list(structure.get_chains())
    chains.sort(key=lambda chain: calculate_chain_angle(chain, centroid))

    if not clockwise:
        chains.reverse()

    # Renumber chains as numbers starting from 1
    for i, chain in enumerate(chains):
        chain.id = str(i + 1)

    # Determine the starting character for chain identifiers
    current_char = ord(start_char)
    # Renumber and rename chains as letters starting with A
    for chain in chains:
        chain.id = chr(current_char)
        current_char += 1

    # Save the modified structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_filename)

    print(f'Renumbered and renamed chains and saved the result to {output_filename}')

def get_sym_vec(pdb_str, order=3, chains=[]):
    protein = from_pdb_str(pdb_str)
    
    if chains == [] or chains == ['']:
        chains = list(protein.chain_id_mapping.keys())
        
    res1_ca = np.stack([protein.atom_positions[:, 1, :][protein.chain_index == i][9] for i in range(order)])
    res2_ca = np.stack([protein.atom_positions[:, 1, :][protein.chain_index == i][19] for i in range(order)])
    
    center1 = np.mean(res1_ca, axis=0)
    center2 = np.mean(res2_ca, axis=0)
    
    sym_vec = center2 - center1
    
    return sym_vec


def get_rot_mat(sym_vec, axis='z', invert=False):
    norm_sym_vec = sym_vec / np.sqrt(np.sum(sym_vec ** 2))
    axis_vec = np.array([0, 0, 0], dtype=np.float32)
    if axis == 'z':
        axis_vec[2] = 1.
    elif axis == 'y':
        axis_vec[1] = 1.
    else:
        axis_vec[0] = 1.
        
    if invert:
        axis_vec *= -1
    
    
    #print(norm_sym_vec, axis_vec)
    rotation_axis = np.cross(norm_sym_vec, axis_vec)
    rotation_ang = np.arccos(np.dot(norm_sym_vec, axis_vec))
    
    #print(rotation_ang, rotation_axis, np.sum(rotation_axis ** 2), np.sqrt(np.sum(rotation_axis ** 2)))
    if np.sum(rotation_axis ** 2) == 0:
        #print("here")
        rot_mat = np.identity(3)
        if invert:
            rot_mat[-1][-1] *= -1
    else:
        rot_mat = R.from_rotvec(rotation_ang * rotation_axis / np.sqrt(np.sum(rotation_axis ** 2))).as_matrix()
    
    return rot_mat


def do_rotation(atom_positions, rot_mat):
    x = atom_positions @ rot_mat.T[None]
    return x


def zalign(pdb_str, order, axis='z', invert=False, sym_chains=[]):

    sym_vec = get_sym_vec(pdb_str, order, chains=sym_chains)
    
    #("Sym vec", sym_vec)
        
    rot_mat = get_rot_mat(sym_vec, axis, invert=invert)
    
    #print("Rot mat", rot_mat)

    protein = from_pdb_str(pdb_str)
    #print(protein.atom_positions)
    
    protein.atom_positions = do_rotation(protein.atom_positions, rot_mat)
    #print(protein.atom_positions)
        
    pdb_str = to_pdb(protein)

    return pdb_str


def calculate_z_range(pdb_str):
    protein = from_pdb_str(pdb_str)
    CA_z = protein.atom_positions[:, 1, 2]
    return min(CA_z), max(CA_z)


def from_pdb_string(pdb_str: str, model_idx: int = 0, chain_id: Optional[str] = None, discard_water: bool = True, mse_to_met: bool = False, ignore_non_std: bool = True) -> Protein:
    """Takes a PDB string and constructs a Protein object.

    WARNING: All non-standard residue types will be converted into UNK. All
        non-standard atoms will be ignored.

    Args:
        pdb_str: The contents of the pdb file
        model_idx: The specific model in the PDB file that will be
            parsed. This is 0-indexed. Default is 0.
        chain_id: If chain_id is specified (e.g. A), then only that chain
            is parsed. Otherwise all chains are parsed.
        discard_water: Boolean specifying whether to ignore water molecules.
            Default is True.
        mse_to_met: Boolean specifying whether to convert MSE residues to MET residues.
            Default is False.
        ignore_non_std: Boolean specifying whether to ignore nonstandard residue types.
            If False, then they will be converted to UNK. Default is True.

    Returns:
        A new `Protein` parsed from the pdb contents.
    """
    pdb_fh = io.StringIO(pdb_str)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('none', pdb_fh)
    models = list(structure.get_models())
    if model_idx is not None and model_idx > len(models) - 1:
        raise ValueError(
            f'Requested model index is out of range. Found {len(models)} models.'
        )
    elif model_idx is not None:
        model = models[model_idx]
    else:
        model = models[0]

    atom_positions = []
    aatype = []
    atom_mask = []
    residue_index = []
    chain_ids = []
    b_factors = []
    insertion_code_offset = 0
    for chain in sorted(model, key=lambda x: x.id):
        if chain_id is not None and chain.id != chain_id:
            continue
        for res in sorted(chain, key=lambda x: x.id[1]):
            # Discard water residues.     
            if discard_water:
                if res.resname == 'HOH':
                    continue
            
            # Convert MSE residues to MET by changing only necessary fields.
            if mse_to_met:
                if res.resname == 'MSE':
                    res.resname = 'MET'
                    for atom in res:
                        if atom.name == 'SE':
                            atom.name = 'SD'
                                    
            # Ignore non-standard residues
            res_shortname = restype_3to1.get(res.resname, 'X')
            if ignore_non_std:
                if res_shortname == 'X':
                    continue
            
            # Increment residue index offset if insertion code is detected.
            if res.id[2] != ' ':
                insertion_code_offset += 1
            
            restype_idx = restype_order.get(
                res_shortname, restype_num)
            pos = np.full((14, 3), fill_value=(0))
            mask = np.zeros((14,))
            res_b_factors = np.zeros((14,))
            for atom in res:
                if atom.name not in restype_name_to_atom14_names[res.resname]:
                    continue
                pos[restype_name_to_atom14_names[res.resname].index(atom.name)] = atom.coord
                mask[restype_name_to_atom14_names[res.resname].index(atom.name)] = 1.
                res_b_factors[restype_name_to_atom14_names[res.resname].index(atom.name)] = atom.bfactor
            if np.sum(mask) < 0.5:
                # If no known atom positions are reported for the residue then skip it.
                continue
            
            # Update protein-level lists
            aatype.append(restype_idx)
            atom_positions.append(pos)
            atom_mask.append(mask)
            residue_index.append(res.id[1] + insertion_code_offset)
            chain_ids.append(chain.id)
            b_factors.append(res_b_factors)

    # Chain IDs are usually characters so map these to ints.
    unique_chain_ids = np.unique(chain_ids)
    chain_id_mapping = {cid: n for n, cid in enumerate(unique_chain_ids)}
    chain_index = np.array([chain_id_mapping[cid] for cid in chain_ids])

    #print(atom_positions)
    return Protein(
            atom_positions=np.array(atom_positions),
            atom_mask=np.array(atom_mask),
            aatype=np.array(aatype),
            residue_index=np.array(residue_index),
            chain_index=chain_index,
            b_factors=np.array(b_factors),
            chain_id_mapping=chain_id_mapping
    )


def from_pdb_file(pdb_file: str, **kwargs) -> Protein:
    # Obtain PDB string from PDB file.
    with open(pdb_file, 'r') as f:
        pdb_str = f.read()
        
    # Parse the string and get Protein.
    return from_pdb_string(pdb_str, **kwargs)


def _chain_end(atom_index, end_resname, chain_name, residue_index) -> str:
    chain_end = 'TER'
    return (f'{chain_end:<6}{atom_index:>5}      {end_resname:>3} '
        f'{chain_name:>1}{residue_index:>4}')


def to_pdb(prot: Protein) -> str:
    """Converts a `Protein` instance to a PDB string.

    Args:
        prot: The protein to convert to PDB.

    Returns:
        PDB string.
    """
    res_1to3 = lambda r: restype_1to3.get(restypes_with_x[r], 'UNK')

    pdb_lines = []

    atom_mask = prot.atom_mask
    aatype = prot.aatype
    atom_positions = prot.atom_positions
    residue_index = prot.residue_index.astype(np.int32)
    chain_index = prot.chain_index.astype(np.int32)
    b_factors = prot.b_factors

    if np.any(aatype > restype_num):
        raise ValueError('Invalid aatypes.')

    # Construct a mapping from chain integer indices to chain ID strings.
    chain_ids = {}
    for i in np.unique(chain_index):  # np.unique gives sorted output.
        if i >= PDB_MAX_CHAINS:
            raise ValueError(
                f'The PDB format supports at most {PDB_MAX_CHAINS} chains.')
        chain_ids[i] = PDB_CHAIN_IDS[i]

    pdb_lines.append('MODEL     1')
    atom_index = 1
    last_chain_index = chain_index[0]
    # Add all atom sites.
    for i in range(aatype.shape[0]):
        # Close the previous chain if in a multichain PDB.
        if last_chain_index != chain_index[i]:
            pdb_lines.append(_chain_end(
                atom_index, res_1to3(aatype[i - 1]), chain_ids[chain_index[i - 1]],
                residue_index[i - 1]))
            last_chain_index = chain_index[i]
            atom_index += 1  # Atom index increases at the TER symbol.

        res_name_3 = res_1to3(aatype[i])
        for atom_name, pos, mask, b_factor in zip(
                restype_name_to_atom14_names[res_name_3], atom_positions[i], atom_mask[i], b_factors[i]):
            if mask < 0.5:
                continue

            record_type = 'ATOM'
            name = atom_name if len(atom_name) == 4 else f' {atom_name}'
            alt_loc = ''
            insertion_code = ''
            occupancy = 1.00
            element = atom_name[0]  # Protein supports only C, N, O, S, this works.
            charge = ''

            # PDB is a columnar format, every space matters here!
            atom_line = (f'{record_type:<6}{atom_index:>5} {name:<4}{alt_loc:>1}'
                         f'{res_name_3:>3} {chain_ids[chain_index[i]]:>1}'
                         f'{residue_index[i]:>4}{insertion_code:>1}   '
                         f'{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}'
                         f'{occupancy:>6.2f}{b_factor:>6.2f}          '
                         f'{element:>2}{charge:>2}')
            pdb_lines.append(atom_line)
            atom_index += 1

    # Close the final chain.
    pdb_lines.append(_chain_end(atom_index, res_1to3(aatype[-1]),
                                chain_ids[chain_index[-1]], residue_index[-1]))
    pdb_lines.append('ENDMDL')
    pdb_lines.append('END')

    # Pad all lines to 80 characters.
    pdb_lines = [line.ljust(80) for line in pdb_lines]
    return '\n'.join(pdb_lines) + '\n'  # Add terminating newline.

def renumber_pdb(pdb_path):

    if os.path.isfile(pdb_path):
    
        protein = from_pdb_file(pdb_path)

        n_chains = max(protein.chain_index)
        for i in range(n_chains + 1):
            n_res = np.sum(protein.chain_index == i)
            protein.residue_index[protein.chain_index == i] = np.arange(1, n_res + 1)

        pdb_str = to_pdb(protein)
        return pdb_str

def renumber_chains_clockwise(pdb_string):
    print("renumbering chains clockwise...")
    
    chains, _, _ = get_coordinates_pdb(pdb_string)
    #print(chains)
    order = len(chains)
    
    #align to z axis
    pdb_str = zalign(pdb_string, order, axis="z", invert=False, sym_chains=chains)
    if os.path.exists("temp1.pdb"):
        os.remove("temp1.pdb")
    with open("temp1.pdb", "w") as f:
        f.write(pdb_str)
    
    #check if correct orientation
    min_z, max_z = calculate_z_range(pdb_str)
    mid_z = (max_z + min_z) / 2
    #print("Zrange", min_z, mid_z, max_z)
    
    chains2, residues2, _ = get_coordinates_pdb(pdb_str)

    #get z coord of first atom of first residue
    res1_z = float(residues2["A1"][0][-1][2])
    if res1_z > mid_z:
        #flip
        print("Flipping")
        pdb_str = zalign(pdb_str, order, axis="z", invert=True, sym_chains=chains2)
    
    if os.path.exists("temp.pdb"):
        os.remove("temp.pdb")
        
    with open("temp.pdb", "w") as f:
        f.write(pdb_str)
    if os.path.exists("temp_clockwise.pdb"):
        os.remove("temp_clockwise.pdb")
    
    #orient chains clockwise
    orient_chains("temp.pdb", "temp_clockwise.pdb", False, "A")   

    #renumber before returning
    return renumber_pdb("temp_clockwise.pdb")

def renumber_chains_clockwise_test(pdb_string, invert = False):
    print("renumbering chains clockwise...")
    
    chains, _, _ = get_coordinates_pdb(pdb_string)
    #print(chains)
    chains = set([x[0] for x in chains])
    #print(chains)    
    order = len(chains)
    print("order", order)
    #align to z axis
    pdb_str = zalign(pdb_string, order, axis="z", invert=False, sym_chains=chains)
    if os.path.exists("temp1.pdb"):
        os.remove("temp1.pdb")
    with open("temp1.pdb", "w") as f:
        f.write(pdb_str)
    
    #check if correct orientation
    min_z, max_z = calculate_z_range(pdb_str)
    mid_z = (max_z + min_z) / 2
    #print("Zrange", min_z, mid_z, max_z)
    
    chains2, residues2, _ = get_coordinates_pdb(pdb_str)

    #get z coord of first atom of first residue
    if invert:
        #flip
        print("Flipping")
        pdb_str = zalign(pdb_str, order, axis="z", invert=True, sym_chains=chains2)
    
    if os.path.exists("temp.pdb"):
        os.remove("temp.pdb")
        
    with open("temp.pdb", "w") as f:
        f.write(pdb_str)
    if os.path.exists("temp_clockwise.pdb"):
        os.remove("temp_clockwise.pdb")
    
    #orient chains clockwise
    orient_chains("temp.pdb", "temp_clockwise.pdb", False, "A")   

    #renumber before returning
    return renumber_pdb("temp_clockwise.pdb")
    
if __name__ == "__main__":
    #path = "./"
    #pdb1 = path + "design.pdb"
    pdb1 = "/work/users/a/m/amritan/evopro_tests/rmsd/temp1.pdb"
    with open(pdb1, "r") as f:
        pdb_str = f.read()
    
    chains2, residues2, _ = get_coordinates_pdb(pdb_str)
    print(chains2)
    zalign(pdb_str, 3, axis='z', invert=True, sym_chains=chains2)
