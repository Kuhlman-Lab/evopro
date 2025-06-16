""" 
Utility functions, includes saving and loading data and others.
"""

# Standard imports.
import copy
import bz2
import pickle
import _pickle as cPickle
import hashlib
from typing import Any, Dict, Union, Tuple
from Bio.PDB import MMCIFParser, PDBIO
from io import StringIO
import numpy as np
import torch
import torch.nn as nn
from omegaconf import OmegaConf

import utils.hbond_constants as rc
import utils.rigid_utils as ru

_worker_rngs = {}
_worker_rng_seed = [120723]

Array = Union[np.ndarray, torch.Tensor]

def check_omegaconf_key_exists(key: str, conf: Any) -> bool:
    UNDEFINED_VALUE = 3430323896892821
    if OmegaConf.select(conf, key, default=UNDEFINED_VALUE) == UNDEFINED_VALUE:
        return False
    return True 

def get_lengths_chains(dsobj, conf):
    lengths = []
    for chains in conf.structure_prediction.structure_pred_chains:
        c = list(chains)
        lengths.append([dsobj.get_lengths(chain) for chain in c])
    
    return lengths

def cif_to_pdb(mmcif_str):
    parser = MMCIFParser()
    cif_fh = StringIO(mmcif_str) 
    structure = parser.get_structure("structure", cif_fh)
    
    # Truncate residue names to 3 letters
    for model in structure:
        for chain in model:
            for residue in chain:
                if len(residue.resname) > 3:
                    residue.resname = residue.resname[:3]
    
    io = PDBIO()
    io.set_structure(structure)
    output = StringIO()
    io.save(output)
    pdb_string = output.getvalue()
    
    return pdb_string

def full_pickle(title: str, data: Any) -> None:
    """
    Saves the 'data' with the 'title' and adds the extension .pkl.
    """
    pikd = open(title + '.pkl', 'wb')
    pickle.dump(data, pikd)
    pikd.close()


def loosen(file: str) -> Any:
    """
    Loads and returns a pickled object.
    """
    pikd = open(file, 'rb')
    data = pickle.load(pikd)
    pikd.close()
    return data


def compressed_pickle(title: str, data: Any) -> None:
    """
    Pickle a file and then compress it into a file with extension .pbz2.
    """

    with bz2.BZ2File(title + '.pbz2', 'w') as f:
        cPickle.dump(data, f)


def decompress_pickle(file: str) -> Any:
    """
    Load any compressed pickle file.
    """
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data


def get_hash(x: str) -> str:
    """
    Looks up and returns a hash of given string.
    """
    return hashlib.sha1(x.encode()).hexdigest()


def print_timing(timing: Dict[str, float]) -> None:
    """
    Prints timing results (stored in dict) in a prettier format.
    """
    for k, v in timing.items():
        print(f'{k} took {v:.2f} sec.')
        
def merge_related_lists(list_of_pairs):
    """
    Takes a list of two-element lists and merges lists that share common elements.
    
    Args:
        list_of_pairs: List of lists, where each inner list contains two elements
        
    Returns:
        List of lists, where related elements are grouped together
    """
    # Initialize result list to store merged groups
    merged_groups = []
    
    # Convert pairs to sets for easier comparison and merging
    sets = [set(pair) for pair in list_of_pairs]
    
    while sets:
        current = sets.pop(0)  # Take the first set
        merged = False
        
        # Check against existing merged groups
        for group in merged_groups:
            if current & group:  # If there's any overlap
                group.update(current)  # Merge current into existing group
                merged = True
                break
                
        if not merged:
            merged_groups.append(current)
            
        # Check if any remaining sets can be merged with existing groups
        i = 0
        while i < len(sets):
            merged_with_existing = False
            for group in merged_groups:
                if sets[i] & group:  # If there's any overlap
                    group.update(sets[i])
                    sets.pop(i)
                    merged_with_existing = True
                    break
            if not merged_with_existing:
                i += 1
    
    # Convert sets back to sorted lists
    return [sorted(list(group)) for group in merged_groups]

def get_worker_rng() -> int:
    worker_info = torch.utils.data.get_worker_info()
    wid = worker_info.id if worker_info is not None else 0
    if wid not in _worker_rngs:
        _worker_rngs[wid] = np.random.RandomState(_worker_rng_seed[0] + wid)
    return _worker_rngs[wid]

def impute_CB(N_xyz: Array, CA_xyz: Array, C_xyz: Array) -> Array:
    """Imputes CB coordinates from N, CA, and C coordinates.

    Args:
        N_xyz (Array): Coordinates of N atom with shape: (..., 3).
        CA_xyz (Array): Coordinates of CA atom with shape: (..., 3).
        C_xyz (Array): Coordinates of C atom with shape: (..., 3).

    Returns:
        Array: Imputed CB coordinates with shape: (..., 3).
    """

    # Make sure N_xyz, CA_xyz, and C_xyz are same class
    assert type(N_xyz) is type(CA_xyz) is type(C_xyz)

    # Calculate a, b, c orientation vectors
    b = CA_xyz - N_xyz
    c = C_xyz - CA_xyz
    if isinstance(N_xyz, np.ndarray):
        a = np.cross(b, c, axis=-1)
    else:
        a = torch.cross(b, c, dim=-1)

    # Calculate CB coordinates
    CB_xyz = -0.58273431 * a + 0.56802827 * b - 0.54067466 * c + CA_xyz

    return CB_xyz

def build_sc_from_chi(
    bb_xyz: Array, aatype: Array, chi_angles: Array, chi_angle_mask: Array
) -> Tuple[Array, Array]:
    """Build side chain atoms from backbone atoms and chi angles.

    Args:
        bb_xyz (Array): 3D coordinates of the backbone atoms, shape: (Nres, 4, 3).
        aatype (Array): Amino acid type, shape (Nres,).
        chi_angles (Array): Chi angles in radians, shape: (Nres, 4).
        chi_angle_mask (Array): Mask of which chi angles are present, shape: (Nres, 4).

    Returns:
        Tuple[Array, Array]: Tuple containing 3D coordinates of each residue's atoms, shape: (Nres, 14, 3), and a mask of which atoms are present, shape: (Nres, 14).
    """
    # Make sure bb_xyz, chi_angles, and chi_angle_mask are same class
    assert type(bb_xyz) is type(aatype) is type(chi_angles) is type(chi_angle_mask)

    # Make sure the shapes are expected
    n_res = bb_xyz.shape[0]
    assert bb_xyz.shape == (n_res, 4, 3)
    assert aatype.shape == (n_res,)
    assert chi_angles.shape == (n_res, 4)
    assert chi_angle_mask.shape == (n_res, 4)

    if isinstance(bb_xyz, np.ndarray):
        # For ease, if using numpy, we convert arrays to tensors and then back.
        is_numpy = True
        bb_xyz = torch.from_numpy(bb_xyz)
        chi_angles = torch.from_numpy(chi_angles)
        chi_angle_mask = torch.from_numpy(chi_angle_mask)
    else:
        is_numpy = False

    # Convert chi_angles to sine and cosine.
    chi_angles = torch.stack(
        [torch.sin(chi_angles), torch.cos(chi_angles)],
        dim=-1,
    )

    # Get the default transformations for the chis
    default_4x4 = torch.from_numpy(rc.restype_rigid_group_default_frame[aatype])[:, -4:]
    default_r = ru.Rigid.from_tensor_4x4(default_4x4)

    # Construct and apply updates to the defaults based on chi values
    chi_rots = torch.zeros(default_r.get_rots().get_rot_mats().shape)
    chi_rots[..., 0, 0] = 1
    chi_rots[..., 1, 1] = chi_angles[..., 1]
    chi_rots[..., 1, 2] = -chi_angles[..., 0]
    chi_rots[..., 2, 1:] = chi_angles
    chi_rots = ru.Rigid(ru.Rotation(rot_mats=chi_rots), None)
    chi_frames = default_r.compose(chi_rots)

    # Build transforms for each chi directly from the bb frame
    chi2_frame_to_frame = chi_frames[:, 1]
    chi3_frame_to_frame = chi_frames[:, 2]
    chi4_frame_to_frame = chi_frames[:, 3]

    chi1_frame_to_bb = chi_frames[:, 0]
    chi2_frame_to_bb = chi1_frame_to_bb.compose(chi2_frame_to_frame)
    chi3_frame_to_bb = chi2_frame_to_bb.compose(chi3_frame_to_frame)
    chi4_frame_to_bb = chi3_frame_to_bb.compose(chi4_frame_to_frame)

    all_frames_to_bb = ru.Rigid.cat(
        [
            chi1_frame_to_bb.unsqueeze(-1),
            chi2_frame_to_bb.unsqueeze(-1),
            chi3_frame_to_bb.unsqueeze(-1),
            chi4_frame_to_bb.unsqueeze(-1),
        ],
        dim=-1,
    )

    # Build backbone frame and transform chi frames to global frame.
    bb_frames = ru.Rigid.from_3_points(bb_xyz[:, 0], bb_xyz[:, 1], bb_xyz[:, 2])
    chi_frames_to_global = bb_frames[..., None].compose(all_frames_to_bb)

    # Construct group mask for assigning atoms to chi groups.
    atom14_group_mask = torch.clamp(
        torch.from_numpy(rc.restype_atom14_to_rigid_group[aatype])[:, 5:] - 4, min=0
    )
    atom14_group_mask_oh = nn.functional.one_hot(atom14_group_mask, num_classes=4)

    # Mask transformations appropriately for each atom.
    atoms_to_global = chi_frames_to_global[:, None] * atom14_group_mask_oh
    atoms_to_global = atoms_to_global.map_tensor_fn(lambda x: torch.sum(x, dim=-1))

    # Get the literature positions for each atom.
    lit_xyz = torch.from_numpy(rc.restype_atom14_rigid_group_positions[aatype])[:, 5:]

    # Apply transformations to lit positions to get final positions.
    xyz = atoms_to_global.apply(lit_xyz)

    # Create an appropriate atom mask.
    atom_mask = copy.deepcopy(torch.from_numpy(rc.restype_atom14_mask[aatype])[:, 5:])
    mask_mask = (chi_angle_mask[:, None] * atom14_group_mask_oh).sum(-1)
    atom_mask = mask_mask * atom_mask

    # Apply mask and construct final coordinates.
    xyz = xyz * atom_mask[..., None]
    xyz = torch.cat(
        [bb_xyz, impute_CB(bb_xyz[:, 0], bb_xyz[:, 1], bb_xyz[:, 2]).unsqueeze(1), xyz],
        dim=1,
    )
    atom_mask = torch.cat(
        [torch.from_numpy(rc.restype_atom14_mask[aatype])[:, :5], atom_mask], dim=1
    )
    xyz = xyz * atom_mask[..., None]

    if is_numpy:
        return xyz.numpy(), atom_mask.numpy()
    else:
        return xyz, atom_mask