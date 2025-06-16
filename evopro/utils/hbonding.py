import numpy as np
import itertools
from scipy.spatial.distance import cdist
from typing import Dict, Union
import torch
import subprocess

import utils.hbond_constants as rc
from utils.parsing_utils import *
from objects.protein import Protein

Array = Union[np.ndarray, torch.Tensor]

REDUCE_EXE = "/proj/kuhl_lab/reduce/reduce_src/reduce"
HET_DICT = "/proj/kuhl_lab/reduce/reduce_wwPDB_het_dict.txt"

def score_hydrogen_bonds(result, reslist1, reslist2, max_don_h_dist = 1.1):
    """returns a list of pairs of residues from reslists that are making hydrogen bonds"""
    pdb_str = result['pdb']
    protein = Protein.from_pdb_string(pdb_str, mse_to_met=True)
    reduced_protein = reduce(protein)
    hbond_pot = HBondPotential(max_don_h_dist=max_don_h_dist)
    hbond_matrix = hbond_pot.score(reduced_protein)
    
    pairs = []
    score = 0
    chains, residues, resindices = get_coordinates_pdb(pdb_str)

    for i in range(len(reslist1)):
        for j in range(len(reslist2)):
            if hbond_matrix[resindices[reslist1[i]], resindices[reslist2[j]]] == 1:
                pairs.append((reslist1[i], reslist2[j]))
                score += 1
    return pairs, score

def reduce(protein, his=True, flip=True):
    # Build the command to run Reduce
    reduce_cmd = [REDUCE_EXE, "-q", "-DB", HET_DICT, "-"]
    if his:
        reduce_cmd += ["-HIS"]
    if flip:
        reduce_cmd += ["-FLIP"]
    pop = subprocess.Popen(
        reduce_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    # Load the Protein into a string
    pdb_str = protein.to_pdb(unk_to_gly=True)
    # Pass the string to Reduce
    out_pdb = pop.communicate(input=pdb_str)[0]
    # Load the output PDB into a new Protein
    new_protein = Protein.from_pdb_string(out_pdb, discard_Hs = False)
    
    return new_protein

class HBondPotential:
    """
    Uses Baker-Hubbard algorithm to detect H-bonds between residues.

    Source:
        E. N. Baker, R. E. Hubbard, “Hydrogen bonding in globular proteins,”
        Progress in Biophysics and Molecular Biology, vol. 44, pp. 97–179,
        January 1984. doi: 10.1016/0079-6107(84)90007-5
    """

    def __init__(self, max_don_h_dist=1.1) -> None:
        # Max D - H bond length for pair assignment
        self.max_don_h_dist = max_don_h_dist

        # Numerical stability epsilon
        self.eps = 1e-8

        # Baker-Hubbard params
        self.max_HA_dist = 2.5
        self.min_DAH_angle = 120.0

    def get_features(self, protein: Protein) -> Dict[str, np.ndarray]:
        """
        Retrieve appropriate atom info for scoring.
        Outputs have shape [N, ...] where N is the number of donor/acceptor groups to check.
        """
        # Compile Acceptor/Donor masks for each residue
        acc_mask = np.zeros((protein.n_res, 14), dtype=bool)
        don_mask = np.zeros((protein.n_res, 14), dtype=bool)
        for res_idx in range(protein.n_res):
            restype = protein.aatype[res_idx]
            acc_mask[res_idx] = np.array(
                rc.restype_hbond_acceptors_atom14[restype]
            ).astype(bool)
            don_mask[res_idx] = np.array(
                rc.restype_hbond_donors_atom14[restype]
            ).astype(bool)

        acc_mask *= protein.atom27_mask[:, :14].astype(bool)
        don_mask *= protein.atom27_mask[:, :14].astype(bool)

        # Collect donor/acceptor xyz and res/atom info
        acc_xyz = protein.atom27_xyz[:, :14][acc_mask]
        acc_res, acc_is_bb = np.where(acc_mask)
        acc_is_bb = acc_is_bb < 4

        don_xyz = protein.atom27_xyz[:, :14][don_mask]  # [N_DON, 3]
        don_res, don_is_bb = np.where(don_mask)
        don_is_bb = don_is_bb < 4
        h_atom_xyz = protein.atom27_xyz[don_res, 14:]  # [N_DON, 13, 3]

        # Find H atom(s) for each donor using distance threshold
        # It would be nice to vectorize this, but it's fast enough already
        new_don_xyz, new_don_h_xyz, new_don_res, new_don_is_bb = [], [], [], []
        for i in range(don_xyz.shape[0]):
            c_don_xyz = don_xyz[i, ...][None, ...]  # [1, 3]
            c_h_atom_xyz = h_atom_xyz[i, ...]  # [13, 3]
            dist = cdist(c_don_xyz, c_h_atom_xyz)  # [1, 13]
            d_mask = (dist < self.max_don_h_dist) * protein.atom27_mask[don_res[i], 14:]
            c_h_atom_xyz = c_h_atom_xyz[np.squeeze(d_mask, axis=0) > 0]
            n_h_atoms = c_h_atom_xyz.shape[0]
            c_don_xyz = np.repeat(c_don_xyz, repeats=n_h_atoms, axis=0)
            new_don_xyz.append(c_don_xyz)
            new_don_h_xyz.append(c_h_atom_xyz)
            new_don_res.extend([don_res[i]] * n_h_atoms)
            new_don_is_bb.extend([don_is_bb[i]] * n_h_atoms)

        don_res = np.array(new_don_res, dtype=np.int32)
        don_is_bb = np.array(new_don_is_bb, dtype=bool)

        don_xyz = np.concatenate(new_don_xyz)
        don_h_xyz = np.concatenate(new_don_h_xyz)
        assert (
            don_xyz.shape[0]
            == don_h_xyz.shape[0]
            == don_res.shape[0]
            == don_is_bb.shape[0]
        )
        assert acc_xyz.shape[0] == acc_res.shape[0] == acc_is_bb.shape[0]

        # Get indices for all-vs-all search
        n_acc, n_don = acc_res.shape[0], don_res.shape[0]
        combos = np.array(
            [
                i
                for i in itertools.product(
                    np.arange(
                        n_acc,
                    ),
                    np.arange(
                        n_don,
                    ),
                )
            ]
        )

        acc_xyz = acc_xyz[combos[:, 0]]
        don_xyz = don_xyz[combos[:, 1]]
        don_h_xyz = don_h_xyz[combos[:, 1]]
        acc_res = acc_res[combos[:, 0]]
        don_res = don_res[combos[:, 1]]
        acc_is_bb = acc_is_bb[combos[:, 0]]
        don_is_bb = don_is_bb[combos[:, 1]]

        features = {
            "acc": acc_xyz,
            "don": don_xyz,
            "don_h": don_h_xyz,
            "acc_res": acc_res,
            "don_res": don_res,
            "acc_is_bb": acc_is_bb,
            "don_is_bb": don_is_bb,
        }
        return features

    def score(
        self, protein: Protein, include_bb_bb: bool = False, include_sc_bb: bool = False
    ) -> np.ndarray:
        def distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
            return np.sqrt(np.sum(np.square(b - a), axis=-1))

        features = self.get_features(protein)

        mask = np.ones_like(features["acc_is_bb"], dtype=bool)

        if not include_bb_bb:
            mask *= ~np.logical_and(features["acc_is_bb"], features["don_is_bb"])

        if not include_sc_bb:
            mask *= ~np.logical_xor(features["acc_is_bb"], features["don_is_bb"])

        for key in features:
            features[key] = features[key][mask]

        # calculate H-A dist
        HA_dist = distance(features["acc"], features["don_h"])
        # Filter DA dist to avoid weird self-bonding behavior
        DA_dist = distance(features["acc"], features["don"])

        # calculate D-H-A angle
        D_to_H = robust_normalize(features["don"] - features["don_h"], eps=self.eps)
        H_to_A = robust_normalize(features["acc"] - features["don_h"], eps=self.eps)
        DHA_angle = np.sum(D_to_H * H_to_A, axis=-1)
        DHA_angle = np.clip(DHA_angle, -1 + self.eps, 1 - self.eps)
        DHA_angle = np.arccos(DHA_angle)
        DHA_angle = np.rad2deg(DHA_angle)

        # Filter bonds
        HA_dist_flag = (
            (HA_dist <= self.max_HA_dist) * (HA_dist > self.eps) * (DA_dist > self.eps)
        )
        DHA_angle_flag = DHA_angle > self.min_DAH_angle
        total_flag = HA_dist_flag * DHA_angle_flag

        # Fill in HBond residue-wise matrix
        # [N, N] with axis 1 for acceptor and axis 2 for donor
        pair_matrix = np.zeros((protein.n_res, protein.n_res), dtype=np.int32)
        acc_resid = features["acc_res"][total_flag]
        don_resid = features["don_res"][total_flag]
        pair_matrix[acc_resid, don_resid] = 1
        return pair_matrix
    
def robust_norm(
    array: Array, axis: int = -1, p: float = 2.0, eps: float = 1e-8
) -> Array:
    """Computes robust p-norm of vectors.

    Args:
        array (Array): Array containing vectors to compute norm.
        axis (int, optional): Axis of array to norm. Defaults to -1.
        p (float, optional): p value for the p-norm to perform. Defaults to 2.
        eps (float, optional): Epsilon for robust norm computation. Defaults to 1e-8.

    Returns:
        Array: Norm of axis of array
    """
    assert p >= 1, "p must be greater than or equal to 1"

    if isinstance(array, np.ndarray):
        return (np.sum(np.abs(array) ** p, axis=axis) + eps) ** (1 / p)
    else:
        return (torch.sum(torch.abs(array) ** p, dim=axis) + eps) ** (1 / p)


def robust_normalize(
    array: Array, axis: int = -1, p: float = 2.0, eps: float = 1e-8
) -> Array:
    """Computes robust p-normalization of vectors.

    Args:
        array (Array): Array containing vectors to normalize.
        axis (int, optional): Axis of array to normalize. Defaults to -1.
        p (float, optional): p value for the p-norm to perform. Defaults to 2.
        eps (float, optional): Epsilon for robust norma computation. Defaults to 1e-8.

    Returns:
        Array: Normalized array
    """
    assert p >= 1, "p must be greater than or equal to 1"

    if isinstance(array, np.ndarray):
        return array / np.expand_dims(
            robust_norm(array, axis=axis, p=p, eps=eps), axis=axis
        )
    else:
        return array / robust_norm(array, axis=axis, p=p, eps=eps).unsqueeze(axis)
    
if __name__ == "__main__":
    print("Testing HBondPotential")
    pdb_file = "ppi.pdb"
    protein = Protein.from_pdb_file(pdb_file, mse_to_met=True)
    reduced_protein = reduce(protein)
    hbond_pot = HBondPotential()
    hbond_matrix = hbond_pot.score(reduced_protein)
    print("Hbond matrix shape:", hbond_matrix.shape)
    # np.set_printoptions(threshold=sys.maxsize)
    # print("Hbond matrix:", hbond_matrix)
    