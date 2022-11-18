from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.utils.write_pdb import PDBio
from folddesign.utils.calc_rmsd import RMSDcalculator
import math
import pickle
import numpy as np

def score_conf_of(ofresults, reslist, resindices):
    plddt = ofresults["confidence"].cpu().data.numpy()
    plddt = list(plddt)

    score = 0
    for res in reslist:
        resid = resindices[res]
        score = score + 100*plddt[resid]

    return score

def score_confidence_pairs_of(ofresults, pairs, resindices):
    """calculates confidence score of all pairwise residue interactions"""

    plddt = ofresults["confidence"].cpu().data.numpy()
    plddt = list(plddt)
    score = 0
    for pair in pairs:
        res1_id = resindices[pair[0]]
        res2_id = resindices[pair[1]]
        score = score + (plddt[res1_id]*plddt[res2_id]/100.0)
    return -(score)

def score_confidence_interface_of(ofresults, reslist1, reslist2, resindices):
    """calculates confidence score of all interactions between two lists of residues"""

    plddt = ofresults["confidence"].cpu().data.numpy()
    plddt = list(plddt)
    score = 0

    for res1 in reslist1:
        res1_id = resindices[res1]
        plddt_res1 = plddt[res1_id]
        for res2 in reslist2:
            res2_id = resindices[res2]
            plddt_res2 = plddt[res2_id]
            score = score + (plddt_res1*plddt_res2/100.0)
    return -(score)
