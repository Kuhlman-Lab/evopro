from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
import os
import subprocess
from alphafold.common import protein
import shutil

def score_backbone_plddt(results):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    score = -confscore2
    return score, (confscore2), pdb, results
