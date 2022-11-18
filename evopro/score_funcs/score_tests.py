from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of
from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
import os
import subprocess
import shutil

def basic_conf_score(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    score = -confscore2
    return score, (confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2, binder_chain="B"):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2)

    return rmsd_binder
