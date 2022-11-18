import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
import os
import subprocess
import shutil

def score_tcr_contacts(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    contacts_list = [34, 38, 42, 61, 62, 77, 82, 113, 116, 182, 183]
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in contacts_list]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]

    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, score_cap = 20)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)

    score = -contactscore*100 + confidencescore/10.0
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_tcr_binder(results, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    confscore2 = confscore2/len(reslist2)
    score = -confscore2
    return 10*score, (confscore2), pdb, results

def score_tcr_complexbinder(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_tcr_contacts(results)
    else:
        return score_tcr_binder(results)

def score_binder_rmsd(pdb1, pdb2, with_linker=False):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False)
    reslist1 = [x for x in residues1.keys() if x.startswith("B")]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2)

    return 100*rmsd
