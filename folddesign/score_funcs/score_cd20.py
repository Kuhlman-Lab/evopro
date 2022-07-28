from score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues
from folddesign.utils.pdb_parser import get_coordinates_pdb
import os
import subprocess
from alphafold.common import protein
import shutil

def score_cd20_1(results):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    reslist1 = [x for x in residues.keys()]
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)

    reslist3 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist3, resindices, fil = False)

    score = -contactscore*100 + confidencescore - confscore2*10
    return score, (contactscore, confidencescore, confscore2), contacts, pdb, results

def score_cd20_2(results):
    pdb = protein.to_pdb(results['unrelaxed_protein'])

if __name__ == "__main__":
    print("no main functionality")
