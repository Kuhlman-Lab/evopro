from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import os
import subprocess
import shutil
import math

## USE SCORE FUNCTIONS THAT AMRITA WROTE FOR AID PROJECT
## SCORE_BINDER_MONOMER = EVALUATE pLDDT FROM ALPHAFOLD PREDICTIONS
## SCORE_BINDER_RMSD_TO_STARTING = EVALUATE FOLDED RMSD TO ORIGINAL INPUT STRUCTURE

def score(results, dsobj, starting_pdb=None, contacts=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    # pLDDT of specific regions
    reslist_A = [x for x in residues.keys() if int(x[1:len(x)]) in range(45,46)]
    plddt_A = score_plddt_confidence(results, reslist_A, resindices, dsobj=dsobj, first_only=False)
    plddt_weight = 5*plddt_A
    # other terms
    rmsd_score = score_binder_rmsd_to_starting(pdb, contacts, starting_pdb=starting_pdb, dsobj=None)
    contact_score = score_contact(results, pdb, dsobj, orient=None) 
    score = -confscore + 0.1*rmsd_score + contact_score - plddt_weight
    #print(score)
    return score, (score, confscore, rmsd_score, contact_score, plddt_weight), pdb, results

def score_contact(results, pdb, dsobj, orient=None):
    #trying to bring termini closer to the main body
    reslist1 = ['A48','A51']
    reslist2 = ['A4','A28','A35'] #[x for x in residues.keys() if x.startswith("B")]
    # set distance cut-off as 10A for now
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dist=10, dsobj=dsobj, first_only=False)
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0
    score = -contactscore + orientation_penalty
    #print(score, (score, len(contacts), contactscore, orientation_penalty))
    return score

def score_binder_rmsd_to_starting(pdb, contacts, starting_pdb=None, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0
    # change from 4A to 2A
    rmsd_cutoff = 2.0

    with open(starting_pdb, 'r') as f:
        pdb_string_starting = f.read()
    reslist1 = contacts

    reslist2 = contacts

    # want to make this specific to the undesigned residues
    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return rmsd_potential*5


### old score function
def score_old(results, dsobj, starting_pdb=None, contacts=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    # pLDDT of specific regions
    reslist_A = [x for x in residues.keys() if int(x[1:len(x)]) in range(45,46)]
    plddt_A = score_plddt_confidence(results, reslist_A, resindices, dsobj=dsobj, first_only=False)
    plddt_weight = 2*plddt_A
    # other terms
    rmsd_score = score_binder_rmsd_to_starting(pdb, contacts, starting_pdb=starting_pdb, dsobj=None)
    contact_score = score_contact(results, pdb, dsobj, orient=None)
    score = -confscore/10 + rmsd_score - contact_score - plddt_weight
    #print(score)
    return score, (score, confscore, rmsd_score, contact_score, plddt_weight), pdb, results


if __name__=="__main__":
    print("no main functionality")
