from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import os
import subprocess
import shutil
import math

def score_binder(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_complex(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0
    # pLDDT score for each chain, total = sum but weigh chain B higher
    reslist_A = [x for x in residues.keys() if x.startswith("A")]
    reslist_B = reslist2
    plddt_A = score_plddt_confidence(results, reslist_A, resindices, dsobj=dsobj, first_only=False)
    plddt_B = score_plddt_confidence(results, reslist_B, resindices, dsobj=dsobj, first_only=False)
    plddt_score = -plddt_A - 2*plddt_B
    # total score has contact score and plddt score
    score = -contactscore + orientation_penalty + plddt_score
    print(score, (score, len(contacts), contactscore, orientation_penalty, plddt_score))
    return score, (score, len(contacts), contactscore, orientation_penalty, plddt_score), contacts, pdb, results

def score_binder_monomer(results, dsobj, starting_pdb=None, contacts=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    rmsd_score = score_binder_rmsd_to_starting(pdb, starting_pdb, dsobj=dsobj, contacts=contacts)
    score = -confscore2/10 + rmsd_score
    print(score)
    return score, (score, confscore2, rmsd_score), pdb, results

def score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)

    return rmsd_binder*5

def score_binder_rmsd_to_starting(pdb, starting_pdb=None, dsobj=None, contacts=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

    with open(starting_pdb, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x in contacts]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys() if x in contacts]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return rmsd_potential*5


if __name__=="__main__":
    print("no main functionality")
