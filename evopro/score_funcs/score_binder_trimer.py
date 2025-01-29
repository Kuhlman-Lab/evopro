import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score, Rg
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient
import os
import subprocess
import shutil
import math
import time

def score_binder(results, dsobj, contacts=None, distance_cutoffs=None):
    from alphafold.common import protein
    #print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex(results, dsobj, contacts, distance_cutoffs)
    else:
        return score_binder_monomer(results, dsobj)
    
def score_binder_3(results, dsobj, contacts=None, distance_cutoffs=None):
    from alphafold.common import protein
    #print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_3(results, dsobj, contacts, distance_cutoffs)
    else:
        return score_binder_monomer_rG(results, dsobj)
    
def score_binder_3_efficient(results, dsobj, contacts=None, distance_cutoffs=None):
    from alphafold.common import protein
    #print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, _, _ = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_3_efficient(results, dsobj, contacts, distance_cutoffs)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_complex(results, dsobj, contacts, distance_cutoffs):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
        for contact in bonus_contacts:
            if (contact[0][0:1] == 'A' or contact[0][0:1] == 'B' or contact[0][0:1] == 'C') and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if (contact[1][0:1] == 'A' or contact[1][0:1] == 'B' or contact[1][0:1] == 'C') and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    print(penalty_resids)
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
        for contact in penalty_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
    
    score = -contactscore + penalty + bonus
    print(score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty))
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results

def score_binder_complex_3(results, dsobj, contacts, distance_cutoffs):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("D") or x.startswith("E") or x.startswith("F")]
    contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0], contact_cap=36*3)
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1], contact_cap=36*3)
        for contact in bonus_contacts:
            if (contact[0][0:1] == 'A' or contact[0][0:1] == 'B' or contact[0][0:1] == 'C') and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if (contact[1][0:1] == 'A' or contact[1][0:1] == 'B' or contact[1][0:1] == 'C') and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    print(penalty_resids)
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2], contact_cap=36*3)
        for contact in penalty_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
    
    score = -contactscore + penalty + bonus
    print(score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty))
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results

def score_binder_complex_3_efficient(results, dsobj, contacts, distance_cutoffs):
    start = time.time()
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    #print(contacts)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("D") or x.startswith("E") or x.startswith("F")]
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0], contact_cap=36*3)
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted_efficient(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1], contact_cap=36*3)
        for contact in bonus_contacts:
            if (contact[0][0:1] == 'A' or contact[0][0:1] == 'B' or contact[0][0:1] == 'C') and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if (contact[1][0:1] == 'A' or contact[1][0:1] == 'B' or contact[1][0:1] == 'C') and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    #print(penalty_resids)
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted_efficient(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2], contact_cap=36*3)
        for contact in penalty_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
    
    score = -contactscore + penalty + bonus
    print("Time to score:", time.time()-start)
    print("Scores:", score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty))
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results

def score_binder_monomer_rG(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    rg = Rg(pdb, chnid="A")
    
    # apply flat-bottom quadratic-shaped potential function
    spring_constant = 10.0 
    rg_potential = 0
    if rg > 12:
        rg_potential = spring_constant*math.pow(rg - 12, 2)
    
    score = -confscore2/10 + rg_potential
    print(score, confscore2, rg, rg_potential)
    return score, (score, confscore2, rg, rg_potential), pdb, results

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    print(score)
    return score, (score, confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2, binder_chain="D", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    print(reslist1, reslist2)
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=None)

    return rmsd_binder*5

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 5.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    print("in rmsd to starting", reslist1, reslist2)

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return rmsd_potential*5

def score_binder_old(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_old(results, dsobj, contacts)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_complex_old(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    pae_score = score_pae_confidence_pairs(results, contacts, resindices, dsobj=dsobj)
    pae_score = pae_score/(len(contacts)*10)
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0
    score = -contactscore + pae_score + orientation_penalty
    return score, (score, contactscore, pae_score, orientation_penalty), contacts, pdb, results

if __name__=="__main__":
    print("no main functionality. just used for testing.")
    #with open("/work/users/a/m/amritan/lpl/angptl3/rfdiff/run1/outputs/_2.pdb", "r") as f:
    #    pdb = f.read()
    
    #print(Rg(pdb, chnid="A"))
    
