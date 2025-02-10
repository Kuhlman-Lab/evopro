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


def score_binder_clashpenalty(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashpenalty(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_clashpenalty_bonusplacement(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashpenalty_bonusplacement(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)


def score_binder_clashpenalty_contactonly_pae(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashpenalty_contactonly_pae(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_clashpenalty_bonusplacement_contactonly_pae(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    print(results)
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashpenalty_bonusplacement_contactonly_pae(results, dsobj, contacts, orient=orient)
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
    score = -contactscore + orientation_penalty
    print(score, (score, len(contacts), contactscore, orientation_penalty))
    return score, (score, len(contacts), contactscore, orientation_penalty), contacts, pdb, results


def score_binder_complex_clashpenalty(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]

    reslist3 = [x for x in residues.keys() if x.startswith("A")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)

    full_contact, contactscore_full = score_contacts_pae_weighted(results, pdb, reslist3, reslist2, dsobj=dsobj, first_only=False)
    num_clashes = 0
    #residue selected on epitope
    clash_resids = {15,54,154,201}
    clash_reslist = ["A"+str(x) for x in clash_resids]
    clash_contacts, clash_contactscore = score_contacts_pae_weighted(results, pdb, clash_reslist, reslist2, dsobj=dsobj, first_only=False, dist=8)
    for contact in clash_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[0]))
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[1]))
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))

    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore_full)/num_contacts)/2

    score = -contactscore + pae_per_contact + orientation_penalty + clash_penalty

    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, pae_per_contact), contacts, pdb, results


def score_binder_complex_clashpenalty_contactonly_pae(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]

    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)

    num_clashes = 0
    #residue selected on epitope
    clash_resids = {15,54,154,201}
    clash_reslist = ["A"+str(x) for x in clash_resids]
    clash_contacts, clash_contactscore = score_contacts_pae_weighted(results, pdb, clash_reslist, reslist2, dsobj=dsobj, first_only=False, dist=8)
    for contact in clash_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[0]))
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[1]))
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))

    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2

    score = -contactscore + orientation_penalty + clash_penalty

    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, pae_per_contact), contacts, pdb, results

def score_binder_complex_clashpenalty_bonusplacement(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    reslist3 = [x for x in residues.keys() if x.startswith("A")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    full_contact, contactscore_full = score_contacts_pae_weighted(results, pdb, reslist3, reslist2, dsobj=dsobj, first_only=False)
    num_clashes = 0
    #new updated on 04012023
    clash_resids = {15,54,154,201}
    clash_reslist = ["A"+str(x) for x in clash_resids]
    clash_contacts, clash_contactscore = score_contacts_pae_weighted(results, pdb, clash_reslist, reslist2, dsobj=dsobj, first_only=False, dist=8)
    for contact in clash_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[0]))
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[1]))
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))

    num_bonus = 0
    bonus_resids = {174,176,182}
    bonus_reslist = ["A"+str(x) for x in clash_resids]
    bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_reslist, reslist2, dsobj=dsobj, first_only=False, dist=4)
    for contact in bonus_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in bonus_resids:
            num_bonus += 1
            
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in bonus_resids:
            num_bonus += 1
            
    # calculate clash penalty
    contact_bonus = num_bonus * -2

    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2

    score = -contactscore + pae_per_contact+orientation_penalty + clash_penalty + contact_bonus

    return score, (score, len(contacts), contactscore, contact_bonus, clash_penalty, pae_per_contact), contacts, pdb, results


def score_binder_complex_clashpenalty_bonusplacement_contactonly_pae(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    num_clashes = 0
    #new updated on 04012023
    clash_resids = {15,54,154,201}
    clash_reslist = ["A"+str(x) for x in clash_resids]
    clash_contacts, clash_contactscore = score_contacts_pae_weighted(results, pdb, clash_reslist, reslist2, dsobj=dsobj, first_only=False, dist=8)
    for contact in clash_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[0]))
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[1]))
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))

    num_bonus = 0
    bonus_resids = {174,176,182}
    bonus_reslist = ["A"+str(x) for x in clash_resids]
    bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_reslist, reslist2, dsobj=dsobj, first_only=False, dist=4)
    for contact in bonus_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in bonus_resids:
            num_bonus += 1
            
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in bonus_resids:
            num_bonus += 1
            
    # calculate clash penalty
    contact_bonus = num_bonus * -2

    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2

    score = -contactscore + orientation_penalty + clash_penalty + contact_bonus

    return score, (score, len(contacts), contactscore, contact_bonus, clash_penalty, pae_per_contact), contacts, pdb, results
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    reslist3 = [x for x in residues.keys() if x.startswith("A")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    full_contact, contactscore_full = score_contacts_pae_weighted(results, pdb, reslist3, reslist2, dsobj=dsobj, first_only=False)
    num_clashes = 0
    #new updated on 04012023
    clash_resids = {15,54,154,201}
    clash_reslist = ["A"+str(x) for x in clash_resids]
    clash_contacts, clash_contactscore = score_contacts_pae_weighted(results, pdb, clash_reslist, reslist2, dsobj=dsobj, first_only=False, dist=8)
    for contact in clash_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[0]))
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact[1]))
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))

    num_bonus = 0
    bonus_resids = {174,176,182}
    bonus_reslist = ["A"+str(x) for x in clash_resids]
    bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_reslist, reslist2, dsobj=dsobj, first_only=False, dist=4)
    for contact in bonus_contacts:
         if contact[0][0:1] == 'A' and int(contact[0][1:]) in bonus_resids:
            num_bonus += 1
            
         if contact[1][0:1] == 'A' and int(contact[1][1:]) in bonus_resids:
            num_bonus += 1
            
    # calculate clash penalty
    contact_bonus = num_bonus * -2

    num_contacts = len(contacts)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2

    score = -contactscore + pae_per_contact+orientation_penalty + clash_penalty + contact_bonus

    return score, (score, len(contacts), contactscore, contact_bonus, clash_penalty, pae_per_contact), contacts, pdb, results

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    print(score)
    return score, (score, confscore2), pdb, results


def score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)

    return rmsd_binder*5


def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

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
    print("no main functionality")
