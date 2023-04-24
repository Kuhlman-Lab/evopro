# Tiger's modified scorefunctions for CD70 system
# changed de novo binder chain to D for CD70 trimer (chains A/B/C)

from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import os
import subprocess
import shutil
import math

path_to_3helixbundle = '/nas/longleaf/home/tigerz/kuhlmanlab/structures/5djt_3-helix-bundle/3helix_bundle_cleaned_better.pdb'
path_to_scaff42 = '/nas/longleaf/home/tigerz/kuhlmanlab/structures/cao_2022_de_novo_binders/scaffolds_amrita/scaffolds/scaffold_42.pdb'

def score_binder(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_clashPenalty(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashPenalty(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_clashPenalty_2(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashPenalty_2(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)        

def score_binder_clashPenaltyImproved(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashes_improved(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_clashPenalty_interface7(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_clashes_interface7(results, dsobj, contacts, penalty_weight = 1, orient=orient)    
    else:
        return score_binder_monomer(results, dsobj)  

def score_binder_clashPenalty_interface8(results, dsobj, contacts, distance_cutoffs=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        # interface 8, clash penalty weight = 3
        return score_binder_complex_clashes_if8(results, dsobj, contacts)    
    else:
        return score_binder_monomer(results, dsobj)               

def score_binder_complex(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)

    orientation_penalty = 0

    num_contacts = len(contacts)    
    score = -contactscore + orientation_penalty
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    #print("line 36: " + str(score))
    #print("line 37: ")
    #print(results)
    #score_list = []
    #score_list.append(score)
    return score, (score, len(contacts), contactscore, orientation_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashPenalty(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts # complex
    reslist2 = [x for x in residues.keys() if x.startswith("D")] # de novo binder
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)

    # calculate contacts with residues within 4 A of glycosylation site
    # but NOT within 4 A of N170-gly (70-72, 127-129, 133, 168-172) or N63-gly (61-64, 74-76, 94, 179)
    # penalty resids (renumbered): 8-11, 17-19, 21-23, 43, 74-76, 80, 115-119, 126 
    num_clashes = 0
    clash_resids = {8, 9, 10, 11, 17, 18, 19, 21, 22, 23, 43, 74, 75, 76, 80, 115, 116, 117, 118, 119, 126}
    print("line 58:")
    print(contacts)
    #print()
    #for contact in contacts: # contacts = list of tuples [(res1, res2), (res1, res5), (res3, res4)]
    # example: [('A107', 'D34'), ('C26', 'D33'), ('C26', 'D34'),... for chains A,B,C,D
    for contact in contacts:
        #print(contact[0])
        #print(contact[0][0:1])
        #print(contact[0][1:])
        if contact[0][0:1] != 'D' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact)) 
        if contact[1][0:1] != 'D' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1   
            print("clash found at: " + str(contact)) 
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))
    print()

    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    num_contacts = len(contacts)    
    score = -contactscore + orientation_penalty + clash_penalty 
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashPenalty_2(results, dsobj, contacts, orient=None):
    # for use with interface 6
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts # complex
    reslist2 = [x for x in residues.keys() if x.startswith("D")] # de novo binder
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)

    # calculate contacts with residues within 4 A of glycosylation site wih some manual edis 
    # penalty resids (renumbered): 10-11, 17-19, 21-23, 43, 76, 80, 115-118, 126 
    num_clashes = 0
    clash_resids = {10, 11, 17, 18, 19, 21, 22, 23, 43, 76, 80, 115, 116, 117, 118, 126}
    
    # check if res1 or res2 contains a clashing residue on CD70 (chains A/B/C) 
    # contacts = list of tuples [(res1, res2), (res1, res5), (res3, res4)...]
    for contact in contacts:
        if contact[0][0:1] != 'D' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact)) 
        if contact[1][0:1] != 'D' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1   
            print("clash found at: " + str(contact)) 
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))
    
    # orientation penalty  
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    # calculate overall score and PAE per contact
    num_contacts = len(contacts)    
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    score = -contactscore + orientation_penalty + clash_penalty

    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashPenalty_singlechain(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts # complex
    reslist2 = [x for x in residues.keys() if x.startswith("B")] # de novo binder
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)

    # penalty residues (next to sugars) : 8-11, 17-19, 21-23, 43, 74-76, 80, 115-119, 126 
    num_clashes = 0
    clash_resids = {8, 9, 10, 11, 17, 18, 19, 21, 22, 23, 43, 74, 75, 76, 80, 115, 116, 117, 118, 119, 126}
    
    # check if res1 or res2 contains a clashing residue on target protein (chain A) 
    # contacts = list of tuples [(res1, res2), (res1, res5), (res3, res4)...]
    for contact in contacts:
        if contact[0][0:1] != 'B' and int(contact[0][1:]) in clash_resids:
            num_clashes += 1
            print("clash found at: " + str(contact)) 
        if contact[1][0:1] != 'B' and int(contact[1][1:]) in clash_resids:
            num_clashes += 1   
            print("clash found at: " + str(contact)) 
    # calculate clash penalty
    clash_penalty = num_clashes*3
    print("num clashes: " +  str(num_clashes))
    
    # orientation penalty  
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    # calculate overall score and PAE per contact
    num_contacts = len(contacts)    
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    score = -contactscore + orientation_penalty + clash_penalty

    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashes_improved(results, dsobj, contacts, orient=None):
    # for if6 
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts # complex
    reslist2 = [x for x in residues.keys() if x.startswith("D")] # de novo binder
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)

    # calculate contacts with residues within 4 A of glycosylation site
    # but NOT within 4 A of N170-gly (70-72, 127-129, 133, 168-172) or N63-gly (61-64, 74-76, 94, 179)
    # penalty resids (renumbered): 8-11, 17-19, 21-23, 43, 74-76, 80, 115-119 

    num_clashes = 0
    clash_resids = {8, 9, 10, 11, 17, 18, 19, 21, 22, 23, 43, 74, 75, 76, 80, 115, 116, 117, 118, 119, 126}
    clash_resid_list = []
    for clash_res in clash_resids:
        clash_resid_list.append('A'+str(clash_res))
        clash_resid_list.append('B'+str(clash_res))
        clash_resid_list.append('C'+str(clash_res))

    #print(clash_resid_list)        
    clashes, num_clashes = score_contacts(pdb, clash_resid_list, reslist2, dsobj=dsobj, first_only=False)
    clash_penalty = num_clashes*3

    print("line 211: clash_penalty = 3*num_clashes ")
    print("num clashes: " +  str(num_clashes))
    #print(contacts)
    print(clashes)
    print()

    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0

    num_contacts = len(contacts)    
    score = -contactscore + orientation_penalty + clash_penalty 
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashes_interface7(results, dsobj, contacts, penalty_weight = 1, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts # complex
    reslist2 = [x for x in residues.keys() if x.startswith("D")] # de novo binder
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)

    # penalty residues (ori numbering): 63-64, 70-72, 74-76, 94, 129, 133, 168-171
    # penalty resids (renumbered):  10-11, 17-19, 21-23, 43, 76, 80, 115-118, 126 

    num_clashes = 0
    clash_resids = {10, 11, 17, 18, 19, 21, 22, 23, 43, 76, 80, 115, 116, 117, 118, 126}
    clash_resid_list = []
    for clash_res in clash_resids:
        clash_resid_list.append('A'+str(clash_res))
        clash_resid_list.append('B'+str(clash_res))
        clash_resid_list.append('C'+str(clash_res))

    clashes, num_clashes = score_contacts(pdb, clash_resid_list, reslist2, dsobj=dsobj, first_only=False)
    clash_penalty = num_clashes*penalty_weight

    print("line 211: clash_penalty = " +str(penalty_weight)+ "*num_clashes ")
    print("num clashes: " +  str(num_clashes))
    #print(contacts)
    #print(clashes)
    print()

    orientation_penalty = 0

    num_contacts = len(contacts)    
    score = -contactscore + orientation_penalty + clash_penalty 
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashes_interface_general(results, dsobj, contacts, penalty_weight, clashResids):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    #reslist1 = contacts # complex
    reslist1 = contacts[0] # first of three args in list
    reslist2 = [x for x in residues.keys() if x.startswith("D")] # de novo binder
    print("line 314")
    print(reslist1)
    print(reslist2)
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
    #print("Binder contacts:")
    #print(contacts)

    # penalty residues (ori numbering): 63-64, 70-72, 74-76, 94, 129, 133, 168-171
    # penalty resids (renumbered):  10-11, 17-19, 21-23, 43, 76, 80, 115-118, 126 

    num_clashes = 0
    clash_resid_list = []
    for clash_res in clashResids:
        clash_resid_list.append('A'+str(clash_res))
        clash_resid_list.append('B'+str(clash_res))
        clash_resid_list.append('C'+str(clash_res))

    clashes, num_clashes = score_contacts(pdb, clash_resid_list, reslist2, dsobj=dsobj, first_only=False)
    clash_penalty = num_clashes*penalty_weight

    print("line 211: clash_penalty = " +str(penalty_weight)+ "*num_clashes ")
    print("num clashes: " +  str(num_clashes))
    #print(contacts)
    #print(clashes)
    print()

    orientation_penalty = 0

    num_contacts = len(contacts)    
    score = -contactscore + orientation_penalty + clash_penalty 
    PAE_per_contact = 0
    if num_contacts > 0:
        PAE_per_contact = (70.0-70.0*contactscore/num_contacts)/2
    return score, (score, len(contacts), contactscore, orientation_penalty, clash_penalty, PAE_per_contact), contacts, pdb, results

def score_binder_complex_clashes_if8(results, dsobj, contacts):  
    #inerface 8
    # penalty_weight = 3
    clashResids = {10, 11, 17, 18, 19, 21, 22, 23, 43, 76, 80, 115, 116, 117, 118}

    return score_binder_complex_clashes_interface_general(results, dsobj, contacts, 3, clashResids)


def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    return score, (score, confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2, binder_chain="D", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    #print(reslist1)
    #print(reslist2)
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    print("RMSD: " + str(rmsd_binder))
    return rmsd_binder*5

def score_binder_rmsd_weight2(pdb1, pdb2, binder_chain="D", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    #print(reslist1)
    #print(reslist2)
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    print("RMSD: " + str(rmsd_binder) + ", Score=RMSDx2: " + str(rmsd_binder*2))
    return rmsd_binder*2

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

def score_binder_rmsd_restraint(pdb1, pdb2, binder_chain="D",  dsobj=None):
    # for CD70 trimer system
    # to keep the de novo binder close in structure to the original binder
    # calculates    1) Ca-only RMSD of de novo binder unbound vs bound to CD70 complex
    #               2) Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential
    # total RMSD score = 1) + 2)
    spring_constant = 10.0 
    rmsd_cutoff = 5.0

    # calculate 1)
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1) # complex
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2) # de novo binder alone
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    #print("RMSD: " + str(rmsd_binder))

    # calculate 2)
    with open(path_to_3helixbundle, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("RMSD scores, RMSD_bound: %.3f, rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_binder, rmsd_to_starting, rmsd_potential) )

    return rmsd_binder*5 + rmsd_potential

def score_binder_rmsd_restr(pdb1, pdb2, spring_constant, rmsd_cutoff, binder_chain="D", dsobj=None):
    # for CD70 trimer system, to keep the de novo binder close in structure to the original binder
    # calculates    1) Ca-only RMSD of de novo binder unbound vs bound to CD70 complex
    #               2) Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential
    # total RMSD score = 1) + 2)

    # calculate 1)
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1) # complex
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2) # de novo binder alone
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True, dsobj=dsobj)

    # calculate 2)
    with open(path_to_3helixbundle, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:     # apply flat-bottom quadratic-shaped potential function
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("RMSD scores, RMSD_bound: %.3f, rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_binder, rmsd_to_starting, rmsd_potential) )

    return rmsd_binder*5 + rmsd_potential

def score_binder_rmsd_to_starting(pdb, path_to_starting, spring_constant, rmsd_cutoff, binder_chain="D", dsobj=None):
    # for CD70 trimer system, to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential

    # import de novo binder alone??
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb) # de novo binder alone
    print(chains1)
    reslist1 = [x for x in residues1.keys()]

    # import starting pdb structure
    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist1, pdb, ca_only=True, dsobj=dsobj)
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:     # apply flat-bottom quadratic-shaped potential function
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_to_starting, rmsd_potential) )

    return rmsd_potential

def score_binder_rmsd_start_5_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 5, 3)

def score_binder_rmsd_start_5_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 5, 5)

def score_binder_rmsd_start_5_7(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 7
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 5, 7)    

def score_binder_rmsd_start_10_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 10
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 10, 3) 

def score_binder_rmsd_start_10_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 10
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 10, 5)  

def score_binder_rmsd_start_20_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 20
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 20, 3) 

def score_binder_rmsd_start_20_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 20
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_starting, 20, 5)     

if __name__=="__main__":
    print("no main functionality")
