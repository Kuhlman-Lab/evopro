from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient
import math

def score_overall(results, dsobj, contacts=None, distance_cutoffs=None,):
    from alphafold.common import protein
    #The results here are a list of results dictionaries
    #one for each AF2 prediction generated for the same design sequence
    #so this number should match the number of items in af2_preds
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    #parsing through each dictionary and scoring them based on which prediction it is
    for result in results:

        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, _, _ = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts, distance_cutoffs))
            score.append(score_complex_confidence(result, dsobj))
        else:
            score.append(score_binder_monomer(result, dsobj))
    
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    score.append(score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    score.append(threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_overall_2(results, dsobj, contacts=None, distance_cutoffs=None,):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:

        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, _, _ = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts, distance_cutoffs))
            score.append(score_complex_confidence(result, dsobj))
        else:
            score.append(score_binder_monomer(result, dsobj))
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[2]['unrelaxed_protein'])
    score.append(score_potential_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_complex(results, dsobj, contacts, distance_cutoffs, binder_chain="B"):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    _, residues, _ = get_coordinates_pdb(pdb)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith(binder_chain)]
    print(reslist1, reslist2)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, _ = score_contacts_pae_weighted_efficient(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
        for contact in bonus_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    if penalty_resids:
        penalty_contacts, _ = score_contacts_pae_weighted_efficient(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
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

def score_complex_confidence(results, dsobj):
    from alphafold.common import protein
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    _, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of complex is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = -confscore2*2 + plddt_potential
    print(score)
    return score, (score, confscore2)

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    _, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of monomer is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = -confscore2*2 + plddt_potential
    print(score)
    return score, (score, confscore2)

def score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    if len(reslist1) == len(reslist2):
        rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    else:
        print("Trying to calculate RMSD between proteins of different length. Setting RMSD to 0.")
        rmsd_binder = 0

    return [-rmsd_binder*30]

def score_potential_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    rmsd_potential = 0
    spring_constant = 10.0
    rmsd_cutoff_max = 10.0
    rmsd_cutoff_min = 5.0
    rmsd_weight = 25
    if rmsd < rmsd_cutoff_min:
        rmsd_potential = rmsd_weight*rmsd
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    elif rmsd <rmsd_cutoff_max:
        rmsd_potential = rmsd_weight*rmsd_cutoff_min
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    else:
        rmsd_potential = rmsd_weight*rmsd_cutoff_min - spring_constant*math.pow(rmsd_cutoff_max - rmsd,2)
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    return [-rmsd_potential]

def threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 10

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]


    rmsd_to_starting = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of 5.0. Calculating RMSD potnetial")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return [rmsd_potential*5]

