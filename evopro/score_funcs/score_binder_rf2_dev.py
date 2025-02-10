from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd
import math

def score_overall(results, dsobj, contacts=None, distance_cutoffs=None, starting_pdb=None):
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = result['pdb']
        pdbs.append(pdb)
        chains, _, _ = get_coordinates_pdb(pdb)
        print(chains, len(chains))
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts, distance_cutoffs, contact_cap=120))
        else:
            score.append(score_binder_monomer(result, dsobj))
            if starting_pdb:
                score.append(score_binder_rmsd_to_starting(pdb, starting_pdb, dsobj=dsobj))
    
    print("SCORE", score)
    score.append(score_binder_rmsd(pdbs[0], pdbs[2], binder_chain="D", dsobj=dsobj, ca_only=True))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_complex(results, dsobj, contacts, distance_cutoffs, contact_cap=36, binder_chain="D"):
    pdb = results['pdb']
    _, residues, _ = get_coordinates_pdb(pdb)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith(binder_chain)]
    print(reslist1, reslist2)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2, contact_cap=contact_cap, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, _ = score_contacts_pae_weighted_efficient(results, pdb, bonus_resids, reslist2, contact_cap=contact_cap, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
        for contact in bonus_contacts:
            if contact[0][0:1] == 'B' and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if contact[1][0:1] == 'B' and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    if penalty_resids:
        penalty_contacts, _ = score_contacts_pae_weighted_efficient(results, pdb, penalty_resids, reslist2, contact_cap=contact_cap, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
        for contact in penalty_contacts:
            if contact[0][0:1] == 'B' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'B' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
    
    score = -contactscore + penalty + bonus
    print(score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty))
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty)

def score_binder_monomer(results, dsobj):
    pdb = results['pdb']
    _, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    print(score)
    return score, (score, confscore2)

def score_binder_rmsd(pdb1, pdb2, binder_chain="A", dsobj=None, ca_only=True):
    _, residues1, _ = get_coordinates_pdb(pdb1)
    _, residues2, _ = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj, ca_only=ca_only)

    return (rmsd_binder*5, rmsd_binder)

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x.startswith("A")]

    _, residues, _ = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return (rmsd_potential*5, rmsd_to_starting, rmsd_potential)

if __name__=="__main__":
    print("no main functionality")
