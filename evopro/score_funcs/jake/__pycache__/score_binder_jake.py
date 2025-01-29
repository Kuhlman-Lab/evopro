from evopro.utils.pdb_parser import get_coordinates_pdb, get_coordinates_pdb_old, get_ca_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd
import math
import numpy as np

def score_overall(results, dsobj, contacts=None, distance_cutoffs=None, starting_pdb=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts, distance_cutoffs))
            #for calculating # of neighbors for a given cys residue
            score.append(score_resi_dist(result, dsobj, contacts, distance_cutoffs))
        else:
            score.append(score_binder_monomer(result, dsobj))
            if starting_pdb:
                score.append(score_binder_rmsd_to_starting(pdb, starting_pdb, dsobj=dsobj))
    
    

    #print(score)
    score.append(score_binder_rmsd(pdbs[0], pdbs[1], dsobj=dsobj))
    print([x[0] for x in score])
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_2targ(results, dsobj, contacts=None, distance_cutoffs=None, starting_pdb=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_2targ(result, dsobj, contacts, distance_cutoffs))
            #for calculating # of neighbors for a given cys residue
            score.append(score_resi_dist(result, dsobj, contacts, distance_cutoffs))
        else:
            score.append(score_binder_monomer(result, dsobj))
            if starting_pdb:
                score.append(score_binder_rmsd_to_starting(pdb, starting_pdb, dsobj=dsobj))
    
    

    #print(score)
    score.append(score_binder_rmsd(pdbs[0], pdbs[1], dsobj=dsobj))
    print([x[0] for x in score])
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_2targ(results, dsobj, contacts, distance_cutoffs):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,15,8)
    reslist1a = [x for x in contacts[0] if x.startswith("A")]
    reslist1b = [x for x in contacts[0] if x.startswith("B")]
    reslist2a = [x for x in residues.keys() if x.startswith("A")]
    reslist2b = [x for x in residues.keys() if x.startswith("B")]
    #print(f'reslist 1a: {reslist1a}')
    #print(f'reslist 1b: {reslist1b}')
    #print(reslist1, reslist2)
    print(f'distance cutoffs: {distance_cutoffs}')
    contact_list_a, contactscore_a = score_contacts_pae_weighted(results, pdb, reslist1a, reslist2b, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0], contact_cap=100)
    contact_list_b, contactscore_b = score_contacts_pae_weighted(results, pdb, reslist1b, reslist2a, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0], contact_cap=100)
    num_contacts_a = len(contact_list_a)
    num_contacts_b = len(contact_list_b)
    contacts_tot = num_contacts_a + num_contacts_b
    contactscore = contactscore_a + contactscore_b
    print(f'contacts a: {num_contacts_a}, contacts b: {num_contacts_b}, score a: {contactscore_a}, score b: {contactscore_b}')
    
    within_contacts = 0
    within_contactscore = 0
    binder_region_resids = contacts[1]
    core_resids = contacts[2]
    if binder_region_resids:
        within_contacts, within_contactscore = score_contacts_pae_weighted(results, pdb, binder_region_resids, core_resids, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
        '''
        for contact in bonus_contacts:
            if contact[0][0:1] == 'B' and int(contact[0][1:]) in bonus_resids:
                within_contacts += 1
                print("bonus found at: " + str(contact[0]))
            if contact[1][0:1] == 'B' and int(contact[1][1:]) in bonus_resids:
                within_contacts += 1
                print("bonus found at: " + str(contact[1]))
        '''
    contactscore += within_contactscore
    print(f'Within Binder Contact: {within_contacts}, Within Score: {within_contactscore}')
    '''
    penalties = 0
    penalty_resids = contacts[2]
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
        for contact in penalty_contacts:
            if contact[0][0:1] == 'B' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'B' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    '''
    #tot_contacts = len(contact_list_a) + len(contact_list_b)
    pae_per_contact_a = 0
    if num_contacts_a > 0:
        pae_per_contact_a = (70.0-(70.0*contactscore_a)/num_contacts_a)/2
        #print("a is working")
    pae_per_contact_b = 0
    if num_contacts_b > 0:
        pae_per_contact_b = (70.0-(70.0*contactscore_b)/num_contacts_b)/2
        #print("b is working")
    
    print(f"pae_a: {pae_per_contact_a}")
    print(f"pae_b: {pae_per_contact_b}")
    pae_tot = (pae_per_contact_a + pae_per_contact_b) / 2
    score = -contactscore #+ bonus #+penalty
    print(score, (score, contacts_tot, contactscore, pae_tot))
    return score, (score, contacts_tot, contactscore, pae_tot), contacts, pdb, results

def score_binder_complex(results, dsobj, contacts, distance_cutoffs):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    print(reslist1, reslist2)
    contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0], contact_cap=100)
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
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
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
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
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results

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
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj, ca_only=True)

    return (rmsd_binder*5, rmsd_binder)

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x.startswith("A")]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return (rmsd_potential*5, rmsd_to_starting, rmsd_potential)

def score_resi_dist(results, dsobj, contacts, distance_cutoffs):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    old_chains, old_residues, old_resindices = get_coordinates_pdb_old(pdb)
    chains, residues, resindices = get_coordinates_pdb(pdb)

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    
    #print(f"old_residue_keys: {old_residues.keys()}")
    #print(f"residue_keys: {residues.keys()}")
    for res in old_residues.keys():
        #print(f'residue: {res}')
        r = res.split('_')
        #print(f'split list = {r}')
        if r[0] =='B' and r[1] == "CYS" and int(r[2]) >=50:
            reslist1 = [str(r[0])+str(r[2])]
            #print(reslist1)
            print('checking cys neighbors')
            break
    reslist2 = [x for x in residues.keys() if x.startswith("A")]
    print(reslist1, reslist2)
    contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=15)
    if len(contact_list) > 0:
        score = -len(contact_list)
    else:
        score = 10
    print(f'cys score: {score}') 
    print(score, (score, len(contact_list), contactscore))
    return score, (score, len(contact_list), contactscore), contacts, pdb, results

def calculate_plane_from_points(p1, p2, p3, p4, p5, p6):
    # Convert points to numpy arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    
    p4 = np.array(p4)
    p5 = np.array(p5)
    p6 = np.array(p6)

    # Calculate the vectors from p1 to p2 and p1 to p3
    v1 = p2 - p1
    v2 = p3 - p1

    v3 = p5 - p4
    v4 = p6 - p4
    
    # Calculate the normal vector as the cross product of v1 and v2
    normal_vector = np.cross(v1, v2)
    
    normal_vector2 = np.cross(v3, v4)

    # Calculate the d constant in the plane equation
    d1 = np.dot(normal_vector, p1)

    d2 = np.dot(normal_vector2, p4)

    #calculate the angle difference between the two planes 
    angle_diff = np.dot(normal_vector, normal_vector2)

    #calculate the distance between the two planes
    d_diff = np.dot(d1, d2)

    return angle_diff, d_diff

def score_pip3(results, dsobj, contacts=None, distance_cutoffs=None, starting_pdb=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    binder = '1h0a'

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        score.append(score_binder_monomer(result, dsobj))
        if starting_pdb:
            score.append(score_binder_rmsd_to_starting(pdb, starting_pdb, dsobj=dsobj))
    
    #select proper residues at the interface
    coords = get_ca_coordinates_pdb(pdb)
    #1h0a has 158 residues, 1mai has 130 residues uncut
    if binder == '1h0a':
        res1 = coords[2]
        res2 = coords[13]
        res3 = coords[68]
        res4 = coords[285]
        res5 = coords[292]
        res6 = coords[306]
    else:
        res1 = coords[32]
        res2 = coords[57]
        res3 = coords[103]
        res4 = coords[257]
        res5 = coords[264]
        res6 = coords[278]


    #determine if the binding interfaces are on the same plane (d1-d2 ~0)
    #determine if the binding interfaces are parallel (normal1 x normal2 ~0)
    angle, plane_dist = calculate_plane_from_points(res1, res2, res3, res4, res5, res6)
    print(f'angle: {angle}, distance: {plane_dist}')

    #print(f'angle: {angle}, distance: {plane_dist}')

    angle_constant = 10.0 
    angle_cutoff = 0.1745
    distance_constant = 10.0
    distance_cutoff = 5.0
    angle_potential = 0
    #exponential potential function for angle and distance between the two planes
    if abs(angle) > angle_cutoff:
        angle_potential = angle_constant*math.pow(angle - angle_cutoff, 2)
    score.append((angle_potential*5, angle_potential))
    
    distance_potential = 0
    if plane_dist > distance_cutoff:
        distance_potential = distance_constant*math.pow(plane_dist - distance_constant, 2)
    score.append((distance_potential*5, distance_potential))

    #print(score)
    #score.append(score_binder_rmsd(pdbs[0], pdbs[1], dsobj=dsobj))
    print([x[0] for x in score])
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results




if __name__=="__main__":
    print("no main functionality")
