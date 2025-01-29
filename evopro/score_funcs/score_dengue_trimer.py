import math
import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd
    

def score_trimer(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb, ca_only=True, dsobj=dsobj)
    else:
        rmsdscore=0
        
    if contacts[0]:
        contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[0], dist=4, contact_cap=36)
    else:
        contact_list = []
        contactscore = 0
    
    if contacts[1]:
        bonus_list, bonusscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[1], dist=4, contact_cap=36)
    else:
        bonus_list = []
        bonusscore = 0
    
    if contacts[2]:
        penalty_list, penaltyscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[1], dist=4, contact_cap=36)
    else:
        penalty_list = []
        penaltyscore = 0
    
    overall_contactscore = -contactscore - len(bonus_list)*5 + len(penalty_list)*10
    
    overall_score = contactscore - confscore/10 + rmsdscore
    score = [(overall_contactscore, -contactscore, len(bonus_list)*5, len(penalty_list)*10),(- confscore/10, confscore),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results

def score_trimer2(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb, ca_only=True, dsobj=dsobj)
    else:
        rmsdscore=0

    if contacts[0]:
        contact_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[0], dist=4, contact_cap=36)
    else:
        contact_list = []
        contactscore = 0

    if contacts[1]:
        bonus_list, bonusscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[1], dist=4, contact_cap=36)
    else:
        bonus_list = []
        bonusscore = 0

    if contacts[2]:
        penalty_list, penaltyscore = score_contacts_pae_weighted(results, pdb, reslist1, contacts[1], dist=4, contact_cap=36)
    else:
        penalty_list = []
        penaltyscore = 0

    overall_contactscore = -contactscore - len(bonus_list)*5 + len(penalty_list)*10

    overall_score = contactscore - confscore/2 + rmsdscore
    score = [(overall_contactscore, -contactscore, len(bonus_list)*5, len(penalty_list)*10),(- confscore/2, confscore),  (rmsdscore,)]

    return overall_score, score, [pdb], results

def score_rmsd_to_starting(pdb, path_to_starting, ca_only=True, dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=ca_only, dsobj=dsobj)

    return rmsd_to_starting
