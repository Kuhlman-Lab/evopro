import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, score_pae_interaction
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient

def score_interchain_pae(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    
    target_length = len(reslist2)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    pae = results['pae_output'][0]
    avg_pae = score_pae_interaction(pae, pdb, target_length=target_length)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore
    return score, (avg_pae, pae_per_contact, confscore, len(contact_list), contactscore), pdb, results

def score_interchain_pae_ang3(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B") or x.startswith("C") or x.startswith("D")]
    
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    
    target_length = len(reslist2)
    pae = results['pae_output'][0]
    avg_pae = score_pae_interaction(pae, pdb, target_length=target_length)
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore
    return score, (avg_pae, pae_per_contact, confscore, len(contact_list), contactscore), pdb, results

def score_interchain_pae_ang3_old(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = [x for x in residues.keys() if x.startswith("A") or x.startswith("B") or x.startswith("C")]
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    
    target_length = len(reslist1)
    pae = results['pae_output'][0]
    avg_pae = score_pae_interaction(pae, pdb, target_length=target_length)
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore
    return score, (avg_pae, pae_per_contact, confscore, len(contact_list), contactscore), pdb, results

def score_interchain_pae_monomer(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore
    return score, (pae_per_contact, confscore, len(contact_list), contactscore), pdb, results

