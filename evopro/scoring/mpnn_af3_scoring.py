from omegaconf import OmegaConf

from standard_score_funcs import calculate_residue_plddt, score_contacts_pae_weighted, get_avg_plddt, score_rmsd
from utils.parsing_utils import get_coordinates_pdb

def score_seq_diff_backbone_plddt_only_af3(result, reslist=None):
    
    #get average plddt of the pdb
    df, avg_plddt = calculate_residue_plddt(result, reslist=reslist)
    
    score = -avg_plddt
    return (score, avg_plddt, 0)

def score_seq_diff_backbone_af3(result, backbone, reslist=None):
    # get rmsd between predicted and true backbone
    pdb = result['pdb']
    rmsd = score_rmsd(pdb, backbone, reslist=reslist)
    
    #get average plddt of the pdb
    df, avg_plddt = calculate_residue_plddt(result, reslist=reslist)
    score = -avg_plddt/10.0 + rmsd
    
    return (score, avg_plddt, rmsd)

def score_seq_diff_backbone_pae_af3(result, backbone, reslist=None):
    conf = OmegaConf.create({"scoring": {"contacts": {"contact_distance": 6, "max_contacts": 50, "score_type": "ca"}}})

    #get average plddt of the pdb
    df, avg_plddt = calculate_residue_plddt(result, reslist=reslist)
    
    #get pae score
    chains, residues, resindices = get_coordinates_pdb(result['pdb'])
    reslist1 = [x for x in residues if x.startswith('A')]
    reslist2 = [x for x in residues if x not in reslist1]
    pairs, pae_contact_score = score_contacts_pae_weighted(result, None, reslist1, reslist2, conf)
    
    score = -avg_plddt/10.0 - pae_contact_score
    
    return (score, avg_plddt, pae_contact_score)

def score_seq_diff_backbone_af2(result, backbone, reslist=None):
    # get rmsd between predicted and true backbone
    pdb = result['pdb']
    rmsd = score_rmsd(pdb, backbone, reslist=reslist)
    
    #get average plddt of the pdb
    avg_plddt = get_avg_plddt(result, reslist=reslist)
    score = -avg_plddt/10.0 + rmsd
    
    return (score, avg_plddt, rmsd)