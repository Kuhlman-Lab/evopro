import torch
import torch.nn as nn
from biopandas.pdb import PandasPdb
import pandas as pd
from typing import List
import numpy as np
from Bio.PDB import PDBParser
import pandas as pd
from io import StringIO
from omegaconf import OmegaConf

from utils.geometry import *
from utils.parsing_utils import *
from utils.multichain_alignment import multichain_permutation_alignment

def score_seq_diff_backbone_plddt_only(result, backbone, reslist=None):
    
    #get average plddt of the pdb
    df, avg_plddt = calculate_residue_plddt(result, reslist=reslist)
    
    score = -avg_plddt
    return (score, avg_plddt, 0)

def score_seq_diff_backbone(result, backbone, reslist=None):
    # get rmsd between predicted and true backbone
    pdb = result['pdb']
    with open("temp1.pdb", "w") as f:
        f.write(pdb)
    with open("temp2.pdb", "w") as f:
        f.write(backbone)
    rmsd = score_rmsd(pdb, backbone, reslist=reslist)
    
    #get average plddt of the pdb
    df, avg_plddt = calculate_residue_plddt(result, reslist=reslist)
    score = -avg_plddt/10.0 + rmsd
    
    return (score, avg_plddt, rmsd)

def score_seq_diff_backbone_pae(result, backbone, reslist=None):
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

def get_avg_plddt(result, reslist=None):
    plddt = result['plddt']
    sum_plddt = 0
    if reslist:
        pdb = result['pdb']
        _, _, resindices = get_coordinates_pdb(pdb)
        for res in reslist:
            resid = resindices[res]
            sum_plddt = sum_plddt + plddt[resid]
        avg_plddt = sum_plddt/len(reslist)
    else:
        avg_plddt = sum(plddt)/len(plddt)
    
    return avg_plddt

def calculate_residue_plddt(result, reslist=None):
    """
    Calculate average pLDDT for each residue, with optional filtering.
    
    Parameters:
    -----------
    pdb_file : str
        Path to the input PDB file
    plddt_list : list or numpy array
        List of pLDDT values corresponding to atoms in the PDB file
    target_residues : list, optional
        List of residues to filter, in format ['A1', 'B2', ...]
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with residue information and average pLDDT
    """
    pdb_str = result['pdb']
    plddt_list = result['atom_plddts']

    # Initialize PDB parser
    parser = PDBParser(QUIET=True)
    # with open("outputs/temp.pdb", "w") as f:
    #     f.write(pdb_str)
    pdb_fh = StringIO(pdb_str)
    structure = parser.get_structure('protein', pdb_fh)
    
    # Prepare data structures
    residue_plddt = {}
    plddt_index = 0
    
    # Iterate through the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # Create a residue identifier string (e.g., 'A1')
                residue_key = f"{chain.id}{residue.id[1]}"
                
                # Check if this residue is in the target list (if specified)
                if reslist is None or residue_key in reslist:
                    # Collect pLDDT values for atoms in this residue
                    residue_atoms_plddt = []
                    
                    for atom in residue:
                        # Check if we have pLDDT for this atom
                        if plddt_index < len(plddt_list):
                            residue_atoms_plddt.append(plddt_list[plddt_index])
                            plddt_index += 1
                    
                    # Calculate average pLDDT for the residue
                    if residue_atoms_plddt:
                        avg_plddt = np.mean(residue_atoms_plddt)
                        
                        # Determine residue type
                        if residue.id[0] == ' ':
                            residue_type = 'Protein'
                        elif residue.id[0][0] == 'H':
                            residue_type = 'Hetero'
                        else:
                            residue_type = 'Other'
                        
                        # Store residue information
                        residue_identifier = (chain.id, residue.id[1], residue.id[2])
                        residue_plddt[residue_identifier] = {
                            'chain': chain.id,
                            'residue_number': residue.id[1],
                            'residue_name': residue.resname,
                            'residue_type': residue_type,
                            'residue_key': residue_key,
                            'avg_plddt': avg_plddt,
                            'atom_count': len(residue_atoms_plddt)
                        }
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(residue_plddt, orient='index')
        
    # Calculate overall average pLDDT for target residues
    if not df.empty:
        overall_avg_plddt = df['avg_plddt'].mean()
        # print(f"\nOverall Average pLDDT for Target Residues: {overall_avg_plddt:.2f}")
    else:
        overall_avg_plddt = 0
        print("WARNING: No pLDDT values found for target residues.\nSetting average pLDDT to 0.")
    
    return df, overall_avg_plddt


def score_contacts_pae_weighted(result, seq, reslist1, reslist2, conf):
    dist = conf.scoring.contacts.contact_distance
    contact_cap = conf.scoring.contacts.max_contacts
    pdb = result['pdb']

    chains, residues, resindices = get_coordinates_pdb_tokens(pdb)
    pae = result['pae']
    
    # print(resindices)
    
    #do this if af2
    prot = ProteinContacts(pdb, top_k=contact_cap)
    _, E_idx, _ = prot._dist()
    neighbors = E_idx.numpy()[0]
    
    reslist1_inds = [resindices[res] for res in reslist1]
    reslist2_inds = [resindices[res] for res in reslist2]
    resindices_rev = dict((v, k) for k, v in resindices.items())

    score = 0
    pairs = []
    for res1_ind in reslist1_inds:
        neighbors_res1 = neighbors[res1_ind]
        #print(neighbors_res1, reslist2_inds)
        hits = list(Intersection(neighbors_res1, reslist2_inds))
        #print(hits)
        for res2_ind in hits:
            contact = 0
            weight = 0
            res1 = resindices_rev[res1_ind]
            res2 = resindices_rev[res2_ind]
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        pair = (res1, res2)
                        pair_rev = (res2, res1)
                        if pair not in pairs and pair_rev not in pairs:
                            if len(pairs)<contact_cap:
                                contact=1
                                pae_contact = pae[res1_ind][res2_ind] + pae[res2_ind][res1_ind]
                                weight = (70-pae_contact)/70
                                pairs.append(pair)

            score = score + contact*weight
    # print(pairs)
    return pairs, score

def gather_edges(edges, neighbor_idx):
    # Features [B,N,N,C] at Neighbor indices [B,N,K] => Neighbor features [B,N,K,C]
    neighbors = neighbor_idx.unsqueeze(-1).expand(-1, -1, -1, edges.size(-1))
    edge_features = torch.gather(edges, 2, neighbors)
    return edge_features

def Intersection(lst1, lst2):
    return set(lst1).intersection(lst2)

class ProteinContacts(nn.Module):
    def __init__(self, pdb_str, top_k=30):
        
        """ Extract protein features """
        
        super(ProteinContacts, self).__init__()

        self.top_k = top_k
        self.X, self.mask = self._get_features(pdb_str)

    def _get_features(self, pdb_str):
        """ Extract protein features
        X: coordinates of Ca atoms
        mask: binary mask of Ca atoms"""
        
        # Get protein features
        X = get_token_coordinates_pdb(pdb_str)
        mask = torch.ones(len(X))
        
        X = torch.tensor(X, dtype=torch.float32)
        #mask = torch.tensor(mask, dtype=torch.float32)
        
        X = X[None]
        mask = mask[None]
        
        return X, mask
    
    def _dist(self, eps=1E-6):
            """ Pairwise euclidean distances """
            # Convolutional network on NCHW
            #print(self.X.shape, self.mask.shape)
            mask_2D = torch.unsqueeze(self.mask,1) * torch.unsqueeze(self.mask,2)
            dX = torch.unsqueeze(self.X,1) - torch.unsqueeze(self.X,2)
            D = mask_2D * torch.sqrt(torch.sum(dX**2, 3) + eps)

            # Identify k nearest neighbors (including self)
            D_max, _ = torch.max(D, -1, keepdim=True)
            D_adjust = D + 2 * (1. - mask_2D) * D_max
            D_neighbors, E_idx = torch.topk(D_adjust, min(self.top_k, self.X.shape[-2]), dim=-1, largest=False)
            mask_neighbors = gather_edges(mask_2D.unsqueeze(-1), E_idx)

            return D_neighbors, E_idx, mask_neighbors
        
        
        
def score_contact(result, reslist1, reslist2, dist=4):
    """returns a list of pairs of residues that are making contacts"""
    pdb = result['pdb']
    chains, residues, resindices = get_coordinates_pdb(pdb)
    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        pair = (res1, res2)
                        pair_rev = (res2, res1)
                        if pair not in pairs and pair_rev not in pairs:
                            pairs.append(pair)
                            contact=1

            score = score+contact

    return pairs, score

        
        
def score_rmsd(pdb1, pdb2, reslist=None):

    with open("temp1.pdb", "w") as f:
        f.write(pdb1)
    pred_struct = PandasPdb().read_pdb("temp1.pdb")

    with open("temp2.pdb", "w") as f:
        f.write(pdb2)
    true_struct = PandasPdb().read_pdb("temp2.pdb")
    
    output_path = None
    
    if not reslist:
        rmsd, mapping, aligned_pdb = multichain_permutation_alignment(
        pred_struct, true_struct, output_path
    )
    
    else:
        pred_struct = filter_residues(pred_struct, reslist)
        true_struct = filter_residues(true_struct, reslist)
        rmsd, mapping, aligned_pdb = multichain_permutation_alignment(
        pred_struct, true_struct, output_path
    )
        
    return rmsd

def filter_residues(pdb_structure, residues_to_keep: List[str]):
    """
    Filter a PandasPdb structure to only include specified residues.
    
    Parameters:
    -----------
    pdb_structure : PandasPdb object
        The PDB structure loaded using PandasPdb().read_pdb()
    residues_to_keep : List[str]
        List of residue identifiers in format "A1", "B3", etc. where the first
        character is the chain ID and the rest is the residue number
        
    Returns:
    --------
    PandasPdb object
        A new PDB structure containing only the specified residues
    """
    # Create a copy of the structure to avoid modifying the original
    filtered_struct = PandasPdb()
    
    # Get the ATOM and HETATM dataframes
    atom_df = pdb_structure.df['ATOM'].copy()
    hetatm_df = pdb_structure.df['HETATM'].copy()
    
    # Ensure residue_number is numeric
    atom_df['residue_number'] = pd.to_numeric(atom_df['residue_number'], errors='coerce')
    hetatm_df['residue_number'] = pd.to_numeric(hetatm_df['residue_number'], errors='coerce')
    
    # Create lists to store matching rows
    atom_matches = []
    hetatm_matches = []
    
    for residue_id in residues_to_keep:
        # Extract chain and residue number using regex
        match = re.match(r'([A-Za-z])(\d+)', residue_id)
        if not match:
            raise ValueError(f"Invalid residue identifier: {residue_id}. Expected format: 'A1', 'B2', etc.")
            
        chain_id, residue_num = match.groups()
        residue_num = int(residue_num)
        
        # Filter ATOM records
        chain_atoms = atom_df[
            (atom_df['chain_id'].astype(str) == chain_id) & 
            (atom_df['residue_number'] == residue_num)
        ]
        if not chain_atoms.empty:
            atom_matches.append(chain_atoms)
        
        # Filter HETATM records
        chain_hetatms = hetatm_df[
            (hetatm_df['chain_id'].astype(str) == chain_id) & 
            (hetatm_df['residue_number'] == residue_num)
        ]
        if not chain_hetatms.empty:
            hetatm_matches.append(chain_hetatms)

    # Combine all matches (only if we found any)
    filtered_struct.df['ATOM'] = pd.concat(atom_matches) if atom_matches else pd.DataFrame(columns=atom_df.columns)
    filtered_struct.df['HETATM'] = pd.concat(hetatm_matches) if hetatm_matches else pd.DataFrame(columns=hetatm_df.columns)
    
    return filtered_struct
