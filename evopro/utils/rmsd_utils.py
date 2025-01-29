import torch
from biopandas.pdb import PandasPdb
from tmtools import tm_align
from itertools import permutations 

def apply_transform(A,R,t):
    A_aligned = torch.bmm(R, A.transpose(1,2)).transpose(1,2) + t
    return A_aligned

def kabsch(A, B):
        a_mean = A.mean(dim=1, keepdims=True).type('torch.DoubleTensor')
        b_mean = B.mean(dim=1, keepdims=True).type('torch.DoubleTensor')
        A_c = A - a_mean
        B_c = B - b_mean
        # Covariance matrix
        H = torch.bmm(A_c.transpose(1,2), B_c)  # [B, 3, 3]
        U, S, V = torch.svd(H)
        # Rotation matrix
        R = torch.bmm(V, U.transpose(1,2))  # [B, 3, 3]
        # Translation vector
        t = b_mean - torch.bmm(R, a_mean.transpose(1,2)).transpose(1,2)
        A_aligned = apply_transform(A,R,t)
        # A_aligned = torch.bmm(R, A.transpose(1,2)).transpose(1,2) + t
        return A_aligned, R, t

def get_chain_ca_coords(df,chain_ids=None):
    if chain_ids!=None:
        df = df[df['chain_id'].isin(chain_ids)]
    else:
        pass 
    df = df[df['atom_name']=='CA']
    X = df['x_coord'].tolist()
    Y = df['y_coord'].tolist()
    Z = df['z_coord'].tolist()
    coords = torch.tensor([X,Y,Z]).transpose(0,1).unsqueeze(0)
    return coords

def crop_strucs(coord1,coord2):
    coord1_length = coord1.shape[1]
    coord2_length = coord2.shape[1]
    if coord1_length < coord2_length:
        min_length = coord1_length
    elif coord2_length < coord1_length:
        min_length = coord2_length
    else:
        min_length = coord1_length
    
    cropped_coord1 = coord1[:,:min_length,:]
    cropped_coord2 = coord2[:,:min_length,:]
    return cropped_coord1, cropped_coord2

def get_coords(df):
    X = df['x_coord'].tolist()
    Y = df['y_coord'].tolist()
    Z = df['z_coord'].tolist()
    coords = torch.tensor([X,Y,Z]).transpose(0,1).unsqueeze(0)
    return coords

def get_transform(A,B):
    a_mean = A.mean(dim=1, keepdims=True).type('torch.DoubleTensor')
    b_mean = B.mean(dim=1, keepdims=True).type('torch.DoubleTensor')
    A_c = A - a_mean
    B_c = B - b_mean
    # Covariance matrix
    H = torch.bmm(A_c.transpose(1,2), B_c)  # [B, 3, 3]
    U, S, V = torch.svd(H)
    # Rotation matrix
    R = torch.bmm(V, U.transpose(1,2))  # [B, 3, 3]
    # Translation vector
    t = b_mean - torch.bmm(R, a_mean.transpose(1,2)).transpose(1,2)
    return R, t

def apply_transform(A,R,t):
    A_aligned = torch.bmm(R, A.transpose(1,2)).transpose(1,2) + t
    return A_aligned

def calc_rmsd(aligned_aa_struc,ref_struc):
    ref_struc_cropped, aligned_aa_struc_cropped = crop_strucs(ref_struc,aligned_aa_struc)
    rmsd = torch.mean(torch.sqrt(torch.sum((aligned_aa_struc_cropped-ref_struc_cropped).pow(2),-1)),-1)
    return rmsd

def align_struc_kabsch(pdb1,pdb2,chain_ids1=None,chain_ids2=None):
    ppdb1 = PandasPdb().read_pdb(pdb1)
    ppdb2 = PandasPdb().read_pdb(pdb2)
    align_region1 = get_chain_ca_coords(ppdb1.df['ATOM'],chain_ids1)
    align_region2 = get_chain_ca_coords(ppdb2.df['ATOM'],chain_ids2)
    cropped_align_region1,cropped_align_region2 = crop_strucs(align_region1,align_region2)
    R, t = get_transform(cropped_align_region2,cropped_align_region1)
    ca_fullstruc_pred = ppdb2.df['ATOM'][ppdb2.df['ATOM']['atom_name']=='CA']
    ca_fullstruc_ref = ppdb1.df['ATOM'][ppdb1.df['ATOM']['atom_name']=='CA']
    aligned_aa_struc = apply_transform(get_coords(ca_fullstruc_pred).type(torch.float64),R,t)
    aa_rmsd = calc_rmsd(aligned_aa_struc,get_coords(ca_fullstruc_ref))
    return aa_rmsd

#function to parse pdb file and get list of chains
#then calculates entity list where an entity is a unique chain sequence and its length is the number of chains with that sequence
def get_entity_list(pdb):
    ppdb = PandasPdb().read_pdb(pdb)
    chains = ppdb.df['ATOM']['chain_id'].unique()
    entity_list = []
    for chain in chains:
        entity_list.append(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id']==chain]['residue_name'].tolist())
    return entity_list
    

def compute_permutation_alignment(pdb1,pdb2,homo1=None,homo2=None):
    
    chains = get_entity_list(pdb1)
    print(chains)
    rmsds = []
    
    ppdb1 = PandasPdb().read_pdb(pdb1)
    ppdb2 = PandasPdb().read_pdb(pdb2)
    
    if homo1:
        num_chains = len(homo1)
        perm = permutations(homo1, num_chains) 
        for p in perm:
            print(p)
        
            
            
    
    # ppdb1 = PandasPdb().read_pdb(pdb1)
    # ppdb2 = PandasPdb().read_pdb(pdb2)
    # align_region1 = get_chain_ca_coords(ppdb1.df['ATOM'],chain_ids1)
    # align_region2 = get_chain_ca_coords(ppdb2.df['ATOM'],chain_ids2)
    # cropped_align_region1,cropped_align_region2 = crop_strucs(align_region1,align_region2)
    # R, t = get_transform(cropped_align_region2,cropped_align_region1)
    # ca_fullstruc_pred = ppdb2.df['ATOM'][ppdb2.df['ATOM']['atom_name']=='CA']
    # ca_fullstruc_ref = ppdb1.df['ATOM'][ppdb1.df['ATOM']['atom_name']=='CA']
    # aligned_aa_struc = apply_transform(get_coords(ca_fullstruc_pred).type(torch.float64),R,t)
    # aa_rmsd = calc_rmsd(aligned_aa_struc,get_coords(ca_fullstruc_ref))
    # return min(rmsds)



if __name__ == "__main__":
    # path = "/work/users/a/m/amritan/evopro_tests/rmsd/for_nikka/test1/"
    # pdb1 = path + "design_56.pdb"
    # pdb2 = path + "seq_0_final_model_1_chainAB.pdb"
    
    path = "/work/users/a/m/amritan/lpl/diff_linker/trimer/binder_2483/run1/evopro/design_8/"
    pdb1 = path + "design_8_split.pdb"
    pdb2 = path + "pdb_renum_flipped.pdb"
    #print(align_struc_kabsch(pdb1,pdb2, chain_ids1=["A","B","C"], chain_ids2=["A","B","C"]))  
    print(compute_permutation_alignment(pdb1,pdb2, homo1=["A","B","C"], homo2=["D","E","F"]))