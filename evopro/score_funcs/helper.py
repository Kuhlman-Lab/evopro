#importing from biopython
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities
from Bio.PDB import NeighborSearch
from Bio.PDB.HSExposure import ExposureCN
from Bio.PDB import Superimposer
from Bio import SeqIO

from Bio.PDB.internal_coords import IC_Chain

#importing other modules
import os
import numpy as np #numpy isn't explicitly used in this code, but biopython runs it under the hood and calls it as np

def parse_file(filename):
    parser = PDBParser(PERMISSIVE=1)
    structure_id = "test"
    structure = parser.get_structure(structure_id, filename)
    model = structure[0]
    chain_A = model["A"]
    residues = unfold_entities(model, 'R')
    atoms = unfold_entities(structure[0], 'A')
    return model, chain_A, residues, atoms

def get_serine(residues):
    """ Scans amino acid sequence for the RRXS phosphorylation motif and returns the coordinate of the serine residue """
    for i in range(3, len(residues)):
        if residues[i].get_resname() == "SER":
            if residues[i - 2].get_resname() == "ARG" and residues[i -3].get_resname() == "ARG":
                return i + 1  #since chain A is 1-indexed but residues are 0-indexed
            

def check_hydrophobic(chain_A, serine_resid):
    """ Returns True if the residue is a canonical hydrophobic amino acid """
    hydrophobic_res = chain_A[serine_resid - 1].get_resname()
    canon_hydrophobic = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]
    # non_hydrophobic = ["ASP", "GLU", "HIS", "ARG", "LYS", "SER", "THR", "ASN", "GLN", "CYS", "SEC", "GLY", "PRO"]

    if hydrophobic_res in canon_hydrophobic:
        return True
    else:
        return False
    

def check_polar_contacts(chain_A, resid, funct_gp, atoms, radius = 4) -> bool:

    """ Checks polar contacts for the functional group of a given residue within a radius of 4 Angstrom (default radius). 
    Returns bool """

    serine_og = chain_A[resid][funct_gp]
    ns = NeighborSearch(atoms)

    close_atoms = ns.search(serine_og.coord, radius) #searches neighborhood of radius 4 (default) based on (x,y,z) coordinates

    neighbor_res = [atom.get_parent() for atom in close_atoms]
    neighbor_res_ids = [res.id[1] for res in neighbor_res]

    non_helix_ids = []
    for res in neighbor_res_ids:
        if abs(res - resid) > 4: #excluding any i+4 hydrogen bonds in the alpha helix
            if res not in non_helix_ids:
                non_helix_ids.append(res)

    #only examining charged (+/-) polar amino acids 
    polar_amino_acids = ["ASP", "GLU", "HIS", "ARG", "LYS"]

    count_polar = 0
    num_non_helix = len(non_helix_ids)

    #need to check this condition, otherwise throws error for NoneType when getting the residue name
    if num_non_helix == 0:
        return False
    else:
        for i in non_helix_ids:
            amino_acid = chain_A[i].get_resname()

            #counting number of polar amino acids within 4 Angstrom of the functional group of a desired residue
            if amino_acid in polar_amino_acids:
                count_polar += 1

    if count_polar >= 1:
        return True
    else:
        return False
    

def check_solvent_exposed(model, chain_A, serine_resid, cutoff = 25) -> bool:

    """ Checks solvent exposure using the half sphere metric from Biopython """

    hse = ExposureCN(model) #not called in function but is used under the hood

    arg1 = chain_A[serine_resid - 2]
    arg2 = chain_A[serine_resid - 3]
        
    score_arg1 = arg1.xtra["EXP_CN"]
    score_arg2 = arg2.xtra["EXP_CN"]
    
    if score_arg1 <= cutoff and score_arg2 <= cutoff:
        return True
    else: 
        return False
    

def check_no_polar_arginines(chain_A, serine_resid, atoms) -> bool:

    """ Scanning each atom of both arginine residues to check if the residues have any polar contacts, returns bool """

    chk1 = check_polar_contacts(chain_A, serine_resid - 2, "NH1", atoms)
    chk2 = check_polar_contacts(chain_A, serine_resid - 2, "NH2", atoms)
    chk3 = check_polar_contacts(chain_A, serine_resid - 2, "NE", atoms)
    chk4 = check_polar_contacts(chain_A, serine_resid - 3, "NH1", atoms)
    chk5 = check_polar_contacts(chain_A, serine_resid - 3, "NH2", atoms)
    chk6 = check_polar_contacts(chain_A, serine_resid - 3, "NE", atoms)

    chk7 = check_polar_contacts(chain_A, serine_resid - 2, "CZ", atoms)
    chk8 = check_polar_contacts(chain_A, serine_resid - 3, "CZ", atoms)

    if chk1 == False and chk2 == False and chk3 == False and chk4 == False and chk5 == False and chk6 == False and chk7 == False and chk8 == False:
        return True
    else:
        return False
    

def determine_rmsd_change(model1, model2, chain1: str = "A", chain2: str = "A") -> float:

    """ Given two input models, returns the RMSD between the pdb structures. The Superimposer function will rotate and translate
    the pdb structures automatically, so any differences in the references frames of the two structures is irrelevant. To determine 
    the C2 symmetry of a homodimer output, use the same model for model1 and model2 and provide chain A and chain B
    instead of the default options. """

    residues = unfold_entities(model1[chain1], 'R') #getting the residues in chain 1
    residues_rotated = unfold_entities(model2[chain2], 'R') #getting the residues in chain 2

    #getting the alpha carbon of each residue -- this assumes no ligands or solvents are present in .pdb
    ref_atoms = [res["CA"] for res in residues] 
    rot_atoms = [res["CA"] for res in residues_rotated]

    #Superimposer requires atom lists of equal length
    if len(ref_atoms) != len(rot_atoms):

        #Pick the smaller atom list to set the length of the two lists
        if len(ref_atoms) > len(rot_atoms):
            slice = len(rot_atoms)
        else:
            slice = len(ref_atoms)
    else:
        slice = len(rot_atoms) #this is arbitrary since the length of the two lists is the same

    #Superimposer rotates and translates the chains for you to minimize rmsd
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms[0:slice], rot_atoms[0:slice])
    super_imposer.apply(rot_atoms)

    # Return RMSD:
    return super_imposer.rms

def get_phi_psi(chain_A):
    chain_A.atom_to_internal_coordinates(verbose=True)
    phi_s = []
    psi_s = []
    for i in range(1, len(chain_A) - 1):
        phi = chain_A[i].internal_coord.get_angle("phi")
        psi = chain_A[i].internal_coord.get_angle("psi")

        #ensures any odd angle values are not included in alpha helices, but that other functions can still run
        if phi is None:
            phi_s.append(-100)
        elif psi is None:
            psi_s.append(-100)
        else:
            phi_s.append(phi)
            psi_s.append(psi)

    return phi_s, psi_s



def get_helix_cutoff(phi_s, psi_s):
    #alpha helices
    alpha_psi = [i for i, x in enumerate(psi_s[1:]) if -48 <= round(x) <= -34]
    alpha_phi = [i for i, x in enumerate(phi_s[1:]) if -71 <= round(x) <= -57]

    alpha_ids = [value for value in alpha_phi if value in alpha_psi]

    #compute consecutive difference between indicies i.e. how large are the jumps between helices
    consecutive_diff = [x - alpha_ids[i - 1] for i, x in enumerate(alpha_ids)][1:]

    #get the fourth largest jump -- theoretically four jumps corresponds to the four loops between helices
    cutoff= sorted(consecutive_diff)[-4]

    return alpha_ids, consecutive_diff, cutoff


def define_helix_boundaries(alpha_ids, consecutive_diff, cutoff):

    start_nterm = alpha_ids[0]
    end_cterm = alpha_ids[-1] + 1
    helices = []
    for i in range(len(consecutive_diff)):
        if consecutive_diff[i] >= cutoff:
            helices.append(alpha_ids[i])

    end_nterm = helices[0] + 1 #because indexing is [start, stop) in python
    last_jump = helices[-1] 

    #get index after last jump, ex. alpha_ids[74] == 93 to alpha_ids[75] == 99 to capture 99 as first index
    start_cterm_alpha = np.where(np.asarray(alpha_ids) == last_jump)[0][0] + 1 # 74 + 1 -> 75 (per example)
    start_cterm = alpha_ids[start_cterm_alpha] #start indexing counts at 99 (per example)

    return start_nterm, end_nterm, start_cterm, end_cterm

def check_helix_contacts(chain_A, helix_ids, funct_gp, atoms, radius = 6) -> bool:

    """ Checks polar contacts for the functional group of a given residue within a radius of 6 Angstrom (default radius). 
    Returns bool """

    count = 0 #counting residue contacts that are not in helix_ids
    non_helix_ids = []

    for resid in helix_ids:
        helix_residue = chain_A[resid + 1][funct_gp]
        ns = NeighborSearch(atoms)

        close_atoms = ns.search(helix_residue.coord, radius) #searches neighborhood of radius 4 (default) based on (x,y,z) coordinates

        neighbor_res = [atom.get_parent() for atom in close_atoms]
        neighbor_res_ids = [res.id[1] for res in neighbor_res]

        for res in neighbor_res_ids:
            if res not in helix_ids:
                if res not in non_helix_ids:
                    count += 1
                    non_helix_ids.append(res)


    return count, non_helix_ids

def get_avg_bfactor(residue_ids, residues):

    plddt_sum = 0
    plddt_count = 0

    for i in residue_ids:
        atoms = residues[i].get_unpacked_list()
        for atom in atoms:
            plddt = atom.get_bfactor()
            plddt_sum += plddt
            plddt_count += 1

    return plddt_sum/plddt_count

def number_contacts(chain_A, resid, funct_gp, atoms, radius = 4):

    """ Checks polar contacts for the functional group of a given residue within a radius of 4 Angstrom (default radius). 
    Returns bool """

    serine_og = chain_A[resid][funct_gp]
    ns = NeighborSearch(atoms)

    close_atoms = ns.search(serine_og.coord, radius) #searches neighborhood of radius 4 (default) based on (x,y,z) coordinates

    neighbor_res = [atom.get_parent() for atom in close_atoms]
    neighbor_res_ids = [res.id[1] for res in neighbor_res]

    unique_res = set(neighbor_res_ids)

    return len(unique_res)

def define_all_helix_boundaries(alpha_ids, consecutive_diff, cutoff):

    helices = []
    for i in range(len(consecutive_diff)):
        if consecutive_diff[i] >= cutoff:
            helices.append(alpha_ids[i])

    start_nterm = alpha_ids[0]
    end_nterm = helices[0]

    start_next_id = np.where(np.asarray(alpha_ids) == end_nterm)[0][0] + 1
    start_next_helix = alpha_ids[start_next_id]
    end_next_helix = helices[1]

    start_next_next_id = np.where(np.asarray(alpha_ids) == end_next_helix)[0][0] + 1
    start_next_next_helix = alpha_ids[start_next_next_id]
    end_next_next_helix = helices[2]

    start_penultimate_id = np.where(np.asarray(alpha_ids) == end_next_next_helix)[0][0] + 1
    start_penultimate = alpha_ids[start_penultimate_id]
    end_penultimate = helices[3]

    start_cterm_alpha = np.where(np.asarray(alpha_ids) == end_penultimate)[0][0] + 1
    start_cterm = alpha_ids[start_cterm_alpha]
    end_cterm = alpha_ids[-1]

    return start_nterm, end_nterm, start_next_helix, end_next_helix, start_next_next_helix, end_next_next_helix, start_penultimate, end_penultimate, start_cterm, end_cterm


def get_fifth_helix(chain_A, residues):
    """ Used to determine the residue ids for the 5th helix in the 5-helical bundle given the chain A and residues
    from the AlphaFold monomer prediction"""
    phi_s, psi_s = get_phi_psi(chain_A)
    alpha_ids, consecutive_diff, cutoff = get_helix_cutoff(phi_s[1:], psi_s[1:])
    # jump1 = sorted(consecutive_diff)[-1]
    # jump2 = sorted(consecutive_diff)[-2]
    # jump3 = sorted(consecutive_diff)[-3]
    # jump4 = sorted(consecutive_diff)[-4]
    # last_diff = jump3 - jump4
    # avg_jump = ((jump1 - jump2) + (jump2 - jump3) + (last_diff))/3
    # if last_diff > avg_jump:
    #     raise Exception("This structure has poorly defined helix boundaries")
    # else:
    start_nterm, end_nterm, start_cterm, end_cterm = define_helix_boundaries(alpha_ids, consecutive_diff, cutoff)
    serine_res = get_serine(residues)
    helix_ids_nterm = [i for i in range(start_nterm, end_nterm + 1)]
    helix_ids_cterm = [i for i in range(start_cterm, end_cterm + 1)]
    if serine_res in helix_ids_nterm:
        return helix_ids_nterm, helix_ids_cterm
    else:
        return helix_ids_cterm, helix_ids_nterm
    
def get_helix_motif(pdb):
    #creating a tempory pdb file to store output
    with open("temp.pdb", "w") as f:
        f.write(pdb)
    
    #changing default max peptide bond, since it is too stringent for design purposes
    IC_Chain.MaxPeptideBond = 4.0
    model, chain_A, residues, atoms = parse_file("temp.pdb")

    chain_A.internal_coord = None  # force re-loading structure data with new cutoff
    chain_A.atom_to_internal_coordinates(verbose=True)
    
    helix_motif, helix_no_motif = get_fifth_helix(chain_A, residues)
    
    os.remove("temp.pdb") #removing the file since only a temporary file is needed
    
    return helix_motif
    
def calc_residue_dist(res, res2) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = res["CA"].coord - res2["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


    
def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    distance_mat = np.zeros((len(chain_one), len(chain_two)))
    for i, res1 in enumerate(chain_one) :
        for j, res2 in enumerate(chain_two) :
            distance_mat[i, j] = calc_residue_dist(res1, res2)
    return distance_mat

def get_serine_dimer_distance(resids, helix_motif, serine_res, distance_matrix_dimer):
    """ Return the average distance between the residues serine is in contact with in the monomer when considering
    those same residues in the homodimer. Ideally, want the distance between the serine and the monomer contacts to
    increase in the dimer, since the helix with the motif will undock.
    Returns float"""
    lst_res = list(resids)
    sum_dist = 0
    counter = 0
    for res in lst_res:
        if res not in helix_motif:
            dist = distance_matrix_dimer[serine_res, res]
            counter += 1
            sum_dist += dist
    if counter != 0:
        return sum_dist/counter
    else:
        return 0
    
def return_contacts(chain_A, resid, funct_gp, atoms, radius = 4):
    serine_og = chain_A[resid][funct_gp]
    ns = NeighborSearch(atoms)
    close_atoms = ns.search(serine_og.coord, radius) #searches neighborhood of radius 4 (default) based on (x,y,z) coordinates
    neighbor_res = [atom.get_parent() for atom in close_atoms]
    neighbor_res_ids = [(res.id[1] - 1) for res in neighbor_res] #gives chain number instead of residue number (i.e. vals need to be converted to zero indexed)
    unique_res = set(neighbor_res_ids)
    return unique_res


def get_weighted_solvent_exposed(pdb_string1, pdb_string2, helix_motif):
    with open("temp.pdb", "w") as f:
        f.write(pdb_string1)
    with open("temp2.pdb", "w") as f:
        f.write(pdb_string2)

    #changing default max peptide bond, since it is too stringent for design purposes
    IC_Chain.MaxPeptideBond = 4.0
    model, chain_A, residues, atoms = parse_file("temp.pdb")

    chain_A.internal_coord = None  # force re-loading structure data with new cutoff
    chain_A.atom_to_internal_coordinates(verbose=True)

    serine_res = get_serine(residues)
    # helix_motif, helix_no_motif = get_fifth_helix(chain_A, residues) #run if helix_motif = None
    model_dimer, chain_A_dimer, residues_dimer, atoms_dimer = parse_file("temp2.pdb")
    distance_matrix_dimer_chainA = calc_dist_matrix(chain_A_dimer, chain_A_dimer)
    resids = return_contacts(chain_A, serine_res, "OG", atoms, radius = 4)
    dist = get_serine_dimer_distance(resids, helix_motif, serine_res, distance_matrix_dimer_chainA)
    resids = return_contacts(chain_A, serine_res, "OG", atoms, radius = 4)
    serine_hse_mono = number_contacts(chain_A, serine_res, "OG", atoms, radius = 8)
    serine_hse_dim = number_contacts(chain_A_dimer, serine_res, "OG", atoms_dimer, radius = 8)
    weight = (serine_hse_mono - serine_hse_dim)/serine_hse_dim
    weighted_dist = weight * dist

    os.remove("temp.pdb")
    os.remove("temp2.pdb")
    
    return weighted_dist


def get_weighted_solvent_exposed_return_ser_contacts(pdb_string1, pdb_string2, helix_motif):
    with open("temp.pdb", "w") as f:
        f.write(pdb_string1)
    with open("temp2.pdb", "w") as f:
        f.write(pdb_string2)

    #changing default max peptide bond, since it is too stringent for design purposes
    IC_Chain.MaxPeptideBond = 4.0
    model, chain_A, residues, atoms = parse_file("temp.pdb")

    chain_A.internal_coord = None  # force re-loading structure data with new cutoff
    chain_A.atom_to_internal_coordinates(verbose=True)

    serine_res = get_serine(residues)
    # helix_motif, helix_no_motif = get_fifth_helix(chain_A, residues) #run if helix_motif = None
    model_dimer, chain_A_dimer, residues_dimer, atoms_dimer = parse_file("temp2.pdb")
    distance_matrix_dimer_chainA = calc_dist_matrix(chain_A_dimer, chain_A_dimer)
    resids = return_contacts(chain_A, serine_res, "OG", atoms, radius = 4)
    dist = get_serine_dimer_distance(resids, helix_motif, serine_res, distance_matrix_dimer_chainA)
    resids = return_contacts(chain_A, serine_res, "OG", atoms, radius = 4)
    serine_hse_mono = number_contacts(chain_A, serine_res, "OG", atoms, radius = 8)
    serine_hse_dim = number_contacts(chain_A_dimer, serine_res, "CA", atoms_dimer, radius = 8)
    weight = (serine_hse_mono - serine_hse_dim)/serine_hse_dim
    weighted_dist = weight * dist

    os.remove("temp.pdb")
    os.remove("temp2.pdb")
    
    return weighted_dist, serine_hse_mono


def rmsd_to_original_pdb(pdb_string1, original_pdb):
    with open("temp.pdb", "w") as f:
        f.write(pdb_string1)

    model, chain_A, residues, atoms = parse_file("temp.pdb")
    model_org, chain_A_org, residues_org, atoms_org = parse_file(original_pdb)

    rmsd = determine_rmsd_change(model, model_org) #default is comparing chain A in model to chain A in model_org

    os.remove("temp.pdb")

    return rmsd



def get_distance_helix3_4(pdb_string1, pdb_string2, helix_motif):
    with open("temp.pdb", "w") as f:
        f.write(pdb_string1)
    with open("temp2.pdb", "w") as f:
        f.write(pdb_string2)

    #changing default max peptide bond, since it is too stringent for design purposes
    IC_Chain.MaxPeptideBond = 4.0
    model, chain_A, residues, atoms = parse_file("temp.pdb")

    chain_A.internal_coord = None  # force re-loading structure data with new cutoff
    chain_A.atom_to_internal_coordinates(verbose=True)

    serine_res = get_serine(residues)
    # helix_motif, helix_no_motif = get_fifth_helix(chain_A, residues) #run if helix_motif = None
    model_dimer, chain_A_dimer, residues_dimer, atoms_dimer = parse_file("temp2.pdb")
    distance_matrix_dimer_chainA = calc_dist_matrix(chain_A_dimer, chain_A_dimer)
    resids = return_contacts(chain_A, serine_res, "OG", atoms, radius = 4)
    dist = get_serine_dimer_distance(resids, helix_motif, serine_res, distance_matrix_dimer_chainA)

    os.remove("temp.pdb")
    os.remove("temp2.pdb")
    
    return dist



def check_serine_polar_contacts(pdb_string1, polar_list = ["ASP", "GLU"]):
    with open("temp3.pdb", "w") as f:
        f.write(pdb_string1)

    model, chain_A, residues, atoms = parse_file("temp3.pdb")

    serine_res = get_serine(residues)

    serine_og = chain_A[serine_res]["OG"]
    ns = NeighborSearch(atoms)

    close_atoms = ns.search(serine_og.coord, 4) #smaller search to increase likelihood of hydrogen bonding

    neighbor_res = [atom.get_parent() for atom in close_atoms]
    neighbor_res_ids = [res.id[1] for res in neighbor_res]

    non_helix_ids = []
    for res in neighbor_res_ids:
        if abs(res - serine_res) > 4:
            if res not in non_helix_ids:
                non_helix_ids.append(res)

    count_polar = 0
    num_non_helix = len(non_helix_ids)
    if num_non_helix == 0:
        return False
    else:
        for i in non_helix_ids:
            amino_acid = chain_A[i].get_resname()
            if amino_acid in polar_list:
                count_polar += 1

    os.remove("temp3.pdb")

    if count_polar >= 1:
        return True
    else:
        return False
    

def load_file(path):
  #loading in the structure
  IC_Chain.MaxPeptideBond = 4.0
  model, chain_A, residues, atoms = parse_file(path)

  chain_A.internal_coord = None  # force re-loading structure data with new cutoff
  chain_A.atom_to_internal_coordinates(verbose=True)

  return model, chain_A, residues, atoms


def check_hydrophobic(residue) -> bool:
    """
    Returns True/False depending on whether a given residue is hydrophobic.

    Args:
        residue (Bio.PDB.Residue.Residue): The residue object from Biopython

    Returns:
        (bool): Returns True if the residue is hydrophobic and false if the residue is not
        hydrophobic
    """
    residue_name = residue.get_resname()

    #defining which amino acids are hydrophobic
    hydrophobic_res_names = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

    if residue_name in hydrophobic_res_names:
        return True
    else:
        return False


def get_hydrophobic_ids(chain_A, residues: list) -> list:
    """
    Returns a list of hydrophobic indices given the path to a pdb file

    Args:
        chain_A: a string pointing to the pdb file
        
    Returns:
        indices (list): a list of the indicies (for chain A) for the 
        hydrophobic residues in a structure
    """

    #getting the residue indicies
    res_ids = [res.id[1] for res in residues]

    #getting the HSE scores for hydrophobic and non-hydrophobic residues
    indices = []

    for i in res_ids:
        res = chain_A[i]
        if check_hydrophobic(res): #checking if the residue is hydrophobic
            indices.append(i) #appending to list if hydrophobic
        else:
            continue

    return indices



from collections import defaultdict

def create_mod_dict(indices: list) -> dict:
    """
    Creates a dictionary of all residues in the same equivalence class. For example, indices 2 and 10 are 
    associated with key 2 and indices 0 and 16 are associated with key 0 (since 0 mod 4 is 0 and 16 mod 4 is 0).
    
    Args:
        indices (list): a list of indices
        
    Returns:
        mod_dict (dict): a dictionary of four keys (0, 1, 2, and 3), where each key contains a list of indices in the 
        same equivalence class
    """

    #initializing the dictionary
    mod_dict = defaultdict(list)

    #adding to dictionary depending on value mod 4, ex. 33 is 1 mod 4 since 8*4 + 1 = 33
    for idx in indices:
        mod_dict[idx % 4].append(idx)

    return mod_dict

def get_hse_dim_vs_mono_score(pdb_string1, pdb_string2):
    """
    Returns the highest difference in HSE scores between dimer and monomer when comparing indices in the same equivalence class (mod 4). This is because
  alpha helices turn every 3.6 residues (~4), so hydrophobic residues 4 apart from each other could create a contact point for dimerization in the absense
  of kinase. By comparing HSE scores, it can be determined if AF2 predicts that those residues indeed dimerize (HSE is higher in dimer than monomer). 
  The highest score across all equivalence classes (residues mod 4) is returned.

  Args:
      monomer_path (str): the path to the monomer pdb file
      dimer_path (str): the path to the dimer pdb file

  Returns:
      highest_score (int): the highest sum of hse scores across all residue equivalence classes
    """
    with open("temp4.pdb", "w") as f:
        f.write(pdb_string1)

    with open("temp5.pdb", "w") as f:
        f.write(pdb_string2)

    #loading in the monomer and dimer files
    model, chain_A, residues, atoms = load_file("temp4.pdb")
    model_dimer, chain_A_dimer, residues_dimer, atoms_dimer = load_file("temp5.pdb")

    #getting the indices for all hydrophobic residues
    hydrophobic_indices = get_hydrophobic_ids(chain_A, residues) #residue numbers are the same for chain A in monomer and chain A in homodimer

    #getting the dictionary of modulo values
    mod_dict = create_mod_dict(hydrophobic_indices)

    #initializing the half solvent exposure models for the monomer and dimer
    hse = ExposureCN(model)
    hse2 = ExposureCN(model_dimer)

    counter = 0
    for key, value in mod_dict.items():
        counter += 1
        print(f"beginning key {key}")
        scores = []
        for i in value:
            score_monomer = chain_A[i].xtra["EXP_CN"]
            score_dimer = chain_A_dimer[i].xtra["EXP_CN"]
            sub = score_dimer - score_monomer
            scores.append(score_dimer - score_monomer)
            print(i, sub) #for debugging purposes

        total_score = sum(scores)
        if counter == 1:
            high_score = total_score
        else:
            if total_score > high_score:
                high_score = total_score
            else:
                high_score = high_score
        print(f"this is the sum of the scores: {total_score}")

    print(f"the high score is {high_score}")
    return high_score



def polar_positive_pentalty(pdb_string1):
    with open("temp6.pdb", "w") as f:
        f.write(pdb_string1)

    model, chain_A, residues, atoms = parse_file("temp6.pdb")

    serine_res = get_serine(residues)

    serine_og = chain_A[serine_res]["OG"]
    ns = NeighborSearch(atoms)

    close_atoms = ns.search(serine_og.coord, 4) #smaller search to increase likelihood of hydrogen bonding

    neighbor_resnames = [atom.get_parent().get_resname() for atom in close_atoms]
    polar_positive = ["ARG", "HIS", "LYS"]

    penalty = 0
    for name in neighbor_resnames:
        if name in polar_positive:
            penalty += 1

    os.remove("temp6.pdb")

    return penalty


def glu_asp_bonus(pdb_string1, helix):
    polar_list = ["ASP", "GLU"]

    with open("temp7.pdb", "w") as f:
        f.write(pdb_string1)

    model, chain_A, residues, atoms = parse_file("temp7.pdb")

    serine_res = get_serine(residues)

    serine_og = chain_A[serine_res]["OG"]
    ns = NeighborSearch(atoms)

    close_atoms = ns.search(serine_og.coord, 4) #smaller search to increase likelihood of hydrogen bonding

    bonus = 0
    for atom in close_atoms:
        atom_name = atom.get_id()
        if "O" in atom_name:
            residue_name = atom.get_parent().get_resname()
            if residue_name in polar_list:
                residue_id = atom.get_parent().id[1]
                if residue_id not in helix:
                    bonus += 1

    os.remove("temp7.pdb")

    return bonus