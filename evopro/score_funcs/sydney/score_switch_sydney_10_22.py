from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import os
from evopro.score_funcs.helper import get_helix_motif, rmsd_to_original_pdb, get_distance_helix3_4, polar_positive_pentalty, glu_asp_bonus
from alphafold.common import protein
import math

import subprocess

#IMPORTANT: DID NOT FINISH IMPLEMENTING

#importing from biopython
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities

from Bio.PDB.internal_coords import IC_Chain

def determine_third_and_fourth_helices(serine_idx, all_helices):

    for helix in all_helices:

        #helix is a tuple where helix[0] is the first index (residue number) of the helix and helix[1] is the last index
        start_helix = helix[0]
        stop_helix = helix[1]

        helix_ids = [*range(start_helix, stop_helix + 1, 1)] #ex. (0,4) becomes the list [0, 1, 2, 3, 4]

        #don't want the fifth helix to be considered, since its average distance to serine will likely be the lowest
        if serine_idx in helix_ids:
            return helix

stride_path = "/nas/longleaf/home/mghvasta/stride/stride"

def get_num_helices(pdb_file):
    """
    Determines the number of alpha helices present in a pdb structure and returns the number.
    
    Args:
        pdb_path (str): path to the pdb file to analyze
        
    Returns:
        number_helices (int): the number of helices found in the structure"""

    stride_command = [stride_path, '-f', pdb_file]
    try:
        stride_output = subprocess.check_output(stride_command, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running STRIDE for {pdb_file}: {e}")

    hel_num = 0
    all_helices = []
    # Process the STRIDE output
    for line in stride_output.splitlines():
        # Check if the line starts with "ASG" (indicating secondary structure assignment)
        if line.startswith("LOC"):
            # Split the line into columns
            columns = line.split()
            # The sixth column contains the secondary structure class
            if columns[1] == "AlphaHelix":
                hel_num += 1
                start_stop = (int(columns[3]), int(columns[6]))
                all_helices.append(start_stop)

    return hel_num, all_helices


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
            

def score_complex_confidence(results, dsobj):
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of complex is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = plddt_potential
    print(f"Final score is {score}")
    return score, (score, confscore2)


def score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dist=4, contact_cap=36, dsobj=None, first_only=False):
    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)

    chains, residues, resindices = get_coordinates_pdb(pdb)
    pae = results['pae_output'][0]

    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            weight = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        pair = (res1, res2)
                        pair_rev = (res2, res1)
                        if pair not in pairs and pair_rev not in pairs:
                            if len(pairs)<contact_cap:
                                contact=1
                                res1_id = resindices[res1]
                                res2_id = resindices[res2]
                                pae_contact = pae[res1_id][res2_id] + pae[res2_id][res1_id]
                                weight = (70-pae_contact)/70
                                pairs.append(pair)

            score = score + contact*weight
    
    return pairs, score



def score_binder_complex_2(results, dsobj, contacts):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = [x for x in contacts[0] if x.startswith("B")] #added for proteinC
    reslist2 = [x for x in residues.keys() if x.startswith("C")] #added for proteinC
    print("reslist1", reslist1) #added for proteinC
    print("reslist2", reslist2) #added for proteinC
    contacts_list, contactscoreBC = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    score = 5*contactscoreBC # + nonpolar_penalty*5
    print(score, (score, len(contacts), contactscoreBC))
    return score, (score, len(contacts), contactscoreBC)



def score_overall(results:list, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    print("Number of predictions being scored:", len(results))
    pdbs = []
    
    #note: results is in the same order as --af2_preds AB,C,CC so in this case results[0] is AB, results[1] is C, and results[2] is CC
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    # ptm_monomer_array = results[1]['ptm']
    # ptm_monomer = ptm_monomer_array.flat[0]
    iptm_phos_dimer_array = results[0]['iptm'] #want confidence in the binder contacts to be high when the 5th helix is disordered
    iptm_phos_dimer = iptm_phos_dimer_array.flat[0]
    iptm_dimer_array = results[2]['iptm'] #want confidence in the binder contacts to be low when the 5th helix is ordered
    iptm_dimer = iptm_dimer_array.flat[0]

    rmsd_to_starting = rmsd_to_original_pdb(pdb2, "original.pdb") #want the fold to be close to the starting structure
    rmsd_cutoff = 4.0
    spring_constant = 10
    if rmsd_to_starting > rmsd_cutoff:
        rmsd = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
    else:
        rmsd = 0

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)

    monomer_file = "temp.pdb"
    with open(monomer_file, "w") as file:
        file.write(pdb1)
    file.close()

    model, chain_A, residues, atoms = parse_file(monomer_file)
    serine_idx = get_serine(residues)
    serine_res = chain_A[serine_idx]
    
    #monomer is the last prediction, we want to use the monomer to get the helix motif
    try:
        hel_num, all_helices = get_num_helices(monomer_file)
        if hel_num != 5:
            return 1000, ((1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results
    except:
        return 1000, ((1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results

    #dict_keys(['pae_output', 'ranking_confidence', 'plddt', 'structure_module', 'ptm', 'iptm', 'unrelaxed_protein'])

    fifth_helix = determine_third_and_fourth_helices(serine_idx, all_helices)

    bonus = glu_asp_bonus(pdb2, fifth_helix) 

    score, score_tuple = score_complex_confidence(results, dsobj)

    #computing the overall score
    overall_score = + iptm_dimer - iptm_phos_dimer - bonus - rmsd + score

    return overall_score, ((ptm_monomer,), (iptm_dimer,), (iptm_phos_dimer,), (penalty,), (bonus,), (rmsd,), score_tuple), pdbs, results