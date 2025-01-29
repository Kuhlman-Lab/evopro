from pyrosetta import *
import numpy as np
import os

from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import os
from evopro.score_funcs.helper import get_helix_motif, get_weighted_solvent_exposed_return_ser_contacts, rmsd_to_original_pdb
from alphafold.common import protein

from pyrosetta.toolbox.mutants import mutate_residue

from pyrosetta.rosetta.core.scoring import get_score_function

init()

def score_overall(results:list, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    
    #note: results is in the same order as --af2_preds AB,A so in this case results[0] is AB and results[1] is A
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])

    rmsd_to_starting = rmsd_to_original_pdb(pdb2, "original.pdb")
    if rmsd_to_starting > 2:
        return 1000, (1000, 1000, (1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results
    

    gyration = Rg(pdb2, chnid="A") #want a limit of 15, per Tomiris's score function
    
    #monomer is the last prediction, we want to use the monomer to get the helix motif
    try:
        hel = get_helix_motif(pdb2)
    except:
        return 1000, (1000, 1000, (1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results

    #hel = get_helix_motif(pdb2) #this was the original code without the try/except block
    
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            # complexscore1, complex_tuple = score_binder_complex(result, dsobj, contacts) #tuple is (score, len(contacts), contactscore)
            complexscore2, confidence_tuple = score_complex_confidence(result, dsobj, helix_motif=hel) #confidence_tuple is (confscore1, confscore2)
        else:
            monomer, monomer_tuple = score_binder_monomer(result, dsobj, helix_motif=hel) #monomer_tuple is (score, confscore2)

    #compare the two monomers in the complex
    #want a low rmsd between the monomers in the homodimer since homodimers tend to be C2 symmetric
    #low rmsd between the monomers is a necessary, though not sufficient, condition for homodimerization
    rmsd_symm = score_rmsd(pdb1, pdb1, chains1="A", chains2="B", dsobj=None) 
    
    #compare one monomer in the dimer to the monomeric prediction
    #negative because we want to maximize this score term
    rmsd_switch = score_rmsd(pdb1, pdb2, chains1="A", chains2="A", dsobj=None)
    #score.append(gyration_score(pdb2, dsobj=None))
    dist, serine_hse_mono = get_weighted_solvent_exposed_return_ser_contacts(pdb2, pdb1, helix_motif = hel)
    print("Dist", dist)


    delta_energy, starting_energy = get_phosphorlyation_energy(pdb2)
    

    # overall_score = complexscore1 + complexscore2 + monomer + rmsd1 + rmsd2 - dist - serine_hse_mono - 1/rmsd3 #negative dist to maximize score term
    overall_score = complexscore2 - monomer + rmsd_symm - rmsd_switch - dist - serine_hse_mono - (1/rmsd_to_starting) + (gyration/15) - delta_energy + starting_energy #negative dist to maximize score term

    return overall_score, (confidence_tuple, monomer_tuple, (rmsd_symm,), (rmsd_switch,), (dist,), (serine_hse_mono,), (rmsd_to_starting,), (gyration,), (delta_energy,), (starting_energy,)), pdbs, results

def score_binder_complex(results, dsobj, contacts):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]

    print("reslist1", reslist1)
    print("reslist2", reslist2)
    
    contacts_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    #test starts here

    score = contactscore
    print(score, (score, len(contacts), contactscore))
    return score, (score, len(contacts), contactscore)

def score_complex_confidence(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
        
    reslist1 = ["A"+str(x+1) for x in helix_motif]
    reslist1 += ["B"+str(x+1) for x in helix_motif]
    reslist2 = [x for x in residues.keys() if x not in reslist1]
    
    confscore1 = score_plddt_confidence(results, reslist1, resindices)

    confscore2 = score_plddt_confidence(results, reslist2, resindices)

    #penalizing a high plddt value and especially penalizing a high plddt for the fifth helix
    score = confscore2/10 + confscore1/5 
    print(score)
    
    #first value is added to overall score, second value is returned for debugging purposes
    return score, (confscore1, confscore2)

def score_binder_monomer(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist = [x for x in residues.keys()]

    confscore = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False) 

    score = confscore/20
    print(score)
    return score, (score, confscore)

def score_rmsd(pdb1, pdb2, chains1="AB", chains2="AB", dsobj=None):
   
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues0.keys() if x[0] in chains1]

    chains, residues, resindices = get_coordinates_pdb(pdb1)
    reslist2 = [x for x in residues.keys() if x[0] in chains2]

    rmsd_to_starting = get_rmsd(reslist1, pdb2, reslist2, pdb1, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting


def get_chain_length(chain_1: int, chain_2: int) -> int:
    """
    Given a chain number (1 corresponds to A, 2 corresponds to chain B,
    etc.), returns the length of the chain. Function is called in find_serine()
    
    Args:
        chain_1 (int): The number corresponding to the first chain
        chain_2 (int): The number corresponding to the second chain (can be equal
            to chain_1)
        
    Returns:
        chain_length (int): The length of the chain
    """

    start = pyrosetta.rosetta.core.pose.conformation().chain_begin(chain_1) - 1 #getting the first index of the chain
    end = pyrosetta.rosetta.core.pose.conformation().chain_end(chain_2) #getting the last index of the chain
    chain_length = end - start + 1 #add 1 to include all indices ex. a chain A with residues 0 to 78 has 79 residues

    return chain_length


def find_serine(pose: pyrosetta.rosetta.core.pose.Pose, chain_1: int, chain_2: int) -> int:
    """
    Given a pose and two chains, returns the index of serine present in the RRXS motif
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose): Pose object to search for motif in
        chain_1 (int): The first chain to begin the search (A is 1, chain B is 2, etc.)
        chain_2 (int): The chain to end the search (A is 1, chain B is 2, etc.)
        
    Returns:
        index (int): The index of the serine corresponding to the RRXS motif.
    """

    chain_length = get_chain_length(chain_1, chain_2)

    for index in range(4, chain_length):
        if pose.residue(index).name()[0:3] == "SER":
            if pose.residue(index - 2).name()[0:3] == "ARG" and pose.residue(index - 3).name()[0:3] == "ARG":
                return index
            

def get_energy(pose:pyrosetta.rosetta.core.pose.Pose) -> float:
    """
    Calculates the energy of a pose using the default PyRosetta energy function described in
    O'Meara et al., 2015
    
    Args:
        pose (pyrosetta.rosetta.core.pose.Pose): The pose to be scored
        
    Returns:
        energy (float): The score"""
    
    #using an all atom score function
    score_fn = get_score_function(True)
    energy = score_fn(pose) #scoring the pose
    return energy


def get_phosphorlyation_energy(pdb_string1):
    with open("temp.pdb", "w") as f:
        f.write(pdb_string1)

    pose = pose_from_file("temp.pdb")
    serine_id = find_serine(pose, 1, 1)

    #saving the pose in its current form before mutating
    original_pose = Pose()
    original_pose.assign(pose)

    mutate_residue(pose, serine_id, "D") #mutating serine residue to aspartate to mimic phosphorylation

    starting_energy = get_energy(original_pose)
    energy_after_mut = get_energy(pose)
    delta_energy = energy_after_mut - starting_energy

    return delta_energy, starting_energy
