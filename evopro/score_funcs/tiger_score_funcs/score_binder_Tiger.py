# Tiger's modified scorefunctions for CD70 system
# changed de novo binder chain to D for CD70 trimer (chains A/B/C)

from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import os
import subprocess
import shutil
import math

path_to_3helixbundle = '/nas/longleaf/home/tigerz/kuhlmanlab/structures/5djt_3-helix-bundle/3helix_bundle_cleaned_better.pdb'

def score_binder(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex(results, dsobj, contacts, orient=orient)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_complex(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    #print("Binder contacts:")
    #print(contacts)
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0
    score = -contactscore + orientation_penalty
    return score, (score, len(contacts), contactscore, orientation_penalty), contacts, pdb, results

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    return score, (score, confscore2), pdb, results

def score_binder_rmsd(pdb1, pdb2, binder_chain="D", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    #print(reslist1)
    #print(reslist2)
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    #print("RMSD: " + str(rmsd_binder))
    return rmsd_binder*5

def score_binder_old(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    if len(chains)>1:
        return score_binder_complex_old(results, dsobj, contacts)
    else:
        return score_binder_monomer(results, dsobj)

def score_binder_complex_old(results, dsobj, contacts, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)
    pae_score = score_pae_confidence_pairs(results, contacts, resindices, dsobj=dsobj)
    pae_score = pae_score/(len(contacts)*10)
    if orient:
        orientation_penalty = orientation_score(orient, pdb, dsobj=dsobj)
    else:
        orientation_penalty = 0
    score = -contactscore + pae_score + orientation_penalty
    return score, (score, contactscore, pae_score, orientation_penalty), contacts, pdb, results

def score_binder_rmsd_restraint(pdb1, pdb2, binder_chain="D",  dsobj=None):
    # for CD70 trimer system
    # to keep the de novo binder close in structure to the original binder
    # calculates    1) Ca-only RMSD of de novo binder unbound vs bound to CD70 complex
    #               2) Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential
    # total RMSD score = 1) + 2)
    spring_constant = 10.0 
    rmsd_cutoff = 5.0

    # calculate 1)
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1) # complex
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2) # de novo binder alone
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    #print("RMSD: " + str(rmsd_binder))

    # calculate 2)
    with open(path_to_3helixbundle, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("RMSD scores, RMSD_bound: %.3f, rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_binder, rmsd_to_starting, rmsd_potential) )

    return rmsd_binder*5 + rmsd_potential

def score_binder_rmsd_restr(pdb1, pdb2, spring_constant, rmsd_cutoff, binder_chain="D", dsobj=None):
    # for CD70 trimer system, to keep the de novo binder close in structure to the original binder
    # calculates    1) Ca-only RMSD of de novo binder unbound vs bound to CD70 complex
    #               2) Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential
    # total RMSD score = 1) + 2)

    # calculate 1)
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1) # complex
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2) # de novo binder alone
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True, dsobj=dsobj)

    # calculate 2)
    with open(path_to_3helixbundle, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist2, pdb2, ca_only=True, dsobj=dsobj)
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:     # apply flat-bottom quadratic-shaped potential function
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("RMSD scores, RMSD_bound: %.3f, rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_binder, rmsd_to_starting, rmsd_potential) )

    return rmsd_binder*5 + rmsd_potential

def score_binder_rmsd_to_starting(pdb, path_to_starting, spring_constant, rmsd_cutoff, binder_chain="D", dsobj=None):
    # for CD70 trimer system, to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original 3hb applied to a flat-bottom quadratic potential

    # import de novo binder alone??
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb) # de novo binder alone
    print(chains1)
    reslist1 = [x for x in residues1.keys()]

    # import starting pdb structure
    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains_0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist0 = [x for x in residues0.keys()]

    rmsd_to_starting = get_rmsd(reslist0, pdb_string_starting, reslist1, pdb, ca_only=True, dsobj=dsobj)
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:     # apply flat-bottom quadratic-shaped potential function
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    print("rmsd_to_starting: %.3f, rmsd_potential: %.3f" %(rmsd_to_starting, rmsd_potential) )

    return rmsd_potential

def score_binder_rmsd_start_5_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 5, 3)

def score_binder_rmsd_start_5_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 5, 5)

def score_binder_rmsd_start_5_7(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 5
    # rmsd cutoff = 7
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 5, 7)    

def score_binder_rmsd_start_10_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 10
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 10, 3) 

def score_binder_rmsd_start_10_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 10
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 10, 5)  

def score_binder_rmsd_start_20_3(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 20
    # rmsd cutoff = 3
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 20, 3) 

def score_binder_rmsd_start_20_5(pdb, path_to_starting, binder_chain="D", dsobj=None):
    # spring constant = 20
    # rmsd cutoff = 5
    return score_binder_rmsd_to_starting(pdb, path_to_3helixbundle, 20, 5)     

if __name__=="__main__":
    print("no main functionality")
