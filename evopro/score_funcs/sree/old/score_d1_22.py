from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import get_seq_indices, score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import math
import time

def score_overall_d1_22(results, dsobj, distance_cutoffs=None, contacts=None, rmsd_pdb=None):
    start = time.time()
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

        score.append(score_binder_d1_22_upweight_linker(result, dsobj, contacts=contacts))


    overall_score = sum([x[0] for x in score])
    score_all = [x[1] for x in score]
    print("Time to score: ", time.time()-start)
    print(overall_score, score_all)
    return overall_score, score_all, pdbs, results

def score_binder_d1_22_upweight_linker(results, dsobj, distance_cutoffs=None, rmsd_pdb=None, contacts=None, orient=None):
    print("contacts",contacts)
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    reslist_designable = [x for x in residues.keys() if x in dsobj._get_designable_positions()]
    reslist_not_designable = [x for x in residues.keys() if x not in dsobj._get_designable_positions()]

    #confscore_designable = score_plddt_confidence(results, reslist_designable, resindices, dsobj=dsobj, first_only=False)
    #confscore_nondesignable = score_plddt_confidence(results, reslist_not_designable, resindices, dsobj=dsobj, first_only=False)
    confscore = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)

    path_to_starting = "/work/users/s/e/sela/fall2023/dengueproj/new_evopro/d1_22/1uzg.pdb"
    rmsd_score = score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=dsobj, resids_orig=contacts[0], resids_new=reslist_not_designable)
    
    #score = -confscore_designable*2 - confscore_nondesignable + rmsd_score*3
    score = -1*confscore + rmsd_score

    print("Overall score", score)
    print("PLDDT",confscore)
    print("RMSD",rmsd_score)
    return score, [confscore,rmsd_score], pdb, results

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None, resids_orig=None, resids_new=None):
    # to keep the designed protein close in structure to the original protein
    # calculates Ca-only RMSD of designed protein vs to original protein applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 10.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x in resids_orig]

   #  with open("/work/users/a/m/amritan/evopro_tests/vary_length/temp.pdb", 'w') as f:
      #  f.write(pdb)
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys() if x in resids_new]
    
    print("Residue lists for RMSD:", reslist1, reslist2)

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, translate=True)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
    print("RMSD to starting ",rmsd_to_starting)
    return rmsd_potential


if __name__=="__main__":
    print("no main functionality")
