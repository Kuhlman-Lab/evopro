from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import get_seq_indices, score_contacts, score_contacts_pae_weighted, score_pae_confidence_pairs, score_pae_confidence_lists, score_plddt_confidence, get_rmsd, orientation_score
import math

def score_overall_1(results, dsobj, contacts=None):
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

        score.append(score_binder_1_d2(result, dsobj, contacts))

    
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results
        
def score_overall_2(results, dsobj, contacts=None):
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

        score.append(score_binder_1(result, dsobj, contacts))

    
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_1(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    #path_to_starting = "/work/users/a/m/amritan/evopro_tests/vary_length/d1truncated.pdb"
    path_to_starting = "/nas/longleaf/home/sela/runevopro/d1/tlvary/d1truncated.pdb"
    rmsd_score = score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=dsobj, resids=contacts)
    print(confscore2, rmsd_score)
    score = -confscore2 + rmsd_score
    print(score)
    return score, (score, confscore2, rmsd_score), pdb, results

def score_binder_2(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    path_to_starting = "/nas/longleaf/home/sela/runevopro/d1/tlvary/d1truncated.pdb"
    rmsd_score = score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=dsobj, resids=contacts)
    print(confscore2, rmsd_score)
    score = -confscore2/10 + rmsd_score
    print(score)
    return score, (score, confscore2, rmsd_score), pdb, results

def score_binder_3(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    path_to_starting = "/nas/longleaf/home/sela/runevopro/d1/tlvary/d1truncated.pdb"
    rmsd_score = score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=dsobj, resids=contacts)
    print(confscore2, rmsd_score)
    score = -confscore2/10 + rmsd_score*2
    print(score)
    return score, (score, confscore2, rmsd_score), pdb, results

def score_binder_1_d2(results, dsobj, contacts=None, orient=None):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    path_to_starting = "/nas/longleaf/home/sela/runevopro/d2/tla/d2af.pdb"
    rmsd_score = score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=dsobj, resids=contacts)
    print(confscore2, rmsd_score)
    score = -confscore2 + rmsd_score
    print(score)
    return score, (score, confscore2, rmsd_score), pdb, results

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None, resids=None):
    # to keep the designed protein close in structure to the original protein
    # calculates Ca-only RMSD of designed protein vs to original protein applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x in resids]

   #  with open("/work/users/a/m/amritan/evopro_tests/vary_length/temp.pdb", 'w') as f:
      #  f.write(pdb)
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys() if x in resids]
    
    print("Residue lists before renumbering:", reslist1, reslist2)
    """if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=True)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=True)

    print("Residue lists after renumbering:", reslist1, reslist2)"""
    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, translate=True)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)

    return rmsd_potential


if __name__=="__main__":
    print("no main functionality")
