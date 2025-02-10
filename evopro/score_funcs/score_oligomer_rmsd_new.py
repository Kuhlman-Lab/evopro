import math
import sys, os
from biopandas.pdb import PandasPdb
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd
# from evopro.utils.renumber_chains_clockwise import renumber_chains_clockwise, renumber_chains_clockwise_test
from evopro.utils.multichain_permutation_alignment import align_structures_and_save
from evopro.utils.multichain_permutation_align_test import multichain_permutation_alignment

from alphafold.common import protein

def score_rmsd_test(pdb_string, path_to_starting, threshold=False, rmsd_cutoff=10, spring_constant=10.0):
    
    with open("temp.pdb", "w") as f:
        f.write(pdb_string)
    rmsd, mapping, aligned_pdb = multichain_permutation_alignment(
            "temp.pdb", path_to_starting)
    
    rmsd_potential = 0
    if threshold:
        if rmsd > rmsd_cutoff:
            # apply flat-bottom quadratic-shaped potential function
            rmsd_potential = spring_constant*math.pow(rmsd - rmsd_cutoff, 2)
            print("Rmsd has reached the threshold of " + str(rmsd_cutoff) + ". Calculating RMSD potential:")
            print(rmsd_potential)

    return rmsd, rmsd_potential


def score_pairs_aa(sequence, pair_aas):
    num_pairs = []
    for aa in pair_aas:
        aa_pair = aa
        aa_pair = aa_pair + aa
        #print(aa_pair)
        num_pairs = num_pairs + [i for i in range(len(sequence)) if sequence.startswith(aa_pair, i)]
    
    return len(num_pairs)
    

def score_seq_diff(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results

def score_seq_diff_2(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results


def score_seq_diff_3(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]

    return overall_score, score, [pdb], results

def score_seq_diff_4(results, dsobj, contacts=None, distance_cutoffs=None,rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsd_pot_score
    score = [(confscore,),  (rmsd_pot_score,)]

    return overall_score, score, [pdb], results

def score_seq_diff_5(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsd_pot_score
    score = [(confscore,),  (rmsd_pot_score,)]

    return overall_score, score, [pdb], results

def score_seq_diff_6(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]
    print(score)

    return overall_score, score, [pdb], results


def score_seq_diff_6new(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=4, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]
    print(score)

    return overall_score, score, [pdb], results


def score_seq_diff_7(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_8(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_9(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_9_2rmsd(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=2, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results


def score_seq_diff_10(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_11(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_12(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_13(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_test(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=8, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore/5 + rmsd_pot_score/5
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results


    
if __name__ == "__main__":
    pdb1_name = "/work/users/a/m/amritan/evopro_tests/rmsd/test/design_manual_renum.pdb"
    pdb2_name = "/work/users/a/m/amritan/evopro_tests/rmsd/test/seq_0_final_model_1_chainABC.pdb"
    
    with open(pdb1_name, 'r') as f:
        pdb1 = f.read()
    
    print(score_rmsd_test(pdb1, pdb2_name))
