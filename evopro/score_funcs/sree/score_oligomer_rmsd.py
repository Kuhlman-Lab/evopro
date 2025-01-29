import math
import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd
from evopro.utils.renumber_chains_clockwise import renumber_chains_clockwise, renumber_chains_clockwise_test

def score_pairs_aa(sequence, pair_aas):
    num_pairs = []
    for aa in pair_aas:
        aa_pair = aa
        aa_pair = aa_pair + aa
        #print(aa_pair)
        num_pairs = num_pairs + [i for i in range(len(sequence)) if sequence.startswith(aa_pair, i)]
    
    return len(num_pairs)
    

def score_seq_diff(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results

def score_seq_diff_2(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results


def score_seq_diff_3(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=False)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]

    return overall_score, score, [pdb], results

def score_seq_diff_4(results, dsobj, contacts=None, distance_cutoffs=None,rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsd_pot_score
    score = [(confscore,),  (rmsd_pot_score,)]

    return overall_score, score, [pdb], results

def score_seq_diff_5(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsd_pot_score
    score = [(confscore,),  (rmsd_pot_score,)]

    return overall_score, score, [pdb], results

def score_seq_diff_6(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]
    print(score)

    return overall_score, score, [pdb], results

def score_seq_diff_7(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_8(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_9(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_10(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_11(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_12(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=10, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

    return overall_score, score, [pdb], results

def score_seq_diff_aa_penalty(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsd_pot_score = threshold_rmsd(pdb, rmsd_pdb)
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
        
    sequence = dsobj.get_sequence_string()
    pair_aas = ["E", "R"]
    num_pairs = score_pairs_aa(sequence, pair_aas)
    pair_aa_penalty = num_pairs*500
    
    overall_score = -confscore + rmsdscore + rmsd_pot_score + pair_aa_penalty
    score = [(confscore,), (rmsdscore, rmsd_pot_score), (pair_aa_penalty,num_pairs)]

    return overall_score, score, [pdb], results

def score_seq_diff_aa_penalty_2(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsd_pot_score = threshold_rmsd(pdb, rmsd_pdb)
        rmsdscore = 0
    else:
        rmsd_pot_score = 0
        rmsdscore=0

    sequence = dsobj.get_sequence_string()
    pair_aas = ["E", "R"]
    num_pairs = score_pairs_aa(sequence, pair_aas)
    pair_aa_penalty = num_pairs

    overall_score = -confscore + rmsdscore + rmsd_pot_score + pair_aa_penalty
    score = [(confscore,), (rmsdscore, rmsd_pot_score), (pair_aa_penalty,num_pairs)]

    return overall_score, score, [pdb], results


def score_rmsd_to_starting(pdb, path_to_starting, ca_only=True, dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting_temp = f.read()
        
    pdb_string_starting = renumber_chains_clockwise(pdb_string_starting_temp)
    
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    pdb_renum = renumber_chains_clockwise(pdb)
    _, residues2, _ = get_coordinates_pdb(pdb_renum)
    reslist3 = [x for x in residues2.keys()]

    rmsd_to_starting_renum = get_rmsd(reslist1, pdb_string_starting, reslist3, pdb_renum, ca_only=ca_only, dsobj=dsobj)

    return rmsd_to_starting_renum

def score_rmsd_to_starting_threshold(pdb, path_to_starting, ca_only=True, dsobj=None, threshold=False, rmsd_cutoff=5, spring_constant=10.0):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting_temp = f.read()
    print("pdb_string_starting_temp",pdb_string_starting_temp)    
    pdb_string_starting = renumber_chains_clockwise_test(pdb_string_starting_temp)
    
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    pdb_renum = renumber_chains_clockwise_test(pdb)
    _, residues2, _ = get_coordinates_pdb(pdb_renum)
    reslist2 = [x for x in residues2.keys()]
    
    rmsd_to_starting_upright = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb_renum, ca_only=ca_only, dsobj=dsobj)
    print("RMSD", rmsd_to_starting_upright)
    
    pdb_renum_flipped = renumber_chains_clockwise_test(pdb, invert = True)
    _, residues3, _ = get_coordinates_pdb(pdb_renum)
    reslist3 = [x for x in residues3.keys()]

    rmsd_to_starting_flipped = get_rmsd(reslist1, pdb_string_starting, reslist3, pdb_renum_flipped, ca_only=ca_only, dsobj=dsobj)
    print("RMSD inverted", rmsd_to_starting_flipped)
    
    rmsd_to_starting = min(rmsd_to_starting_upright, rmsd_to_starting_flipped)
    
    rmsd_potential = 0
    if threshold:
        if rmsd_to_starting > rmsd_cutoff:
            # apply flat-bottom quadratic-shaped potential function
            rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
            print("Rmsd has reached the threshold of " + str(rmsd_cutoff) + ". Calculating RMSD potential:")
            print(rmsd_potential)

    return rmsd_to_starting, rmsd_potential


def threshold_rmsd(pdb, path_to_starting, dsobj=None, ca_only=True):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    #rmsd_cutoff = 5
    rmsd_cutoff = 10

    with open(path_to_starting, 'r') as f:
        pdb_string_starting_temp = f.read()
        
    pdb_string_starting = renumber_chains_clockwise(pdb_string_starting_temp)
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    
    reslist1 = [x for x in residues0.keys()]
    
    pdb_renum = renumber_chains_clockwise(pdb)
    _, residues2, _ = get_coordinates_pdb(pdb_renum)
    reslist3 = [x for x in residues2.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist3, pdb_renum, ca_only=ca_only, dsobj=dsobj)


    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of " + str(rmsd_cutoff) + ". Calculating RMSD potential:")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return rmsd_potential*5

def threshold_rmsd_small(pdb, path_to_starting, dsobj=None, ca_only=True):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0
    rmsd_cutoff = 5
    #rmsd_cutoff = 10

    with open(path_to_starting, 'r') as f:
        pdb_string_starting_temp = f.read()
        
    pdb_string_starting = renumber_chains_clockwise(pdb_string_starting_temp)
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    
    reslist1 = [x for x in residues0.keys()]
    
    pdb_renum = renumber_chains_clockwise(pdb)
    _, residues2, _ = get_coordinates_pdb(pdb_renum)
    reslist3 = [x for x in residues2.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist3, pdb_renum, ca_only=ca_only, dsobj=dsobj)
    print("RMSD to starting", rmsd_to_starting)
    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of " + str(rmsd_cutoff) + ". Calculating RMSD potential:")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return rmsd_potential*5

if __name__ == "__main__":
    pdb1_name = "/work/users/a/m/amritan/evopro_tests/rmsd/design_manual_renum.pdb"
    pdb2_name = "/work/users/a/m/amritan/evopro_tests/rmsd/seq_0_final_model_1_chainABC.pdb"
    
    with open(pdb1_name, 'r') as f:
        pdb1 = f.read()
    
    with open(pdb2_name, 'r') as f:
        pdb2 = f.read()

    #print(threshold_rmsd_small_test(pdb2, pdb1_name, dsobj=None, ca_only=True))
