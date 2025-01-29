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
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=5, spring_constant=10.0)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

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
        rmsdscore, rmsd_pot_score = score_rmsd_to_starting_threshold(pdb, rmsd_pdb, threshold=True, rmsd_cutoff=1.5, spring_constant=10.0)
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
    overall_score = -confscore + rmsdscore + rmsd_pot_score
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


def score_rmsd_to_starting_threshold(pdb, path_to_starting, ca_only=True, dsobj=None, threshold=False, rmsd_cutoff=5, spring_constant=10.0):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting_temp = f.read()
        
    pdb_string_starting = renumber_chains_clockwise_test(pdb_string_starting_temp)
    
    _, residues0, _ = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    pdb_renum = renumber_chains_clockwise_test(pdb)
    _, residues2, _ = get_coordinates_pdb(pdb_renum)
    reslist2 = [x for x in residues2.keys()]
    
    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb_renum, ca_only=ca_only, dsobj=dsobj)
    print("RMSD", rmsd_to_starting)
    
    rmsd_potential = 0
    if threshold:
        if rmsd_to_starting > rmsd_cutoff:
            # apply flat-bottom quadratic-shaped potential function
            rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
            print("Rmsd has reached the threshold of " + str(rmsd_cutoff) + ". Calculating RMSD potential:")
            print(rmsd_potential)

    return rmsd_to_starting, rmsd_potential
