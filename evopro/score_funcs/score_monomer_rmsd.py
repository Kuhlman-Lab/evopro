import math
import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd

def score_pairs_aa(sequence, pair_aas):
    num_pairs = []
    for aa in pair_aas:
        aa_pair = aa
        aa_pair = aa_pair + aa
        #print(aa_pair)
        num_pairs = num_pairs + [i for i in range(len(sequence)) if sequence.startswith(aa_pair, i)]
    
    return len(num_pairs)
    

def score_seq_diff(results, dsobj, contacts=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results

def score_seq_diff_2(results, dsobj, contacts=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
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
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
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
        rmsdscore = threshold_rmsd(pdb, rmsd_pdb)
    else:
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]

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
        rmsd_pot_score = threshold_rmsd(pdb, rmsd_pdb)
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore + rmsd_pot_score
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
        rmsd_pot_score = threshold_rmsd(pdb, rmsd_pdb)
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
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
        rmsd_pot_score = threshold_rmsd_small(pdb, rmsd_pdb)
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
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
        rmsd_pot_score = threshold_rmsd(pdb, rmsd_pdb)
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb)
    else:
        rmsd_pot_score = 0
        rmsdscore=0
    overall_score = -confscore/2 + rmsdscore + rmsd_pot_score
    score = [(confscore,),  (rmsdscore, rmsd_pot_score)]

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
        rmsdscore = threshold_rmsd_small(pdb, rmsd_pdb)
    else:
        rmsdscore=0
    overall_score = -confscore/5 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]

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
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=ca_only, dsobj=dsobj)

    return rmsd_to_starting


def threshold_rmsd(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    #rmsd_cutoff = 5
    rmsd_cutoff = 10

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]


    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of 5.0. Calculating RMSD potential")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return rmsd_potential*5

def threshold_rmsd_small(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0
    rmsd_cutoff = 5
    #rmsd_cutoff = 10

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]


    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of 5.0. Calculating RMSD potential")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return rmsd_potential*5


if __name__ == "__main__":
    sequence = "RVLVIIIVEEEEIRERVIEVVEEREVIVIEIEEVEEEVVRILEELEIEVRIVVILIREEERERVERLIIELIEEEELEVVIVRVREVEVVERLVEEILER"
    pair_aas = ["E", "R"]
    print(score_pairs_aa(sequence, pair_aas))
