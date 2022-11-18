import sys
sys.path.append('/proj/kuhl_lab/alphafold/alphafold')
from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of
from folddesign.utils.pdb_parser import get_coordinates_pdb
import os
import subprocess
import shutil
#from alphafold.common import protein
from score_sec_struc_stride import score_sec_struc
    
def score_sec_struc_debug(path_to_target, pdb):    
    
    import subprocess
    from needlemanwunsch import nw
    import numpy as np
    import os

    target_formatted = stride_formatted(path_to_target)
    af2_formatted = stride_formatted(pdb)
    target_aligned, af2_aligned = align_structures(target_formatted, af2_formatted)
    target_gaps = insert_gaps(target_formatted, target_aligned)
    print(f'target aligned: {target_aligned}')
    print(f'af2 aligned: {af2_aligned}')
    af2_gaps = insert_gaps(af2_formatted, af2_aligned)
    structure_score, length_score = determine_score(target_gaps, af2_gaps)

    return structure_score, length_score

def tim_barrel_score1_debug(pdb, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(pdb, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc_debug(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 1*length_score

    return score, (conf_score, -structure_score, length_score), pdb, results

x, y = tim_barrel_score1_debug(/pine/scr/s/m/smyersn/folddesign/test_scaffolds/score_function4_testing/tim_barrel/fn1/seq_0_final_model_1_complex.pdb)
print(x, y)
