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

def tim_barrel_score1(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 1*length_score

    return score, (conf_score, -structure_score, length_score), pdb, results

def tim_barrel_score2(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 0.5*conf_score - 1*structure_score + 1*length_score

    return score, (0.5*conf_score, -structure_score, length_score), pdb, results

def tim_barrel_score3(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 0.5*structure_score + 1*length_score

    return score, (conf_score, -0.5*structure_score, length_score), pdb, results

def tim_barrel_score4(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 0.5*length_score

    return score, (conf_score, -structure_score, 0.5*length_score), pdb, results

def tim_barrel_score5(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 0.5*structure_score + 0.5*length_score

    return score, (conf_score, -0.5*structure_score, 0.5*length_score), pdb, results

def tim_barrel_score6(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 1*length_score

    return score, (0.5*conf_score, -0.5*structure_score, length_score), pdb, results

def tim_barrel_score7(results, orient=None):

    target = '/pine/scr/s/m/smyersn/folddesign/test_scaffolds/test_scaffolds/ferr_ems_00120/ferr_ems_00120.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 0.5*conf_score - 1*structure_score + 0.5*length_score

    return score, (conf_score, -structure_score, length_score), pdb, results

#above are score functions for test scaffolds.  Below are score functions for tim barrel.

def tim_barrel_score8(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 1*length_score

    return score, (conf_score, -structure_score, length_score), pdb, results

def tim_barrel_score9(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 0.5*conf_score - 1*structure_score + 1*length_score

    return score, (0.5*conf_score, -structure_score, length_score), pdb, results

def tim_barrel_score10(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 0.5*structure_score + 1*length_score

    return score, (conf_score, -0.5*structure_score, length_score), pdb, results

def tim_barrel_score11(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 0.5*length_score

    return score, (conf_score, -structure_score, 0.5*length_score), pdb, results

def tim_barrel_score12(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 0.5*structure_score + 0.5*length_score

    return score, (conf_score, -0.5*structure_score, 0.5*length_score), pdb, results

def tim_barrel_score13(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 1.0*conf_score - 1*structure_score + 1*length_score

    return score, (0.5*conf_score, -0.5*structure_score, length_score), pdb, results

def tim_barrel_score14(results, orient=None):

    target = '/pine/scr/s/m/smyersn/proteinmpnn/tim_barrel/t_1/3d5h.pdb'
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist2, resindices, fil = False)
    conf_score = -confscore2/len(reslist2)
    structure_score, length_score = score_sec_struc(target, pdb)
    score = 0.5*conf_score - 1*structure_score + 0.5*length_score

    return score, (conf_score, -structure_score, length_score), pdb, results
