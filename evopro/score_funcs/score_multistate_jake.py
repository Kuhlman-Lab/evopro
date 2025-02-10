from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient
"""
import sys
sys.path.append("/proj/kuhl_lab/ThermoMPNN")
from analysis.thermompnn_benchmarking import compute_centrality
"""
import math, time
from Bio import SeqIO
from Bio.PDB import PDBParser

import pandas as pd
import numpy as np
import math 
import sys
import os
import shutil
import torch

CODON_TABLE = pd.DataFrame({
    'AA': np.array(['*', '*', '*', 'A', 'A', 'A', 'A', 'C', 'C', 'D', 'D', 'E', 'E',
       'F', 'F', 'G', 'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'K', 'K',
       'L', 'L', 'L', 'L', 'L', 'L', 'M', 'N', 'N', 'P', 'P', 'P', 'P',
       'Q', 'Q', 'R', 'R', 'R', 'R', 'R', 'R', 'S', 'S', 'S', 'S', 'S',
       'S', 'T', 'T', 'T', 'T', 'V', 'V', 'V', 'V', 'W', 'Y', 'Y']), 
    'Codon': np.array(['TAA', 'TAG', 'TGA', 'GCA', 'GCC', 'GCG', 'GCT', 'TGC', 'TGT',
       'GAC', 'GAT', 'GAA', 'GAG', 'TTC', 'TTT', 'GGA', 'GGC', 'GGG',
       'GGT', 'CAC', 'CAT', 'ATA', 'ATC', 'ATT', 'AAA', 'AAG', 'CTA',
       'CTC', 'CTG', 'CTT', 'TTA', 'TTG', 'ATG', 'AAC', 'AAT', 'CCA',
       'CCC', 'CCG', 'CCT', 'CAA', 'CAG', 'AGA', 'AGG', 'CGA', 'CGC',
       'CGG', 'CGT', 'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT', 'ACA',
       'ACC', 'ACG', 'ACT', 'GTA', 'GTC', 'GTG', 'GTT', 'TGG', 'TAC',
       'TAT']), 
    'Frequency': np.array([0.3 , 0.24, 0.47, 0.23, 0.4 , 0.11, 0.27, 0.54, 0.46, 0.54, 0.46,
       0.42, 0.58, 0.54, 0.46, 0.25, 0.34, 0.25, 0.16, 0.58, 0.42, 0.17,
       0.47, 0.36, 0.43, 0.57, 0.07, 0.2 , 0.4 , 0.13, 0.08, 0.13, 1.  ,
       0.53, 0.47, 0.28, 0.32, 0.11, 0.29, 0.27, 0.73, 0.21, 0.21, 0.11,
       0.18, 0.2 , 0.08, 0.24, 0.15, 0.15, 0.22, 0.05, 0.19, 0.28, 0.36,
       0.11, 0.25, 0.12, 0.24, 0.46, 0.18, 1.  , 0.56, 0.44])
})

BIDIR_MATRIX = np.array([[0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1.,
        0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1.,
        1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.],
       [0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.],
       [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
        0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1.,
        1., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
        1., 0., 0., 0., 0.],
       [0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1.,
        0., 0., 0., 0., 0.],
       [0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
        0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 1., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0.]])

ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'

def score_bidir_seqs(seq1, seq2):
    """Calculate % of bidirectionally compliant residues in two sequences"""
    assert len(seq1) == len(seq2)
    score = 0
    score_all = []
    for pos_a, pos_b in zip(seq1, seq2):
        idx_a, idx_b = ALPHABET.index(pos_a), ALPHABET.index(pos_b)

        # load compatibility table
        # print(pos_a, pos_b, '\tReverse complementary:\t', bool(BIDIR_MATRIX[idx_a, idx_b]))
        score += BIDIR_MATRIX[idx_a, idx_b]
        score_all.append(BIDIR_MATRIX[idx_a, idx_b])

    score = score / len(seq1)
    # print('BIDIRECTIONALITY SCORE (percent of codons matched):', round(score * 100, 2), '%')
    return score

def score_overall(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer(result, dsobj))
    
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_complex(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist1 = contacts
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False)

    score = -contactscore
    print(score, (score, len(contacts), contactscore))
    return score, (score, len(contacts), contactscore)

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -confscore2/10
    print(score)
    return score, (score, confscore2)


def score_overall_bidir(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer(result, dsobj))
    
    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_overall_bidir_gyration(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score = []
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer(result, dsobj))

    for e in pdbs:
        score.append(((-100 * (1/Rg(e))), ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_overall_bidir_gyration_mult(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score = []
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        score.append(((-100 * (1/Rg(e))), ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = 1
    for a in score:
        overall_score = a[0] * overall_score
    return overall_score, score, pdbs, results

def score_pLDDT_bidir_gyration(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score = []
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        score.append(((-100 * (1/Rg(e))), ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_bidirectional(dsobj):
    """Calculates bidirectional coding match % for scoring purposes"""
    seq1, seq2 = dsobj.get_sequence_string().split(',')
    # check reverse complements
    raw_score = score_bidir_seqs(seq1, seq2[::-1])
    score = -raw_score * 10  # scale to between -10 (best) and 0 (worst)
    return score, (score, raw_score)

def score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    if len(reslist1) == len(reslist2):
        rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    else:
        print("Trying to calculate RMSD between proteins of different length. Setting RMSD to 0.")
        rmsd_binder = 0

    return -rmsd_binder

def score_monomer_combined(results, path_to_starting, dsobj):
    """Calculates both RMSD and pLDDT and combines them"""
    from alphafold.common import protein
    # af2 function
    pdbs, scores = [], []
    for result in results:
        # convert dictionary of protein info to scoring-compatible PDB format
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        # do RMSD calculation
        rmsd_scores = score_binder_rmsd_to_starting(pdb, path_to_starting)
        scores.append(rmsd_scores)
        # do pLDDT calculation
        pLDDT_scores = score_binder_monomer(result, dsobj)
        scores.append(pLDDT_scores)

    overall_score = sum(x[0] for x in scores)
    return overall_score, scores, pdbs, results

def score_rmsd_multi_wrapper(results, path_to_starting, dsobj):
    """Run RMSD-to-starting for multiple states"""
    import os
    from alphafold.common import protein

    # need a directory passed as starting structure
    assert os.path.isdir(path_to_starting)
    pdbs_available = sorted([s for s in os.listdir(path_to_starting) if s.endswith('.pdb')])
    # NOTE that pdbs available must be in alphabetical order that matches chain order
    # need same # of reference PDBs and predicted PDBs to compare
    assert len(pdbs_available) == len(results)
    
    score, pdbs = [], []
    for ref, pred in zip(pdbs_available, results):
        pdb = protein.to_pdb(pred['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)

        # currently only does one chain calculation
        if len(chains) > 1:
            raise ValueError('RMSD Wrapper not configured for >1 chain')
        else:
            print(path_to_starting)
            score.append(score_binder_rmsd_to_starting(pdb, os.path.join(path_to_starting, ref)))

    # take total score of all states combined
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results
    
    return

def score_rmsd_wrapper(results, path_to_starting, dsobj):
    """Use this in EvoPro flags file"""
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        
        # currently only does one chain calculation
        if len(chains) > 1:
            raise ValueError('RMSD Wrapper not configured for >1 chain')
        else:
            print(path_to_starting)
            score.append(score_binder_rmsd_to_starting(pdb, path_to_starting))
    
    # overall score is currently just RMSD term here
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_rmsd_to_starting(pdb, path_to_starting, dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 4.0

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
    return rmsd_to_starting, (rmsd_to_starting, rmsd_potential)

if __name__=="__main__":
    print("no main functionality")


def score_rmsd_multi_bidir_pLDDT(results, path_to_starting, dsobj, contacts = None):
    """Run RMSD-to-starting for multiple states"""
    import os
    from alphafold.common import protein

    # need a directory passed as starting structure
    assert os.path.isdir(path_to_starting)
    pdbs_available = sorted([s for s in os.listdir(path_to_starting) if s.endswith('.pdb')])
    # NOTE that pdbs available must be in alphabetical order that matches chain order
    # need same # of reference PDBs and predicted PDBs to compare
    assert len(pdbs_available) == len(results)
    
    score, pdbs = [], []
    for ref, pred in zip(pdbs_available, results):
        pdb = protein.to_pdb(pred['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)

        # currently only does one chain calculation
        if len(chains) > 1:
            raise ValueError('RMSD Wrapper not configured for >1 chain')
        else:
            print(path_to_starting)
            score.append(score_binder_rmsd_to_starting(pdb, os.path.join(path_to_starting, ref)))
    
    #pLDDT Score Calculation
    for result in results:
            pdb = protein.to_pdb(result['unrelaxed_protein'])
            pdbs.append(pdb)
            chains, residues, resindices = get_coordinates_pdb(pdb)
            print(chains, '***')
            if len(chains)>1:
                score.append(score_binder_complex(result, dsobj, contacts))
            else:
                score.append(score_binder_monomer(result, dsobj))

    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)

    #bidirectional score calculation
    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))

    
    # take total score of all states combined
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_pLDDT_bidir_gyration_gradient(results, dsobj, curr_iter, num_iters, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score = []
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        score.append(((-100 * (1/Rg(e))), ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional_Weighted(dsobj, curr_iter, num_iters))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_pLDDT_bidir_gyration_sqrt(results, dsobj, curr_iter, num_iters, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score = []
    pdbs = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        score.append(((-100 * (1/Rg(e))), ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional_sqrt(dsobj, curr_iter, num_iters))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_pLDDT_bidir_gyration_cutoff(results, dsobj, contacts=None):
    """
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    spring_const = 10
    rg_cutoff = 10
    score = []
    pdbs = []

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 12. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

#Score Calculators Jake

def score_bidirectional_Weighted(dsobj, curr_iter, num_iters):
    """Calculates bidirectional coding match % for scoring purposes"""
    seq1, seq2 = dsobj.get_sequence_string().split(',')
    # check reverse complements
    raw_score = score_bidir_seqs(seq1, seq2[::-1])
    if (curr_iter*2)/num_iters >= 1:
        score = -raw_score * 10  # scale to between -10 (best) and 0 (worst)
    else:
        score = -raw_score * 10 * ((curr_iter*2)/num_iters) # scale to between -10 (best) and 0 (worst) on a gradient
    return score, (score, raw_score)

def score_bidirectional_sqrt(dsobj, curr_iter, num_iters):
    """Calculates bidirectional coding match % for scoring purposes"""
    seq1, seq2 = dsobj.get_sequence_string().split(',')
    # check reverse complements
    raw_score = score_bidir_seqs(seq1, seq2[::-1])
    score = -raw_score * 10 * (math.sqrt(curr_iter/num_iters)) # scale to between -10 (best) and 0 (worst) on a gradient
    return score, (score, raw_score)

def score_binder_monomer_heavy(results, dsobj):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    score = -1.2 * confscore2/10
    print(score)
    return score, (score, confscore2)

def compute_centrality(xyz, basis_atom: str = "CA", radius: float = 10.0, core_threshold: int = 20, surface_threshold: int = 15, backup_atom: str = "C", chain: str = 'A') -> torch.Tensor:

    coords = xyz[basis_atom + f'_chain_{chain}']
    coords = torch.tensor(coords)
    # Compute distances and number of neighbors.
    pairwise_dists = torch.cdist(coords, coords)
    pairwise_dists = torch.nan_to_num(pairwise_dists, nan=2 * radius)
    num_neighbors = torch.sum(pairwise_dists < radius, dim=-1) - 1
    return num_neighbors

#RosettaFold2 Score Functions


def score_overall_rf2(results, dsobj, contacts=None, distance_cutoffs=None):
    
    print("Number of predictions being scored:", len(results))
    rg_cutoff = 12
    spring_const = 10
    score=[]
    pdbs = []
    for result in results:

        pdb = result['pdb']
        pdbs.append(pdb)
        chains, _, _ = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex_rf2(result, dsobj, contacts, distance_cutoffs))
            score.append(score_complex_confidence_rf2(result, dsobj))
        else:
            score.append(score_binder_monomer_rf2(result, dsobj))
    
    """
    pdb1 = results[0]['pdb']
    pdb2 = results[1]['pdb']
    score.append(score_binder_rmsd_rf2(pdb1, pdb2, binder_chain="B", dsobj=None))
    score.append(threshold_rmsd_rf2(pdb1, pdb2, binder_chain="B", dsobj=None))
    """

    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 12. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
    
    #bidirectional score calculation
    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))

    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_monomer_rf2(results, dsobj):
    spring_constant = 10.0
    plddt_cutoff = 70.0
    pdb = pdb = results['pdb']
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of monomer is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)

    score = -1.1*confscore2/10 + plddt_potential
    print(score)
    return score, (score, confscore2)

def score_binder_complex_rf2(results, dsobj, contacts, distance_cutoffs):
    start = time.time()
    pdb = results['pdb']
    _, residues, _ = get_coordinates_pdb(pdb)
    print(f"contacts === {contacts}")

    if not contacts:
        contacts=(None,None,None)
    if not distance_cutoffs:
        distance_cutoffs=(4,4,8)
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("D")]
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[0])
    
    bonuses = 0
    bonus_resids = contacts[1]
    if bonus_resids:
        bonus_contacts, bonus_contactscore = score_contacts_pae_weighted_efficient(results, pdb, bonus_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[1])
        for contact in bonus_contacts:
            if (contact[0][0:1] == 'A' or contact[0][0:1] == 'B' or contact[0][0:1] == 'C') and int(contact[0][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[0]))
            if (contact[1][0:1] == 'A' or contact[1][0:1] == 'B' or contact[1][0:1] == 'C') and int(contact[1][1:]) in bonus_resids:
                bonuses += 1
                print("bonus found at: " + str(contact[1]))
        
    bonus = -bonuses * 3
    
    penalties = 0
    penalty_resids = contacts[2]
    print(penalty_resids)
    if penalty_resids:
        penalty_contacts, penalty_contactscore = score_contacts_pae_weighted_efficient(results, pdb, penalty_resids, reslist2, dsobj=dsobj, first_only=False, dist=distance_cutoffs[2])
        for contact in penalty_contacts:
            if contact[0][0:1] == 'A' and int(contact[0][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[0]))
            if contact[1][0:1] == 'A' and int(contact[1][1:]) in penalty_resids:
                penalties += 1
                print("penalty found at: " + str(contact[1]))
        
    penalty = penalties * 3
    
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
    
    score = -contactscore + penalty + bonus
    print("Time to score:", time.time()-start)
    print(score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty))
    return score, (score, len(contact_list), contactscore, pae_per_contact, bonus, penalty), contacts, pdb, results

def score_complex_confidence_rf2(results, dsobj):
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = pdb = results['pdb']
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of complex is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = -confscore2*2 + plddt_potential
    print(score)
    return score, (score, confscore2)

def threshold_rmsd_rf2(pdb1, pdb2, binder_chain="B", dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 10

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]


    rmsd_to_starting = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of 5.0. Calculating RMSD potnetial")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return [rmsd_potential*5]

def score_binder_rmsd_rf2(pdb1, pdb2, binder_chain="B", dsobj=None):
    _, residues1, _ = get_coordinates_pdb(pdb1)
    _, residues2, _ = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    if len(reslist1) == len(reslist2):
        rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    else:
        print("Trying to calculate RMSD between proteins of different length. Setting RMSD to 0.")
        rmsd_binder = 0

    return [-rmsd_binder*20]

#Current Functions

def score_23_12_18(results, dsobj, contacts=None):
    """
    Composed of weighted pLDDT, ROG with Cutoff and bidirectionality
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    spring_const = 5
    rg_cutoff = 9
    score = []
    pdbs = []

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 10. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
        #score.append((rge, ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_24_05_10_longer(results, dsobj, contacts=None):
    """
    Composed of weighted pLDDT, ROG with Cutoff and bidirectionality
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    spring_const = 5
    rg_cutoff = 9
    score = []
    pdbs = []

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 10. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
        #score.append((rge, ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_neighbors(results, dsobj, contacts=None):
    """
    Composed of weighted pLDDT, number of neighbors and bidirectionality
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    import os
    print("Number of predictions being scored:", len(results))

    #spring_const = 10
    #rg_cutoff = 9
    score = []
    pdbs = []
    s = []

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        with open(os.path.join(os.getcwd(), "temp.pdb"), "w") as pdbf:
            pdbf.write(str(pdb))
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))
    
    temp_pdb = PDBParser().get_structure("temp", os.path.join(os.getcwd(), "temp.pdb"))
    #calculating number of neighbors 
    for model in temp_pdb:
        cords = []
        for chain in model:
            for residue in chain:
                cords.append(residue['CA'].get_vector())
        with open(os.path.join(os.getcwd(), "cords.txt"), "w") as cord:
            cord.write(str(cords))
        score.append((compute_centrality(cords),))

    """
    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 10. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
        #score.append((rge, ))
    """
    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results



def score_neighbors_diff(results, dsobj, contacts=None):
    """
    Composed of weighted pLDDT, number of neighbors and bidirectionality
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    spring_const = 10
    rg_cutoff = 9
    score = []
    pdbs = []
    cords = []
    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    #calculating number of neighbors
    for d in pdbs:
        s = PDBParser.get_structure("a", d)      
    
    for model in s:
        for chain in model:
            for residue in chain:      
                for atom in residue:
                    cords.append(atom.get_vector())

    neighs = compute_centrality(cords)
    for i in (len(neighs)/2) - 1:
        l = len(neighs)-i-1
        diff += abs(neighs[i]-neighs[l])
    avg_diff = diff/(len(neighs)/2)
    score.append((avg_diff,))

    #Radius of Gyration Score
    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 10. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
        #score.append((rge, ))

    #Bidirectional Score
    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results



def score_len90(results, dsobj, contacts=None):
    """
    Composed of weighted pLDDT, ROG with Cutoff and bidirectionality
    General purpose fxn for handling both complex and monomer scoring based on # of chains
    Also handles MULTI STATE predictions by taking the sum of each state score (each "result" score)
    Also calculates a bidirectional encoding score (% of residues paired) and adds this to the score
    Also calculates a radius of gyration score to prevent long helicies from being preffered structures
    """
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    spring_const = 5
    rg_cutoff = 11
    score = []
    pdbs = []

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        print(chains, '***')
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
        else:
            score.append(score_binder_monomer_heavy(result, dsobj))

    for e in pdbs:
        rge = Rg(e)

        if rge > rg_cutoff:
            rog_potential = spring_const*math.pow(rge - rg_cutoff, 2)
            print("ROG has reached the threshold of 11. Calculating ROG potnetial")
            print(rog_potential)
        else:
            rog_potential = 0
        score.append((rog_potential, ))
        #score.append((rge, ))

    if len(pdbs) != 2:
        raise AssertionError("Cannot do bidirectional scoring for >2 PDBs at once!")
    else:
        score.append(score_bidirectional(dsobj))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results
