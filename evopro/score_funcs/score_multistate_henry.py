from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd
import math

from Bio import SeqIO

import pandas as pd
import numpy as np

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
    return rmsd_to_starting * 5, (rmsd_to_starting, rmsd_potential * 5)

if __name__=="__main__":
    print("no main functionality")
