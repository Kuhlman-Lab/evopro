from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.utils.write_pdb import PDBio
from folddesign.utils.calc_rmsd import RMSDcalculator
import math
import pickle
import numpy as np

def distance(p1, p2):
    """returns the distance between two 3D points represented as tuples"""

    dist = math.sqrt((float(p2[0])-float(p1[0]))**2+(float(p2[1])-(float(p1[1])))**2+(float(p2[2])-float(p1[2]))**2)
    return dist

def score_contacts(pdbfile, reslist1, reslist2, fil = True, orient=None, dist=4, of=False):
    if fil:
        chains, residues, resindices = get_coordinates_pdb(pdbfile, of=of)
    else:
        chains, residues, resindices = get_coordinates_pdb(pdbfile, fil = False, of=of)
    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        contact=1
                        pair = (res1, res2)
                        if pair not in pairs:
                            pairs.append(pair)
            score = score+contact
    if score>36:
        score=36

    if orient:
        helres = orient[0]
        pockres = orient[1]
        correct = False
        for res in pockres:
            for atom1 in residues[res]:
                for atom2 in residues[helres]:
                    if distance(atom1[2], atom2[2])<=20:
                        correct = True
        if not correct:
            score = score - 100

    return pairs, score

def score_confidence_pairs(resultsfile, reslist1, reslist2, resindices, fil = True):
    """calculates confidence score of all pairwise residue interactions"""
    score = 0
    if fil:
        with open(resultsfile,'rb') as f:
            p = pickle.load(f)
            pae = p['pae_output'][0]
            for res1,res2 in zip(reslist1,reslist2):
                res1_id = resindices[res1]
                res2_id = resindices[res2]
                score = score + pae[res1_id][res2_id]
    else:
        pae = resultsfile['pae_output'][0]
        score = 0
        for res1,res2 in zip(reslist1,reslist2):
            res1_id = resindices[res1]
            res2_id = resindices[res2]
            score = score + pae[res1_id][res2_id]
    return score

def score_confidence_lists(resultsfile, reslist1, reslist2, resindices, fil = True):
    """calculates confidence score of all interactions between two lists of residues"""
    pae = {}
    score = 0
    if fil:
        with open(resultsfile, 'rb') as f:
            p = pickle.load(f)
            pae = p['pae_output'][0]
    else:
        pae = resultsfile['pae_output'][0]

    for res1 in reslist1:
        res1_id = resindices[res1]
        for res2 in reslist2:
            res2_id = resindices[res2]
            score = score + pae[res1_id][res2_id]
    return score

def score_confidence_residues(resultsfile, reslist, resindices, fil = False):
    score = 0
    if fil:
        with open(resultsfile, 'rb') as f:
            p = pickle.load(f)
            plddt = p['plddt']
    else:
        plddt = resultsfile['plddt']

    for res in reslist:
        resid = resindices[res]
        score = score + plddt[resid]

    return score

def get_rmsd(reslist1, pdb1, reslist2, pdb2):

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil=False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil=False)
    A = []
    for res in reslist1:
        for atom in residues1[res]:
            A.append(list(atom[-1]))
    B = []
    for res in reslist2:
        for atom in residues2[res]:
            B.append(list(atom[-1]))
    A = np.array(A)
    B = np.array(B)
    A = A.astype(float)
    B = B.astype(float)
    from calculate_rmsd import kabsch_rmsd
    rmsd = kabsch_rmsd(A, B)
    return rmsd

def write_raw_plddt(results, filename):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist = [x for x in residues.keys()]
    plddt = results['plddt']
    with open(filename, "w") as opf:
        opf.write("plddt: per residue confidences\n")
        for res in reslist:
            resid = resindices[res]
            opf.write(str(res) + "\t" + str(plddt[resid]) + "\n")

def write_pairwise_scores(pairs, results, filename):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    pae = results['pae_output'][0]
    with open(filename, "w") as opf:
        opf.write("pae: pairwise confidence errors\n")
        for pair in pairs:
            res1_id = resindices[pair[0]] 
            res2_id = resindices[pair[1]] 
            opf.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(pae[res1_id][res2_id]) + "\n")
    

if __name__ == "__main__":
    #f1 = "/pine/scr/a/m/amritan/kuhlmanlab/folddesign/folddesign/data/A1_CD20_helix_design.pdb"
    pdb1 = "sequences_0_model_1_multimer_v2_0_unrelaxed.pdb"
    pdb2 = "sequences_0_model_1_multimer_v2_0_unrelaxed.pdb"

    rmsd = get_rmsd(reslist1, pdb1, reslits2, pdb2)
    print(rmsd)
