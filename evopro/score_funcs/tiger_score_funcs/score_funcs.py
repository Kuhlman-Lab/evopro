import sys
#sys.path.append("/proj/kuhl_lab/evopro/")
sys.path.append("/nas/longleaf/home/amritan/Desktop/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.utils.write_pdb import PDBio
from evopro.utils.calc_rmsd import RMSDcalculator
from evopro.score_funcs.calculate_rmsd import kabsch_rmsd
import math
import pickle
import numpy as np
import re

def distance(p1, p2):
    """returns the distance between two 3D points represented as tuples"""

    dist = math.sqrt((float(p2[0])-float(p1[0]))**2+(float(p2[1])-(float(p1[1])))**2+(float(p2[2])-float(p1[2]))**2)
    return dist

def get_seq_indices(dsobj, reslist, first_only=True):
    """
    returns a list with updated residue indices for this sequence (after insertions/deletions)
    if first_only is True, only the first instance of the resid will be returned 
        (insertions at that position are ignored)
    """
    new_reslist = []
    for resid_orig in reslist:
        chain = re.split('(\d+)', resid_orig)[0]
        numbering = dsobj.numbering[chain]
        if first_only:
            try:
                new_resind = numbering.index(resid_orig)
                resid_new = chain + str(new_resind+1)
                new_reslist.append(resid_new)
            except:
                print(resid_orig + " has been deleted")
        else:
            new_resinds = [i for i, x in enumerate(numbering) if x == resid_orig]
            resids_new = [chain + str(x+1) for x in new_resinds]
            new_reslist = new_reslist + resids_new

    return new_reslist

def score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dist=4, contact_cap=36, dsobj=None, first_only=False):
    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)

    chains, residues, resindices = get_coordinates_pdb(pdb)
    pae = results['pae_output'][0]

    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            weight = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        pair = (res1, res2)
                        pair_rev = (res2, res1)
                        if pair not in pairs and pair_rev not in pairs:
                            if len(pairs)<contact_cap:
                                contact=1
                                res1_id = resindices[res1]
                                res2_id = resindices[res2]
                                pae_contact = pae[res1_id][res2_id] + pae[res2_id][res1_id]
                                weight = (70-pae_contact)/70
                                pairs.append(pair)

            score = score + contact*weight
    print("Contact residues:")
    print(pairs)
    print("Contact_score: %.3f, num_contacts: %.3f" %(score, len(pairs)))
    return pairs, score

def score_contacts(pdbfile, reslist1, reslist2, dist=4, score_cap=36, dsobj=None, first_only=False):
    """returns a list of pairs of residues that are making contacts, and the contact score"""
    chains, residues, resindices = get_coordinates_pdb(pdbfile)
    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)
    score = 0
    pairs = []
    for res1 in reslist1:
        for res2 in reslist2:
            contact = 0
            for atom1 in residues[res1]:
                for atom2 in residues[res2]:
                    if distance(atom1[2], atom2[2])<=dist:
                        pair = (res1, res2)
                        pair_rev = (res2, res1)
                        if pair not in pairs and pair_rev not in pairs:
                            pairs.append(pair)
                            contact=1

            score = score+contact
    if score>score_cap:
        score=score_cap

    return pairs, score

def orientation_score(pdb, pairs, orient_dist = 10, penalty = 10, dsobj=None, first_only=True):
    corrects = []
    chains, residues, resindices = get_coordinates_pdb(pdb)
    for tup in pairs:
        correct = 0
        if dsobj:
            tup = get_seq_indices(dsobj, [tup[0],tup[1]], first_only=first_only)
        for atom1 in residues[tup[0]]:
            for atom2 in residues[tup[1]]:
                if distance(atom1[2], atom2[2])<=orient_dist:
                    correct = 1
        
        corrects.append(correct)

    orientation_score = sum([penalty for x in corrects if x==0])

    return orientation_score, corrects

def score_pae_confidence_pairs(resultsfile, pairs, resindices, fil = False, dsobj=None, first_only=True):
    """calculates confidence score of all pairwise residue interactions"""
    score = 0
    reslist1 = []
    reslist2 = []
    for pair in pairs:
        reslist1.append(pair[0])
        reslist2.append(pair[1])

    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)
    if fil:
        with open(resultsfile,'rb') as f:
            p = pickle.load(f)
            pae = p['pae_output'][0]
            for res1,res2 in zip(reslist1,reslist2):
                res1_id = resindices[res1]
                res2_id = resindices[res2]
                score = score + pae[res1_id][res2_id] + pae[res2_id][res1_id]
    else:
        pae = resultsfile['pae_output'][0]
        score = 0
        for res1,res2 in zip(reslist1,reslist2):
            res1_id = resindices[res1]
            res2_id = resindices[res2]
            score = score + pae[res1_id][res2_id] + pae[res2_id][res1_id]
    return score

def score_pae_confidence_lists(resultsfile, reslist1, reslist2, resindices, fil = False, dsobj=None, first_only=False):
    """calculates confidence score of all permutations of pairwise interactions between two lists of residues"""
    pae = {}
    score = 0
    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)
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

def score_plddt_confidence(resultsfile, reslist, resindices, fil = False, dsobj=None, first_only=False):
    score = 0
    if dsobj:
        reslist = get_seq_indices(dsobj, reslist, first_only=first_only)

    if fil:
        with open(resultsfile, 'rb') as f:
            p = pickle.load(f)
            plddt = p['plddt']
    else:
        plddt = resultsfile['plddt']

    for res in reslist:
        resid = resindices[res]
        score = score + plddt[resid]

    return score/len(reslist)

def get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=False, translate=True, dsobj=None, first_only=True):
    if dsobj:
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    A = []
    for res in reslist1:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A.append(list(atom[-1]))
            else:
                A.append(list(atom[-1]))
    B = []
    for res in reslist2:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B.append(list(atom[-1]))
            else:
                B.append(list(atom[-1]))
    A = np.array(A)
    B = np.array(B)
    A = A.astype(float)
    B = B.astype(float)
    rmsd = kabsch_rmsd(A, B, translate=translate)
    return rmsd

def write_raw_plddt(results, filename):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
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
    chains, residues, resindices = get_coordinates_pdb(pdb)
    pae = results['pae_output'][0]
    with open(filename, "w") as opf:
        opf.write("pae: pairwise confidence errors\n")
        for pair in pairs:
            res1_id = resindices[pair[0]] 
            res2_id = resindices[pair[1]] 
            opf.write(str(pair[0]) + "\t" + str(pair[1]) + "\t" + str(pae[res1_id][res2_id]) + "\t" + str(pae[res2_id][res1_id]) + "\n")
    
if __name__ == "__main__":
    #f1 = "/pine/scr/a/m/amritan/kuhlmanlab/folddesign/folddesign/data/A1_CD20_helix_design.pdb"
    path = "/nas/longleaf/home/amritan/Desktop/evopro/evopro/user_inputs/"
    pdb1 = path + "bad_model.pdb"
    with open(pdb1, "r") as f:
        pdb_string = f.read()

    chains, residues, resindices = get_coordinates_pdb(pdb1, fil=True)
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    print(score_contacts(pdb_string, reslist1, reslist2))
    #rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True)
    #print(rmsd)
