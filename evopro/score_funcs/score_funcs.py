import sys
#sys.path.append("/nas/longleaf/home/amritan/Desktop/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.calculate_rmsd import kabsch_rmsd, kabsch_rmsd_superimposeall
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
    #print("Renumbering...", dsobj, dsobj.numbering, dsobj.jsondata)
    new_reslist = []
    for resid_orig in reslist:
        #print("Original resid: " + resid_orig)
        chain = re.split('(\d+)', resid_orig)[0]
        numbering = dsobj.numbering[chain]
        #print("Numbering for chain", chain, numbering)
        if first_only:
            try:
                new_resind = numbering.index(resid_orig)
                resid_new = chain + str(new_resind+1)
                new_reslist.append(resid_new)
                #print("New resid: " + resid_new)
            except:
                print("Error....perhaps", resid_orig, "has been deleted")
        else:
            new_resinds = [i for i, x in enumerate(numbering, start=1) if x == resid_orig]
            #print(new_resinds)
            resids_new = [chain + str(x) for x in new_resinds]
            new_reslist = new_reslist + resids_new

    #print(reslist)
    #print(new_reslist)
    return new_reslist


def score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dist=4, contact_cap=36, dsobj=None, first_only=False):
    if dsobj:
        reslist1 = get_seq_indices(dsobj, reslist1, first_only=first_only)
        reslist2 = get_seq_indices(dsobj, reslist2, first_only=first_only)

    chains, residues, resindices = get_coordinates_pdb(pdb)
    pae = results['pae_output'][0]
    
    #print(pae, pae.shape)

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

def score_pae_interaction(pae, pdb, target_length=113):
    _, residues, resindices = get_coordinates_pdb(pdb)

    full_length = len(residues)
    #print(full_length)
    binder_length = full_length - target_length
    #print(binder_length)
    pae_interaction1 = np.mean( pae[:binder_length,binder_length:] )
    pae_interaction2 = np.mean( pae[binder_length:,:binder_length] )

    return ( pae_interaction1 + pae_interaction2 ) / 2

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

def score_plddt_confidence(resultsfile, reslist, resindices, fil = False, dsobj=None, first_only=False, rf2 = False):
    score = 0
    if dsobj:
        reslist = get_seq_indices(dsobj, reslist, first_only=first_only)

    if fil:
        with open(resultsfile, 'rb') as f:
            p = pickle.load(f)
            if rf2:
                plddt = p['lddt']
            else:
                plddt = p['plddt']
    else:
        if rf2:
            plddt = resultsfile['lddt']
        else:
            plddt = resultsfile['plddt']

    for res in reslist:
        resid = resindices[res]
        score = score + plddt[resid]

    if len(reslist)<1:
        return 100
    else:
        return score/len(reslist)

def score_plddt_confidence_overall(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    plddt = results['plddt']
    confscore = sum(plddt)/len(plddt)
    return confscore, (confscore, ), pdb, results

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

def get_rmsd_superimposeall(reslist1, reslist1_2, pdb1, reslist2, reslist2_2, pdb2, ca_only=False, translate=True, dsobj=None, first_only=True):

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
    A2 = []      
    for res in reslist1_2:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A2.append(list(atom[-1]))
            else:
                A2.append(list(atom[-1]))
                
    B = []
    for res in reslist2:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B.append(list(atom[-1]))
            else:
                B.append(list(atom[-1]))
    
    B2 = []
    for res in reslist2_2:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B2.append(list(atom[-1]))
            else:
                B2.append(list(atom[-1]))
                
                
    A = np.array(A)
    B = np.array(B)
    A = A.astype(float)
    B = B.astype(float)
    A2 = np.array(A2)
    B2 = np.array(B2)
    A2 = A2.astype(float)
    B2 = B2.astype(float)
    
    rmsd = kabsch_rmsd_superimposeall(A, B, A2, B2, translate=translate)
    return rmsd


def Rg(pdb, chnid="A"):
    '''
    Calculates the Radius of Gyration (Rg) of a protein given its .pdb 
    structure file. Returns the Rg integer value in Angstrom.
    '''
    coord = list()
    mass = list()
    #Structure = open(filename, 'r')
    pdb = pdb.split("\n")
    for line in pdb:
        line = line.strip().split()
        #print(line)
        if len(line) > 2 and line[4] == chnid:
            try:
                x = float(line[6])
                y = float(line[7])
                z = float(line[8])
                coord.append([x, y, z])
                if line[2] == 'C' or line[2] == 'CA':
                    mass.append(12.0107)
                elif line[2] == 'O':
                    mass.append(15.9994)
                elif line[2] == 'N':
                    mass.append(14.0067)
                elif line[2] == 'S':
                    mass.append(32.065)
            except:
                pass
    #print(coord)
    #print(mass)
    #print(mass, coord)
    xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
    #print(xm)
    tmass = sum(mass)
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
    mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
    rg = math.sqrt(rr / tmass-mm)
    return round(rg, 3)

def radius_of_gyration(pdb, reslist=None):
    coord = list()
    mass = list()
    chains, residues, resindices = get_coordinates_pdb(pdb)
    #print(residues)
    pdb = pdb.split("\n")
    print("Number of pdb lines", len(pdb))
    if not reslist:
        reslist = [x for x in residues.keys()]
    
    for line in pdb:
        #print(line)
        try:
            line = line.split()
            #print("checkpoint1")
            x = float(line[6])
            y = float(line[7])
            z = float(line[8])
            #print("checkpoint2")
            coord.append([x, y, z])
            #print("checkpoint3")
            if line[-1] == 'C':
                mass.append(12.0107)
            elif line[-1] == 'O':
                mass.append(15.9994)
            elif line[-1] == 'N':
                mass.append(14.0067)
                
            elif line[-1] == 'S':
                mass.append(32.065)
            #print("checkpoint4")
        except:
            #print("passing")
            pass
    print(mass, coord)
    xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
    tmass = sum(mass)
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
    mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
    print(xm, tmass, mm, rr)
    if tmass-mm == 0:
        rg = 0
    else:
        rg = math.sqrt(rr / tmass-mm)
    return(round(rg, 3))

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
    
def score_rmsd_to_starting(pdb, path_to_starting, dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

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
    #f1 = "/pine/scr/a/m/amritan/kuhlmanlab/folddesign/folddesign/data/A1_CD20_helix_design.pdb"
    
    path = "/work/users/a/m/amritan/dv_evopro/test_rmsd_calc/"
    pdb1 = path + "design_renumbered_swap.pdb"
    with open(pdb1, "r") as f:
        pdb1_string = f.read()
        
    pdb2 = path + "seq_0_model_1.pdb"
    with open(pdb2, "r") as f:
        pdb2_string = f.read()

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1_string)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2_string)
    reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys()]
    #print(score_contacts(pdb_string, reslist1, reslist2))
    rmsd = get_rmsd_sorted(reslist1, pdb1_string, reslist2, pdb2_string, ca_only=True)
    print(rmsd)
    #print(radius_of_gyration(pdb_string))
    
    """for i in range(100):
        path = "/work/users/a/m/amritan/cd19/mutcd19/mpnn/"
        pdb1 = path + "diff_3.pdb"
        pdb2 = path + "outputs/seq_"+str(i)+"_model_1.pdb"
        with open(pdb1, "r") as f:
            pdb1_string = f.read()
        with open(pdb2, "r") as f:
            pdb2_string = f.read()
        chains1, residues1, resindices1 = get_coordinates_pdb(pdb1_string)
        chains2, residues2, resindices2 = get_coordinates_pdb(pdb2_string)
        #print(residues1, residues2)
        reslist1 = [x for x in residues1.keys()]
        reslist2 = [x for x in residues2.keys()]
        reslist1_2 = [x for x in residues1.keys() if x.startswith("B")]
        reslist2_2 = [x for x in residues2.keys() if x.startswith("B")]
        #print(reslist1, reslist2, reslist1_2, reslist2_2)
        print(i, get_rmsd_superimposeall(reslist1, reslist1_2, pdb1_string, reslist2, reslist2_2, pdb2_string, ca_only=True, translate=True))"""
