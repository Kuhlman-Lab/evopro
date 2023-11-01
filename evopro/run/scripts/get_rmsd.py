from calculate_rmsd import kabsch_rmsd, kabsch_rmsd_superimposeall
import numpy as np

def get_coordinates_pdb(pdb, fil=False):
    lines = []
    chains = []
    residues = {}
    residueindices = {}
    if fil:
        with open(pdb,"r") as f:
            i=0
            for lin in f:
                x = lin[30:38].strip(' ')
                y = lin[38:46].strip(' ')
                z = lin[46:54].strip(' ')
                l = lin.strip().split()
                if l:
                    if 'ATOM' in l[0] or 'HETATM' in l[0]:
                        resid = l[4]+l[5]
                        atominfo = (l[1], l[2], (x, y, z))
                        if l[4] not in chains:
                            chains.append(l[4])
                        if resid not in residues:
                            residues[resid] = [atominfo]
                        else:
                            residues[resid].append(atominfo)
                        if resid not in residueindices:
                            residueindices[resid] = i
                            i = i+1

    else:
        pdb_split = pdb.split("\n")
        i=0
        pdb_split = [x for x in pdb_split if x]
        for lin in pdb_split:
            x = lin[30:38].strip(' ')
            y = lin[38:46].strip(' ')
            z = lin[46:54].strip(' ')
            l = lin.strip().split()
            if 'ATOM' in l[0] or 'HETATM' in l[0]:
                resid = l[4]+l[5]
                atominfo = (l[1], l[2], (x, y, z))
                if l[4] not in chains:
                    chains.append(l[4])
                if resid not in residues:
                    residues[resid] = [atominfo]
                else:
                    residues[resid].append(atominfo)
                if resid not in residueindices:
                    residueindices[resid] = i
                    i = i+1

    return chains, residues, residueindices

def get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=False, translate=True):

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

def get_rmsd_superimposeall(reslist1, pdb1, reslist2, pdb2, ca_only=False, translate=True):

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    
    reslistA = [x for x in residues1.keys()]
    reslistB = [x for x in residues2.keys()]
    for res in reslistA:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A.append(list(atom[-1]))
            else:
                A.append(list(atom[-1]))
    B = []
    for res in reslistB:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B.append(list(atom[-1]))
            else:
                B.append(list(atom[-1]))
    
    A2 = []
    for res in reslist1:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A.append(list(atom[-1]))
            else:
                A.append(list(atom[-1]))
    B2 = []
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
    rmsd = kabsch_rmsd_superimposeall(A, B, A2, B2, translate=translate)
    return rmsd

def get_rmsd_superimposeselected(reslist1, reslist1_all, pdb1, reslist2, reslist2_all, pdb2, ca_only=False, translate=True):

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    
    A=[]
    for res in reslist1_all:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A.append(list(atom[-1]))
            else:
                A.append(list(atom[-1]))
    B = []
    for res in reslist2_all:
        for atom in residues2[res]:
            if ca_only:
                if atom[1] == 'CA':
                    B.append(list(atom[-1]))
            else:
                B.append(list(atom[-1]))
    
    A2 = []
    for res in reslist1:
        for atom in residues1[res]:
            if ca_only:
                if atom[1] == 'CA':
                    A2.append(list(atom[-1]))
            else:
                A2.append(list(atom[-1]))
    B2 = []
    for res in reslist2:
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
    
    #print(A.shape, B.shape, A2.shape, B2.shape)
    
    rmsd = kabsch_rmsd_superimposeall(A, B, A2, B2, translate=translate)
    return rmsd

def score_rmsd(pdb1, pdb2):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2)

    return rmsd

