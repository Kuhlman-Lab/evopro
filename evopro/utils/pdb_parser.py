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

def get_coordinates_pdb_old(pdb, fil = False):
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
                        resid = l[4]+'_'+l[3]+'_'+l[5]
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
                resid = l[4]+'_'+l[3]+'_'+l[5]
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

if __name__ == "__main__":
	c, r, ri = get_coordinates_pdb("test.pdb", fil=True)
	print(c, r, ri)
