def get_coordinates_pdb(pdb, fil = False, of=False):
    lines = []
    chains = []
    residues = {}
    residueindices = {}
    if fil:
        with open(pdb,"r") as f:
            i=0
            for lin in f:
                l = lin.strip().split()
                if l:
                    if 'ATOM' in l[0] or 'HETATM' in l[0]:
                        if of:
                            resid = l[4]+'_'+l[3]+'_'+str(int(l[5])+1)
                        else:
                            resid = l[4]+'_'+l[3]+'_'+l[5]
                        atominfo = (l[1], l[2], (l[6], l[7], l[8]))
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
            l = lin.strip().split()
            if 'ATOM' in l[0] or 'HETATM' in l[0]:
                if of:
                    resid = l[4]+'_'+l[3]+'_'+str(int(l[5])+1)
                else:
                    resid = l[4]+'_'+l[3]+'_'+l[5]
                atominfo = (l[1], l[2], (l[6], l[7], l[8]))
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
	c, r, ri = get_coordinates_pdb("/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/user_inputs/sequences_0_model_1_multimer_v2_0_unrelaxed.pdb")
	print(c, r, ri)
