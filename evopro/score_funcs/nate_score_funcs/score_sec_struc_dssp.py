path_to_pdb = "/pine/scr/s/m/smyersn/af2/test_scaffolds/g_insertions/kggnaggd/outputs/sequences_0_model_1_ptm_0_unrelaxed.pdb"
target =  'LEEEEEELLLLLLLHHHHHHHHHHHLLEEEEELLLEEEEELLLLHHHHHHHHHLLLLEEEEEL'

def sec_struc_score(target, pdb_string):
    
    from Bio.PDB.DSSP import dssp_dict_from_pdb_file
    import os

    with open('temp.pdb', 'w') as f:
        f.write(pdb_string)

    dssp_tuple = dssp_dict_from_pdb_file('temp.pdb')
    dssp_dict = dssp_tuple[0]

    struc = []
    for i in dssp_dict:
        res = dssp_dict[i][1]
        if res == 'E':
            struc.append(res)
        elif res == 'H':
            struc.append(res)
        else:
            struc.append('L')

    seq=''.join(struc)

    points = 0

    for i in range(len(target)):
        if target[i] == seq[i]:
            points += 1
        else:
            points += 0

    score = (points/len(target)) * 100

    if os.path.exists('temp.pdb'):
        os.remove('temp.pdb')
    
    return(score)

if __name__ == "__main__":
    sec_struc_score(target, pdb_string)
