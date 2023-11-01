import os
from Bio import SeqIO
from get_rmsd import get_rmsd_superimposeselected, get_coordinates_pdb

import warnings
warnings.filterwarnings("ignore")

def get_seq_from_pdb(pdbfilename):
    seq = []
    with open(pdbfilename, "r") as f:
        for record in SeqIO.parse(f, 'pdb-atom'):
            seq.append(record.seq)
    
    return str("".join(seq[-1]))

rf2_path = "/work/users/a/m/amritan/lpl/angptl3/evopro/rf2_evopro/large_run1/all_pdbs/complex/"
af2_path = "/work/users/a/m/amritan/lpl/angptl3/evopro/rf2_evopro/run_af2_on_outputs/large_run1/dist_run/outputs/"

rf2_files = []
for i in range(1, 21):
    for j in range(1, 6):
        for k in range(0, 10):
            rf2_files.append(rf2_path + "pair" + str(i) + "_run" + str(j) + "_seq_" + str(k) + "_final_model_1_complex.pdb")
        
pdb_dict_by_seq = {}
for file in rf2_files:
    seq = get_seq_from_pdb(file)
    with open(file, "r") as f:
        pdb_str = f.read()
    pdb_dict_by_seq[seq] = (pdb_str, file)
    
af2_files = [os.path.join(af2_path, "seq_" + str(x) + "_model_1.pdb") for x in range(950)]
#rf2_files = [os.path.join(rf2_path + x, "S_00_pred.pdb") for x in os.listdir(rf2_path) if os.path.isdir(os.path.join(rf2_path + x)) and os.path.isfile(os.path.join(rf2_path + x, "S_00_pred.pdb"))]

#print(len(pdb_dict_by_seq), len(rf2_files))
data = []

for file in af2_files:
    seq = get_seq_from_pdb(file)
    #print(seq)
    rf2_pdb, rf2_pdb_path = pdb_dict_by_seq[seq]
    with open(file, "r") as f:
        af2_pdb = f.read()
    
    chains1, residues1, resindices1 = get_coordinates_pdb(af2_pdb)
    chains2, residues2, resindices2 = get_coordinates_pdb(rf2_pdb)
    reslist1_all = [x for x in residues2.keys()]
    reslist2_all = [x for x in residues2.keys()]
    reslist1 = [x for x in residues1.keys() if x.startswith("D") or x.startswith("E") or x.startswith("F")]
    reslist2 = [x for x in residues2.keys() if x.startswith("D") or x.startswith("E") or x.startswith("F")]

    #print(reslist1_all, reslist2_all)
    #print(reslist1, reslist2)
    rmsd = get_rmsd_superimposeselected(reslist1, reslist1_all, af2_pdb, reslist2, reslist2_all, rf2_pdb, ca_only=True)
    data.append((rf2_pdb_path, file, seq, rmsd))
    print(rf2_pdb_path, file, seq, rmsd)

data.sort(key=lambda x: x[3])
with open("af2_rmsd.txt", "w") as f:
    for d in data:
        f.write("\t".join([str(x) for x in d]) + "\n")
