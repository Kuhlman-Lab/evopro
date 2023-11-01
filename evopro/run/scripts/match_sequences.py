import os, sys
import argparse
from Bio import SeqIO
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import get_rmsd

import warnings
warnings.filterwarnings("ignore")

def get_seq_from_pdb(pdbfilename):
    seq = []
    with open(pdbfilename, "r") as f:
        for record in SeqIO.parse(f, 'pdb-atom'):
            seq.append(record.seq)
    
    return str("".join(seq[-1]))

def score_rmsd(pdb1, pdb2):
    _, residues1, _ = get_coordinates_pdb(pdb1)
    _, residues2, _ = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2, ca_only=True)

    return rmsd

def get_average_plddt(pdb_str):

    plddts = []

    pdb_split = pdb_str.split("\n")

    pdb_split = [x for x in pdb_split if x]
    for lin in pdb_split:
        b = lin[60:66].strip(' ')
        l = lin.strip().split()
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            if l[2] == "CA":

                plddts.append(float(b))

    return sum(plddts)/len(plddts)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--rf2_dir", type=str, default="./rf2_pdbs/")
    parser.add_argument("--af2_dir", type=str, default='./outputs/')
    parser.add_argument("--output_file", type=str, default='af2_rf2_rmsd.csv')
    
    parser.add_argument("--get_rmsd", action='store_true', help='Default is False.')
    parser.add_argument("--get_rf2_plddt", action='store_true', help='Default is False.')
    
    args = parser.parse_args()
    
    af2_pdbs = [os.path.join(args.af2_dir, x) for x in os.listdir(args.af2_dir) if x.endswith(".pdb")]
    rf2_pdbs = [os.path.join(args.rf2_dir, x) for x in os.listdir(args.rf2_dir) if x.endswith(".pdb")]
    
    pdb_dict_by_seq = {}
    for file in rf2_pdbs:
        seq = get_seq_from_pdb(file)
        with open(file, "r") as f:
            pdb_str = f.read()
        if args.get_rf2_plddt:
            plddt = get_average_plddt(pdb_str)
        else:
            plddt = None
        pdb_dict_by_seq[seq] = (pdb_str, file, plddt)
        
    for file in af2_pdbs:
        seq = get_seq_from_pdb(file)
        print(seq)
        rf2_pdb, rf2_pdb_path, rf2_plddt = pdb_dict_by_seq[seq]
        with open(file, "r") as f:
            af2_pdb = f.read()
        
        if args.get_rmsd:
            rmsd = score_rmsd(af2_pdb, rf2_pdb)
        else:
            rmsd = None
        
        af2_file_name = file.split("/")[-1].split(".")[0]
        rf2_file_name = rf2_pdb_path.split("/")[-1].split(".")[0]
        index = int(af2_file_name.split("_")[1])
        
        pdb_dict_by_seq[seq] = [index, rf2_file_name, af2_file_name, rmsd, rf2_plddt, rf2_pdb_path, file, rf2_pdb, af2_pdb]
    
    with open(args.output_file, "w") as opf:
        for seq in pdb_dict_by_seq:
            opf.write(seq + "," + ",".join([str(x) for x in pdb_dict_by_seq[seq][:-4]]) + "\n")