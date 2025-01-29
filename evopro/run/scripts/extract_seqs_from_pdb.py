"""Script that parses a directory of multi-chain PDB files and writes out an AF2-compatible sequence file."""

import os, sys
from Bio import SeqIO   
import argparse

import warnings
warnings.filterwarnings("ignore") 

def get_seq_from_pdb(pdbfilename):
    seq = ""
    with open(pdbfilename, "r") as f:
        for record in SeqIO.parse(f, 'pdb-atom'):
            seq += "," + record.seq
        print(seq)
    
    return seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, default='./all_pdbs/complex/')
    parser.add_argument("--new_filename", type=str, default='sequences.csv')

    args = parser.parse_args()
    
    pdbdir = args.pdb_dir

    seqs = []
    pdbs = [os.path.join(pdbdir, x) for x in os.listdir(pdbdir) if x.endswith(".pdb")]
    for pdb in pdbs:
        seq = get_seq_from_pdb(pdb)
        seqs.append(seq)
    
    with open(args.new_filename, "w") as f:
        for seq in seqs:
            f.write(str(seq) + "\n")
    
    
