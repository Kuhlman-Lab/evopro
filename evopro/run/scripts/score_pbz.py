"""Script that parses a directory of PBZ2 files and scores them based on a custom scoring function. (By default, returns average plDDT score)
Generates a sorted file containing the filename and score for each compressed file."""

import argparse
import os, sys
import importlib
from Bio import SeqIO
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall
from alphafold.common import protein

import bz2
import _pickle as cPickle
from typing import Any

def get_seq_from_pdb(pdbfilename, binder_chain_index=-1):
    seq = []
    with open(pdbfilename, "r") as f:
        for record in SeqIO.parse(f, 'pdb-atom'):
            seq.append(record.seq)

    return str("".join(seq[binder_chain_index]))

def decompress_pickle(file: str) -> Any:
    """
    Load any compressed pickle file.
    """
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pbz_dir", type=str, default='outputs/')
    parser.add_argument('--custom_score_file',
                        default=None,
                        type=str,
                        help='File containing score function to use to generate scores for each prediction (optional).')
    parser.add_argument('--custom_score_func',
                        default=None,
                        type=str,
                        help='function to use from score file above.')
    parser.add_argument("--file_prefix", type=str, default='seq_')
    parser.add_argument("--new_filename", type=str, default='seqs_and_scores.csv')
    parser.add_argument("--binder_chain_index", type=int, default=-1, help='Default is -1 (binder is last chain). This value is zero-indexed.')


    args = parser.parse_args()
    
    pbzdir = args.pbz_dir
    
    pbzs = [os.path.join(pbzdir, f) for f in os.listdir(pbzdir) if f.endswith(".pbz2") and f.startswith(args.file_prefix)]
    #print(pbzs)
    
    seq_to_data = {}
    for pbz in pbzs:
        data = decompress_pickle(pbz)
        while type(data)==list:
            data = data[0]
        pdb = protein.to_pdb(data['unrelaxed_protein'])
        with open("temp.pdb", "w") as f:
            f.write(pdb)
        seq = get_seq_from_pdb("temp.pdb", binder_chain_index=args.binder_chain_index)
        os.remove("temp.pdb")
        #i = int(pbz.split("/")[-1].split("_")[1])
        seq_to_data[seq] = data
        
    score_func = None
    if args.custom_score_file:
        file = args.custom_score_file
        function = args.custom_score_func
        
        try:
            scorefile = file.rsplit("/", 1)
            scorepath = scorefile[0]
            scorefilename = scorefile[1].split(".")[0]

            sys.path.append(scorepath)
            mod = importlib.import_module(scorefilename)
            score_func = getattr(mod, function)
        except:
            raise ValueError("Invalid score function. Please provide a valid python file and function name within that file, separated by a space.")
    else:
        score_func = score_plddt_confidence_overall
    
    scores = {}
    for val in seq_to_data:
        print(val)
        score = score_func(seq_to_data[val])
        scores[val] = score
        
    #myKeys = list(scores.keys())
    #myKeys.sort()
    #sorted_scores = {i: scores[i] for i in myKeys}
    #sorted_scores = sorted(scores)
    with open(args.new_filename, "w") as f:
        for val in scores:
            score_extended = ",".join([str(x) for x in list(scores[val][1])])
            f.write(str(val) + "," + str(scores[val][0]) + "," + score_extended + "\n")