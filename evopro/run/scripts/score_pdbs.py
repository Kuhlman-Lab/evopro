"""Script that parses a directory of PDB files and scores them based on a custom scoring function. (By default, returns average plDDT score)
Generates a sorted file containing the filename and score for each PDB."""


import argparse
import os, sys
import importlib
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, default='./')
    parser.add_argument('--custom_score_file',
                        default=None,
                        type=str,
                        help='File containing score function to use to generate scores for each prediction (optional).')
    parser.add_argument('--custom_score_func',
                        default=None,
                        type=str,
                        help='function to use from score file above.')
    parser.add_argument("--new_filename", type=str, default='seqs_and_scores.txt')

    args = parser.parse_args()
    
    pdbdir = args.pdb_dir
    
    pdbs = [os.path.join(pdbdir, f) for f in os.listdir(pdbdir) if f.endswith(".pdb")]
    
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
        
    scores = []
    for pdb in pdbs:
        score = score_func(pdb)
        scores.append(score)
        
    with open(args.new_filename, "w") as f:
        for pdb, score in zip(pdbs, scores):
            f.write("{} {}\n".format(pdb, score))
