import argparse
import os, sys
import importlib
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, default='outputs/')
    parser.add_argument('--custom_score',
                        default=None,
                        type=str,
                        help='Score function to use to generate scores for each prediction (optional).')
    parser.add_argument("--new_filename", type=str, default='seqs_and_scores.txt')

    args = parser.parse_args()
    
    pdbdir = args.pdb_dir
    
    pdbs = [os.path.join(pdbdir, f) for f in os.listdir(pdbdir) if f.endswith(".pdb")]
    
    score_func = None
    if args.custom_score:
        file = args.custom_score.split(" ")[0]
        function = args.custom_score.split(" ")[1]
        
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
        
    for pdb in pdbs:
        score = score_func(pdb)
