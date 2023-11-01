import os, sys
import importlib

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser

def score_diff_structures(diff_backbones, score_func, outfile="scores.txt", diff_chain="B", sort_by_score=0):
    structs_and_scores = []
    for pdbfile in diff_backbones:
        with open(pdbfile, "r") as fil:
            pdb = fil.read()
        
        score = score_func(pdb, diff_chain=diff_chain)
        structs_and_scores.append((pdbfile, score)) # (pdb name, score_tuple)
        
    #sort by first (zero-th) score column unless otherwise specified
    structs_and_scores.sort(key=lambda x: x[1][sort_by_score])
    
    with open(outfile, "w") as outf:
        for s in structs_and_scores:
            outf.write(s[0] + "\t" + str(s[1]) + "\n")
            
def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--pdb_dir',
                        default='./',
                        type=str,
                        help='Path to and name of directory with input pdb files.')
    
    parser.add_argument('--score_file',
                        default=None,
                        type=str,
                        help='Name of PDB file (if only one) to generate sequences using MPNN and check RMSD against.')
    
    parser.add_argument('--score_func',
                        default=None,
                        type=str,
                        help='Name of PDB file (if only one) to generate sequences using MPNN and check RMSD against.')
    
    return parser
        
if __name__=="__main__":
    
    parser = getFlagParser()
    args = parser.parse_args(sys.argv[1:])
    
    #pdb_dir = "/work/users/a/m/amritan/lpl/angptl3/rfdiff/run1/outputs/"
    #score_file = "/proj/kuhl_lab/evopro/evopro/score_funcs/diff_score_backbone_binder.py"
    #score_func_name = "score_diff_backbone_rG"
    
    try:
        scorefile = args.score_file.rsplit("/", 1)
        scorepath = scorefile[0]
        scorefilename = scorefile[1].split(".")[0]

        sys.path.append(scorepath)
        mod = importlib.import_module(scorefilename)
        scorefunc = getattr(mod, args.score_func)
    except:
        raise ValueError("Invalid score function")
    
    diff_backbones = []
    for pdb in os.listdir(args.pdb_dir):
        if pdb.endswith(".pdb"):
            diff_backbones.append(os.path.join(args.pdb_dir, pdb))
    
    score_diff_structures(diff_backbones, scorefunc)
    
    