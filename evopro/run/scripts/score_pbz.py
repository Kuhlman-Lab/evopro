import argparse
import os, sys
import importlib
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

import bz2
import _pickle as cPickle
from typing import Any

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
    parser.add_argument('--custom_score',
                        default=None,
                        type=str,
                        help='Score function to use to generate scores for each zip file (optional).')
    parser.add_argument("--file_prefix", type=str, default='seq_')
    parser.add_argument("--new_filename", type=str, default='seqs_and_scores.csv')

    args = parser.parse_args()
    
    pbzdir = args.pbz_dir
    
    pbzs = [os.path.join(pbzdir, f) for f in os.listdir(pbzdir) if f.endswith(".pbz2") and f.startswith(args.file_prefix)]
    #print(pbzs)
    
    seq_to_data = {}
    for pbz in pbzs:
        data = decompress_pickle(pbz)
        i = int(pbz.split("/")[-1].split("_")[1])
        seq_to_data[i] = data
        
    score_func = None
    if args.custom_score:
        file = args.custom_score.split(",")[0]
        function = args.custom_score.split(",")[1]
        
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
        
    myKeys = list(scores.keys())
    myKeys.sort()
    sorted_scores = {i: scores[i] for i in myKeys}
    with open(args.new_filename, "w") as f:
        for val in sorted_scores:
            score_extended = ",".join([str(x) for x in list(scores[val][1])])
            f.write(str(val) + "," + str(scores[val][0]) + "," + score_extended + "\n")