"""Script that runs AF2 (distributor, runtime-optimized) on a list of sequences and returns the scores and PDB files.
By default, the PDBs are scored for average plDDT score and sorted by that score. A custom scoring function can be provided instead."""

import sys, os
import importlib

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

import pandas as pd

def list_of_strings(arg):
    return arg.split(',')

def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--af2_flags_file',
                        default='./af2.flags',
                        type=str,
                        help='Path to and name of af2.flags file.')
    

    parser.add_argument('--seq_file_list', type=list_of_strings, help='List of files paths for csv files containing sequences')
    
    parser.add_argument('--sequence_file',
                        default="sequences.csv",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--output_dir',
                        default="outputs",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--max_chains_length',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, overall max length of each chain separated by commas.')
    
    parser.add_argument('--custom_score',
                        default=None,
                        type=str,
                        help='Score function to use to generate scores for each prediction (optional).')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for AlphaFold2. Default is 1.')

    return parser


def run_af2_dist(args, score_func=None):

    file_list = args.seq_file_list

    counter = 0
    seq_dict = {}
    for file in file_list:
        with open(file, "r") as f:
            seq_num = 0
            if counter == 0:
                for line in f:
                    seq = [x for x in line.strip("\n").split(",") if x]
                    seq_dict[seq_num] = seq
                    seq_num += 1
            else:
                for line in f:
                    seq = [x for x in line.strip("\n").split(",") if x]
                    current_seq_list = seq_dict[seq_num]
                    current_seq_list.append(seq)
                    seq_num += 1
        counter += 1

    seqs_list = list(seq_dict.values())
    print(seqs_list)

    lengths = []

    for file in file_list:
        df = pd.read_csv(file, header=None)
        col_num = df.shape[1]
        for i in range(0, len(col_num) - 1):
            i += 1
            max_val = df.iloc[:,i].str.len().max()
            lengths.append(max_val)

    print("Compiling AF2 models for lengths:", lengths)
    dist = Distributor(args.n_workers, af2_init, args.af2_flags_file, lengths)

    results = dist.churn(seqs_list)

    print("done churning")
    dist.spin_down()

    scores = []
    for result in results:
        scores.append(score_func(result))


    print("scores before normalizing")
    for score in scores:
        print(score[1])
        print(f"overall score is {score[0]}")
    print("done writing old scores")


    scores_to_norm = []
    scores_no_norm = []
    for score in scores:
        if score[0] == 1000:
            scores_no_norm.append(score)
        else:
            scores_to_norm.append(score)

    if len(scores_to_norm) != 0:

        #here 0 is arbitrary since the number of score terms is the same across all scores
        num_score_terms = len(scores_to_norm[0][1]) 

        #each tuple in score_cols is a column of scores ex. monomer prediction for (design 1, design 2), pae for (design 1, design 2)
        score_cols = []
        for j in range(num_score_terms):
            #where i in the index of that design, 1 holds the score, j is the position of the score term in the score tuple, and each tuple begins with the actual score to return (other values in the tuple should be used for debugging, not scoring)
            col = tuple(scores_to_norm[i][1][j][0] for i in range(len(scores_to_norm)))
            score_cols.append(col)


        min_vals = [min(col) for col in score_cols]
        max_vals = [max(col) for col in score_cols]

        for index, tup in enumerate(scores_to_norm):
            scores = tup[1]
            new_score = ()
            overall_score = 0
            for i in range(len(scores)):
                value = scores[i][0]

                if max_vals[i] != min_vals[i]:
                    scaled_val = (value - min_vals[i])/(max_vals[i] - min_vals[i])
                else:
                    scaled_val = 0 #if min = max then all values = min = max so scaled_val is zero

                if len(scores[i]) > 1:
                    new_tuple = ((scaled_val,) + tuple(scores[i][1:]),) #keeping the rest of the tuple as is
                else:
                    new_tuple = ((scaled_val,),) #making a nested tuple to ensure the format is the same as the initial scores

                new_score += new_tuple #creating a new scoring tuple to replace the current tuple
                overall_score += scaled_val

            scores_to_norm[index] = [overall_score, new_score, tup[2], tup[3]]

    final_scores = scores_to_norm + scores_no_norm

    print("normalizing scores")
    for scores in final_scores:
        print(scores[1])
        print(f"overall score is {scores[0]}")
    print("done normalizing")