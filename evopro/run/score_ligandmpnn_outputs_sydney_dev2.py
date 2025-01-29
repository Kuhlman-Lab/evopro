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

    parser.add_argument('--seq_file_list', type=list_of_strings, help='List of files paths for csv files containing sequences')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for AlphaFold2. Default is 1.')
    
    parser.add_argument('--af2_flags_file',
                        default='./af2.flags',
                        type=str,
                        help='Path to and name of af2.flags file.')

    return parser


def run_af2_dist(args, score_func=None):

    file_list = args.seq_file_list
    print(file_list)

    seqs_list = []

    counter = 0
    for file in file_list:
        with open(file, "r") as f:
            seq_num = 0
            if counter == 0:
                for line in f:
                    seq = [x for x in line.strip("\n").split(",") if x]
                    print(seq)
                    print("\n")
                    seqs_list.append([seq])
            else:
                for line in f:
                    seq = [x for x in line.strip("\n").split(",") if x]
                    seqs_list[seq_num].append(seq)
                    seq_num += 1
        counter += 1

    print(seqs_list)

    lengths = []

    for file in file_list:
        df = pd.read_csv(file, header=None)
        col_num = df.shape[1]
        for i in range(0, col_num - 1):
            i += 1
            max_val = df.iloc[:,i].str.len().max()
            lengths.append(max_val)


    print("Compiling AF2 models for lengths:", lengths)
    dist = Distributor(args.n_workers, af2_init, args.af2_flags_file, lengths)

    results = dist.churn(seqs_list)

    print("done churning")
    dist.spin_down()

    print("results begin here")
    print(results)


if __name__ == "__main__":
    parser = getFlagParser()
    args = parser.parse_args(sys.argv[1:])
    
    score_func = None
        
    run_af2_dist(args, score_func=score_func)