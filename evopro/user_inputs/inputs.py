import argparse
import sys
from typing import Optional, Sequence, Tuple, Dict

class FileArgumentParser(argparse.ArgumentParser):
    """Overwrites default ArgumentParser to better handle flag files."""

    def convert_arg_line_to_args(self, arg_line: str) -> Optional[Sequence[str]]:
        """ Read from files where each line contains a flag and its value, e.g.
        '--flag value'. Also safely ignores comments denoted with '#' and
        empty lines.
        """

        # Remove any comments from the line
        arg_line = arg_line.split('#')[0]

        # Escape if the line is empty
        if not arg_line:
            return None

        # Separate flag and values
        split_line = arg_line.strip().split(' ')

        # If there is actually a value, return the flag value pair
        if len(split_line) > 1:
            return [split_line[0], ' '.join(split_line[1:])]
        # Return just flag if there is no value
        else:
            return split_line

def getEvoProParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run EvoPro."""

    parser = FileArgumentParser(description='EvoPro runner script that can take '
                                'command-line arguments or read from an '
                                'argument flag file.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--input_dir',
                        default='.',
                        type=str,
                        help='Path to directory that contains input files. Default is ./')

    parser.add_argument('--num_iter',
                        default='50',
                        type=int,
                        help='Number of iterations of genetic algorithm. Default is 50.')

    parser.add_argument('--pool_size',
                        default='20',
                        type=int,
                        help='Size of "genetic pool", or the number of sequences evaluated per '
                        'iteration. Default is 20.')

    parser.add_argument('--pool_size_variable',
                        action='store_true',
                        help='Specify a file pool_sizes.txt with pool sizes for every iteration (or until pool size stops changing). '
                        'Defaults to False and uses constant pool size.')

    parser.add_argument('--num_gpus',
                        default='1',
                        type=int,
                        help='Number of gpus available. Default is 1.')

    parser.add_argument('--score_file',
                        default='',
                        type=str,
                        help='Path and file name of python script containing the score'
                        ' function used to evaluate fitness of the alphafold predictions. Required.')

    parser.add_argument('--score_func',
                        default='',
                        type=str,
                        help='Name of the score function in the score file'
                        ' used to evaluate fitness of the alphafold predictions. Required.')

    parser.add_argument('--rmsd_func',
                        default=None,
                        type=str,
                        help='Name of the rmsd function in the score file'
                        ' used to evaluate fitness of the alphafold predictions. Optional, requires stabilize_binder=True.')

    parser.add_argument('--rmsd_to_starting',
                        default=None,
                        type=str,
                        help='Name of the rmsd function in the score file and the path/filename to pdb for RMSD.'
                        ' used to evaluate rmsd to the starting scaffold using a U-shaped potential. ')
    parser.add_argument('--path_to_starting',
                        default=None,
                        type=str,
                        help='path/filename to pdb to pass to scoring function for RMSD to starting.')

    parser.add_argument('--define_contact_area',
                        default=None,
                        type=str,
                        help='Defining residues on target interface to be targeted for contacts. Default is None.')

    parser.add_argument('--bonus_contacts',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be given a bonus for making contacts, followed by the'
                        'distance cutoff. Default is None and 4A.')

    parser.add_argument('--penalize_contacts',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be given a penalty for making contacts, followed by the'
                        'distance cutoff. Default is None and 4A.')
    
    parser.add_argument('--no_repeat_af2',
                         action='store_false',
                         help='Use this flag to specify if you do not want AF2 to be run multiple times on the same sequence. Default is False.')

    parser.add_argument('--dont_write_compressed_data',
                         action='store_false',
                         help='Default is True.')

    parser.add_argument('--write_pdbs',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--mpnn_iters',
                        default=None,
                        type=str,
                        help='Iteration numbers at which MPNN is used to refill the pool. Defaults to mpnn_freq')

    parser.add_argument('--mpnn_freq',
                        default='10',
                        type=int,
                        help='Protein MPNN is used to refill the pool once every _ iterations.'
                        'Default is 10.')

    parser.add_argument('--skip_mpnn',
                        default=None,
                        type=str,
                        help='Skip MPNN refilling in these iterations. Default is None.')
    
    parser.add_argument('--mpnn_temp',
                        default='0.1',
                        type=str,
                        help='Protein MPNN is used to refill the pool at this sampling temperature.'
                        'Default is 0.1.')
    
    parser.add_argument('--mpnn_version',
                        default="s_48_020",
                        type=str,
                        help='Model version used to run MPNN. Default is s_48_020 (soluble).')

    parser.add_argument('--plot_confidences',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--plot_scores_avg',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--plot_scores_median',
                         action='store_true',
                         help='Default is False.')
    
    parser.add_argument('--plot_scores_top',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--crossover_percent',
                        default='0.2',
                        type=float,
                        help='Fraction of pool refilled by crossover.'
                        'Default is 0.2.')

    parser.add_argument('--mutation_percents',
                        default='0.125',
                        type=str,
                        help='Number of mutations made as a percentage of sequence length.'
                        'Default is 0.125 for every iteration. If more than one value is provided, number of iterations will be split evenly and assigned.')
    
    parser.add_argument('--af2_preds_extra',
                        default=None,
                        type=str,
                        help='Chain ID permutations to run through individual AF2 runs, separated by commas. Default is None.')

    return parser

if __name__ == "__main__":
    parser = getEvoProParser()
    args = parser.parse_args(sys.argv[1:])
    print(args)
