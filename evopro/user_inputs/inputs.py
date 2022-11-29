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
                        default='40',
                        type=int,
                        help='Size of "genetic pool", or the number of sequences evaluated per '
                        'iteration. Default is 40.')

    parser.add_argument('--num_gpus',
                        default='4',
                        type=int,
                        help='Number of gpus available.'
                        'Default is 4.')

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

    parser.add_argument('--vary_length',
                        default='0',
                        type=int,
                        help='How much the length is allowed to vary. Default is 0.')

    parser.add_argument('--define_contact_area',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be targeted for contacts. Default is None.')

    parser.add_argument('--write_plddt',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--write_pae',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--stabilize_monomer',
                        default=None,
                        type=str,
                        help='Chain ID to run through AF2 individually (optional). Default is none.')

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

    parser.add_argument('--len_binder',
                        default='0',
                        type=int,
                        help='Number of residues in binder, used to find binder when linker is included. Default is 0.')

    parser.add_argument('--use_omegafold',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--pae_output',
                         action='store_true',
                         help='Default is False. THIS DOES NOT WORK YET')

    parser.add_argument('--plot_scores',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--crossover_percent',
                        default='0.2',
                        type=float,
                        help='Fraction of pool refilled by crossover.'
                        'Default is 0.2.')

    parser.add_argument('--substitution_insertion_deletion_weights',
                        default=None,
                        type=str,
                        help='Specify probability of substitutions, insertions, and deletions (in that order) during mutation. Default is 0.8,0.1,0.1')

    parser.add_argument('--mutation_percents',
                        default='0.125',
                        type=str,
                        help='Number of mutations made as a percentage of sequence length.'
                        'Default is 0.125 for every iteration. If more than one value is provided, number of iterations will be split evenly and assigned.')

    return parser

if __name__ == "__main__":
    parser = getEvoProParser()
    args = parser.parse_args(sys.argv[1:])
    print(args)
