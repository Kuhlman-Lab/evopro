import argparse
import sys
from typing import Sequence

class FileArgumentParser(argparse.ArgumentParser):
    """Overwrites default ArgumentParser to better handle flag files."""

    def convert_arg_line_to_args(self, arg_line: str) -> Sequence[str]:
        # Read from files where each line contains a flag and its value, e.g.
        # '--flag value'.
        split_line = arg_line.split(' ')

        if len(split_line) > 1:
            return [split_line[0], ' '.join(split_line[1:])]
        else:
            return split_line


def getFDDParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run FDD."""

    parser = FileArgumentParser(description='FDD runner script that can take '
                                'command-line arguments or read from an '
                                'argument flag file.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--input_dir',
                        default='',
                        type=str,
                        help='Path to directory that contains input files.')

    parser.add_argument('--num_iter',
                        default='100',
                        type=int,
                        help='Number of iterations of genetic algorithm. Default is 100.')

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

    parser.add_argument('--multimer',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--write_plddt',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--write_pae',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--stabilize_binder',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--write_pdbs',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--mpnn_freq',
                        default='10',
                        type=int,
                        help='Protein MPNN is used to refill the pool once every _ iterations.'
                        'Default is 10.')

    parser.add_argument('--symmetry',
                        default='',
                        type=str,
                        help='PDB chains that should be symmetric, separated by commas')

    parser.add_argument('--len_binder',
                        default='0',
                        type=int,
                        help='Number of residues in binder, used to find binder when linker is included. Default is 0.')

    return parser

if __name__ == "__main__":
    parser = getFDDParser()
    args = parser.parse_args(sys.argv[1:])
    print(args)
