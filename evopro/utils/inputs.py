import json
import os
import argparse, math
from typing import Optional, Sequence

from utils.parsing_utils import *

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

def parse_arguments(conf):
    
    if not conf.flags.run_dir.endswith("/"):
        conf.flags.run_dir = conf.flags.run_dir + "/"
    
    #TODO: cd to run dir?
    
    if not os.path.exists(conf.structure_prediction.structure_pred_flags_file):
        raise ValueError("Flags file for structure prediction runs not found.")
    
    if not os.path.exists(conf.flags.residue_specs_file):
        raise ValueError("Residue specifications file not found.")
    
    if not conf.flags.starting_seqs_file or not os.path.exists(conf.flags.starting_seqs_file):
            print("No starting sequences file provided. Random mutations will be generated for starting pool.")

    mpnn_iters = []
    if conf.sequence_prediction.mpnn_freq:
        freq = int(conf.sequence_prediction.mpnn_freq)
        for i in range(1, conf.flags.num_iter+1):
            if i % freq == 0:
                mpnn_iters.append(i)
    
    elif conf.sequence_prediction.mpnn_iters:
        mpnn_iters = [int(x.strip()) for x in conf.sequence_prediction.mpnn_iters.split(",")]
    
    if conf.sequence_prediction.skip_mpnn:
        mpnn_skips = [int(x.strip()) for x in conf.sequence_prediction.skip_mpnn.split(",")]
        for elem in mpnn_skips:
            if elem in mpnn_iters:
                mpnn_iters.remove(elem)
    
    conf.sequence_prediction.mpnn_iters = mpnn_iters

    plot_styles=[]
    if conf.flags.plot_scores:
        styles = conf.flags.plot_scores.split(",")
        for style in styles:
            if style in ["average", "top", "median"]:
                plot_styles.append(style)
    
    conf.flags.plot_scores = plot_styles

    seeds = []

    if conf.structure_prediction.seed:
        seed_list = str(conf.structure_prediction.seed)
        seed_list = seed_list.split(",") 
        if len(seed_list)>1:
            mult = math.ceil(conf.flags.num_iter/len(seed_list))
        else:
            mult = conf.flags.num_iter 
    else:
        seed_list = ["random"]
        mult = conf.flags.num_iter
    for s in seed_list:
        seeds = seeds + [s for _ in range(mult)]
    
    seeds = seeds[:conf.flags.num_iter]

    conf.structure_prediction.seed = seeds

    #TODO: pool sizes without file
    pool_sizes = []
    if conf.flags.pool_sizes_file and os.path.exists(conf.flags.pool_sizes_file):
        with open(conf.flags.pool_sizes_file, "r") as pf:
            for lin in pf:
                pool_sizes.append(int(lin.strip()))
    else:
        for i in range(conf.flags.num_iter):
            pool_sizes.append(conf.flags.pool_size)
    
    while len(pool_sizes) < conf.flags.num_iter:
        pool_sizes.append(pool_sizes[-1])
    conf.flags.pool_sizes = pool_sizes
    
    conf.structure_prediction.structure_pred_chains = conf.structure_prediction.structure_pred_chains.strip().split(",")
    if not conf.sequence_prediction.sequence_pred_chains:
        conf.sequence_prediction.sequence_pred_chains = conf.structure_prediction.structure_pred_chains
    else:
        conf.sequence_prediction.sequence_pred_chains = conf.sequence_prediction.sequence_pred_chains.strip().split(",")

    bias_by_res_dict = None
    if conf.sequence_prediction.bias_by_res:
        if os.path.isfile(conf.sequence_prediction.bias_by_res):
            with open(conf.sequence_prediction.bias_by_res, 'r') as json_file:
                bias_by_res_dict = json.load(json_file)
    
    conf.sequence_prediction.bias_by_res = bias_by_res_dict
    
    mut_percents = []
    mult = conf.flags.num_iter
    if conf.flags.mutation_percents:
        mutation_percents = str(conf.flags.mutation_percents)
        mutation_percents = mutation_percents.split(",")
        if len(mutation_percents)>1:
            mult = math.ceil(conf.flags.num_iter/len(mutation_percents))
    else:
        mutation_percents = [0.125]
        
    for mut_percent in mutation_percents:
        mut_percents = mut_percents + [float(mut_percent) for x in range(mult)]

    mut_percents = mut_percents[:conf.flags.num_iter]

    conf.flags.mutation_percents = mut_percents

    # print("Writing compressed data:", arguments["write_compressed_data"])
    # print("Writing pdb files every iteration:", arguments["write_pdbs"])
    # print("Repeating AF2:", arguments["no_repeat_af2"])
    return conf