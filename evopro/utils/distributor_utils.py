import sys

from utils.distributor import Distributor
# sys.path.append('/proj/kuhl_lab/alphafold/run')

def get_lengths(dsobj, conf):
    lengths = []
    for chains in conf.structure_prediction.structure_pred_chains:
        c = list(chains)
        lengths.append([dsobj.get_lengths(chain) for chain in c])
    
    return lengths

def initialize_distributor(conf, lengths=None, pre_func=False):
    mode = conf.structure_prediction.structure_prediction_tool
    
    sys.path.append(conf.structure_prediction.structure_prediction_tool_location)
    
    if mode == "af2":
        from run_af2 import af2_init
        dist = Distributor(n_workers=conf.flags.num_gpus, f_init=af2_init, arg_file=conf.structure_prediction.structure_pred_flags_file, lengths=lengths, pre_func=pre_func)
    
    elif mode == "rf2":
        from run_rf2 import rf2_init
        dist = Distributor(n_workers=conf.flags.num_gpus, f_init=rf2_init, arg_file=conf.structure_prediction.structure_pred_flags_file, lengths=lengths, pre_func=pre_func)
    
    else:
        raise ValueError("Invalid structure prediction tool")
    
    return dist