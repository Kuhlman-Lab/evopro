import json
import string
import sys, random
import omegaconf
import sys

import random

from utils.parsing_utils import get_coordinates_pdb
from utils.pdb_utils import change_chainid_pdb, append_pdbs

chain_names = list(string.ascii_uppercase)
all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def create_new_seq_mpnn(conf, pool):
    sys.path.append(conf.sequence_prediction.sequence_prediction_tool_location)
    from run_mpnn import main as run_mpnn
    
    mpnn_chains = conf.sequence_prediction.sequence_pred_chains
    if not mpnn_chains:
        mpnn_chains = conf.structure_prediction.structure_pred_chains[0]
    try:
        mpnn_conf = omegaconf.OmegaConf.load(conf.sequence_prediction.sequence_prediction_tool_conf)
    except:
        mpnn_conf = conf.sequence_prediction.sequence_prediction_tool_conf
    
    while True:
        dsobj = random.choice(pool.pool)
        mpnn_pdb = None
        
        #check if dsobj has just been added after latest prediction/scoring round
        if len(dsobj.score) == 0:
            continue
                
        #concat all preds into one pdb
        for pdb, pred in zip(dsobj.score[0][1], conf.structure_prediction.structure_pred_chains):
            if pred in mpnn_chains:

                if mpnn_pdb is None:
                    mpnn_pdb = pdb
                else:
                    chains, _, _ = get_coordinates_pdb(mpnn_pdb)
                    new_pdb = pdb
                    chains_new, _, _ = get_coordinates_pdb(pdb)
                    chains_mod = []
                    chain_ind = 0
                    while len(chains_mod)<len(chains_new):
                        if chain_names[chain_ind] not in chains:
                            chains_mod.append(chain_names[chain_ind])
                        chain_ind+=1
                    for cm, cn in zip(chains_mod, chains_new):
                        new_pdb = change_chainid_pdb(new_pdb, old_chain=cn, new_chain=cm)
                    
                    mpnn_pdb = append_pdbs(pdb, new_pdb)
        
        with open("mpnn.pdb", "w") as mpnnf:
            mpnnf.write(str(mpnn_pdb))
        jsondata = dsobj._get_json()
        
        new_seq = run_mpnn(mpnn_conf, design_run=True, json_data=jsondata, pdb_paths=["mpnn.pdb"])
                
        if type(new_seq) == str:
            return ",".join(new_seq.split("/"))
        else:
            return ",".join(new_seq)