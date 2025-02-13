from utils.parsing_utils import *

def example_custom_score(result, conf):
    #get inputs from yaml file
    pdb = result['pdb']
    print(result.keys())
    
    #get arguments from yaml file
    if conf.scoring.custom.custom_score_args:
        args = conf.scoring.custom.custom_score_args.split(",")
    
    #get residues to score
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    #calculate score
    score = 1000
    
    return score
