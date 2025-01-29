import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import Rg


def score_diff_backbone_rG(pdbfile, diff_chain="A"):
    
    with open(pdbfile, "r") as f:
        pdb = f.read()
    
    rG = Rg(pdb, chnid=diff_chain)
    
    return rG
    