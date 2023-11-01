import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import Rg


def score_diff_backbone_rG(pdb, diff_chain="B"):

    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist = [x for x in residues.keys() if x.startswith(diff_chain)]
    
    #doesnt work!
    #print(reslist)
    #rG = Rg(pdb, chnid=diff_chain)
    #print(rG)
    
    return (len(reslist), 0)
    