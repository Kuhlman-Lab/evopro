import math
import sys
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd

def score_seq_diff(results, dsobj, contacts=None, rmsd_pdb=None):
    from alphafold.common import protein
    result = results[0]
    pdb = protein.to_pdb(result['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(result, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting(pdb, rmsd_pdb, dsobj=dsobj)
    else:
        rmsdscore=0
    overall_score = -confscore/10 + rmsdscore
    score = [(confscore,),  (rmsdscore,)]
    
    return overall_score, score, [pdb], results

def score_rmsd_to_starting(pdb, path_to_starting, dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting