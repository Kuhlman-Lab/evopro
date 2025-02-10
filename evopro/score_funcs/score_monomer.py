from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd

def score_plddt_only(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    results = results[0]
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    if rmsd_pdb:
        rmsdscore = score_rmsd_to_starting_selection(pdb, rmsd_pdb, chains1="A", chains2="A")
    else:
        rmsdscore=0
    score = -confscore/10 + rmsdscore
    
    return score, ((confscore,), (rmsdscore,)), pdb, results

def score_rmsd_to_starting_selection(pdb, path_to_starting, chains1="AB", chains2="AB", dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x[0] in chains1]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys() if x[0] in chains2]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting