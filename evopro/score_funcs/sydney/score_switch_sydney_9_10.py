from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import os
from evopro.score_funcs.helper import get_helix_motif, rmsd_to_original_pdb, get_distance_helix3_4, polar_positive_pentalty, glu_asp_bonus
from alphafold.common import protein

def score_overall(results:list, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    print("Number of predictions being scored:", len(results))
    pdbs = []
    
    #note: results is in the same order as --af2_preds AB,C,CC so in this case results[0] is AB, results[1] is C, and results[2] is CC
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    ptm_monomer_array = results[1]['ptm']
    ptm_monomer = ptm_monomer_array.flat[0]
    iptm_phos_dimer_array = results[0]['iptm'] #want confidence in the binder contacts to be high when the 5th helix is disordered
    iptm_phos_dimer = iptm_phos_dimer_array.flat[0]
    iptm_dimer_array = results[2]['iptm'] #want confidence in the binder contacts to be low when the 5th helix is ordered
    iptm_dimer = iptm_dimer_array.flat[0]

    rmsd_to_starting = rmsd_to_original_pdb(pdb2, "original.pdb") #want the fold to be close to the starting structure

    for result in results:
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
    
    #monomer is the last prediction, we want to use the monomer to get the helix motif
    try:
        hel = get_helix_motif(pdb2)
    except:
        return 1000, ((1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results

    #dict_keys(['pae_output', 'ranking_confidence', 'plddt', 'structure_module', 'ptm', 'iptm', 'unrelaxed_protein'])

    dist = get_distance_helix3_4(pdb2, pdb1, hel)
    print("Dist", dist)

    penalty = polar_positive_pentalty(pdb2)
    bonus = glu_asp_bonus(pdb2, hel)

    #computing the overall score
    overall_score = -ptm_monomer + iptm_dimer - iptm_phos_dimer + dist + penalty - bonus + rmsd_to_starting

    return overall_score, ((ptm_monomer,), (iptm_dimer,), (iptm_phos_dimer,), (dist,), (penalty,), (bonus,), (rmsd_to_starting,)), pdbs, results