from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of

def score_cd70_dimer_interface(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])<286]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>291]
    reslist3 = [x for x in residues.keys() if x.startswith("C")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results
