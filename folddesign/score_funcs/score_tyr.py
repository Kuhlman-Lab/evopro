from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues
from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of

def score_sh3_tyr(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])<286]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>291]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results

def score_sh3_tyr_of(results, pdbstr):    
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>26 and int(x.split("_")[2])<46]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>135 and int(x.split("_")[2])<143]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>298 and int(x.split("_")[2])<303]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>310 and int(x.split("_")[2])<320]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>322 and int(x.split("_")[2])<334]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results
