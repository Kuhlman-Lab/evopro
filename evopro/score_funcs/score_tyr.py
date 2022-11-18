from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of

def score_sh3_redirect(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    if len(chains)>1:
        return score_sh3_tyr_contacts_design4(results)
    else:
        return score_sh3_conf(results)

def score_sh3_conf(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist, resindices, fil = False)
    score = -confscore2
    return score, (confscore2, ), pdb, results

def score_sh3_rmsd(pdb1, pdb2, with_linker=False):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False)
    if len(chains1)>1:
        reslist1 = [x for x in residues1.keys() if x.startswith("B")]
    else:
        reslist1 = [x for x in residues1.keys() if int(x.split("_")[2])>135]
    reslist2 = [x for x in residues2.keys()]
    rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2)

    rmsd = rmsd_binder
    if with_linker:
        with open(pd1, "r") as pd1f:
            pd1_pdb = pd1f.read()
        chains3, residues3, resindices3 = get_coordinates_pdb(pd1_pdb, fil = False)
        reslist3 = [x for x in residues3.keys()]
        reslist4 = [x for x in residues1.keys() if int(x.split("_")[2])<120]
        rmsd_pd1 = get_rmsd(reslist3, pd1_pdb, reslist4, pdb1)
        rmsd = rmsd + rmsd_pd1
    return rmsd

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

def score_sh3_tyr_contacts(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    pdb_crystal = "/nas/longleaf/home/mulikova/test/AF_custom/templates/4p6s.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171 and x.startswith("A")]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198 and x.startswith("A")]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>9 and int(x.split("_")[2])<12 and x.startswith("B")]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>14 and int(x.split("_")[2])<23 and x.startswith("B")]
    reslist6 = [x for x in residues.keys() if x.startswith("B")]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 
    reslist7 = [x for x in residues1.keys() if x.startswith("A")]
    print(reslist_ty, reslist_sh3, reslist6, reslist7)
    rmsd = get_rmsd(reslist7, pdbstring, reslist6, pdb, ca_only=True)
    contacts, contactscore = score_contacts(pdb, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_lists(results, reslist_ty, reslist_sh3, resindices, fil = False)
    confscore2 = score_confidence_residues(results, reslist6, resindices, fil = False)
    score = -contactscore*1000000000 + confidencescore/100
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_sh3_tyr_contacts_design1(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    pdb_crystal = "/nas/longleaf/home/mulikova/test/AF_custom/templates/4p6s.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171 and x.startswith("A")]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>2 and int(x.split("_")[2])<6 and x.startswith("B")]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>14 and int(x.split("_")[2])<16 and x.startswith("B")]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>16 and int(x.split("_")[2])<20 and x.startswith("B")]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>20 and int(x.split("_")[2])<30 and x.startswith("B")]
    reslist8 = [x for x in residues.keys() if int(x.split("_")[2])>35 and int(x.split("_")[2])<50 and x.startswith("B")]
    reslist9 = [x for x in residues.keys() if int(x.split("_")[2])>53 and int(x.split("_")[2])<58 and x.startswith("B")]
    reslist6 = [x for x in residues.keys() if x.startswith("B")]
    reslist_ty = reslist1 
    reslist_sh3 = reslist3 + reslist4 + reslist2 + reslist5 + reslist8 + reslist9
    reslist7 = [x for x in residues1.keys() if x.startswith("A")]
    rmsd = get_rmsd(reslist7, pdbstring, reslist6, pdb, ca_only=True)
    orient = [['A_GLN_95'], ['B_GLU_32']]
    contacts, contactscore = score_contacts(pdb, reslist_ty, reslist_sh3, fil = False, orient=orient, orient_dist=15)
    reslist10 = [x[0] for x in contacts]
    reslist11 = [x[1] for x in contacts]
    contact_confidence = score_confidence_pairs(results, reslist10, reslist11, resindices, fil = False)
    confidencescore = score_confidence_lists(results, reslist_ty, reslist_sh3, resindices, fil = False)
    confscore2 = score_confidence_residues(results, reslist6, resindices, fil = False)
    if len(contacts) == 0:
        score = -contactscore*10000 + rmsd*100
    else:
        score = -contactscore*10000 + contact_confidence*10000/len(contacts) + rmsd*1000
    return score, (contactscore, contact_confidence, confscore2, rmsd), contacts, pdb, results

def score_sh3_tyr_contacts_design4(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    pdb_crystal = "/nas/longleaf/home/mulikova/test/AF_custom/templates/4p6s.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171 and x.startswith("A")]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>2 and int(x.split("_")[2])<6 and x.startswith("B")]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>14 and int(x.split("_")[2])<16 and x.startswith("B")]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>16 and int(x.split("_")[2])<20 and x.startswith("B")]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>20 and int(x.split("_")[2])<30 and x.startswith("B")]
    reslist8 = [x for x in residues.keys() if int(x.split("_")[2])>35 and int(x.split("_")[2])<50 and x.startswith("B")]
    reslist9 = [x for x in residues.keys() if int(x.split("_")[2])>53 and int(x.split("_")[2])<58 and x.startswith("B")]
    reslist6 = [x for x in residues.keys() if x.startswith("B")]
    reslist_ty = reslist1
    reslist_sh3 = reslist3 + reslist4 + reslist2 + reslist5 + reslist8 + reslist9
    reslist7 = [x for x in residues1.keys() if x.startswith("A")]
    rmsd = get_rmsd(reslist7, pdbstring, reslist6, pdb, ca_only=True)
    orient = [['A_GLN_95'], ['B_GLU_32']]
    contacts, contactscore = score_contacts(pdb, reslist_ty, reslist_sh3, fil = False, orient=orient, orient_dist=15)
    confidencescore = score_confidence_lists(results, reslist_ty, reslist_sh3, resindices, fil = False)
    confscore2 = score_confidence_residues(results, reslist6, resindices, fil = False)
    score = -contactscore*10000 + confidencescore + rmsd*10000
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_sh3_tyr_contacts_2(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    pdb_crystal = "/nas/longleaf/home/mulikova/test/AF_custom/templates/4p6s.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171 and x.startswith("A")]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198 and x.startswith("A")]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>9 and int(x.split("_")[2])<12 and x.startswith("B")]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>14 and int(x.split("_")[2])<23 and x.startswith("B")]
    reslist6 = [x for x in residues.keys() if x.startswith("B")]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4
    reslist7 = [x for x in residues1.keys() if x.startswith("A")]
    print(reslist_ty, reslist_sh3, reslist6, reslist7)
    rmsd = get_rmsd(reslist7, pdbstring, reslist6, pdb, ca_only=True)
    orient = [['A_VAL_164'], ['B_PHE_19']]
    contacts, contactscore = score_contacts(pdb, reslist_ty, reslist_sh3, fil = False, orient=orient, orient_dist=15)
    confidencescore = score_confidence_lists(results, reslist_ty, reslist_sh3, resindices, fil = False)
    confscore2 = score_confidence_residues(results, reslist6, resindices, fil = False)
    score = -contactscore*1000 + confidencescore/100
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_sh3_tyr_contacts_withoutrmsd(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>26 and int(x.split("_")[2])<40]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>135 and int(x.split("_")[2])<143]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>317 and int(x.split("_")[2])<322]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>329 and int(x.split("_")[2])<339]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>341 and int(x.split("_")[2])<353]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdb, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_lists(results, reslist_ty, reslist_sh3, resindices, fil = False)
    score = -contactscore*100 + confidencescore*10
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

def score_sh3_tyr_of_longlinker(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>26 and int(x.split("_")[2])<40]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>135 and int(x.split("_")[2])<143]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>315 and int(x.split("_")[2])<320]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>330 and int(x.split("_")[2])<333]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>338 and int(x.split("_")[2])<353]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def score_sh3_tyr_of_longlinker_1(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>159 and int(x.split("_")[2])<170]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>329 and int(x.split("_")[2])<332]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>334 and int(x.split("_")[2])<336]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>338 and int(x.split("_")[2])<343]
    reslist_ty = reslist1
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def score_sh3_tyr_of_longlinker_3(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>329 and int(x.split("_")[2])<332]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>334 and int(x.split("_")[2])<336]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>338 and int(x.split("_")[2])<343]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*10000 + confidencescore/10
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def score_sh3_tyr_shortlinker(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>315 and int(x.split("_")[2])<318]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>320 and int(x.split("_")[2])<322]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>324 and int(x.split("_")[2])<329]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def score_sh3_tyr_10g(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>320 and int(x.split("_")[2])<323]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>325 and int(x.split("_")[2])<327]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>329 and int(x.split("_")[2])<334]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*10000 + confidencescore/10
    return score, (contactscore, confidencescore), contacts, pdbstr, results


def score_sh3_tyr_13g(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>308 and int(x.split("_")[2])<311]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>312 and int(x.split("_")[2])<322]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4
    print(reslist_ty, reslist_sh3)
    orient = [['A_VAL_163'], ['A_PHE_323']]
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False, orient=orient, orient_dist=15)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100000 + confidencescore/10
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def run2(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>315 and int(x.split("_")[2])<318]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>320 and int(x.split("_")[2])<322]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>324 and int(x.split("_")[2])<329]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*100 - confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def run3(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>315 and int(x.split("_")[2])<318]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>320 and int(x.split("_")[2])<322]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>324 and int(x.split("_")[2])<329]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*1000 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results

def run4(results, pdbstr):
    chains, residues, resindices = get_coordinates_pdb(pdbstr, fil = False)
    reslist1 = [x for x in residues.keys() if int(x.split("_")[2])>158 and int(x.split("_")[2])<171]
    reslist2 = [x for x in residues.keys() if int(x.split("_")[2])>184 and int(x.split("_")[2])<198]
    reslist3 = [x for x in residues.keys() if int(x.split("_")[2])>315 and int(x.split("_")[2])<318]
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>320 and int(x.split("_")[2])<322]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>324 and int(x.split("_")[2])<329]
    reslist_ty = reslist1 + reslist2
    reslist_sh3 = reslist3 + reslist4 + reslist5
    contacts, contactscore = score_contacts(pdbstr, reslist_ty, reslist_sh3, fil = False)
    confidencescore = score_confidence_interface_of(results, reslist_ty, reslist_sh3, resindices)
    score = -contactscore*1000 - confidencescore
    return score, (contactscore, confidencescore), contacts, pdbstr, results
