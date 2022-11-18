from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
import os
import subprocess
import shutil

def score_cd20_basic(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    reslist1 = [x for x in residues.keys()]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)

    reslist3 = [x for x in residues.keys()]
    confscore2 = score_confidence_residues(results, reslist3, resindices, fil = False)

    score = -contactscore*100 + confidencescore - confscore2*10
    return score, (contactscore, confidencescore, confscore2), contacts, pdb, results

def score_cd20_rmsd(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    #use specific contacts for dimer here
    contacts_list = [6, 9, 10, 13, 16, 17, 20, 21, 24, 27, 144, 147, 148, 151, 155]
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in contacts_list]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in contacts_list]

    reslist3 = [x for x in residues.keys()]

    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)

    confscore2 = score_confidence_residues(results, reslist3, resindices, fil = False)
    #create reslist6 here for all loop residues that we are calculating rmsd for
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>21 and int(x.split("_")[2])<44]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>91 and int(x.split("_")[2])<143]
    reslist6 = reslist4 + reslist5

    pdb_crystal = "/proj/kuhl_lab/folddesign/folddesign/score_funcs/zoey_score_funcs/cd20_dimer_rmsd.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)

    #create reslist5 for crystal strucrture comparison of the loop region
    reslist7 = [x for x in residues1.keys() if int(x.split("_")[2])>66 and int(x.split("_")[2])<89]
    reslist8 = [x for x in residues1.keys() if int(x.split("_")[2])>136 and int(x.split("_")[2])<188]
    reslist9 = reslist7 + reslist8

    rmsd = get_rmsd(reslist9, pdbstring, reslist6, pdb, ca_only=True)
    score = -contactscore*100 + confidencescore - confscore2*10 + rmsd*100
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_cd20_interface(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    #use specific contacts for dimer here
    contacts_list = [6, 9, 10, 13, 16, 17, 20, 21, 24, 27, 144, 147, 148, 151, 155]
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in contacts_list]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in contacts_list]

    reslist3 = [x for x in residues.keys()]

    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    pairs1 = [x[0] for x in contacts]
    pairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence_pairs(results, pairs1, pairs2, resindices, fil = False)

    if len(contacts) < 1:
        confidencescore = 0
    else:
        confidencescore = confidencescore/len(contacts)

    confscore2 = score_confidence_residues(results, reslist3, resindices, fil = False)
    #create reslist6 here for all loop residues that we are calculating rmsd for
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>21 and int(x.split("_")[2])<44]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>91 and int(x.split("_")[2])<143]
    reslist6 = reslist4 + reslist5

    pdb_crystal = "/proj/kuhl_lab/folddesign/folddesign/score_funcs/zoey_score_funcs/cd20_dimer_rmsd.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)

    #create reslist5 for crystal structure comparison of the loop region
    reslist7 = [x for x in residues1.keys() if int(x.split("_")[2])>66 and int(x.split("_")[2])<89]
    reslist8 = [x for x in residues1.keys() if int(x.split("_")[2])>136 and int(x.split("_")[2])<188]
    reslist9 = reslist7 + reslist8

    rmsd = get_rmsd(reslist9, pdbstring, reslist6, pdb, ca_only=True)
    score = -contactscore*100 + confidencescore*10 - confscore2*10 + rmsd*100
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_cd20_interface_new(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    #use specific contacts for dimer here
    contacts_list = [6, 9, 10, 13, 16, 17, 20, 21, 24, 27, 144, 147, 148, 151, 155]
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in contacts_list]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in contacts_list]

    reslist3 = [x for x in residues.keys()]

    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False)
    pairs1 = [x[0] for x in contacts]
    pairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence_pairs(results, pairs1, pairs2, resindices, fil = False)

    if len(contacts) < 1:
        confidencescore = 0
    else:
        confidencescore = confidencescore/len(contacts)

    confscore2 = score_confidence_residues(results, reslist3, resindices, fil = False)
    #create reslist6 here for all loop residues that we are calculating rmsd for
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>21 and int(x.split("_")[2])<44]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>91 and int(x.split("_")[2])<143]
    reslist6 = reslist4 + reslist5

    pdb_crystal = "/proj/kuhl_lab/folddesign/folddesign/score_funcs/zoey_score_funcs/cd20_dimer_rmsd.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)

    #create reslist5 for crystal structure comparison of the loop region
    reslist7 = [x for x in residues1.keys() if int(x.split("_")[2])>66 and int(x.split("_")[2])<89]
    reslist8 = [x for x in residues1.keys() if int(x.split("_")[2])>136 and int(x.split("_")[2])<188]
    reslist9 = reslist7 + reslist8

    rmsd = get_rmsd(reslist9, pdbstring, reslist6, pdb, ca_only=True)
    score = -contactscore*100 + confidencescore - confscore2*10 + rmsd*10000
    #score = -contactscore*100 + confidencescore - confscore2*10 + rmsd*1000
    #score = -contactscore*100 + confidencescore - confscore2*10 + rmsd*100
    #score = -contactscore*10000 + confidencescore*10 - confscore2*10 + rmsd*100
    return score, (contactscore, confidencescore, confscore2, rmsd), contacts, pdb, results

def score_cd20_interface_noloopplddt(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    #use specific contacts for dimer here
    contacts_list = [6, 9, 10, 13, 16, 17, 20, 21, 24, 27, 144, 147, 148, 151, 155]
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in contacts_list]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in contacts_list]

    contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, dist=8)
    #contacts, contactscore = score_contacts(pdb, reslist1, reslist2, fil = False, dist=4)
    pairs1 = [x[0] for x in contacts]
    pairs2 = [x[1] for x in contacts]
    confidencescore = score_confidence_pairs(results, pairs1, pairs2, resindices, fil = False)

    if len(contacts) < 1:
        confidencescore = 0
    else:
        confidencescore = confidencescore/len(contacts)

    #create reslist6 here for all loop residues that we are calculating rmsd for
    reslist4 = [x for x in residues.keys() if int(x.split("_")[2])>21 and int(x.split("_")[2])<44]
    reslist5 = [x for x in residues.keys() if int(x.split("_")[2])>91 and int(x.split("_")[2])<143]
    reslist6 = reslist4 + reslist5
    reslist3 = [x for x in residues.keys() if x not in reslist6]
    confscore3 = score_confidence_residues(results, reslist3, resindices, fil = False)

    pdb_crystal = "/proj/kuhl_lab/folddesign/folddesign/score_funcs/zoey_score_funcs/cd20_dimer_rmsd.pdb"
    with open(pdb_crystal, "r") as f:
        pdbstring = f.read()
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbstring, fil = False)

    #create reslist5 for crystal structure comparison of the loop region
    reslist7 = [x for x in residues1.keys() if int(x.split("_")[2])>66 and int(x.split("_")[2])<89]
    reslist8 = [x for x in residues1.keys() if int(x.split("_")[2])>136 and int(x.split("_")[2])<188]
    reslist9 = reslist7 + reslist8

    rmsd = get_rmsd(reslist9, pdbstring, reslist6, pdb, ca_only=True)
    score = -contactscore*100 + confidencescore - confscore3*10 + rmsd*10000
    return score, (contactscore, confidencescore, confscore3, rmsd), contacts, pdb, results

if __name__ == "__main__":
    print("no main functionality")
