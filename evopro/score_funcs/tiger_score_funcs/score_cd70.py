import sys
sys.path.append('/proj/kuhl_lab/alphafold/alphafold/')
sys.path.append('/proj/kuhl_lab/folddesign/folddesign/score_funcs/')
from folddesign.score_funcs.score_funcs import score_contacts, score_confidence_pairs, score_confidence_lists, score_confidence_residues, get_rmsd
from folddesign.utils.pdb_parser import get_coordinates_pdb
from folddesign.score_funcs.score_funcs_of import score_conf_of, score_confidence_pairs_of, score_confidence_interface_of

# Tiger's scorefunctions: /pine/scr/t/i/tigerz/kuhlman_lab/folddesign
# Copy of Tiger's scorefunctions for running:  /proj/kuhl_lab/folddesign/folddesign/score_funcs/tiger_score_funcs/

def score_cd70_dimer_interface(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    # residues contacting CD27 in either CD70 chain and some adjacent residues
    # note need to renumber from 1 (substract 53)
    # original numbering: 61, 76, 79-83, 111-116, 124-125, 127, 129, 132, 135, 137-142, 144, 146, 148, 151, 160, 162, 165, 170, 172, 178-183

    # using
    # 61, 76-83, 111-116, 124-151, 160-165, 170-172, 178-183
    # renumbered as
    # 8, 23-30, 58-63, 71-98, 107-112, 117-119, 125-130
    CD70_contacts = {23,24,25,26,27,28,29,30,58,59,60,61,62,63,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,107,108,109,110,111,112,117,118,119,125,126,127,128,129,130}
    
    # create list of residue to look at contacts
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in CD70_contacts]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in CD70_contacts]
    reslist3 = [x for x in residues.keys() if x.startswith("C")]
    contacts_CD70_A, contactscore_CD70_A = score_contacts(pdb, reslist1, reslist2, fil = False)
    contacts_CD70_B, contactscore_CD70_B = score_contacts(pdb, reslist1, reslist2, fil = False)
    confidencescore_A = score_confidence_lists(results, reslist1, reslist2, resindices, fil = False)
    score = -contactscore_CD70_A*50 - contactscore_CD70_B*50 + confidencescore
    return score, (contactscore, confidencescore), contacts, pdb, results
# end funcion

def score_cd70_trimer_interface(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)

    #print(residues)
    # residues contacting CD27 in either CD70 chain and some adjacent residues
    # calculated using PyMol interface residues: https://pymolwiki.org/index.php/InterfaceResidues
    # note need to renumber from 1 (substract 53)
    # original numbering: 61, 76, 79-83, 111-116, 124-125, 127, 129, 132, 135, 137-142, 144, 146, 148, 151, 160, 162, 165, 170, 172, 178-183

    # using: 61, 76-83, 111-116, 124-151, 160-165, 170-172, 178-183
    # renumbered as
    # 8, 23-30, 58-63, 71-98, 107-112, 117-119, 125-130
    CD70_contacts = {23,24,25,26,27,28,29,30,58,59,60,61,62,63,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,107,108,109,110,111,112,117,118,119,125,126,127,128,129,130}
    
    # analyzing binding interface contacts
    # residue keys are in format chain_number_residue, like B_Tyr_50
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in CD70_contacts]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in CD70_contacts]
    reslist3 = [x for x in residues.keys() if x.startswith("C") and int(x.split("_")[2]) in CD70_contacts]
    reslist4 = [x for x in residues.keys() if x.startswith("D")]
    # score_contacts returns list of tuples (contacts) and an integer score (contactscore)
    contacts_CD70_A, contactscore_CD70_A = score_contacts(pdb, reslist1, reslist4, fil = False)
    contacts_CD70_B, contactscore_CD70_B = score_contacts(pdb, reslist2, reslist4, fil = False)
    contacts_CD70_C, contactscore_CD70_C = score_contacts(pdb, reslist3, reslist4, fil = False)
    # want to only count contacts between 2 CD70 moleulces forming binding interface
    # want to maximize contactscore
    contactscore_CD70 = contactscore_CD70_A + contactscore_CD70_B + contactscore_CD70_C - min(contactscore_CD70_A, contactscore_CD70_B, contactscore_CD70_C)

    # figure out which two chains are contacting the new protein 
    minContactScore = min(contactscore_CD70_A, contactscore_CD70_B, contactscore_CD70_C)
    if minContactScore == contactscore_CD70_A:
        contacts_CD70_chain1 = contacts_CD70_B
        contacts_CD70_chain2 = contacts_CD70_C
    elif minContactScore == contactscore_CD70_B:
        contacts_CD70_chain1 = contacts_CD70_A
        contacts_CD70_chain2 = contacts_CD70_C
    elif minContactScore == contactscore_CD70_C:
        contacts_CD70_chain1 = contacts_CD70_A
        contacts_CD70_chain2 = contacts_CD70_B
    contacts_CD70_all = contacts_CD70_chain1 + contacts_CD70_chain2 # merge two lists, each is list of tuples [(), (), (), ()...]

    # extract list of residues in CD70 chains and in novel binder
    reslist_contacts_CD70_extracted = []
    reslist_contacts_novel_binder_extracted = []
    for i in range(len(contacts_CD70_all)):
        reslist_contacts_CD70_extracted.append(i[0])
        reslist_contacts_novel_binder_extracted.append(i[1])

    # calculates confidence score of all interactions between 
    confidencescore = score_confidence_lists(results, reslist_contacts_CD70_extracted, reslist_contacts_novel_binder_extracted, resindices, fil = False) # check this

    # now calculate overall score
    # want to minimize overall score
    score = -contactscore_CD70*100 + confidencescore
    return score, (contactscore_CD70, confidencescore), contacts_CD70_all, pdb, results
# end function

def score_cd70_trimer_interface_contacts(results):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb, fil = False)
    # residues contacting CD27 in either CD70 chain and some adjacent residues
    # calculated using PyMol interface residues: https://pymolwiki.org/index.php/InterfaceResidues
    # note need to renumber from 1 (substract 53)
    # original numbering: 61, 76, 79-83, 111-116, 124-125, 127, 129, 132, 135, 137-142, 144, 146, 148, 151, 160, 162, 165, 170, 172, 178-183

    # using
    # 61, 76-83, 111-116, 124-151, 160-165, 170-172, 178-183
    # renumbered as
    # 8, 23-30, 58-63, 71-98, 107-112, 117-119, 125-130
    CD70_contacts = {23,24,25,26,27,28,29,30,58,59,60,61,62,63,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,107,108,109,110,111,112,117,118,119,125,126,127,128,129,130}
    
    # analyzing binding interface contacts
    # residue keys are in format chain_number_residue, like B_Tyr_50
    reslist1 = [x for x in residues.keys() if x.startswith("A") and int(x.split("_")[2]) in CD70_contacts]
    reslist2 = [x for x in residues.keys() if x.startswith("B") and int(x.split("_")[2]) in CD70_contacts]
    reslist3 = [x for x in residues.keys() if x.startswith("C") and int(x.split("_")[2]) in CD70_contacts]
    reslist4 = [x for x in residues.keys() if x.startswith("D")]
    # score_contacts returns list of tuples (contacts) and an integer score (contactscore)
    contacts_CD70_A, contactscore_CD70_A = score_contacts(pdb, reslist1, reslist4, fil = False)
    contacts_CD70_B, contactscore_CD70_B = score_contacts(pdb, reslist2, reslist4, fil = False)
    contacts_CD70_C, contactscore_CD70_C = score_contacts(pdb, reslist3, reslist4, fil = False)
    # want to only count contacts between 2 CD70 moleulces forming binding interface
    # want to maximize contactscore
    contactscore_CD70 = contactscore_CD70_A + contactscore_CD70_B + contactscore_CD70_C - min(contactscore_CD70_A, contactscore_CD70_B, contactscore_CD70_C)

    # figure out which two chains are contacting the new protein 
    minContactScore = min(contactscore_CD70_A, contactscore_CD70_B, contactscore_CD70_C)
    if minContactScore == contactscore_CD70_A:
        contacts_CD70_chain1 = contacts_CD70_B
        contacts_CD70_chain2 = contacts_CD70_C
    elif minContactScore == contactscore_CD70_B:
        contacts_CD70_chain1 = contacts_CD70_A
        contacts_CD70_chain2 = contacts_CD70_C
    elif minContactScore == contactscore_CD70_C:
        contacts_CD70_chain1 = contacts_CD70_A
        contacts_CD70_chain2 = contacts_CD70_B
    contacts_CD70_all = contacts_CD70_chain1 + contacts_CD70_chain2 # merge two lists, each is list of tuples [(), (), (), ()...]

    # now calculate overall score, want to minimize overall score
    score = -contactscore_CD70*100
    return score, (contactscore_CD70), contacts_CD70_all, pdb, results
# end function

def score_binder_rmsd(pdb1, pdb2, with_linker=False):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1, fil = False) # complex
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2, fil = False) # monomer
    if len(chains1)>1:
        print()
        print("Line 148, Length of chains is: " + str(chains1))
        print()
        reslist1 = [x for x in residues1.keys() if x.startswith("D")]
    else:
        reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2)
    return rmsd
