from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd, get_rmsd_superimposeall
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient

def score_seq_diff_monomers(results, diff_backbone):
    from alphafold.common import protein
    
    #print("Results scoring", len(results))
    #0 = complex, 1 = monomer A, 2 = monomer B
    pdbs = [protein.to_pdb(results[i]['unrelaxed_protein']) for i in range(len(results))]
    
    #first, score the complex
    chains0, residues0, resindices0 = get_coordinates_pdb(pdbs[0])
    reslist0 = [x for x in residues0.keys()]
    reslist0_A = [x for x in residues0.keys() if x.startswith("A")]
    reslist0_B = [x for x in residues0.keys() if x.startswith("B")]
    confscore_complex = score_plddt_confidence(results[0], reslist0, resindices0)
    rmsd_todiff_score = score_rmsd_to_starting(pdbs[0], diff_backbone)
    
    #then, score monomer A
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbs[1])
    reslist1 = [x for x in residues1.keys()]
    confscore_chainA = score_plddt_confidence(results[1], reslist1, resindices1)
    rmsd_score_chainA = get_rmsd(reslist1, pdbs[1], reslist0_A, pdbs[0], ca_only=True)
    
    #finally, score monomer B
    chains2, residues2, resindices2 = get_coordinates_pdb(pdbs[2])
    reslist2 = [x for x in residues2.keys()]
    confscore_chainB = score_plddt_confidence(results[2], reslist2, resindices2)
    rmsd_score_chainB = get_rmsd(reslist2, pdbs[2], reslist0_B, pdbs[0], ca_only=True)
    
    avg_confscore = (confscore_complex + confscore_chainA + confscore_chainB)/3
    score = -avg_confscore/10 + rmsd_todiff_score - rmsd_score_chainA - rmsd_score_chainB
    
    return score, (avg_confscore, confscore_complex, confscore_chainA, confscore_chainB, rmsd_todiff_score, rmsd_score_chainA, rmsd_score_chainB), pdbs, results

def score_seq_mutation_scanning(results, diff_backbone):
    from alphafold.common import protein

    #print("Results scoring", len(results))
    #0 = complex, 1 = monomer A, 2 = monomer B
    pdbs = [protein.to_pdb(results[i]['unrelaxed_protein']) for i in range(len(results))]

    #first, score the complex
    chains0, residues0, resindices0 = get_coordinates_pdb(pdbs[0])
    reslist0 = [x for x in residues0.keys()]
    reslist0_A = [x for x in residues0.keys() if x.startswith("A")]
    reslist0_B = [x for x in residues0.keys() if x.startswith("B")]
    confscore_complex = score_plddt_confidence(results[0], reslist0, resindices0)
    rmsd_todiff_score = score_rmsd_to_starting(pdbs[0], diff_backbone)

    #then, score monomer A
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbs[1])
    reslist1 = [x for x in residues1.keys()]
    confscore_chainA = score_plddt_confidence(results[1], reslist1, resindices1)
    #rmsd_score_chainA = get_rmsd(reslist1, pdbs[1], reslist0_A, pdbs[0], ca_only=True)

    #finally, score monomer B
    chains2, residues2, resindices2 = get_coordinates_pdb(pdbs[2])
    reslist2 = [x for x in residues2.keys()]
    confscore_chainB = score_plddt_confidence(results[2], reslist2, resindices2)
    rmsd_score_chainB = get_rmsd(reslist2, pdbs[2], reslist0_B, pdbs[0], ca_only=True)

    avg_confscore = (confscore_complex + confscore_chainA + confscore_chainB)/3
    score = -avg_confscore/10 + rmsd_todiff_score - rmsd_score_chainB

    return score, (avg_confscore, confscore_complex, confscore_chainA, confscore_chainB, rmsd_todiff_score, rmsd_score_chainB), pdbs, results

def score_seq_diff_dimer(results, diff_backbone):
    from alphafold.common import protein

    #print("Results scoring", len(results))
    #0 = complex, 1 = monomer A, 2 = monomer B
    pdbs = [protein.to_pdb(results[i]['unrelaxed_protein']) for i in range(len(results))]

    #first, score the complex
    chains0, residues0, resindices0 = get_coordinates_pdb(pdbs[0])
    reslist0 = [x for x in residues0.keys()]
    reslist0_A = [x for x in residues0.keys() if x.startswith("A")]
    reslist0_B = [x for x in residues0.keys() if x.startswith("B")]
    confscore_complex = score_plddt_confidence(results[0], reslist0, resindices0)
    rmsd_todiff_score = score_rmsd_to_starting(pdbs[0], diff_backbone)

    #then, score monomer A
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbs[1])
    reslist1 = [x for x in residues1.keys()]
    confscore_chainA = score_plddt_confidence(results[1], reslist1, resindices1)
    rmsd_score_chainA = get_rmsd(reslist1, pdbs[1], reslist0_A, pdbs[0], ca_only=True)

    #finally, score monomer B
    chains2, residues2, resindices2 = get_coordinates_pdb(pdbs[2])
    reslist2 = [x for x in residues2.keys()]
    confscore_chainB = score_plddt_confidence(results[2], reslist2, resindices2)
    rmsd_score_chainB = get_rmsd(reslist2, pdbs[2], reslist0_B, pdbs[0], ca_only=True)

    avg_confscore = (confscore_complex + confscore_chainA + confscore_chainB)/3
    score = -avg_confscore/10 + rmsd_todiff_score + rmsd_score_chainA + rmsd_score_chainB

    return score, (avg_confscore, confscore_complex, confscore_chainA, confscore_chainB, rmsd_todiff_score, rmsd_score_chainA, rmsd_score_chainB), pdbs, results

def score_seq_diff_monomers_3(results, diff_backbone):
    from alphafold.common import protein

    #print("Results scoring", len(results))
    #0 = complex, 1 = monomer A, 2 = monomer B
    pdbs = [protein.to_pdb(results[i]['unrelaxed_protein']) for i in range(len(results))]

    #first, score the complex
    chains0, residues0, resindices0 = get_coordinates_pdb(pdbs[0])
    reslist0 = [x for x in residues0.keys()]
    reslist0_A = [x for x in residues0.keys() if x.startswith("A")]
    reslist0_B = [x for x in residues0.keys() if x.startswith("B")]
    reslist0_C = [x for x in residues0.keys() if x.startswith("C")]
    confscore_complex = score_plddt_confidence(results[0], reslist0, resindices0)
    rmsd_todiff_score = score_rmsd_to_starting(pdbs[0], diff_backbone)

    #then, score monomer A
    chains1, residues1, resindices1 = get_coordinates_pdb(pdbs[1])
    reslist1 = [x for x in residues1.keys()]
    confscore_chainA = score_plddt_confidence(results[1], reslist1, resindices1)
    rmsd_score_chainA = get_rmsd(reslist1, pdbs[1], reslist0_A, pdbs[0], ca_only=True)

    #finally, score monomer B
    chains2, residues2, resindices2 = get_coordinates_pdb(pdbs[2])
    reslist2 = [x for x in residues2.keys()]
    confscore_chainB = score_plddt_confidence(results[2], reslist2, resindices2)
    rmsd_score_chainB = get_rmsd(reslist2, pdbs[2], reslist0_B, pdbs[0], ca_only=True)

    #finally, score monomer B
    chains3, residues3, resindices3 = get_coordinates_pdb(pdbs[3])
    reslist3 = [x for x in residues3.keys()]
    confscore_chainC = score_plddt_confidence(results[3], reslist3, resindices3)
    rmsd_score_chainC = get_rmsd(reslist3, pdbs[3], reslist0_C, pdbs[0], ca_only=True)

    avg_confscore = (confscore_complex + confscore_chainA + confscore_chainB + confscore_chainC)/3
    score = -avg_confscore/10 + rmsd_todiff_score - rmsd_score_chainC

    return score, (avg_confscore, confscore_complex, confscore_chainA, confscore_chainB, rmsd_todiff_score, rmsd_score_chainA, rmsd_score_chainB), pdbs, results


def score_seq_diff_pae(results, diff_backbone, binder_chain="A"):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x not in reslist1]
    #print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting(pdb, diff_backbone)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2

    score = -pae_per_contact -confscore/10 + rmsdscore

    return score, (pae_per_contact, confscore, rmsdscore), pdb, results

def score_seq_diff(results, diff_backbone):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting_selection(pdb, diff_backbone, chains1="AB", chains2="AB")
    score = -confscore/10 + rmsdscore
    
    return score, (confscore, rmsdscore), pdb, results

def score_seq_diff_binderB(results, diff_backbone):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys() if x.startswith("B")]
    print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting_binder(pdb, diff_backbone, binder_chain="B")
    score = -confscore/10 + rmsdscore
    
    return score, (confscore, rmsdscore), pdb, results

def score_seq_diff_binderA(results, diff_backbone):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting_binder(pdb, diff_backbone, binder_chain="A")
    score = -confscore/10 + rmsdscore
    
    return score, (confscore, rmsdscore), pdb, results

def score_seq_diff_binder_pae(results, diff_backbone, binder_chain="A"):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print("Result keys: ")
    print(results.keys())
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x not in reslist1]
    print("Reslist1")
    print(reslist1)
    print("Reslist2")
    print(reslist2)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting_binder(pdb, diff_backbone, binder_chain=binder_chain)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -pae_per_contact -confscore/10 + rmsdscore
    print(f"{score}, {pae_per_contact}, {confscore}, {rmsdscore}")
    return score, (pae_per_contact, confscore, rmsdscore), pdb, results

def score_rmsd_to_starting_binder(pdb, path_to_starting, dsobj=None, binder_chain="A"):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]
    reslist1_2 = [x for x in residues0.keys() if x.startswith(binder_chain)]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    reslist2_2 = [x for x in residues.keys() if x.startswith(binder_chain)]

    #rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)
    rmsd_to_starting = get_rmsd_superimposeall(reslist1, reslist1_2, pdb_string_starting, reslist2, reslist2_2, pdb, ca_only=True, translate=True)

    return rmsd_to_starting

def score_rmsd_to_starting(pdb, path_to_starting, dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys()]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting

def score_rmsd_to_starting_selection(pdb, path_to_starting, chains1="AB", chains2="AB", dsobj=None):

    with open(path_to_starting, 'r') as f:
        pdb_string_starting = f.read()
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb_string_starting)
    reslist1 = [x for x in residues0.keys() if x[0] in chains1]

    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys() if x[0] in chains2]

    rmsd_to_starting = get_rmsd(reslist1, pdb_string_starting, reslist2, pdb, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting

def score_seq_diff_binder_contacts(results, diff_backbone, binder_chain="A"):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    #print("Result keys: ")
    #print(results.keys())

    reslist1 = [x for x in residues.keys() if x.startswith("A")] # chain A binders
    reslist2 = [x for x in residues.keys() if x not in reslist1] # other chains

    confscore = score_plddt_confidence(results, reslist1, resindices) # pLDDT chain A
    confscore2 = score_plddt_confidence(results, reslist2, resindices) # pLDDT other chains (chain B)
    rmsdscore = score_rmsd_to_starting_binder(pdb, diff_backbone, binder_chain=binder_chain)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore -confscore/10 + rmsdscore

    ptm = results['ptm']
    ptm = float(ptm)

    #print(results['iptm'])
    if "iptm" in results:
        iptm = results['iptm']
        iptm = float(iptm)
    else:
        iptm = None

    print("Score, contactscore, pLDDT, rmsd, pae/contact, num_contacts, ptm, iptm, B_pLDDT")
    if iptm is not None:
        print(f"{score:.3f}, {contactscore:.3f}, {confscore:.3f}, {rmsdscore:.3f}, {pae_per_contact:.3f}, {num_contacts}, {ptm:.2f}, {iptm:.3f}, {confscore2:.3f}")
    else:
        print(f"{score:.2f}, {contactscore:.2f}, {confscore:.2f}, {rmsdscore:.2f}, {pae_per_contact:.3f}, {num_contacts}, {ptm:.2f}, {iptm}, {confscore2:.3f}, ")
    return score, (contactscore, confscore, rmsdscore, pae_per_contact, num_contacts, ptm, iptm, confscore2), pdb, results

def score_seq_diff_binder_contacts_pS(results, diff_backbone, binder_chain="A"):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys() if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x not in reslist1]

    confscore = score_plddt_confidence(results, reslist1, resindices)
    confscore2 = score_plddt_confidence(results, reslist2, resindices)
    rmsdscore = score_rmsd_to_starting_binder(pdb, diff_backbone, binder_chain=binder_chain)
    contact_list, contactscore = score_contacts_pae_weighted_efficient(results, pdb, reslist1, reslist2)
    num_contacts = len(contact_list)
    pae_per_contact = 0
    if num_contacts > 0:
        pae_per_contact = (70.0-(70.0*contactscore)/num_contacts)/2
        
    score = -contactscore -confscore/10 - confscore2/10 + rmsdscore

    #print(results['iptm'])
    if "iptm" in results:
        iptm = results['iptm']
        iptm = float(iptm)
    else:
        iptm = None
    
    ptm = results['ptm']
    ptm = float(ptm)

    print("Score, contactscore, pLDDT, rmsd, pae/contact, num_contacts, ptm, iptm, B_pLDDT")
    if iptm is not None:
        print(f"{score:.2f}, {contactscore:.2f}, {confscore:.2f}, {rmsdscore:.2f}, {pae_per_contact:.3f}, {num_contacts}, {ptm:.2f}, {iptm:.2f}, {confscore2:.3f}, ")
    else:
        print(f"{score:.2f}, {contactscore:.2f}, {confscore:.2f}, {rmsdscore:.2f}, {pae_per_contact:.3f}, {num_contacts}, {ptm:.2f}, {iptm}, {confscore2:.3f}, ")
    return score, (contactscore, confscore, rmsdscore, pae_per_contact, num_contacts, ptm, iptm, confscore2), pdb, results


if __name__=="__main__":
    print("no main functionality")
