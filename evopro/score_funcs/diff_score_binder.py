from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_plddt_confidence, get_rmsd, get_rmsd_superimposeall

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
    rmsd_score_chainA = get_rmsd(reslist1, pdbs[1], reslist0_A, pdbs[0])
    
    #finally, score monomer B
    chains2, residues2, resindices2 = get_coordinates_pdb(pdbs[2])
    reslist2 = [x for x in residues2.keys()]
    confscore_chainB = score_plddt_confidence(results[2], reslist2, resindices2)
    rmsd_score_chainB = get_rmsd(reslist2, pdbs[2], reslist0_B, pdbs[0])
    
    avg_confscore = (confscore_complex + confscore_chainA + confscore_chainB)/3
    score = -avg_confscore/10 + rmsd_todiff_score - rmsd_score_chainA - rmsd_score_chainB
    
    return score, (avg_confscore, confscore_complex, confscore_chainA, confscore_chainB, rmsd_todiff_score, rmsd_score_chainA, rmsd_score_chainB), pdbs, results

def score_seq_diff(results, diff_backbone):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    reslist1 = [x for x in residues.keys()]
    print(reslist1)
    confscore = score_plddt_confidence(results, reslist1, resindices)
    rmsdscore = score_rmsd_to_starting(pdb, diff_backbone)
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

if __name__=="__main__":
    print("no main functionality")
