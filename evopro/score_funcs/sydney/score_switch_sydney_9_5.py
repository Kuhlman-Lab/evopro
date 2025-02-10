from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import os
from evopro.score_funcs.helper import get_helix_motif, get_weighted_solvent_exposed_return_ser_contacts, rmsd_to_original_pdb, check_serine_polar_contacts, get_hse_dim_vs_mono_score
from alphafold.common import protein

def score_overall(results:list, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    
    #note: results is in the same order as --af2_preds AB,C,CC so in this case results[0] is AB, results[1] is C, and results[2] is CC
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    pdb3 = protein.to_pdb(results[2]['unrelaxed_protein'])

    gyration = Rg(pdb2, chnid="A") 
    gyration2 = Rg(pdb1, chnid="A")
    gyration3 = Rg(pdb3, chnid="A")
    
    #monomer is the last prediction, we want to use the monomer to get the helix motif
    try:
        hel = get_helix_motif(pdb2)
    except:
        return 1000, ((1000, 1000, 1000), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,), (1000,)), pdbs, results
    
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            # complexscore1, complex_tuple = score_binder_complex(result, dsobj, contacts) #tuple is (score, len(contacts), contactscore)
            complexscore2, confidence_tuple = score_complex_confidence(result, dsobj, helix_motif=hel) #confidence_tuple is (confscore1, confscore2)
        else:
            monomer, monomer_tuple = score_binder_monomer(result, dsobj, helix_motif=hel) #monomer_tuple is (score, confscore2)
            print(f"The confidence score of the monomer is: {monomer}")
    
    #compare one monomer in the dimer to the monomeric prediction
    rmsd_switch = score_rmsd(pdb1, pdb2, chains1="A", chains2="A", dsobj=None)

    dist, serine_hse_mono = get_weighted_solvent_exposed_return_ser_contacts(pdb2, pdb1, helix_motif = hel)
    print(f"The distance is {dist}")

    hse_score =get_hse_dim_vs_mono_score(pdb2, pdb3)

    if check_serine_polar_contacts(pdb2):
        serine_hse_mono *= 5 #include bonus if serine is close to glutamate or aspartate
        print(f"Serine is in contact with glu or asp and its new score is {serine_hse_mono}")
    else:
        serine_hse_mono = serine_hse_mono

    overall_score = - monomer - rmsd_switch - (dist/10) - serine_hse_mono + (gyration) + (gyration2/2) + (gyration3/2) + hse_score #negative dist to maximize score term

    return overall_score, (monomer_tuple, (rmsd_switch,), (dist,), (serine_hse_mono,), (gyration,), (gyration2,), (gyration3,), (hse_score,)), pdbs, results

def score_complex_confidence(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
        
    reslist1 = ["A"+str(x+1) for x in helix_motif]
    reslist1 += ["B"+str(x+1) for x in helix_motif]
    reslist2 = [x for x in residues.keys() if x not in reslist1]
    
    confscore1 = score_plddt_confidence(results, reslist1, resindices)

    confscore2 = score_plddt_confidence(results, reslist2, resindices)

    #penalizing a high plddt value and especially penalizing a high plddt for the fifth helix
    score = confscore2 + (2*confscore1)
    print(score)
    
    #first value is added to overall score, second value is returned for debugging purposes
    return score, (confscore1, confscore2)

def score_binder_monomer(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)

    reslist = [x for x in residues.keys()]

    confscore = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False) 

    score = confscore
    print(score)
    return score, (score, confscore)

def score_rmsd(pdb1, pdb2, chains1="AB", chains2="AB", dsobj=None):
   
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues0.keys() if x[0] in chains1]

    chains, residues, resindices = get_coordinates_pdb(pdb1)
    reslist2 = [x for x in residues.keys() if x[0] in chains2]

    rmsd_to_starting = get_rmsd(reslist1, pdb2, reslist2, pdb1, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting