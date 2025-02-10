from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import os
from evopro.score_funcs.helper import get_helix_motif, get_weighted_solvent_exposed
from alphafold.common import protein

def score_overall(results:list, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    
    #note: results is in the same order as --af2_preds AB,A so in this case results[0] is AB and results[1] is A
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    
    #monomer is the last prediction, we want to use the monomer to get the helix motif
    try:
        hel = get_helix_motif(pdb2)
    except:
        return 0, ((0, 0, 0), (0, 0), (0, 0), (0,), (0,), (0,) ), pdbs, results

    #hel = get_helix_motif(pdb2) #this was the original code without the try/except block
    
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            complexscore1, complex_tuple = score_binder_complex(result, dsobj, contacts) #tuple is (score, len(contacts), contactscore)
            complexscore2, confidence_tuple = score_complex_confidence(result, dsobj, helix_motif=hel) #confidence_tuple is (confscore1, confscore2)
        else:
            monomer, monomer_tuple = score_binder_monomer(result, dsobj, helix_motif=hel) #monomer_tuple is (score, confscore2)

    #compare the two monomers in the complex
    rmsd1 = score_rmsd(pdb1, pdb1, chains1="A", chains2="B", dsobj=None)
    
    #compare one monomer in the dimer to the monomeric prediction
    #negative because we want to maximize this score term
    rmsd2 = -score_rmsd(pdb1, pdb2, chains1="A", chains2="A", dsobj=None)
    #score.append(gyration_score(pdb2, dsobj=None))
    dist = get_weighted_solvent_exposed(pdb2, pdb1, helix_motif = hel)
    print("Dist", dist)
    
    overall_score = complexscore1 + complexscore2 + monomer + rmsd1 + rmsd2 - dist #negative dist to maximize score term
    # return overall_score, (complexscore1, complexscore2, monomer, rmsd1, rmsd2, dist), pdbs, results
    return overall_score, (complex_tuple, confidence_tuple, monomer_tuple, (rmsd1,), (rmsd2,), (dist,) ), pdbs, results

def score_binder_complex(results, dsobj, contacts):
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = contacts[0]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]

    print("reslist1", reslist1)
    print("reslist2", reslist2)
    
    contacts_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    #test starts here

    score = -contactscore  #+ nonpolar_penalty*10
    print(score, (score, len(contacts), contactscore))
    return score, (score, len(contacts), contactscore)

def score_complex_confidence(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
        
    reslist1 = ["A"+str(x+1) for x in helix_motif]
    reslist1 += ["B"+str(x+1) for x in helix_motif]
    reslist2 = [x for x in residues.keys() if x not in reslist1]
    
    confscore1 = score_plddt_confidence(results, reslist1, resindices)

    confscore2 = score_plddt_confidence(results, reslist2, resindices)
    
    score = -confscore2/10 + confscore1/10
    print(score)
    
    
    
    #first value is added to overall score, second value is returned for debugging purposes
    return score, (confscore1, confscore2)

def score_binder_monomer(results, dsobj, helix_motif=None):

    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    
    
        
    reslist1 = ["A"+str(x+1) for x in helix_motif]

    reslist2 = [x for x in residues.keys() if x not in reslist1]
    
    confscore1 = score_plddt_confidence(results, reslist1, resindices)
    
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    

    score = -confscore2/10 + confscore1/10
    print(score)
    return score, (score, confscore2)

"""def threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 5

    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]


    rmsd_to_starting = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)

    # apply flat-bottom quadratic-shaped potential function
    rmsd_potential = 0
    if rmsd_to_starting > rmsd_cutoff:
        rmsd_potential = spring_constant*math.pow(rmsd_to_starting - rmsd_cutoff, 2)
        print("Rmsd has reached the threshold of 5.0. Calculating RMSD potnetial")
        print(rmsd_potential)
    #add this as penalty on top of usual rmsd
    return [rmsd_potential*5]"""

def score_rmsd(pdb1, pdb2, chains1="AB", chains2="AB", dsobj=None):

   
    chains0, residues0, resindices0 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues0.keys() if x[0] in chains1]

    chains, residues, resindices = get_coordinates_pdb(pdb1)
    reslist2 = [x for x in residues.keys() if x[0] in chains2]

    rmsd_to_starting = get_rmsd(reslist1, pdb2, reslist2, pdb1, ca_only=True, dsobj=dsobj)

    return rmsd_to_starting

