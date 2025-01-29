from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd, Rg
import math

def score_overall(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
            score.append(score_complex_confidence(result, dsobj))
        else:
            score.append(score_binder_monomer(result, dsobj))
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    score.append(score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    score.append(threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    #score.append(gyration_score(pdb2, dsobj=None))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_overall_3(results, dsobj, contacts=None, distance_cutoffs=None, rmsd_pdb=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>2:
            score.append(score_binder_complex_3(result, dsobj, contacts)) #this is test
            score.append(score_complex_confidence(result, dsobj)) #this is test
        elif len(chains)>1:
            score.append(score_binder_complex_2(result, dsobj, contacts)) #this is test
            score.append(score_complex_confidence(result,dsobj)) #this is test
        else:
            score.append(score_binder_monomer(result, dsobj))
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    pdb3 = protein.to_pdb(results[2]['unrelaxed_protein'])
    pdb4 = protein.to_pdb(results[3]['unrelaxed_protein'])
    score.append(score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    score.append(threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    score.append(score_binder_rmsd_BC(pdb2, pdb4, binder_chain="B", dsobj=None))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_overall_2(results, dsobj, contacts=None):
    from alphafold.common import protein
    print("Number of predictions being scored:", len(results))

    score=[]
    pdbs = []
    for result in results:
        #print(len(result))
        #print(result)
        pdb = protein.to_pdb(result['unrelaxed_protein'])
        pdbs.append(pdb)
        chains, residues, resindices = get_coordinates_pdb(pdb)
        if len(chains)>1:
            score.append(score_binder_complex(result, dsobj, contacts))
            score.append(score_complex_confidence(result, dsobj))
        else:
            score.append(score_binder_monomer(result, dsobj))
    pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    pdb2 = protein.to_pdb(results[2]['unrelaxed_protein'])
    score.append(score_potential_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
    overall_score = sum([x[0] for x in score])
    return overall_score, score, pdbs, results

def score_binder_complex(results, dsobj, contacts):
    from alphafold.common import protein
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
    # seqA = []
    #seqB = []
    
    #print("contacts are:", contacts)
    #for contact in contacts_list:
       # if contact[0].startswith("A"):
        #    seqA.append(contact[0])
       # else:
        #    seqB.append(contact[0])

    #for contact in contacts_list:
        #if contact[1].startswith("A"):
         #   seqA.append(contact[1])
       # else:
        #    seqB.append(contact[1])

    #seqA = list(set(seqA))
    #seqB = list(set(seqB))

    #print("Contacts on A chain:", seqA, "contacts on B chain:", seqB)


    #aalist_A = []
    #aalist_B = []

    #for aatype in seqA:
     #       aalist_A.append(dsobj.sequence[aatype])

    #for aatype in seqB:
     #       aalist_B.append(dsobj.sequence[aatype])

    #print("aalist_A:", aalist_A, "aalist_B:", aalist_B)

    #hydrophobic_residues = ["A", "V", "I", "L", "M", "F", "Y", "W"]

    #hydrophobic_number = 0

    #for resA in aalist_A:
     #   if resA in hydrophobic_residues:
      #      hydrophobic_number = hydrophobic_number + 1
       #     print(resA, hydrophobic_number)
    #for resB in aalist_B:
     #   if resB in hydrophobic_residues:
      #      hydrophobic_number = hydrophobic_number + 1
       #     print(resB, hydrophobic_number)

    #hydrophobic_score = hydrophobic_number/2

   # nonpolar_penalty = 0
    
   # if hydrophobic_score > 3:
    #   nonpolar_penalty = hydrophobic_score
     #  print("too many hydrophobics")

    score = -contactscore  #+ nonpolar_penalty*10
    print(score, (score, len(contacts), contactscore))
    return score, (score, len(contacts), contactscore)

def score_binder_complex_3(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = [x for x in contacts[0] if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    reslist3 = [x for x in contacts[0] if x.startswith("B")] #added for proteinC
    reslist4 = [x for x in residues.keys() if x.startswith("C")] #added for proteinC
    

    print("reslist1", reslist1)
    print("reslist2", reslist2)
    print("reslist3", reslist3) #added for proteinC
    print("reslist4", reslist4) #added for proteinC

    contacts_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    contacts_list2, contactscore2 = score_contacts_pae_weighted(results, pdb, reslist3, reslist4, dsobj=None, first_only=False)
    score = -contactscore - 5*contactscore2 # + nonpolar_penalty*5
    print(score, (score, len(contacts), contactscore, contactscore2))
    return score, (score, len(contacts), contactscore, contactscore2)

def score_binder_complex_2(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = [x for x in contacts[0] if x.startswith("B")] #added for proteinC
    reslist2 = [x for x in residues.keys() if x.startswith("C")] #added for proteinC
    print("reslist1", reslist1) #added for proteinC
    print("reslist2", reslist2) #added for proteinC
    contacts_list, contactscoreBC = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    score = 5*contactscoreBC # + nonpolar_penalty*5
    print(score, (score, len(contacts), contactscoreBC))
    return score, (score, len(contacts), contactscoreBC)

def score_complex_confidence(results, dsobj):
    from alphafold.common import protein
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of complex is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = -confscore2*2 + plddt_potential
    print(score)
    return score, (score, confscore2)

def score_binder_monomer(results, dsobj):
    from alphafold.common import protein
    spring_constant = 10.0
    plddt_cutoff = 80.0
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    reslist2 = [x for x in residues.keys()]
    confscore2 = score_plddt_confidence(results, reslist2, resindices, dsobj=dsobj, first_only=False)
    plddt_potential = 0
    if confscore2 < plddt_cutoff:
        plddt_potential = spring_constant*math.pow(plddt_cutoff - confscore2, 2)
        print("Confidence of monomer is lower than 80.0. Calculating Confidence potnetial")
        print(plddt_potential)
    score = -confscore2*2 + plddt_potential
    print(score)
    return score, (score, confscore2)

def score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    if len(reslist1) == len(reslist2):
        rmsd_binder = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    else:
        print("Trying to calculate RMSD between proteins of different length. Setting RMSD to 0.")
        rmsd_binder = 0

    return [-rmsd_binder*25]

def score_binder_rmsd_BC(pdb2, pdb4, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb2)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb4)
    reslist1 = [x for x in residues1.keys()]
    reslist2 = [x for x in residues2.keys() if x.startswith(binder_chain)]
    if len(reslist1) == len(reslist2):
        rmsd_binder = get_rmsd(reslist1, pdb2, reslist2, pdb4, dsobj=dsobj)
    else:
        print("Trying to calculate RMSD between proteins of different length. Setting RMSD to 0.")
        rmsd_binder = 0

    return [rmsd_binder*25]

def score_potential_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    chains1, residues1, resindices1 = get_coordinates_pdb(pdb1)
    chains2, residues2, resindices2 = get_coordinates_pdb(pdb2)
    reslist1 = [x for x in residues1.keys() if x.startswith(binder_chain)]
    reslist2 = [x for x in residues2.keys()]
    rmsd = get_rmsd(reslist1, pdb1, reslist2, pdb2, dsobj=dsobj)
    rmsd_potential = 0
    spring_constant = 10.0
    rmsd_cutoff_max = 10.0
    rmsd_cutoff_min = 5.0
    rmsd_weight = 25
    if rmsd < rmsd_cutoff_min:
        rmsd_potential = rmsd_weight*rmsd
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    elif rmsd <rmsd_cutoff_max:
        rmsd_potential = rmsd_weight*rmsd_cutoff_min
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    else:
        rmsd_potential = rmsd_weight*rmsd_cutoff_min - spring_constant*math.pow(rmsd_cutoff_max - rmsd,2)
        print("RMSD is ", rmsd, "rmsd_potential is ", rmsd_potential)
    return [-rmsd_potential]

def threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None):
    # to keep the de novo binder close in structure to the original binder
    # calculates Ca-only RMSD of de novo binder unbound vs to original scaffold applied to a flat-bottom quadratic potential
    spring_constant = 10.0 
    rmsd_cutoff = 25

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
    return [rmsd_potential*5]

def gyration_score(pdb2, dsobj=None):
    gyration = Rg(pdb2, chnid="A")
    gyration_limit = 15
    gyration_penalty = 0
    if gyration > gyration_limit:
        gyration_penalty = gyration*20

    return [gyration_penalty]

if __name__=="__main__":
    print("no main functionality")
