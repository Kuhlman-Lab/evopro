from evopro.utils.pdb_parser import get_coordinates_pdb
from evopro.score_funcs.score_funcs import score_contacts_pae_weighted, score_plddt_confidence, get_rmsd
import math
from evopro.utils.pdb_parser import get_seq_from_pdb

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
   #pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
    #pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
    #score.append(score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
   # score.append(threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
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
            seqs = get_seq_from_pdb(pdb)
            if seqs[0] == seqs[1]:
                homotrimer_score = -score_binder_complex_3(result, dsobj, contacts)[0]
                score.append((homotrimer_score,))
                print('Sequences are the same')
            else:
                 score.append(score_binder_complex_3(result, dsobj, contacts)) #this is test
                 score.append(score_complex_confidence(result, dsobj)) #this is test
        elif len(chains)>1:
            score.append(score_binder_complex_2(result, dsobj, contacts)) #this is test
            score.append(score_complex_confidence(result,dsobj)) #this is test
        else:
            score.append(score_binder_monomer(result, dsobj))
   # pdb1 = protein.to_pdb(results[0]['unrelaxed_protein'])
   # pdb2 = protein.to_pdb(results[1]['unrelaxed_protein'])
   # pdb3 = protein.to_pdb(results[2]['unrelaxed_protein'])
   # pdb4 = protein.to_pdb(results[3]['unrelaxed_protein'])
   # score.append(score_binder_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
   # score.append(threshold_rmsd(pdb1, pdb2, binder_chain="B", dsobj=None))
   # score.append(score_binder_rmsd_BC(pdb2, pdb4, binder_chain="B", dsobj=None))
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
    #seqA = []
    #seqB = []
    
   # print("contacts are:", contacts)
    #for contact in contacts_list:
       # if contact[0].startswith("A"):
           # seqA.append(contact[0])
       # else:
           # seqB.append(contact[0])

    #for contact in contacts_list:
       # if contact[1].startswith("A"):
           # seqA.append(contact[1])
        #else:
           # seqB.append(contact[1])

   # seqA = list(set(seqA))
   # seqB = list(set(seqB))

    #print("Contacts on A chain:", seqA, "contacts on B chain:", seqB)


    #aalist_A = []
    #aalist_B = []

   # for aatype in seqA:
           # aalist_A.append(dsobj.sequence[aatype])

    #for aatype in seqB:
           # aalist_B.append(dsobj.sequence[aatype])

   # print("aalist_A:", aalist_A, "aalist_B:", aalist_B)

    #hydrophobic_residues = ["A", "V", "I", "L", "M", "F", "Y", "W"]

   # hydrophobic_number = 0

   # for resA in aalist_A:
       # if resA in hydrophobic_residues:
           # hydrophobic_number = hydrophobic_number + 1
           # print(resA, hydrophobic_number)
    #for resB in aalist_B:
       # if resB in hydrophobic_residues:
           # hydrophobic_number = hydrophobic_number + 1
           # print(resB, hydrophobic_number)

   # hydrophobic_score = hydrophobic_number/2

    #nonpolar_penalty = 0
    
   # if hydrophobic_score > 3:
      # nonpolar_penalty = hydrophobic_score
      # print("too many hydrophobics")

    score = -contactscore # + nonpolar_penalty*10
    print(score, (score, len(contacts), contactscore))
    return score, (score, len(contacts), contactscore)

def score_binder_complex_2(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = [x for x in contacts if x.startswith("B")]#added for proteinC
    reslist2 = [x for x in contacts if x.startswith("C")]
    #reslist2 = [x for x in residues.keys() if x.startswith("C")] #added for proteinC
    print("reslist1", reslist1) #added for proteinC
    print("reslist2", reslist2) #added for proteinC
    contacts_list, contactscoreBC = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    score = -contactscoreBC # + nonpolar_penalty*5, + to penalize, - to minimize energy
    print(score, (score, len(contacts), contactscoreBC))
    return score, (score, len(contacts), contactscoreBC)

def score_binder_complex_3(results, dsobj, contacts):
    from alphafold.common import protein
    pdb = protein.to_pdb(results['unrelaxed_protein'])
    chains, residues, resindices = get_coordinates_pdb(pdb)
    print(contacts)
    print(contacts[0])
    reslist1 = [x for x in contacts[0] if x.startswith("A")]
    reslist2 = [x for x in residues.keys() if x.startswith("B")]
    reslist3 = [x for x in contacts[0] if x.startswith("B")]
    reslist4 = [x for x in residues.keys() if x.startswith("C")]
    reslist5 = [x for x in contacts[0] if x.startswith("C")]
    reslist6 = [x for x in residues.keys() if x.startswith("A")]


    print("reslist1", reslist1)
    print("reslist2", reslist2)
    print("reslist3", reslist3)
    print("reslist4", reslist4)
    print("reslist4", reslist5)
    print("reslist4", reslist6)

    contacts_list, contactscore = score_contacts_pae_weighted(results, pdb, reslist1, reslist2, dsobj=None, first_only=False)
    contacts_list2, contactscore2 = score_contacts_pae_weighted(results, pdb, reslist3, reslist4, dsobj=None, first_only=False)
    contacts_list3, contactscore3 = score_contacts_pae_weighted(results, pdb, reslist5, reslist6, dsobj=None, first_only=False)
   

   #Penalizing contacts for dimer
   # penalty = 0
   # contact_cutoff=#
   # if contactscore2 > contact_cutoff:
       # penalty = math.pow(contactscore2 - contact_cutoff, 2)
       # print("There aren't enough contacts, calculating penalty")
       # print(penalty)
   
    score = -contactscore - contactscore2 -contactscore3 #+ penalty 
    print(score, (score, len(contacts), contactscore, contactscore2, contactscore3))
    return score, (score, len(contacts), contactscore, contactscore2, contactscore3)

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
    score = -confscore2*10 + plddt_potential #the bigger the difference betwenn plddt_cuttof and confscore the bigger the penalty
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
        print("Confidence of monomer is lower than 80.0. Calculating Confidence potential")
        print(plddt_potential)
    score = -confscore2*10 + plddt_potential
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
    rmsd_cutoff = 10

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

if __name__=="__main__":
    print("no main functionality")
