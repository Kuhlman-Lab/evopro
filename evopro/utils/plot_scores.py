import numpy as np 
import pandas as pd

def plot_scores_stabilize_monomer(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores over iterations and a .csv with the data
    # The score for each sequence is broken into 4 parts: 
    # 1) overall_score, 2) contact_score, 3) binder_only_confid (pLDDT), 4) RMSD_score   
    # where the overall score is the sum of 2), 3), and 4)
    s1 = []
    s2 = []
    s3 = []
    if rmsd:
        s4 = []

    for tuplist in seqs_and_scores:
        cutoff = int(len(tuplist)/2)
        v1=0
        v2=0
        v3=0
        if rmsd:
            v4=0
        for tup in tuplist[:cutoff]:
            seq = tup[0]
            score = tup[1]
            v1 += score[0]
            v2 += score[1][0]
            v3 += score[2][0]
            if rmsd:
                v4 += score[3][0]
        if cutoff > 0:
            # v1 = overall score, v2 = contact score, v3 = confidence score, v4 = RMSD
            # averaged using top half of scores
            v1 = v1/cutoff
            v2 = v2/cutoff
            v3 = v3/cutoff
            if rmsd:
                v4 = v4/cutoff
            s1.append(v1)
            s2.append(v2)
            s3.append(v3)
            if rmsd:
                s4.append(v4)

    # calculate axes sizes
    if rmsd:
        min_score = min(s1 + s2 + s3 + s4)
        max_score = max(s1 + s2 + s3 + s4)
    else:
        min_score = min(s1 + s2 + s3)
        max_score = max(s1 + s2 + s3)
    # generate plot
    import matplotlib.pyplot as plt
    num_iters = len(s1)
    plt.rcParams.update({'font.size': 25})
    plt.figure(figsize=(15, 12), dpi=1200)
    plt.plot(s1[:num_iters], color='k', label="Overall")
    plt.plot(s2[:num_iters], color='g', label= "Complex contacts" )
    plt.plot(s3[:num_iters], color='b', label="Binder (pLDDT)" )
    if rmsd:
        plt.plot(s4[:num_iters], color='r', label="RMSD" )
    plt.xlabel("Iteration number")
    plt.ylabel("Raw score value")
    plt.ylim([int(min_score - abs(min_score*0.1)), int(max_score + abs(max_score*0.1))])
    plt.legend()
    plt.savefig(opdir + "scores.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores.csv")
