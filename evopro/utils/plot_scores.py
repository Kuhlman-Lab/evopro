import numpy as np 
import pandas as pd
import statistics
import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser

def plot_scores_stabilize_monomer_avg(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE AVERAGE OF TOP 50% over iterations and a .csv with the data
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
            v2 += score[1]
            v3 += score[2]
            if rmsd:
                v4 += score[3]
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
    plt.savefig(opdir + "scores_avg.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_avg.csv")

def plot_scores_stabilize_monomer_median(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE TOP SCORE over iterations and a .csv with the data
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
        v1=[]
        v2=[]
        v3=[]
        if rmsd:
            v4=[]
        for tup in tuplist[:cutoff]:
            seq = tup[0]
            score = tup[1]
            v1.append(score[0])
            v2.append(score[1])
            v3.append(score[2])
            if rmsd:
                v4.append(score[3])

        s1.append(statistics.median(v1))
        s2.append(statistics.median(v2))
        s3.append(statistics.median(v3))
        if rmsd:
            s4.append(statistics.median(v4))

    # calculate axes sizes
    if rmsd:
        min_score = min(s1 + s2 + s3 + s4)
        max_score = max(s1 + s2 + s3 + s4)
    else:
        min_score = min(s1 + s2 + s3)
        max_score = max(s1 + s2 + s3)
    # generate plot
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
    plt.savefig(opdir + "scores_median.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_median.csv")
    
def plot_scores_avg(seqs_and_scores, opdir, rmsd=True):
    pass

def plot_scores_stabilize_monomer_top(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE TOP SCORE over iterations and a .csv with the data
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
        tup = tuplist[0]
        score = tup[1]
        s1.append(score[0])
        s2.append(score[1])
        s3.append(score[2])
        if rmsd:
            s4.append(score[3])

    # calculate axes sizes
    if rmsd:
        min_score = min(s1 + s2 + s3 + s4)
        max_score = max(s1 + s2 + s3 + s4)
    else:
        min_score = min(s1 + s2 + s3)
        max_score = max(s1 + s2 + s3)
    # generate plot
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
    plt.savefig(opdir + "scores_top.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_top.csv")

def plot_scores_stabilize_monomer_avg_old(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE AVERAGE OF TOP 50% over iterations and a .csv with the data
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
    plt.savefig(opdir + "scores_avg.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_avg.csv")

def plot_scores_stabilize_monomer_median_old(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE TOP SCORE over iterations and a .csv with the data
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
        v1=[]
        v2=[]
        v3=[]
        if rmsd:
            v4=[]
        for tup in tuplist[:cutoff]:
            seq = tup[0]
            score = tup[1]
            v1.append(score[0])
            v2.append(score[1][0])
            v3.append(score[2][0])
            if rmsd:
                v4.append(score[3][0])

        s1.append(statistics.median(v1))
        s2.append(statistics.median(v2))
        s3.append(statistics.median(v3))
        if rmsd:
            s4.append(statistics.median(v4))

    # calculate axes sizes
    if rmsd:
        min_score = min(s1 + s2 + s3 + s4)
        max_score = max(s1 + s2 + s3 + s4)
    else:
        min_score = min(s1 + s2 + s3)
        max_score = max(s1 + s2 + s3)
    # generate plot
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
    plt.savefig(opdir + "scores_median.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_median.csv")

def plot_scores_stabilize_monomer_top_old(seqs_and_scores, opdir, rmsd=True):
    # outputs a .png file graphing EvoPro scorefunction scores OF THE TOP SCORE over iterations and a .csv with the data
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
        tup = tuplist[0]
        score = tup[1]
        s1.append(score[0])
        s2.append(score[1][0])
        s3.append(score[2][0])
        if rmsd:
            s4.append(score[3][0])

    # calculate axes sizes
    if rmsd:
        min_score = min(s1 + s2 + s3 + s4)
        max_score = max(s1 + s2 + s3 + s4)
    else:
        min_score = min(s1 + s2 + s3)
        max_score = max(s1 + s2 + s3)
    # generate plot
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
    plt.savefig(opdir + "scores_top.png")
    # create .csv file output
    iter_list = list(range(1, num_iters+1))
    if rmsd:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3, 'RMSD':s4}
    else:
        d = {'iter':iter_list, 'overall':s1, 'complex_contacts':s2, 'binder_pLDDT':s3}

    df = pd.DataFrame(data = d)
    df.to_csv(opdir + "scores_top.csv")
    
def plot_scores_general(plot, af2_preds, seqs_per_iteration, output_dir, rmsd=False):
    print(len(seqs_per_iteration), seqs_per_iteration[0])
    num_lists = len(seqs_per_iteration[0][0])
    print(num_lists)
    
def analyze_log(path_to_file, filename, prefix): 
    iters = []
    i=-1
    with open(path_to_file+filename,"r") as f:
        for lin in f:
            if "iter" in lin:
                i+=1
                iters.append([])
            else:
                l = lin.strip().split()
                iters[i].append(l)
    #print(len(iters))

    print(iters)
    score_cols = []

    for l, i in zip(iters, range(len(iters))):
        # l = list of lists for each iteration
        vals = []
        for seq in l:
            vals_seq = seq[1:] # strip out sequences from each line
            vals.append([float(x.strip('(),[]')) for x in vals_seq if x.strip('(),[]')]) # lists of list
            print(vals)

        cutoff = int(len(vals)/2) # num seqs per iter/2 

        num_vals = int(len(vals[0])/2)
        # create list of average values for that iter
        avg_vals = [0] * num_vals 
        # create list of lists, each list stores scores over iterations
        if len(score_cols) ==0:
            for j in range(num_vals):
                score_cols.append([]) 

        # loop over values in a scorefunction      
        for j in range(num_vals):
            # loop over sequences in an iter
            for v in vals[:cutoff]:
                print(v, j, num_vals)
                avg_vals[j] += v[j]
        
        if cutoff > 0:
            for j in range(num_vals):
                avg_vals[j] = avg_vals[j]/cutoff
                score_cols[j].append(avg_vals[j])
        else:
            pass
        #print()

    # create .csv file output
    # import columns in a pandas dataframe
    num_iters = len(score_cols[0])
    iter_list = list(range(1, num_iters+1))

    d = {'iter':iter_list}
    for j in range(num_vals):
        d["val" + str(j+1)] = score_cols[j]
    df = pd.DataFrame(data = d)

    # can customize column names
    index = path_to_file.index("/")
    df.to_csv(prefix  + "_parsed.csv")
    print(df)
    
def create_plot(plot_name, vals):
    flat_vals = [item for sublist in vals for item in sublist]
    plt.rcParams.update({'font.size': 25})
    plt.figure(figsize=(15, 12), dpi=1200)
    for l, i in zip(vals, range(len(vals))):
        plt.plot(l, label="Score "+str(i))
    plt.xlabel("Iteration number")
    plt.ylabel("Raw score value")
    plt.ylim([int(min(flat_vals) - abs(min(flat_vals)*0.1)), int(max(flat_vals) + abs(max(flat_vals)*0.1))])
    plt.legend()
    plt.savefig(plot_name)
    
def save_csv(csv_name, vals):
    iter_list = list(range(1, len(vals[0])+1))
    d = {'iter':iter_list}
    for l, i in zip(vals, range(len(vals))):
        d['score '+str(i)] = l
    print(d)
    df = pd.DataFrame(data = d)
    df.to_csv(csv_name)
    
def plot_scores_general_dev(plot, seqs_per_iteration, output_dir):
    return_list = []
    print("All scores per iteration", seqs_per_iteration)
    #print(len(seqs_per_iteration), seqs_per_iteration[0])
    num_scores = len(seqs_per_iteration[0][0][1])
    score_cols = []
    #for each score type, create a list
    for j in range(num_scores):
        score_cols.append([])
    print(len(score_cols), num_scores)
    
    top_seqs_per_iter = [x[:int(len(x)/2)] for x in seqs_per_iteration]
    top_scores_per_iter = []
    for l in top_seqs_per_iter:
        top_scores_per_iter.append([x[1] for x in l])
    print("Scores of top sequences per iteration", top_scores_per_iter, len(top_scores_per_iter), len(top_scores_per_iter[0]))
    
    for score_num in range(num_scores):
        for l in top_scores_per_iter:
            score_cols[score_num].append([])
            for score_list in l:
                score_cols[score_num][-1].append(score_list[score_num])
                
    print("all scores", score_cols)
    
    """for l, i in zip(top_scores_per_iter, range(len(top_scores_per_iter))):
        for scores_list in l:
            for j in range(num_scores):
                score_cols[j][-1].append(scores_list[j])
    print("all scores", score_cols)
    
    for l, i in zip(top_scores_per_iter, range(len(top_scores_per_iter))):
        score_cols[i].append([])
        for tup in l:
            print("tuple", tup)
            for j in range(num_scores):
                score_cols[i][-1].append(tup[1][j])
            print("all scores", score_cols)"""
            
    if "avg" in plot:
        score_avgs = []
        for score_lists in score_cols:
            score_avgs.append([sum(l)/len(l) for l in score_lists])
        print(score_avgs)
        create_plot(output_dir + "/score_avgs.png", score_avgs)
        save_csv(output_dir + "/score_avgs.csv", score_avgs)
        return_list.append("score_avgs.png")
    if "median" in plot:
        score_medians = []
        for score_lists in score_cols:
            score_medians.append([statistics.median(l) for l in score_lists])
        print(score_medians)
        create_plot(output_dir + "/score_medians.png", score_medians)
        save_csv(output_dir + "/score_medians.csv", score_medians)
        return_list.append("score_medians.png")
    if "top" in plot:
        score_tops = []
        for score_lists in score_cols:
            score_tops.append([l[0] for l in score_lists])
        print(score_tops)
        create_plot(output_dir + "/score_top.png", score_tops)
        save_csv(output_dir + "/score_top.csv", score_tops)
        return_list.append("score_top.png")
    
    return return_list

    
def plot_scores_main(plot, af2_preds, seqs_per_iteration, output_dir, rmsd=False):
    #plot_path = output_dir + "/" + plot + ".png"
    plot_path = None
    print(len(seqs_per_iteration), seqs_per_iteration[0])
    num_lists = len(seqs_per_iteration[0][0])
    print(num_lists)
    
    if "avg" in plot:
        if len(af2_preds)>1:
            if not rmsd:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir, rmsd=False)
                plot_path = output_dir + "/scores_avg.png"
            else:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir)
                plot_path = output_dir + "/scores_avg.png"
    if "top" in plot:
        if len(af2_preds)>1:
            if not rmsd:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir, rmsd=False)
                plot_path = output_dir + "/scores_top.png"
            else:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir)
                plot_path = output_dir + "/scores_top.png"

    if "median" in plot:
        if len(af2_preds)>1:
            if not rmsd:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir, rmsd=False)
                plot_path = output_dir + "/scores_median.png"
            else:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir)
                plot_path = output_dir + "/scores_median.png"
                
    return plot_path

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--input_dir',
                        type=str, default='./',
                        help='Directory which contains runtime log files. Default is "./".')
    parser.add_argument('--output_dir',
                        type=str,
                        default='.',
                        help='Directory where figures will be outputted. '
                        'Default is "./".')
    args = parser.parse_args()

    onlyfiles = [f for f in os.listdir(args.input_dir) if os.path.isfile(
        os.path.join(args.input_dir, f))]

    files = {}
    for filename in onlyfiles:
        extension = filename.split('.')[-1]
    logfiles = [f for f in os.listdir(args.input_dir) if f.endswith('.log')]
    print(logfiles)
    
    for logfile in logfiles:
        prefix = logfile.split('.')[0]
        analyze_log(args.input_dir, logfile, prefix)
        
    
    
