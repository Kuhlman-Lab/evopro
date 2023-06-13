# script located in: /Users/tiger/TIGER/Medicine/UNC/Kuhlman_Lab/scripts/plot_scores_evopro_general_clean.py
# and in: /nas/longleaf/home/tigerz/kuhlmanlab/scripts/plot_scores_evopro_general_clean.py
# run using: py /Users/tiger/TIGER/Medicine/UNC/Kuhlman_Lab/scripts/plot_scores_evopro_general_clean.py
#        or: python /nas/longleaf/home/tigerz/kuhlmanlab/scripts/plot_scores_evopro_general_clean.py

# INPUTS: path to an evopro output directory
#   ex. .../run8/outputs/, which contains the file .../run8/outputs/runtime_seqs_and_scores.log

# OUTPUTS: .csv file with score values over iterations
#   run using the function: analyze_log(path_to_file)

import numpy as np 
import pandas as pd

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
# end function

if __name__ == "__main__":
    # analyze evopro runs
    path = ""
    jobid = ""
    analyze_log(path)

    #for i in range(1, 9):
    #    print("Analyzing run: " + str(i))
    #    path = "run" + str(i) + "/outputs/"
    #    analyze_log(path)

    #print()