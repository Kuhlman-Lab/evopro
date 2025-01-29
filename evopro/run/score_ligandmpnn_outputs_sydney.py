"""Script that runs AF2 (distributor, runtime-optimized) on a list of sequences and returns the scores and PDB files.
By default, the PDBs are scored for average plDDT score and sorted by that score. A custom scoring function can be provided instead."""

import sys, os
import importlib

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--af2_flags_file',
                        default='./af2.flags',
                        type=str,
                        help='Path to and name of af2.flags file.')
    
    parser.add_argument('--sequence_file',
                        default="sequences.csv",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--output_dir',
                        default="outputs",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--max_chains_length',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, overall max length of each chain separated by commas.')
    
    parser.add_argument('--custom_score',
                        default=None,
                        type=str,
                        help='Score function to use to generate scores for each prediction (optional).')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for AlphaFold2. Default is 1.')

    return parser


def run_af2_dist(args, score_func=None):
    
    # Load sequences
    seqs_list = []
    with open(args.sequence_file, "r") as sf:
        for line in sf:
            seq = [x for x in line.strip("\n").split(",") if x]
            seqs_list.append([seq])
    
    print(seqs_list)
    if args.max_chains_length:
        lengths = [[int(x) for x in args.max_chains_length.split(",")]]
    else:
        lengths = [[len(x) for x in seqs_list[0].split(",")]]
        
    if not os.path.isfile(args.af2_flags_file):
        raise ValueError("Invalid path to af2 flags file. Please provide a valid af2 flags file.")
    
    if not score_func:
        score_func = score_plddt_confidence_overall
    
    print("Compiling AF2 models for lengths:", lengths)
    dist = Distributor(args.n_workers, af2_init, args.af2_flags_file, lengths)
    
    results = dist.churn(seqs_list)
    
    print("done churning")
    dist.spin_down()
    
    seqs_and_scores = {}
    data = {}
    scores = []
    for seq, result in zip(seqs_list, results):
        s = ",".join(seq[0])
        score = score_func(result[0])
        scores.append(score)
        
        seqs_and_scores[s] = (score[0])
        data[s] = (score[1], score[-1], score[-2])
        
        #print(s, score[0])
    
    sorted_seqs_and_scores = sorted(seqs_and_scores.items(), key=lambda x:x[1])
    
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    #print(sorted_seqs_and_scores)
    for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
        pdb = data[seq[0]][-1]
        r = data[seq[0]][-2]
        with open(os.path.join(args.output_dir,"seq_" + str(i) + "_model_1.pdb"), "w") as pdbf:
            pdbf.write(str(pdb))
        compressed_pickle(os.path.join(args.output_dir, "seq_" + str(i) + "_result"), r)
        #print(seq[0], data[seq[0]][0])
    
    with open(os.path.join(args.output_dir, "seqs_and_scores.csv"), "w") as opf:
        for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
            opf.write(str(seq[0]) + "\t" + str(seqs_and_scores[seq[0]]) + "\t" + str(data[seq[0]][0]) + "\n")
            
               
if __name__ == "__main__":
    parser = getFlagParser()
    args = parser.parse_args(sys.argv[1:])
    
    score_func = None
    if args.custom_score:
        file = args.custom_score.split(" ")[0]
        function = args.custom_score.split(" ")[1]
        
        try:
            scorefile = file.rsplit("/", 1)
            scorepath = scorefile[0]
            scorefilename = scorefile[1].split(".")[0]

            sys.path.append(scorepath)
            mod = importlib.import_module(scorefilename)
            score_func = getattr(mod, function)
        except:
            raise ValueError("Invalid score function. Please provide a valid python file and function name within that file, separated by a space.")
        
    run_af2_dist(args, score_func=score_func)








import multiprocessing as mp
import importlib
import math
import sys, os
import pickle
import json
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.utils.plot_scores import plot_scores_general_dev
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import read_starting_seqs, create_new_seqs, create_new_seqs_mpnn
from evopro.user_inputs.inputs import getEvoProParser
from evopro.utils.plots import get_chain_lengths, plot_pae, plot_plddt
from evopro.utils.utils import compressed_pickle

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def run_genetic_alg_multistate(run_dir, af2_flags_file, score_func, startingseqs, poolsizes = [],
                               num_iter = 50, n_workers=1, mut_percents=None, single_mut_only=False, contacts=None, distance_cutoffs=None,
                               mpnn_temp="0.1", mpnn_version="s_48_020", skip_mpnn=[], mpnn_iters=None, 
                               mpnn_chains=None, bias_AA_dict=None, bias_by_res_dict=None, rmsd_pdb=None,
                               repeat_af2=True, af2_preds=[], crossover_percent=0.2, vary_length=0, sid_weights=[0.8,0.1,0.1],
                               write_pdbs=False, plot=[], conf_plot=False, write_compressed_data=True):

    num_af2=0
    
    lengths = []
    num_preds = len(af2_preds)
    print(af2_preds, num_preds)
    for chains in af2_preds:
        c = list(chains)
        lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])][0] for chain in c])

    print("Compiling AF2 models for lengths:", lengths)

    print("Initializing distributor")
    dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)

    scored_seqs = {}
    if repeat_af2:
        repeat_af2_seqs = {}
    curr_iter = 1
    seqs_per_iteration = []
    newpool = startingseqs

    while len(poolsizes) < num_iter:
        poolsizes.append(poolsizes[-1])

    output_dir = run_dir + "outputs/"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if write_pdbs:
        pdb_folder = output_dir + "pbz2_per_iter/"
        if not os.path.isdir(pdb_folder):
            os.makedirs(pdb_folder)
            
    if not mpnn_chains:
        mpnn_chains = [af2_preds[0]]

    #start genetic algorithm iteration
    while curr_iter <= num_iter:

        pool = newpool
        print("Current pool", pool)

        if curr_iter == 1:
            print("\nIteration 1: Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs(pool, 
                                   poolsizes[curr_iter-1], 
                                   crossover_percent=0, 
                                   mut_percent=mut_percents[curr_iter-1], 
                                   all_seqs = list(scored_seqs.keys()), 
                                   vary_length=vary_length,
                                   sid_weights=sid_weights,
                                   single_mut_only=single_mut_only)

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters and curr_iter not in skip_mpnn:
            print("\nIteration " + str(curr_iter) + ": refilling with ProteinMPNN.")

            pool = create_new_seqs_mpnn(pool, 
                                        scored_seqs, 
                                        poolsizes[curr_iter-1], 
                                        run_dir, 
                                        curr_iter, 
                                        all_seqs = list(scored_seqs.keys()), 
                                        af2_preds=af2_preds,
                                        mpnn_temp=mpnn_temp, 
                                        mpnn_version=mpnn_version, 
                                        mpnn_chains=mpnn_chains,
                                        bias_AA_dict=bias_AA_dict, 
                                        bias_by_res_dict=bias_by_res_dict)

        #otherwise refilling pool with just mutations and crossovers
        else:
            print("\nIteration " + str(curr_iter) + ": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs(pool, 
                                   poolsizes[curr_iter-1], 
                                   crossover_percent=crossover_percent, 
                                   mut_percent=mut_percents[curr_iter-1], 
                                   all_seqs = list(scored_seqs.keys()), 
                                   vary_length=vary_length,
                                   sid_weights=sid_weights,
                                   single_mut_only=single_mut_only)

        if repeat_af2:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in repeat_af2_seqs]
        else:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in scored_seqs]
        
        work_list_all = []
        for p in scoring_pool:
            for c in af2_preds:
                work_list_all.append([[p.jsondata["sequence"][chain] for chain in c]])
        
        print("work list", work_list_all)
        num_af2 += len(work_list_all)

        results = dist.churn(work_list_all)
        results_all = []
        for result in results:
            while type(result) == list:
                result = result[0]
            results_all.append(result)
            
        print("done churning")
        
        results_packed = [results_all[i:i+num_preds] for i in range(0, len(results_all), num_preds)]

        print("These should be equal", len(work_list_all), len(results_all), num_preds*len(results_packed))
        
        scores = []
        if not contacts:
            contacts=(None, None, None)
        if not distance_cutoffs:
            distance_cutoffs=(4, 4, 8)
        for results, dsobj in zip(results_packed, scoring_pool):
            if rmsd_pdb:
                scores.append(score_func(results, dsobj, contacts=contacts, distance_cutoffs=distance_cutoffs, rmsd_pdb=rmsd_pdb))
            else:
                scores.append(score_func(results, dsobj, contacts=contacts, distance_cutoffs=distance_cutoffs))


        print("scores before normalizing")
        for score in scores:
            print(score[1])
            print(f"overall score is {score[0]}")
        print("done writing old scores")


        scores_to_norm = []
        scores_no_norm = []
        for score in scores:
            if score[0] == 1000:
                scores_no_norm.append(score)
            else:
                scores_to_norm.append(score)

        if len(scores_to_norm) != 0:

            #here 0 is arbitrary since the number of score terms is the same across all scores
            num_score_terms = len(scores_to_norm[0][1]) 

            #each tuple in score_cols is a column of scores ex. monomer prediction for (design 1, design 2), pae for (design 1, design 2)
            score_cols = []
            for j in range(num_score_terms):
                #where i in the index of that design, 1 holds the score, j is the position of the score term in the score tuple, and each tuple begins with the actual score to return (other values in the tuple should be used for debugging, not scoring)
                col = tuple(scores_to_norm[i][1][j][0] for i in range(len(scores_to_norm)))
                score_cols.append(col)


            min_vals = [min(col) for col in score_cols]
            max_vals = [max(col) for col in score_cols]

            #TODO: add a flag (?) to ensure the weights can be specified
            # overall_score = -ptm_monomer + iptm_dimer - iptm_phos_dimer + penalty - bonus - rmsd
            weights = [-1, +4, -3, +1, -2, -5]

            for index, tup in enumerate(scores_to_norm):
                scores = tup[1]
                new_score = ()
                overall_score = 0
                for i in range(len(scores)):
                    value = scores[i][0]

                    if max_vals[i] != min_vals[i]:
                        scaled_val = (value - min_vals[i])/(max_vals[i] - min_vals[i])
                    else:
                        scaled_val = 0 #if min = max then all values = min = max so scaled_val is zero

                    if len(scores[i]) > 1:
                        new_tuple = ((scaled_val,) + tuple(scores[i][1:]),) #keeping the rest of the tuple as is
                    else:
                        new_tuple = ((scaled_val,),) #making a nested tuple to ensure the format is the same as the initial scores

                    new_score += new_tuple #creating a new scoring tuple to replace the current tuple
                    overall_score += scaled_val * weights[i]

                scores_to_norm[index] = [overall_score, new_score, tup[2], tup[3]]

        final_scores = scores_to_norm + scores_no_norm

        print("normalizing scores")
        for scores in final_scores:
            print(scores[1])
            print(f"overall score is {scores[0]}")
        print("done normalizing")
        
        #adding sequences and scores into the dictionary
        for score, results, dsobj in zip(final_scores, results_packed, scoring_pool):
            key_seq = dsobj.get_sequence_string()

            overall_scores = [score[0]]
            overall_scores = overall_scores + [x[0] for x in score[1]]
            score_split = score[1]
            pdbs = score[-2]
            score_all = (overall_scores, score_split)
            print(key_seq, overall_scores, score_split)
            
            if key_seq in scored_seqs and repeat_af2:
                scored_seqs[key_seq]["data"].append({"score": score_all, "pdb": pdbs, "result": results})

                sum_score = [0 for _ in overall_scores]
                for elem in scored_seqs[key_seq]["data"]:
                    sum_score = [x+y for x,y in zip(sum_score, elem["score"][0])]
                
                avg_score = [x/len(scored_seqs[key_seq]["data"]) for x in sum_score]
                scored_seqs[key_seq]["average"] = avg_score
            
            else:
                scored_seqs[key_seq] = {"data": [{"score": score_all, "pdb": pdbs, "result": results}]}
                scored_seqs[key_seq]["dsobj"] =  dsobj
                
                #set average score here
                scored_seqs[key_seq]["average"] = overall_scores
            
            #print(scored_seqs[key_seq]["average"], len(scored_seqs[key_seq]["data"]), scored_seqs[key_seq]["data"][-1])
       
        if repeat_af2:
            for key_seq in scored_seqs:
                if len(scored_seqs[key_seq]["data"]) >=5:
                    repeat_af2_seqs[key_seq] = scored_seqs[key_seq]["average"]
            print("repeat_af2_seqs", repeat_af2_seqs)
        
        #creating sorted list version of sequences and scores in the pool
        sorted_scored_pool = []
        #print(pool)
        for dsobj, j in zip(pool, range(len(pool))):
            key_seq = dsobj.get_sequence_string()
            print(key_seq, scored_seqs[key_seq]["average"])
            if repeat_af2:
                if key_seq in repeat_af2_seqs:
                    sorted_scored_pool.append((key_seq, repeat_af2_seqs[key_seq]))
                    #print(sorted_scored_pool)
                else:
                    sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["average"]))
                    #print(sorted_scored_pool)
            else:
                sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["average"]))            
            
            # if write_pdbs:
            #     print("Writing pdbs...")
            #     pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
            #     if type(pdbs) == str:
            #         pdbs = [pdbs]
            #     for pdb, chains in zip(pdbs, af2_preds):
            #         with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_chain"+str(chains)+".pdb", "w") as pdbf:
            #                     pdbf.write(str(pdb))


            if write_pdbs:
                print("Writing pdbs...")
                for i in range(len(pool)):
                    results = [scored_seqs[key_seq]["data"][0]["result"][x] for x in range(len(scored_seqs[key_seq]["data"][0]["result"]))]
                    for r, k in zip(results, range(len(results))):
                        compressed_pickle(os.path.join(pdb_folder, "seq_" + str(i) + "_iter_" + str(curr_iter) + "_result_" + str(k)), r)
        
        print("before sorting", sorted_scored_pool)
        sorted_scored_pool.sort(key = lambda x: x[1][0])
        print("after sorting", sorted_scored_pool)
        seqs_per_iteration.append(sorted_scored_pool)

        #writing log
        with open(output_dir+"runtime_seqs_and_scores.log", "a") as logf:
            logf.write("starting iteration " + str(curr_iter)+ " log\n")
            for elem in sorted_scored_pool:
                logf.write(str(elem[0]) + "\t" + str(elem[1]) + "\n")
        print("done writing runtime results")

        #create a new pool of only 50% top scoring sequences for the next iteration
        newpool_size = round(len(sorted_scored_pool)/2)
        if curr_iter < len(poolsizes):
            newpool_size = round(poolsizes[curr_iter]/2)
        else:
            pass
        newpool_seqs = []
        for sp in sorted_scored_pool[:newpool_size]:
            newpool_seqs.append(sp[0])
        print("newpool", newpool_seqs)

        newpool = []
        #pulling back DS objects for each sequence in new pool
        for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
            newpool.append(scored_seqs[key_seq]["dsobj"])
            #pdbs = scored_seqs[key_seq]["data"][0]["pdb"]

        curr_iter+=1

    with open(output_dir+"seqs_and_scores.log", "w") as logf:
        for tuplist, i in zip(seqs_per_iteration, range(len(seqs_per_iteration))):
            logf.write(">iteration"+str(i+1)+"\n")
            for val in tuplist:
                logf.write(str(val[0])+"\t"+str(val[1])+"\n")

    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        seq = key_seq.split(" ")
        pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
        if type(pdbs) == str:
            pdbs = [pdbs]
        result = scored_seqs[key_seq]["data"][0]["result"][0]
        results = [scored_seqs[key_seq]["data"][0]["result"][x] for x in range(len(scored_seqs[key_seq]["data"][0]["result"]))]
        
        if write_compressed_data:
            for r, k in zip(results, range(len(results))):
                compressed_pickle(os.path.join(output_dir, "seq_" + str(j) + "_result_" + str(k)), r)

        if conf_plot:
            # Plot confidence.
            Ls = get_chain_lengths(result)

        if conf_plot:
            # Plot confidence.
            Ls = get_chain_lengths(result)

            # Plot pAEs.
            ptm, iptm = None, None
            if 'ptm' in result:
                ptm = result['ptm']
            if 'iptm' in result:
                iptm = result['iptm']
            pae_fig = plot_pae(result['pae_output'][0], Ls, ptm=ptm, iptm=iptm)
            pae_fig.savefig(output_dir + "seq_" + str(j) + "_final_model_1_pae.png")

            # Plot pLDDTs.
            plddt_fig = plot_plddt(result['plddt'], Ls)
            plddt_fig.savefig(output_dir + "seq_" + str(j) + "_final_model_1_plddt.png")
        
        for pdb, chains in zip(pdbs, af2_preds):
            with open(output_dir + "seq_" + str(j) + "_final_model_1_chain"+str(chains)+".pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
    
    print("Number of AlphaFold2 predictions: ", num_af2)

    try:
        plotting = plot_scores_general_dev(plot, seqs_per_iteration, output_dir)
        print("plots created at", str(plotting))
    except:
        print("plotting failed")
    
    dist.spin_down()

if __name__ == "__main__":
    parser = getEvoProParser()
    args = parser.parse_args(sys.argv[1:])

    input_dir = args.input_dir
    if not input_dir.endswith("/"):
        input_dir = input_dir + "/"

    if input_dir not in ['', None]:
            onlyfiles = [f for f in os.listdir(input_dir) if os.path.isfile(
                         os.path.join(input_dir, f))]

            flagsfile=None
            resfile=None
            starting_seqs=None
            for filename in onlyfiles:
                if "flag" in filename and "af2" in filename:
                    flagsfile = filename
                    break
                if "flag" in filename and "of" in filename:
                    flagsfile = filename
                    break

            for filename in onlyfiles:
                if "residue" in filename and "spec" in filename:
                    resfile = filename
                    break

            if flagsfile is None:
                raise ValueError("Flags file for Alphafold runs not provided.")
            if resfile is None:
                raise ValueError("Please provide a residue specifications file.")

            #print(input_dir + resfile)
            dsobj1 = DesignSeq(jsonfile=input_dir + resfile)

            for filename in onlyfiles:
                if "starting" in filename and "seqs" in filename:
                    if "mut" in filename:
                        starting_seqs = read_starting_seqs(input_dir + filename, dsobj1, mut_only=True)
                        starting_seqs.append(dsobj1)
                    else:
                        starting_seqs = read_starting_seqs(input_dir + filename, dsobj1)
                        starting_seqs.append(dsobj1)
                    break
            if not starting_seqs:
                print("No starting sequences file provided. Random sequences will be generated.")

                starting_seqs = [dsobj1]
                starting_seqs = create_new_seqs(starting_seqs, args.pool_size, crossover_percent=0)

            poolfile=None
            if args.pool_size_variable:
                for filename in onlyfiles:
                    if "pool" in filename and "size" in filename:
                        poolfile = filename
                        break


    #get score function from flags file
    try:
        scorefile = args.score_file.rsplit("/", 1)
        scorepath = scorefile[0]
        scorefilename = scorefile[1].split(".")[0]

        sys.path.append(scorepath)
        mod = importlib.import_module(scorefilename)
        scorefunc = getattr(mod, args.score_func)
    except:
        raise ValueError("Invalid score function")

    if args.rmsd_func:
        rmsdfunc = getattr(mod, args.rmsd_func)
    else:
        rmsdfunc = None
    
    if args.rmsd_to_starting:
        try:
            func_name, path_to_starting = args.rmsd_to_starting.split(" ")
            rmsd_to_starting_func = getattr(mod, func_name)
        except:
            raise ValueError("Please provide name of function and path to pdb.")
    else:
        rmsd_to_starting_func=None
        path_to_starting=None

    mult = args.num_iter
    mut_percents = []
    if args.mutation_percents:
        mutation_percents = args.mutation_percents.split(",")
        if len(mutation_percents)>1:
            mult = math.ceil(args.num_iter/len(mutation_percents))
    else:
        mutation_percents = [0.125]
        
    for mut_percent in mutation_percents:
        mut_percents = mut_percents + [float(mut_percent) for x in range(mult)]
    
    mut_percents = mut_percents[:args.num_iter]

    if args.mpnn_iters:
        mpnn_iters = [int(x.strip()) for x in args.mpnn_iters.split(",")]
    else:
        mpnn_iters = []
        freq = int(args.mpnn_freq)
        for i in range(1, args.num_iter+1):
            if i % freq == 0:
                mpnn_iters.append(i)

    contacts=[None, None, None]
    distance_cutoffs = [4, 4, 8]
    if args.define_contact_area:
        c = args.define_contact_area.split(" ")
        contacts[0] = parse_mutres_input(c[0])
        if len(c)>1:
            distance_cutoffs[0] = int(c[1])
    if args.bonus_contacts:
        b = args.bonus_contacts.split(" ")
        contacts[1] = parse_mutres_input(b[0])
        if len(b)>1:
            distance_cutoffs[1] = int(b[1])
    if args.penalize_contacts:
        p = args.penalize_contacts.split(" ")
        contacts[2] = parse_mutres_input(p[0])
        if len(p)>1:
            distance_cutoffs[2] = int(p[1])

    mpnn_skips = []
    if args.skip_mpnn:
        mpnn_skips_temp = args.skip_mpnn.split(",")
        for elem in mpnn_skips_temp:
            if "-" in elem:
                r = elem.split("-")
                start = int(r[0])
                end = int(r[-1])
                for val in range(start, end+1):
                    mpnn_skips.append(val)
            else:
                mpnn_skips.append(elem)

    plot_style=[]
    if args.plot_scores_avg:
        plot_style.append("avg")
    if args.plot_scores_top:
        plot_style.append("top")
    if args.plot_scores_median:
        plot_style.append("median")
    pool_sizes = []
    if poolfile:
        with open(poolfile, "r") as pf:
            for lin in pf:
                pool_sizes.append(int(lin.strip()))
    else:
        for i in range(args.num_iter):
            pool_sizes.append(args.pool_size)

    print("Writing compressed data:", not args.dont_write_compressed_data)
    print("Writing pdb files every iteration:", args.write_pdbs)
    print("Varying protein length:", args.vary_length>0)
    print("Repeating AF2:", not args.no_repeat_af2)
    
    af2_preds = args.af2_preds.strip().split(",")
    if args.mpnn_chains:
        mpnn_chains = args.mpnn_chains.strip().split(",")
    else:
        mpnn_chains = None
        
    sid_weights = [0.8, 0.1, 0.1]
    if args.vary_length>0:
        if args.substitution_insertion_deletion_weights:
            sid_weights = [float(x) for x in args.substitution_insertion_deletion_weights.split(",")]
        
    #TODO
    bias_AA_dict = None
    """if os.path.isfile(args.mpnn_bias_AA):
        with open(args.mpnn_bias_AA, 'r') as json_file:
            json_list = list(json_file)
        for json_str in json_list:
            bias_AA_dict = json.loads(json_str)"""
        
    bias_by_res_dict = None
    if args.mpnn_bias_by_res:
        if os.path.isfile(args.mpnn_bias_by_res):
            with open(args.mpnn_bias_by_res, 'r') as json_file:
                bias_by_res_dict = json.load(json_file)
                
    #print(bias_by_res_dict)
                                   
    run_genetic_alg_multistate(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, mut_percents=mut_percents, single_mut_only=args.single_mut_only, contacts=contacts, distance_cutoffs=distance_cutoffs,
        mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, skip_mpnn=mpnn_skips, mpnn_iters=mpnn_iters, 
        mpnn_chains=mpnn_chains, bias_AA_dict=bias_AA_dict, bias_by_res_dict=bias_by_res_dict, rmsd_pdb=args.path_to_starting,
        repeat_af2=not args.no_repeat_af2, af2_preds = af2_preds, crossover_percent=args.crossover_percent, vary_length=args.vary_length, sid_weights=sid_weights,
        write_pdbs=args.write_pdbs, plot=plot_style, conf_plot=args.plot_confidences, write_compressed_data=not args.dont_write_compressed_data)
