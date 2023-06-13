import multiprocessing as mp
import importlib
import math
import sys, os
from typing import Sequence, Union
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq, DesignSeqMSD
from evopro.utils.distributor import Distributor
from evopro.utils.plot_scores import plot_scores_general
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import read_starting_seqs, create_new_seqs, create_new_seqs_mpnn, create_new_seqs_henry, create_new_seqs_mpnn_henry
from evopro.user_inputs.inputs import getEvoProParser
from evopro.utils.plots import get_chain_lengths, plot_pae, plot_plddt
from evopro.utils.utils import compressed_pickle
from evopro.utils.pdb_parser import change_chainid_pdb, append_pdbs

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def run_genetic_alg_multistate(run_dir, af2_flags_file, score_func, startingseqs, poolsizes = [], 
                               num_iter = 50, n_workers=1, mut_percents=None, contacts=None, 
                               mpnn_temp="0.1", mpnn_version="s_48_020", skip_mpnn=[], mpnn_iters=None, mpnn_chains=None,
                               repeat_af2=True, af2_preds=[], crossover_percent=0.2, vary_length=0, 
                               write_pdbs=False, plot=[], conf_plot=False, write_compressed_data=True, bidirectional=False):

    num_af2=0
    lengths = []
    num_preds = len(af2_preds)
    print(af2_preds, num_preds, mpnn_chains)
    # af2_preds tells EvoPro which chains to send to separate AF2 runs
    for chains in af2_preds:
        c = list(chains)
        lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])][0] for chain in c])

    print("Compiling AF2 models for lengths:", lengths)

    print("Initializing distributor")
    dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)

    scored_seqs = {}  # collection of ALL seqs that have been tested so far
    if repeat_af2:
        repeat_af2_seqs = {}
    curr_iter = 1
    seqs_per_iteration = []
    newpool = startingseqs
    all_results = {}

    while len(poolsizes) < num_iter:
        poolsizes.append(poolsizes[-1])
    #print(len(poolsizes))

    output_dir = run_dir + "outputs/"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    if write_pdbs:
        pdb_folder = output_dir + "pdbs_per_iter/"
        if not os.path.isdir(pdb_folder):
            os.makedirs(pdb_folder)

    #start genetic algorithm iteration
    while curr_iter <= num_iter:
        all_seqs = []
        all_scores = []
        pool = newpool  # collection of seqs in the current pool
        if curr_iter == 1:
            print("\nIteration 1: Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs_henry(pool, 
                                   poolsizes[curr_iter-1], 
                                   crossover_percent=0, 
                                   all_seqs = list(scored_seqs.keys()), 
                                   vary_length=vary_length)

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters and curr_iter not in skip_mpnn:
            print("\nIteration " + str(curr_iter) + ": refilling with ProteinMPNN.")
            pool = create_new_seqs_mpnn_henry(pool, 
                                        scored_seqs, 
                                        poolsizes[curr_iter-1], 
                                        run_dir, 
                                        curr_iter, 
                                        all_seqs = list(scored_seqs.keys()), 
                                        af2_preds=af2_preds,
                                        mpnn_temp=mpnn_temp, 
                                        mpnn_version=mpnn_version, 
                                        mpnn_chains=mpnn_chains, bidir=bidirectional)

        #otherwise refilling pool with just mutations and crossovers
        else:
            print("\nIteration " + str(curr_iter) + ": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs_henry(pool, 
                                   poolsizes[curr_iter-1], 
                                   crossover_percent=crossover_percent, 
                                   mut_percent=mut_percents[curr_iter-1], 
                                   all_seqs = list(scored_seqs.keys()), 
                                   vary_length=vary_length)

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
        
        seqs_packed = [work_list_all[i:i+num_preds] for i in range(0, len(work_list_all), num_preds)]
        results_packed = [results_all[i:i+num_preds] for i in range(0, len(results_all), num_preds)]

        print("These should be equal", len(work_list_all), len(results_all), num_preds*len(results_packed))
        
        scores = []
        if not contacts:
            contacts=(None, None, None)
        for seqs, results, dsobj in zip(seqs_packed, results_packed, scoring_pool):
            scores.append(score_func(results, dsobj, contacts=contacts))
        
        #adding sequences and scores into the dictionary
        for score, seqs, results, dsobj in zip(scores, seqs_packed, results_packed, scoring_pool):
            key_seq = dsobj.get_sequence_string()
            print(seqs)
            #if key_seq != ",".join(seqs[0]):
            #        print(key_seq, str(seqs[0]))
            #        raise AssertionError("Sequence does not match DSobj sequence")
            
            overall_scores = [score[0]]
            overall_scores = overall_scores + [x[0] for x in score[1]]
            score_split = score[1]
            pdbs = score[-2]
            score_all = (overall_scores, score_split)
            print(key_seq, overall_scores, score_split)
            
            if key_seq in scored_seqs and repeat_af2:
                scored_seqs[key_seq]["data"].append({"score": score_all, "pdb": pdbs, "result": result})
                #don't need to sort when using all 5 for average score
                #scored_seqs[key_seq]["data"].sort(key=lambda x: x["score"][0][0])
                
                sum_score = [0 for _ in overall_scores]
                print("Before", sum_score)
                for elem in scored_seqs[key_seq]["data"]:
                    sum_score = [x+y for x,y in zip(sum_score, elem["score"][0])]
                avg_score = [x/len(scored_seqs[key_seq]["data"]) for x in sum_score]
                print("After", sum_score, avg_score)
                scored_seqs[key_seq]["average"] = avg_score
            else:
                # AF2 structures compiled here
                scored_seqs[key_seq] = {"data": [{"score": score, "pdb": pdbs, "result": result}]}
                scored_seqs[key_seq]["dsobj"] =  dsobj
                
                #set average score here
                scored_seqs[key_seq]["average"] = overall_scores
       
        if repeat_af2:
            for key_seq in scored_seqs:
                if len(scored_seqs[key_seq]["data"]) >=5:
                    repeat_af2_seqs[key_seq] = scored_seqs[key_seq]["average"]
            print("repeat_af2_seqs", repeat_af2_seqs)
        
        #creating sorted list version of sequences and scores in the pool
        sorted_scored_pool = []
        for dsobj, j in zip(pool, range(len(pool))):
            key_seq = dsobj.get_sequence_string()
            print(key_seq, scored_seqs[key_seq]["average"])
            if repeat_af2:
                if key_seq in repeat_af2_seqs:
                    sorted_scored_pool.append((key_seq, repeat_af2_seqs[key_seq]))
            else:
                sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["average"]))
            
            pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
            if write_pdbs:
                for pdb, chains in zip(pdbs, af2_preds):
                    with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_chain"+str(chains)+".pdb", "w") as pdbf:
                                pdbf.write(str(pdb))
        
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
        result = scored_seqs[key_seq]["data"][0]["result"]
        results = [scored_seqs[key_seq]["data"][x]["result"] for x in range(len(scored_seqs[key_seq]["data"]))]
        
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
            with open(output_dir + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_chain"+str(chains)+".pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
    
    #plotting
    plotting = plot_scores_general(plot, af2_preds, seqs_per_iteration, output_dir)
    print("plot created at", plotting)
    
    print("Number of AlphaFold2 predictions: ", num_af2)
    dist.spin_down()

if __name__ == "__main__":
    parser = getEvoProParser()
    parser.add_argument('--bidirectional', action='store_true', help='Enable bidirectional coding for any 2-part MSD constraints.')

    args = parser.parse_args(sys.argv[1:])

    input_dir = args.input_dir
    if not input_dir.endswith("/"):
        input_dir = input_dir + "/"

    if input_dir not in ['', None]:
            onlyfiles = [f for f in os.listdir(input_dir) if os.path.isfile(
                         os.path.join(input_dir, f))]
            print('Only Files:', onlyfiles)

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

            dsobj1 = DesignSeqMSD(jsonfile=input_dir + resfile)

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
                starting_seqs = create_new_seqs_henry(starting_seqs, args.pool_size, crossover_percent=0)

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

    contacts=None
    if args.define_contact_area:
        contacts = parse_mutres_input(args.define_contact_area)
    else:
        contacts=None

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

    print(args.dont_write_compressed_data)
    
    
    af2_preds = args.af2_preds.strip().split(",")
    if args.mpnn_chains is not None:
        mpnn_chains = args.mpnn_chains.strip().split(",")
    else:
        mpnn_chains = None
                                   
    run_genetic_alg_multistate(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, mut_percents=mut_percents, contacts=contacts, 
        mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, skip_mpnn=mpnn_skips, mpnn_iters=mpnn_iters, mpnn_chains=mpnn_chains,
        repeat_af2=not args.no_repeat_af2, af2_preds = af2_preds, crossover_percent=args.crossover_percent, vary_length=args.vary_length, 
        write_pdbs=args.write_pdbs, plot=plot_style, conf_plot=args.plot_confidences, write_compressed_data=not args.dont_write_compressed_data, 
        bidirectional=args.bidirectional)
        
        
        
