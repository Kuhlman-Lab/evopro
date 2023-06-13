import multiprocessing as mp
import importlib
import collections
import random
import time
import math
import sys, os
from typing import Sequence, Union
sys.path.append("/proj/kuhl_lab/evopro/")
#sys.path.append("/nas/longleaf/home/amritan/Desktop/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.utils.plot_scores import plot_scores_stabilize_monomer_top, plot_scores_stabilize_monomer_avg, plot_scores_stabilize_monomer_median
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import generate_random_seqs, read_starting_seqs, create_new_seqs, create_new_seqs_mpnn_old, create_new_seqs_mpnn
from evopro.user_inputs.inputs import getEvoProParser
from evopro.score_funcs.score_funcs import write_raw_plddt, write_pairwise_scores

sys.path.append('/proj/kuhl_lab/alphafold/run')
from functools import partial
import numpy as np

def run_genetic_alg_gpus(run_dir, af2_flags_file, score_func, startingseqs, poolsizes = [], num_iter = 50, 
    n_workers=2, write_raw_plddts=False, write_pair_confs=False, rmsd_func=None, 
    rmsd_to_starting_func=None, rmsd_to_starting_pdb=None, starting_pdb=None,
    write_pdbs=False, mpnn_iters=None, len_binder=0, use_of=False, crossover_percent=0.2, vary_length=0, 
    mut_percents=None, stabilize_monomer=None, contacts=None, plot=[], conf_plot=False, mpnn_temp="0.1", 
    skip_mpnn=[], repeat_af2=False):

    num_af2=0
    lengths = [[x + vary_length for x in startingseqs[0].get_lengths()]]
    if stabilize_monomer:
        for chain in stabilize_monomer:
            lengths.append([x + vary_length for x in startingseqs[0].get_lengths([chain])])

    print(lengths)

    if use_of:
        from evopro.genetic_alg.geneticalg_helpers import of_init
        sys.path.append('/proj/kuhl_lab/OmegaFold/')
        dist = Distributor(n_workers, of_init, af2_flags_file, lengths)
    
    else:
        #from evopro.genetic_alg.geneticalg_helpers import af2_init
        from run_af2 import af2, af2_init
        print("initializing distributor")
        dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)

    scored_seqs = {}
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
        pool = newpool

        if curr_iter == 1:
            print("Iteration 1: Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs(pool, poolsizes[curr_iter-1], crossover_percent=0, all_seqs = list(scored_seqs.keys()), vary_length=vary_length)

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters and curr_iter not in skip_mpnn:
            print("Iteration "+str(curr_iter)+": refilling with ProteinMPNN.")
            pool = create_new_seqs_mpnn_old(pool, scored_seqs, poolsizes[curr_iter-1], run_dir, curr_iter, all_seqs = list(scored_seqs.keys()), mpnn_temp=mpnn_temp)

        #otherwise refilling pool with just mutations and crossovers
        else:
            print("Iteration "+str(curr_iter)+": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs(pool, poolsizes[curr_iter-1], crossover_percent=crossover_percent, mut_percent=mut_percents[curr_iter-1], all_seqs = list(scored_seqs.keys()), vary_length=vary_length)

        #print("Starting iteration " + str(curr_iter))
        
        if repeat_af2:
            scoring_pool = [p for p in pool]
        else:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in scored_seqs]
        work_list = [[[dsobj.jsondata["sequence"][chain] for chain in dsobj.jsondata["sequence"]]] for dsobj in scoring_pool]

        #if needed, also running af2 on the binder alone
        work_list_all = []
        work_list_2 = []
        if stabilize_monomer:
            for chain in stabilize_monomer:
                work_list_2 = work_list_2 + [[[dsobj.jsondata["sequence"][chain]]] for dsobj in scoring_pool]

            work_list_all = work_list + work_list_2
        else:
            work_list_all = work_list

        print("work list", work_list_all)
        num_af2 += len(work_list_all)

        results = dist.churn(work_list_all)
        print("done churning")

        if stabilize_monomer:
            num_lists = 1 + len(stabilize_monomer)
        else:
            num_lists = 1

        scoring_dsobjs = []
        for n in range(num_lists):
            scoring_dsobjs = scoring_dsobjs + scoring_pool

        all_scores = []
        print("These should be equal", len(work_list_all), len(results), len(scoring_dsobjs))
        for seq, result, dsobj in zip(work_list_all, results, scoring_dsobjs):
            while type(result) is list:
                print(result)
                result = result[0]
            print(seq, result, dsobj, contacts)
            if contacts is not None:
                if starting_pdb:
                    all_scores.append(score_func(result, dsobj, starting_pdb=starting_pdb, contacts=contacts))

                else:
                    all_scores.append(score_func(result, dsobj, contacts=contacts))
            else:
                all_scores.append(score_func(result, dsobj))
            print("All scores length", len(all_scores))
            #all_scores.append(score_func(result, dsobj, contacts=contacts))

        #separating complex and binder sequences, if needed
        complex_seqs = []
        binder_seqs_all = []
        for seq, i in zip(work_list_all, range(len(work_list_all))):
            if stabilize_monomer:
                if i<len(work_list):
                    complex_seqs.append(seq)
                else:
                    binder_seqs_all.append(seq)
            else:
                complex_seqs.append(seq)

        if stabilize_monomer:
            binder_seqs_split = []
            binder_seqs = []
        
            for i in range(0, len(binder_seqs_all), len(work_list)):
                binder_seqs_split.append(binder_seqs_all[i:i + len(work_list)])

            for tup in zip(*binder_seqs_split):
                binder_seqs.append(list(tup))

        #separating complex and binder scores, if needed
        complex_scores = []
        binder_scores_all = []
        for score, i in zip(all_scores, range(len(all_scores))):
            if stabilize_monomer:
                if i<len(work_list):
                    complex_scores.append(score)
                else:
                    binder_scores_all.append(score)
            else:
                complex_scores.append(score)

        if stabilize_monomer:
            binder_scores_split = []
            binder_scores = []

            for i in range(0, len(binder_scores_all), len(work_list)):
                binder_scores_split.append(binder_scores_all[i:i + len(work_list)])

            for tup in zip(*binder_scores_split):
                binder_scores.append(list(tup))

        #adding sequences and scores into the dictionary
        if stabilize_monomer:
            print("adding sequences and scores into the dictionary")
            for dsobj, seq, cscore, bscores in zip(scoring_pool, complex_seqs, complex_scores, binder_scores):
                key_seq = dsobj.get_sequence_string()
                if key_seq != ",".join(seq[0]):
                    print(key_seq, str(seq[0]))
                    raise ValueError("Sequence does not match DSobj sequence")
                rmsd_score = 0
                rmsd_score_list = []
                if rmsd_func:
                    for bscore,chain in zip(bscores, stabilize_monomer):
                        rmsd_score_list.append(rmsd_func(cscore[-2], bscore[-2], binder_chain=chain, dsobj=dsobj))
                if rmsd_to_starting_func:
                    for bscore,chain in zip(bscores, stabilize_monomer):
                        rmsd_score_list.append(rmsd_to_starting_func(bscore[-2], rmsd_to_starting_pdb, dsobj=dsobj))
                
                rmsd_score = sum(rmsd_score_list)
                bscore = sum([b[0] for b in bscores])
                score = (cscore[0] + bscore + rmsd_score, cscore[1], (bscore, [bscor[1] for bscor in bscores]), (rmsd_score, rmsd_score_list))
                pdb = (cscore[-2], [bscore[-2] for bscore in bscores])
                result = (cscore[-1], [bscore[-1] for bscore in bscores])
                print(key_seq, dsobj, score)
                scored_seqs[key_seq] = {"dsobj": dsobj, "score": score, "pdb": pdb, "result": result}
        else:
            print("adding sequences and scores into the dictionary, no stabilize_monomer")
            for dsobj, seq, cscore in zip(scoring_pool, complex_seqs, complex_scores):
                key_seq = dsobj.get_sequence_string()
                score = (cscore[0], cscore[1])
                pdb = (cscore[-2], None)
                result = (cscore[-1], )
                scored_seqs[key_seq] = {"dsobj": dsobj, "score": score, "pdb": pdb, "result": result}

        print("scored sequences dictionary", list(scored_seqs.keys()))
        
        #creating sorted list version of sequences and scores in the pool
        sorted_scored_pool = []
        for dsobj, j in zip(pool, range(len(pool))):
            key_seq = dsobj.get_sequence_string()
            sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["score"]))
            pdbs = scored_seqs[key_seq]["pdb"]
            if write_pdbs:
                if stabilize_monomer:
                    for pdb, k in zip(pdbs, range(len(pdbs))):
                        if k==0:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdb))
                        else:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_binderonly_chain" + stabilize_monomer[k-1] + ".pdb", "w") as pdbf:
                                pdbf.write(str(pdb[0]))
                else:
                    with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdbs[0]))

        sorted_scored_pool.sort(key = lambda x: x[1][0])

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
            pdbs = scored_seqs[key_seq]["pdb"]

        curr_iter+=1

    with open(output_dir+"seqs_and_scores.log", "w") as logf:
        for tuplist, i in zip(seqs_per_iteration, range(len(seqs_per_iteration))):
            logf.write(">iteration"+str(i+1)+"\n")
            for val in tuplist:
                logf.write(str(val[0])+"\t"+str(val[1])+"\n")

    if "avg" in plot:
        if stabilize_monomer:
            if not rmsd_func:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir)
    if "top" in plot:
        if stabilize_monomer:
            if not rmsd_func:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir)

    if "median" in plot:
        if stabilize_monomer:
            if not rmsd_func:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir)

    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        seq = key_seq.split(" ")
        pdbs = scored_seqs[key_seq]["pdb"]
        if conf_plot:

            from evopro.utils.plots import get_chain_lengths, plot_pae, plot_plddt, plot_ticks
            result = scored_seqs[key_seq]["result"][0]
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

        if stabilize_monomer:
            for pdb, k in zip(pdbs, range(len(pdbs))):
                if k==0:
                    with open(output_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
                else:
                    with open(output_dir + "seq_" + str(j) + "_final_model_1_binderonly_chain" + stabilize_monomer[k-1] + ".pdb", "w") as pdbf:
                        pdbf.write(str(pdb[0]))
        else:
            with open(output_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
                        pdbf.write(str(pdbs[0]))

    print("Number of AlphaFold2 predictions: ", num_af2)
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
                raise ValueError("Flags file for Alphafold/OmegaFold runs not provided. (if omegafold, please provide an empty flags file for now)")
            if resfile is None:
                raise ValueError("Please provide a residue specifications file.")

            print(input_dir + resfile)
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

    if args.len_binder:
        len_binder=args.len_binder
    else:
        len_binder=0

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

    stabilize=None

    if args.path_to_starting:
        rmsd_starting = args.path_to_starting
    else:
        rmsd_starting=None

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

    print(pool_sizes)

    run_genetic_alg_gpus(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, write_raw_plddts=args.write_plddt, 
        write_pair_confs=args.write_pae, stabilize_monomer=stabilize, rmsd_func=rmsdfunc, starting_pdb=rmsd_starting,
        write_pdbs=args.write_pdbs, mpnn_iters=mpnn_iters, use_of=args.use_omegafold, 
        crossover_percent=args.crossover_percent, vary_length=args.vary_length, mut_percents=mut_percents, 
        contacts=contacts, plot=plot_style, conf_plot=args.plot_confidences, mpnn_temp=args.mpnn_temp, 
        skip_mpnn=mpnn_skips, repeat_af2=args.repeat_af2)
