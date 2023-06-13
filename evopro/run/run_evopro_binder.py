import multiprocessing as mp
import importlib
import math
import sys, os
from typing import Sequence, Union
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.utils.plot_scores import plot_scores_stabilize_monomer_top, plot_scores_stabilize_monomer_avg, plot_scores_stabilize_monomer_median
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import read_starting_seqs, create_new_seqs, create_new_seqs_mpnn
from evopro.user_inputs.inputs import getEvoProParser
from evopro.utils.plots import get_chain_lengths, plot_pae, plot_plddt
from evopro.utils.utils import compressed_pickle

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init
from functools import partial
import numpy as np

def run_genetic_alg_multistate(run_dir, af2_flags_file, score_func, startingseqs, poolsizes = [], 
                               num_iter = 50, n_workers=1, mut_percents=None, contacts=None, distance_cutoffs=None,
                               rmsd_func=None, rmsd_to_starting_func=None, rmsd_to_starting_pdb=None,
                               mpnn_temp="0.1", mpnn_version="s_48_020", skip_mpnn=[], mpnn_iters=None, 
                               repeat_af2=True, af2_preds_extra=[], crossover_percent=0.2, vary_length=0, 
                               write_pdbs=False, plot=[], conf_plot=False, write_compressed_data=True):

    num_af2=0
    
    print("Repeating AF2", repeat_af2)
    
    lengths = [[x + vary_length for x in startingseqs[0].get_lengths()]]
    for chains in af2_preds_extra:
        c = list(chains)
        lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])] for chain in c][0])

    print(lengths)

    print("initializing distributor")
    dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)

    scored_seqs = {}
    if repeat_af2:
        repeat_af2_seqs = {}
    curr_iter = 1
    seqs_per_iteration = []
    newpool = startingseqs
    all_results = {}

    while len(poolsizes) < num_iter:
        poolsizes.append(poolsizes[-1])

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
            print("Iteration " + str(curr_iter) + ": refilling with ProteinMPNN.")
            pool = create_new_seqs_mpnn(pool, scored_seqs, poolsizes[curr_iter-1], run_dir, curr_iter, all_seqs = list(scored_seqs.keys()), mpnn_temp=mpnn_temp, mpnn_version=mpnn_version)

        #otherwise refilling pool with just mutations and crossovers
        else:
            print("Iteration " + str(curr_iter) + ": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs(pool, poolsizes[curr_iter-1], crossover_percent=crossover_percent, mut_percent=mut_percents[curr_iter-1], all_seqs = list(scored_seqs.keys()), vary_length=vary_length)
        
        if repeat_af2:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in repeat_af2_seqs]
        else:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in scored_seqs]
        
        work_list = [[[dsobj.jsondata["sequence"][chain] for chain in dsobj.jsondata["sequence"]]] for dsobj in scoring_pool]

        #if needed, also running af2 on the binder alone
        work_list_all = []
        work_list_2 = []
        if af2_preds_extra:
            for chains in af2_preds_extra:
                c = list(chains)
                work_list_2 = [[[dsobj.jsondata["sequence"][chain] for chain in c]] for dsobj in scoring_pool]
                work_list_all = work_list + work_list_2
        else:
            work_list_all = work_list

        print("work list", work_list_all)
        num_af2 += len(work_list_all)

        results = dist.churn(work_list_all)
        print("done churning")

        if af2_preds_extra:
            num_lists = 1 + len(af2_preds_extra)
        else:
            num_lists = 1

        scoring_dsobjs = []
        for n in range(num_lists):
            scoring_dsobjs = scoring_dsobjs + scoring_pool

        all_scores = []
        #score the af2 results
        for seq, result, dsobj in zip(work_list_all, results, scoring_dsobjs):
            while type(result) is list:
                result = result[0]
            if contacts is not None:
                all_scores.append(score_func(result, dsobj, contacts=contacts, distance_cutoffs=distance_cutoffs))
            else:
                all_scores.append(score_func(result, dsobj))

        #separating complex and binder sequences, if needed
        complex_seqs = []
        binder_seqs_all = []
        for seq, i in zip(work_list_all, range(len(work_list_all))):
            if af2_preds_extra:
                if i<len(work_list):
                    complex_seqs.append(seq)
                else:
                    binder_seqs_all.append(seq)
            else:
                complex_seqs.append(seq)

        if af2_preds_extra:
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
            if af2_preds_extra:
                if i<len(work_list):
                    complex_scores.append(score)
                else:
                    binder_scores_all.append(score)
            else:
                complex_scores.append(score)

        if af2_preds_extra:
            binder_scores_split = []
            binder_scores = []

            for i in range(0, len(binder_scores_all), len(work_list)):
                binder_scores_split.append(binder_scores_all[i:i + len(work_list)])

            for tup in zip(*binder_scores_split):
                binder_scores.append(list(tup))

        #adding sequences and scores into the dictionary
        if af2_preds_extra:
            print("adding sequences and scores into the dictionary, with af2_preds_extra")
            for dsobj, seq, cscore, bscores in zip(scoring_pool, complex_seqs, complex_scores, binder_scores):
                key_seq = dsobj.get_sequence_string()
                print(key_seq)
                if key_seq != ",".join(seq[0]):
                    print(key_seq, str(seq[0]))
                    raise AssertionError("Sequence does not match DSobj sequence")
                rmsd_score = 0
                rmsd_score_list = []
                if rmsd_func:
                    print("calculating rmsd score")
                    for bscore,chain in zip(bscores, af2_preds_extra):
                        rmsd_score_list.append(rmsd_func(cscore[-2], bscore[-2], binder_chain=chain, dsobj=dsobj))
                if rmsd_to_starting_func:
                    #print("calculating rmsd to starting score")
                    for bscore,chain in zip(bscores, af2_preds_extra):
                        rmsd_score_list.append(rmsd_to_starting_func(bscore[-2], rmsd_to_starting_pdb, dsobj=dsobj))
                
                rmsd_score = sum(rmsd_score_list)
                bscore = sum([b[0] for b in bscores])
                overall_score = cscore[0] + bscore + rmsd_score
                score = ((overall_score, cscore[0], bscore, rmsd_score), cscore[1], (bscore, [bscor[1] for bscor in bscores]), (rmsd_score, rmsd_score_list))
                print(cscore[0], bscore, rmsd_score, score)
                pdb = (cscore[-2], [bscore[-2] for bscore in bscores])
                result = (cscore[-1], [bscore[-1] for bscore in bscores])
                if repeat_af2:
                    if key_seq in scored_seqs:
                        scored_seqs[key_seq]["data"].append({"score": score, "pdb": pdb, "result": result})
                        scored_seqs[key_seq]["data"].sort(key=lambda x: x["score"][0][0])
                        
                        #recalculate average score here
                        avg_score = [0 for x in score[0]]
                        for elem in scored_seqs[key_seq]["data"][:3]:
                            avg_score = [x+y for x,y in zip(avg_score, elem["score"][0])]
                        avg_score = [x/len(scored_seqs[key_seq]["data"][:3]) for x in avg_score]
                        scored_seqs[key_seq]["average"] = avg_score
                    else:
                        scored_seqs[key_seq] = {"data": [{"score": score, "pdb": pdb, "result": result}]}
                        scored_seqs[key_seq]["dsobj"] =  dsobj
                        
                        #set average score here
                        scored_seqs[key_seq]["average"] = score[0]
                else:
                    scored_seqs[key_seq] = {"dsobj": dsobj, "data": [{"score": score, "pdb": pdb, "result": result}]}
                    scored_seqs[key_seq]["average"] = score[0]
        else:
            print("adding sequences and scores into the dictionary, no extra af2 preds")
            for dsobj, seq, cscore in zip(scoring_pool, complex_seqs, complex_scores):
                key_seq = dsobj.get_sequence_string()
                print(key_seq)
                if rmsd_to_starting_func:
                    rmsd_score = rmsd_to_starting_func(bscore[-2], rmsd_to_starting_pdb, dsobj=dsobj)
                overall_score = cscore[0] + rmsd_score
                score = ((overall_score, cscore[0], rmsd_score), cscore[1], rmsd_score)
                pdb = (cscore[-2], None)
                result = (cscore[-1], )
                if repeat_af2:
                    if key_seq in scored_seqs:
                        scored_seqs[key_seq]["data"].append({"score": score, "pdb": pdb, "result": result})
                        print("before sorting")
                        for d in scored_seqs[key_seq]["data"]:
                            print(d["score"])
                        scored_seqs[key_seq]["data"].sort(key=lambda x: x["score"][0][0])
                        print("after sorting")
                        for d in scored_seqs[key_seq]["data"]:
                            print(d["score"])
                        #recalculate average score here
                        avg_score = [0 for x in score[0]]
                        for elem in scored_seqs[key_seq]["data"][:3]:
                            avg_score = [x+y for x,y in zip(avg_score, elem["score"][0])]
                        avg_score = [x/len(scored_seqs[key_seq]["data"][:3]) for x in avg_score]
                        scored_seqs[key_seq]["average"] = avg_score
                    else:
                        scored_seqs[key_seq] = {"data": [{"score": score, "pdb": pdb, "result": result}]}
                        scored_seqs[key_seq]["dsobj"] =  dsobj
                        #set average score here
                        scored_seqs[key_seq]["average"] = score[0]
                else:
                    scored_seqs[key_seq] = {"dsobj": dsobj, "data": [{"score": score, "pdb": pdb, "result": result}]}
                    scored_seqs[key_seq]["average"] = score[0]

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
            else:
                sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["average"]))
            
            pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
            if write_pdbs:
                if af2_preds_extra:
                    for pdb, k in zip(pdbs, range(len(pdbs))):
                        if k==0:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdb))
                        else:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_binderonly_chain" + af2_preds_extra[k-1] + ".pdb", "w") as pdbf:
                                pdbf.write(str(pdb[0]))
                else:
                    with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdbs[0]))
        
        print("before sorting", sorted_scored_pool)
        sorted_scored_pool.sort(key = lambda x: x[1][0])
        print("after sorting", sorted_scored_pool)
        seqs_per_iteration.append(sorted_scored_pool)

        #writing log
        
        with open(output_dir+"runtime_seqs_and_scores.log", "a") as logf:
            logf.write("starting iteration " + str(curr_iter)+ " log\n")
            for elem in sorted_scored_pool:
                key_seq = elem[0]
                logf.write(str(key_seq) + "\t" + str(elem[1]) + "\t")
                for s in scored_seqs[key_seq]["data"]:
                    logf.write(str(s["score"]) + "\t")
                logf.write("\n")
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

    if "avg" in plot:
        if af2_preds_extra:
            if not rmsd_func:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_avg(seqs_per_iteration, output_dir)
    if "top" in plot:
        if af2_preds_extra:
            if not rmsd_func:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_top(seqs_per_iteration, output_dir)

    if "median" in plot:
        if af2_preds_extra:
            if not rmsd_func:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir, rmsd=False)
            else:
                plot_scores_stabilize_monomer_median(seqs_per_iteration, output_dir)

    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        seq = key_seq.split(" ")
        pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
        result = scored_seqs[key_seq]["data"][0]["result"][0]

        if write_compressed_data:
            print("writing compressed data")
            compressed_pickle(os.path.join(output_dir, "seq_" + str(j) + "_result"), result)

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

        if af2_preds_extra:
            for pdb, k in zip(pdbs, range(len(pdbs))):
                if k==0:
                    with open(output_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
                else:
                    with open(output_dir + "seq_" + str(j) + "_final_model_1_binderonly_chain" + af2_preds_extra[k-1] + ".pdb", "w") as pdbf:
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

    print("Writing compressed data", (not args.dont_write_compressed_data))
    
    if args.af2_preds_extra:
        af2_preds_extra = args.af2_preds_extra.strip().split(",")
                                   
    run_genetic_alg_multistate(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, mut_percents=mut_percents, contacts=contacts, distance_cutoffs=distance_cutoffs,
        rmsd_func=rmsdfunc, rmsd_to_starting_func=rmsd_to_starting_func, rmsd_to_starting_pdb=path_to_starting,
        mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, skip_mpnn=mpnn_skips, mpnn_iters=mpnn_iters, 
        repeat_af2=not args.no_repeat_af2, af2_preds_extra = af2_preds_extra, crossover_percent=args.crossover_percent, vary_length=args.vary_length, 
        write_pdbs=args.write_pdbs, plot=plot_style, conf_plot=args.plot_confidences, write_compressed_data=not args.dont_write_compressed_data)
        
        
        
