import multiprocessing as mp
import importlib
import collections
import random
import time
import math
import sys, os
from DesignSeq import DesignSeq
from typing import Sequence, Union
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.utils.distributor import Distributor
from geneticalg_helpers import generate_random_seqs, read_starting_seqs, create_new_seqs, create_new_seqs_mpnn
from evopro.user_inputs.inputs import getEvoProParser
from evopro.score_funcs.score_funcs import write_raw_plddt, write_pairwise_scores

sys.path.append('/proj/kuhl_lab/alphafold/run')
from functools import partial
import numpy as np

def run_genetic_alg_gpus(run_dir, af2_flags_file, score_func, startingseqs, poolsize = 40, num_iter = 50, n_workers=2, write_raw_plddts=False, write_pair_confs=False, stabilize_binder=None, rmsd_func=None ,write_pdbs=False, mpnn_iters=None, len_binder=0, use_of = False, crossover_percent = 0.2, vary_length = 0, mut_percents=None, stabilize_monomer=None):

    lengths = [[x + vary_length for x in startingseqs[0].get_lengths()]]
    if stabilize_monomer:
        for chain in stabilize_monomer:
            lengths.append([x + vary_length for x in startingseqs[0].get_lengths([chain])])

    if use_of:
        from geneticalg_helpers import of_init
        sys.path.append('/proj/kuhl_lab/OmegaFold/')
        dist = Distributor(n_workers, of_init, af2_flags_file, lengths, score_func)
    
    else:
        from geneticalg_helpers import af2_init
        from run_af2 import af2
        dist = Distributor(n_workers, af2_init, af2_flags_file, lengths, score_func)

    scored_seqs = {}
    curr_iter = 1
    seqs_per_iteration = []
    newpool = startingseqs

    if write_pdbs:
        pdb_folder = run_dir + "pdbs_per_iter/"
        if not os.path.isdir(pdb_folder):
            os.makedirs(pdb_folder)

    #start genetic algorithm iteration
    while curr_iter <= num_iter:
        all_seqs = []
        all_scores = []
        pool = newpool

        if curr_iter == 1:
            print("Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs(pool, poolsize, crossover_percent=0, all_seqs = list(scored_seqs.keys()))

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters:
            print("Iteration "+str(curr_iter)+": refilling with ProteinMPNN.")
            pool = create_new_seqs_mpnn(pool, scored_seqs, poolsize, run_dir, curr_iter, all_seqs = list(scored_seqs.keys()))

        #otherwise refilling pool with just mutations and crossovers
        else:
            print("Iteration "+str(curr_iter)+": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            print(curr_iter)
            pool = create_new_seqs(pool, poolsize, crossover_percent=crossover_percent, mut_percent=mut_percents[curr_iter-1], all_seqs = list(scored_seqs.keys()))

        #print("Starting iteration " + str(curr_iter))

        scoring_pool = [p for p in pool if p.get_sequence_string() not in scored_seqs]
        work_list = [[[dsobj.jsondata["sequence"][chain] for chain in dsobj.jsondata["sequence"]]] for dsobj in scoring_pool]

        #if needed, also running af2 on the binder alone
        work_list_all = []
        work_list_2 = []
        if stabilize_binder:
            for chain in stabilize_binder:
                work_list_2.append([dsobj.jsondata["sequence"][chain] for dsobj in scoring_pool])

            work_list_all = work_list + work_list_2
        else:
            work_list_all = work_list

        results = dist.churn(work_list_all)

        #separating complex and binder sequences, if needed
        complex_seqs = []
        binder_seqs_all = []
        for seq, i in zip(work_list_all, range(len(work_list_all))):
            if stabilize_binder:
                if i<len(work_list):
                    complex_seqs.append(seq)
                else:
                    binder_seqs_all.append(seq)
            else:
                complex_seqs.append(seq)

        if stabilize_binder:
            binder_seqs_split = []
            binder_seqs = []
        
            for i in range(0, len(binder_seqs_all), len(work_list)):
                binder_seqs_split.append(binder_seqs_all[i:i + len(work_list)])

            for tup in zip(*binder_seqs_split):
                binder_seqs.append(list(tup))

        #separating complex and binder scores, if needed
        complex_scores = []
        binder_scores_all = []
        for score, i in zip(results, range(len(results))):
            if stabilize_binder:
                if i<len(work_list):
                    complex_scores.append(score)
                else:
                    binder_scores_all.append(score)
            else:
                complex_scores.append(score)

        if stabilize_binder:
            binder_scores_split = []
            binder_scores = []

            for i in range(0, len(binder_scores_all), len(work_list)):
                binder_scores_split.append(binder_scores_all[i:i + len(work_list)])

            for tup in zip(*binder_scores_split):
                binder_scores.append(list(tup))

        #adding sequences and scores into the dictionary
        if stabilize_binder:
            for dsobj, seq, cscore, bscores in zip(scoring_pool, complex_seqs, complex_scores, binder_scores):
                key_seq = dsobj.get_sequence_string()
                if key_seq == ",".join(seq[0]):
                    raise ValueError("Sequence does not match DSobj sequence")
                rmsd_score = 0
                rmsd_score_list = []
                if rmsd_func:
                    for bscore in bscores:
                        rmsd_score_list.append(rmsd_func(cscore[0][-2], bscore[0][-2]))
                    rmsd_score = sum(rmsd_score_list)
                score = (cscore[0][0] + bscore[0][0] + rmsd_score, cscore[0][1], bscore[0][1], (rmsd_score, rmsd_score_list))
                pdb = (cscore[0][-2], [bscore[0][-2] for bscore in bscores])
                scored_seqs[key_seq] = {"dsobj": dsobj, "score": score, "pdb": pdb}
        else:
            for dsobj, seq, cscore in zip(scoring_pool, complex_seqs, complex_scores):
                key_seq = dsobj.get_sequence_string()
                score = (cscore[0][0], cscore[0][1])
                pdb = (cscore[0][-2], None)
                scored_seqs[key_seq] = {"dsobj": dsobj, "score": score, "pdb": pdb}

        
        #creating sorted list version of sequences and scores in the pool
        sorted_scored_pool = []
        for dsobj in pool:
            key_seq = dsobj.get_sequence_string()
            sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["score"]))
        sorted_scored_pool.sort(key = lambda x: x[1][0])

        seqs_per_iteration.append(sorted_scored_pool)

        #writing log
        with open(run_dir+"runtime_seqs_and_scores.log", "a") as logf:
            logf.write("starting iteration " + str(curr_iter)+ " log\n")
            for elem in sorted_scored_pool:
                logf.write(str(elem[0]) + "\t" + str(elem[1]) + "\n")
        print("done writing runtime results")

        #create a new pool of only 50% top scoring sequences for the next iteration
        newpool_seqs = []
        for sp in sorted_scored_pool[:round(len(sorted_scored_pool)/2)]:
            newpool_seqs.append(sp[0])
        print("newpool", newpool_seqs)

        newpool = []
        #pulling back DS objects for each sequence in new pool
        for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
            newpool.append(scored_seqs[key_seq]["dsobj"])
            pdbs = scored_seqs[key_seq]["pdb"]
            if write_pdbs:
                if stabilize_binder:
                    for pdb, k in zip(pdbs, range(len(pdbs))):
                        if k==0:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdb))
                        else:
                            with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_binderonly_chain" + stabilize_binder[k-1] + ".pdb", "w") as pdbf:
                                pdbf.write(str(pdb))
                else:
                    with open(pdb_folder + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                                pdbf.write(str(pdbs[0]))

        curr_iter+=1

    with open(run_dir+"seqs_and_scores.log", "w") as logf:
        for tuplist, i in zip(seqs_per_iteration, range(len(seqs_per_iteration))):
            logf.write(">iteration"+str(i+1)+"\n")
            for val in tuplist:
                logf.write(str(val[0])+"\t"+str(val[1])+"\n")

    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        seq = key_seq.split(" ")
        pdbs = scored_seqs[key_seq]["pdb"]
        if stabilize_binder:
            for pdb, k in zip(pdbs, range(len(pdbs))):
                if k==0:
                    with open(run_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
                else:
                    with open(run_dir + "seq_" + str(j) + "_final_model_1_binderonly_chain" + stabilize_binder[k-1] + ".pdb", "w") as pdbf:
                        pdbf.write(str(pdb))
        else:
            with open(run_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
                        pdbf.write(str(pdbs[0]))

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
    print(mut_percents)
    
    mut_percents = mut_percents[:args.num_iter]

    stabilize = [x.strip() for x in args.stabilize_monomer.split(",")]

    if args.mpnn_iters:
        mpnn_iters = [int(x.strip()) for x in args.mpnn_iters.split(",")]
    else:
        mpnn_iters = []
        freq = int(args.mpnn_freq)
        for i in range(1, args.num_iter):
            if i % freq == 0:
                mpnn_iters.append(i)

    run_genetic_alg_gpus(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsize = args.pool_size, num_iter = args.num_iter, n_workers=args.num_gpus, write_raw_plddts=args.write_plddt, write_pair_confs=args.write_pae, stabilize_monomer=stabilize, rmsd_func=rmsdfunc, write_pdbs=args.write_pdbs, mpnn_iters=mpnn_iters, use_of=args.use_omegafold, crossover_percent=args.crossover_percent, vary_length=args.vary_length, mut_percents=mut_percents)
