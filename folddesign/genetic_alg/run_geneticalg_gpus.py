import multiprocessing as mp
import importlib
import collections
import random
import time
import sys, os
from DesignSeq import DesignSeq
from typing import Sequence, Union
from folddesign.utils.distributor import Distributor
from geneticalg_helpers import af2_init_2models, generate_random_seqs, read_starting_seqs, create_new_seqs, create_new_seqs_mpnn
from folddesign.user_inputs.inputs import getFDDParser
from folddesign.score_funcs.score_funcs import write_raw_plddt, write_pairwise_scores
from folddesign.score_funcs.score_pd1_gpus import score_binder_rmsd

sys.path.append('/proj/kuhl_lab/alphafold/run')

from run_af2 import af2
from functools import partial
import numpy as np

def run_genetic_alg_gpus(run_dir, af2_flags_file, score_func, startingseqs, poolsize = 50, num_iter = 20, n_workers=2, write_raw_plddts=False, write_pair_confs=False, stabilize_binder=False, write_pdbs=False, mpnn_freq=10, symmetry=None, len_binder=0):

    with_linker=False
    if stabilize_binder:
        if len(startingseqs[0].seq.values())>1:
            lengths = [[len(x) for x in startingseqs[0].seq.values()], [len(x) for x in startingseqs[0].seq.values()][-1]]
        else:
            lengths = [[len(x) for x in startingseqs[0].seq.values()], len_binder]
            with_linker=True

    else:
        lengths = [[len(x) for x in startingseqs[0].seq.values()]]
    
    print(lengths)
    dist = Distributor(n_workers, af2_init_2models, af2_flags_file, lengths, score_func)

    scored_seqs = {}
    curr_iter = 1
    seqs_per_iteration = []
    newpool = startingseqs

    while curr_iter <= num_iter:
        all_seqs = []
        all_scores = []
        pool = newpool

        if mpnn_freq<=0:
            print("Not a valid frequency for Protein MPNN usage.")

        if curr_iter == 1:
            pool = create_new_seqs(pool, poolsize)

        #using protein mpnn to refill pool every _ steps
        elif curr_iter % mpnn_freq == 0:
            pool = create_new_seqs_mpnn(pool, scored_seqs, poolsize, run_dir, curr_iter)
            #pool = create_new_seqs(pool, poolsize)

        #otherwise refilling pool with just mutations and crossovers
        else:
            pool = create_new_seqs(pool, poolsize, symmetry=symmetry)

        print("starting iteration " + str(curr_iter))

        scoring_pool = [p for p in pool if p not in scored_seqs.keys()]
        work_list = [[[dsobj.seq[chain] for chain in dsobj.seq]] for dsobj in scoring_pool]

        #if needed, also running af2 on the binder alone
        work_list_2 = []
        if stabilize_binder:
            if with_linker:
                work_list_2 = [[[dsobj.seq[chain] for chain in dsobj.seq][0][-len_binder:]] for dsobj in scoring_pool]
            else:
                work_list_2 = [[[dsobj.seq[chain] for chain in dsobj.seq][-1]] for dsobj in scoring_pool]

        work_list_all = work_list + work_list_2

        results = dist.churn(work_list_all)

        #separating complex and binder sequences, if needed
        all_seqs = []
        complex_seqs = []
        binder_seqs = []
        for seq, i in zip(work_list_all, range(len(work_list_all))):
            if stabilize_binder:
                if i<len(work_list):
                    complex_seqs.append(seq)
                else:
                    binder_seqs.append(seq)
            else:
                all_seqs.append(seq)

        #separating complex and binder scores, if needed
        complex_scores = []
        binder_scores = []
        all_scores = []
        for score, i in zip(results, range(len(results))):
            if stabilize_binder:
                if i<len(work_list):
                    complex_scores.append(score)
                else:
                    binder_scores.append(score)
            else:
                all_scores.append(score)

        #adding sequences and scores into the dictionary
        if stabilize_binder:
            for comp, bind, cscore, bscore in zip(complex_seqs, binder_seqs, complex_scores, binder_scores):
                key_seq = " ".join(comp[0])
                seq = (" ".join(comp[0]), bind[0])
                rmsd_score = score_binder_rmsd(cscore[0][-2], bscore[0][-2])
                score = ((cscore[0][0]/10.0 + 5*bscore[0][0] + 100*rmsd_score, cscore[0][0], bscore[0][0], rmsd_score), (cscore[0][-2], bscore[0][-2]))
                scored_seqs[key_seq] = (seq, score)
        else:
            for seq, score in zip(all_seqs, all_scores):
                key_seq = " ".join(seq[0])
                seq = (" ".join(seq[0]))
                score = ((score[0][0],), (score[0][-2],))
                scored_seqs[key_seq] = (seq, score)

        #writing log
        with open(run_dir+"runtime_seqs_and_scores_unsorted.log", "a") as logf:
            logf.write("starting iteration " + str(curr_iter)+ " log\n")
            for key_seq in scored_seqs:
                logf.write(str(scored_seqs[key_seq][0]) + "\t" + str(scored_seqs[key_seq][1][0]) + "\n")
        print("done writing runtime results")
        
        #creating sorted list version of sequences and scores in the pool
        sorted_scored_pool = []
        for dsobj in pool:
            key_seq = " ".join([dsobj.seq[chain] for chain in dsobj.seq])
            sorted_scored_pool.append((key_seq, scored_seqs[key_seq][1][0]))
        print(sorted_scored_pool)
        sorted_scored_pool.sort(key = lambda x: x[1][0])

        seqs_per_iteration.append(sorted_scored_pool)

        #create a new pool of only 50% top scoring sequences for the next iteration
        newpool_seqs = []
        for sp in sorted_scored_pool[:round(len(sorted_scored_pool)/2)]:
            newpool_seqs.append(sp[0])
        print("newpool", newpool_seqs)

        newpool = []
        for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
            seq = key_seq.split(" ")
            seq_dict = {}
            for chain, s in zip(pool[0].seq, seq):
                seq_dict[chain] = s
            dsobj = DesignSeq(seq=seq_dict, mutable=pool[0].mutable)
            newpool.append(dsobj)
            pdb = scored_seqs[key_seq][1][1][0]
            if write_pdbs:
                with open(run_dir + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_complex.pdb", "w") as pdbf:
                    pdbf.write(str(pdb))
                if stabilize_binder:
                    binder_pdb = scored_seqs[key_seq][1][1][1]
                    with open(run_dir + "seq_" + str(j) + "_iter_" + str(curr_iter) + "_model_1_binder-only.pdb", "w") as bpdbf:
                        bpdbf.write(str(binder_pdb))

        curr_iter+=1

    with open(run_dir+"seqs_and_scores.log", "w") as logf:
        for tuplist, i in zip(seqs_per_iteration, range(len(seqs_per_iteration))):
            logf.write(">iteration"+str(i+1)+"\n")
            for val in tuplist:
                logf.write(str(val[0])+"\t"+str(val[1])+"\n")

    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        seq = key_seq.split(" ")
        pdb = scored_seqs[key_seq][1][1][0]
        with open(run_dir + "seq_" + str(j) + "_final_model_1_complex.pdb", "w") as pdbf:
            pdbf.write(str(pdb))
        if stabilize_binder:
            binder_pdb = scored_seqs[key_seq][1][1][1]
            with open(run_dir + "seq_" + str(j) + "_final_model_1_binder-only.pdb", "w") as bpdbf:
                bpdbf.write(str(binder_pdb))

    dist.spin_down()

if __name__ == "__main__":
    parser = getFDDParser()
    args = parser.parse_args(sys.argv[1:])

    input_dir = args.input_dir 
    if not input_dir.endswith("/"):
        input_dir = input_dir + "/"

    sym_list=None
    if args.symmetry:
        sym_list = args.symmetry.strip().split(",")

    if input_dir not in ['', None]:
            onlyfiles = [f for f in os.listdir(input_dir) if os.path.isfile(
                         os.path.join(input_dir, f))]

            flagsfile=None
            resfile=None
            starting_seqs=None
            for filename in onlyfiles:
                if "flag" in filename and "af2" in filename:
                    flagsfile = filename
            
            for filename in onlyfiles:
                if "residue" in filename and "spec" in filename:
                    resfile = filename

            if flagsfile is None:
                raise ValueError("Flags file for Alphafold runs not provided.")
            if resfile is None:
                raise ValueError("Please provide a residue specifications file.")    

            print(input_dir + resfile)
            dsobj1 = DesignSeq(jsonfile=input_dir + resfile)

            for filename in onlyfiles:
                if "starting" in filename:
                    starting_seqs = read_starting_seqs(input_dir + filename, dsobj1)
                    starting_seqs.append(dsobj1)
            if not starting_seqs:
                print("No starting sequences file provided. Random sequences will be generated.")
                
                starting_seqs = [dsobj1]
                starting_seqs = create_new_seqs(starting_seqs, args.pool_size, symmetry=sym_list)


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

    sym_list=None
    if args.symmetry:
        sym_list = args.symmetry.strip().split(",")

    if args.len_binder:
        len_binder=args.len_binder
    else:
        len_binder=0

    run_genetic_alg_gpus(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsize = args.pool_size, num_iter = args.num_iter, n_workers=args.num_gpus, write_raw_plddts=args.write_plddt, write_pair_confs=args.write_pae, stabilize_binder=args.stabilize_binder, write_pdbs=args.write_pdbs, mpnn_freq=args.mpnn_freq, symmetry=sym_list, len_binder=len_binder)

    
