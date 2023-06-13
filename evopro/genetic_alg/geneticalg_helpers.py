#import standard packages
import multiprocessing as mp
import collections
import random
import string
import time
import sys, os
import re
import json
from typing import Sequence, Union

#import custom packages
from evopro.utils.distributor import Distributor
from evopro.genetic_alg.DesignSeq import DesignSeq, DesignSeqMSD
from evopro.utils.pdb_parser import get_coordinates_pdb, change_chainid_pdb, append_pdbs

sys.path.append('/proj/kuhl_lab/alphafold/run')
from functools import partial
import numpy as np

chain_names = list(string.ascii_uppercase)
all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def of_init(proc_id: int, arg_file: str, lengths: Sequence[Union[str, Sequence[str]]]):
    sys.path.append('/proj/kuhl_lab/OmegaFold/')
    from omegafold.__main__ import main 
    from omegafold import pipeline
    import omegafold as of
    print('initialization of process', proc_id)

    args, state_dict, forward_config = pipeline.get_args(arg_file, proc_id)

    # get the model
    print("constructing omegafold...")
    model = of.OmegaFold(of.make_config())
    if "model" in state_dict:
        state_dict = state_dict.pop("model")
    model.load_state_dict(state_dict)
    model.eval()
    model.to(args.device)

    of_partial = partial(main, arg_file=arg_file, proc_id=proc_id, model=model)

    return of_partial

def af2_init(proc_id: int, arg_file: str, lengths: Sequence[Union[str, Sequence[str]]]):
    from run_af2 import af2
    print('initialization of process', proc_id)

    os.environ['TF_FORCE_UNITED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '2.0'
    os.environ['TF_XLA_FLAGS'] = '--tf_xla_cpu_global_jit'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(proc_id)

    import jax
    from features import (getRawInputs, getChainFeatures, getInputFeatures)
    from run.setup import (getAF2Parser, QueryManager, getOutputDir)
    from model import (getModelNames, getModelRunner, predictStructure, getRandomSeeds)
    from utils.query_utils import generate_random_sequences

    parser = getAF2Parser()
    args = parser.parse_args([f'@{arg_file}'])

    output_dir = getOutputDir(out_dir=args.output_dir, suffix="_"+str(proc_id))

    print(f'lengths: {lengths}')
    sequences = []
    #sequences = generate_random_sequences(lengths, 1, aalist=all_aas)[0]

    for length in lengths:
        if type(length) is list:
            sequence = []
            for chain in length:
                sequence.append("".join(random.choices(all_aas, k=chain)))
        else:
            sequence = "".join(random.choices(all_aas, k=length))
        sequences.append(sequence)

    #sequences = [sequences]

    print(f'determined sequences: {sequences}')

    qm = QueryManager(
        sequences=sequences,
        min_length=args.min_length,
        max_length=args.max_length,
        max_multimer_length=args.max_multimer_length)
    qm.parse_files()
    qm.parse_sequences()
    queries = qm.queries
    del qm
    print("finished query parsing", queries)

    raw_inputs = getRawInputs(
        queries=queries,
        msa_mode=args.msa_mode,
        use_templates=args.use_templates,
        output_dir=output_dir,
        proc_id=proc_id)

    print('finished generating raw inputs')
    print(raw_inputs)

    model_names = getModelNames(
        first_n_seqs=len(queries[0][1]),
        last_n_seqs=len(queries[-1][1]),
        use_ptm=args.use_ptm, num_models=args.num_models,
        use_multimer=not args.no_multimer_models,
        use_multimer_v1=args.use_multimer_v1,
        use_multimer_v2=args.use_multimer_v2)
    print('model names', model_names)

    query_features = []
    for query in queries:
        sequences = query[1]

        features_for_chain = getChainFeatures(
            sequences=sequences,
            raw_inputs=raw_inputs,
            use_templates=args.use_templates,
            proc_id=proc_id)

        input_features = getInputFeatures(
            sequences=sequences,
            chain_features=features_for_chain)

        query_features.append( (sequences, input_features) )

    results_list = []
    models = []
    for model_name in model_names:
        model_runner = getModelRunner(
            model_name=model_name,
            num_ensemble=args.num_ensemble,
            is_training=args.is_training,
            num_recycle=args.max_recycle,
            recycle_tol=args.recycle_tol,
            params_dir=args.params_dir)

        run_multimer = False
        if 'multimer' in model_name:
            run_multimer = True

        for query in query_features:
            sequences = query[0]

            if len(sequences) > 1 and not run_multimer:
                continue
            elif len(sequences) == 1 and run_multimer:
                continue

            input_features = query[1]

            t = time.time()
            result = predictStructure(
                model_name=model_name,
                model_runner=model_runner,
                feature_dict=input_features,
                run_multimer=run_multimer,
                use_templates=args.use_templates)
            print(f'Model {model_name} took {time.time()-t} sec on GPU {proc_id}.')

        models.append((model_name, model_runner))

    af2_partial = partial(af2, arg_file=arg_file, proc_id=proc_id, compiled_runners=models)

    return af2_partial

#NEEDS REWRITE
def generate_random_seqs(num_seqs, lengths):
    oplist = []
    for i in range(num_seqs):
        sequences = []
        for chain, j in zip(lengths, range(len(lengths))):
            sequences.append("".join(random.choices(all_aas, k=chain)))
        newseq = {chain_names[k]:seq for k,seq in zip(range(len(sequences)), sequences)}
        oplist.append(DesignSeq(seq=newseq))

    return oplist

def read_starting_seqs(seqfile, dsobj):
    seqs = []
    newseqs = []
    with open(seqfile, "r") as inpf:
        for l in inpf:
            seqs.append(l.strip().split(","))
    for seq in seqs:
        newseq = "".join(seq)
        newseqobj = DesignSeq(seq=newseq, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
        newseqs.append(newseqobj)

    return newseqs

def create_new_seqs(startseqs, num_seqs, crossover_percent = 0.2, vary_length=0, sid_weights=[0.8, 0.1, 0.1], mut_percent=0.125, all_seqs = []):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""

    if not all_seqs:
        all_seqs = []
    pool = startseqs.copy()
    max_allowed_length = len(pool[0].sequence.keys()) + vary_length
    min_allowed_length = len(pool[0].sequence.keys()) - vary_length
    #print(mut_percent, crossover_percent, max_allowed_length, min_allowed_length, vary_length)

    #filling pool by mutation
    while len(pool)<num_seqs*(1-crossover_percent):
        obj = random.choice(startseqs)
        newseq = obj.mutate(var=vary_length, var_weights = sid_weights, mut_percent=mut_percent)
        len_new = len("".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]))
        newseq_sequence = ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]])
        if len_new<=max_allowed_length and len_new >= min_allowed_length:
            if newseq not in pool:
                if newseq_sequence not in all_seqs:
                    pool.append(newseq)

    #filling rest of pool by crossover
    crossover_loop_count = 0
    while len(pool)<num_seqs:
        if len(startseqs) > 1:
            oldseqs = random.sample(startseqs, 2)
        else:
            oldseqs = random.sample(pool, 2)
        newseq = oldseqs[0].crossover(oldseqs[1])
        crossover_loop_count+=1
        if crossover_loop_count>100:
            print("crossovers are failing to create new sequences. defaulting to mutation.")
            newseq = oldseqs[0].mutate(var=vary_length, var_weights = sid_weights, mut_percent=mut_percent)
        len_new = len("".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]))
        newseq_sequence = ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]])
        if len_new<=max_allowed_length and len_new >= min_allowed_length and newseq not in pool and newseq_sequence not in all_seqs:
            pool.append(newseq)

    return pool

def create_new_seqs_henry(startseqs, num_seqs, crossover_percent = 0.2, vary_length=0, sid_weights=[0.8, 0.1, 0.1], mut_percent=0.125, all_seqs = []):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""
    print('**creating new seqs by mutation/crossover')

    if not all_seqs:
        all_seqs = []
    pool = startseqs.copy()
    max_allowed_length = len(pool[0].sequence.keys()) + vary_length
    min_allowed_length = len(pool[0].sequence.keys()) - vary_length
    #print(mut_percent, crossover_percent, max_allowed_length, min_allowed_length, vary_length)
    print('**length of sequences being generated:', max_allowed_length)

    #filling pool by mutation
    while len(pool)<num_seqs*(1-crossover_percent):
        obj = random.choice(startseqs)
        # Note: mutation already preserves symmetry by default
        newseq = obj.mutate(var=vary_length, var_weights=sid_weights, mut_percent=mut_percent)
        len_new = len("".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]))
        newseq_sequence = ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]])
        if len_new<=max_allowed_length and len_new >= min_allowed_length:
            if newseq not in pool:
                if newseq_sequence not in all_seqs:
                    pool.append(newseq)

    #filling rest of pool by crossover
    crossover_loop_count = 0
    while len(pool)<num_seqs:
        if len(startseqs) > 1:
            oldseqs = random.sample(startseqs, 2)
        else:
            oldseqs = random.sample(pool, 2)
        # Note: crossover should also preserve symmetry by default
        newseq = oldseqs[0].crossover(oldseqs[1])
        crossover_loop_count+=1
        if crossover_loop_count>100:
            print("crossovers are failing to create new sequences. defaulting to mutation.")
            newseq = oldseqs[0].mutate(var=vary_length, var_weights = sid_weights, mut_percent=mut_percent)
        len_new = len("".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]))
        newseq_sequence = ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]])
        if len_new<=max_allowed_length and len_new >= min_allowed_length and newseq not in pool and newseq_sequence not in all_seqs:
            pool.append(newseq)

    return pool


def mutate_by_protein_mpnn(pdb_dir, dsobj, mpnn_temp, mpnn_version="s_48_020", bidir=False):
    sys.path.append('/proj/kuhl_lab/proteinmpnn/run/')
    from run_protein_mpnn import run_protein_mpnn_func
    #print(pdb_dir, dsobj.jsondata, mpnn_temp, mpnn_version, bidir)
    results = run_protein_mpnn_func(pdb_dir, json.dumps(dsobj.jsondata), sampling_temp=mpnn_temp, model_name=mpnn_version, bidir=bidir)

    return results

def create_new_seqs_mpnn(startseqs, scored_seqs, num_seqs, run_dir, iter_num, all_seqs = [], af2_preds=["AB", "B"],
                         mpnn_temp="0.1", mpnn_version="s_48_020", mpnn_chains=None):
    pool = startseqs.copy()
    pdb_dirs = []
    #create a directory within running dir to run protein mpnn
    output_folder = run_dir + "MPNN_" + str(iter_num) + "/"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    for dsobj, j in zip(startseqs, range(len(startseqs))):
        key_seq = dsobj.get_sequence_string()
        #print(key_seq)
        if mpnn_chains:
            pdb = None
            for i, pred in zip(range(len(af2_preds)), af2_preds):
                if pred in mpnn_chains:
                    if not pdb:
                        pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                    else:
                        chains, residues, resindices = get_coordinates_pdb(pdb)
                        new_pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                        chains_new, residues_new, resindices_new = get_coordinates_pdb(new_pdb)
                        chains_mod = []
                        chain_ind = 0
                        while len(chains_mod)<len(chains_new):
                            if chain_names[chain_ind] not in chains:
                                    chains_mod.append(chain_names[chain_ind])
                            chain_ind+=1
                        #print(chains, chains_new, chains_mod)
                        for cm, cn in zip(chains_mod, chains_new):
                            new_pdb = change_chainid_pdb(new_pdb, old_chain=cn, new_chain=cm)
                        
                        pdb = append_pdbs(pdb, new_pdb)
                            
        else:
            pdb = scored_seqs[key_seq]["data"][0]["pdb"][0]
            
        #pdb = scored_seqs[key_seq]["data"][0]["pdb"][0]
        pdb_dir = output_folder + "seq_" + str(j) + "/"
        if not os.path.isdir(pdb_dir):
            os.makedirs(pdb_dir)
        with open(output_folder + "seq_" + str(j) + "/seq_" + str(j) + ".pdb", "w") as pdbf:
            pdbf.write(str(pdb))
        pdb_dirs.append(pdb_dir)

    k=0
    inf = 0
    while len(pool) < num_seqs:
        results = mutate_by_protein_mpnn(pdb_dirs[k], startseqs[k], mpnn_temp, mpnn_version=mpnn_version)
        dsobj = startseqs[k]
        for result in results:
            inf +=1
            #print(pool, num_seqs)
            seq = result[-1][-1].strip().split("/")
            newseq_sequence = "".join(seq)
            newseq_sequence_check = ",".join(seq)
            newseqobj = DesignSeq(seq=newseq_sequence, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
            #print(newseq_sequence_check, all_seqs)
            if newseq_sequence_check not in all_seqs:
                pool.append(newseqobj)
            if inf>=50:
                print("too many mpnn runs without generating a new sequence, using random mutation")
                newseq = dsobj.mutate()
                if ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]) not in all_seqs:
                    pool.append(newseq)
        k+=1
        if k>=len(startseqs):
            k=0

    return pool

def create_new_seqs_mpnn_henry(startseqs, scored_seqs, num_seqs, run_dir, iter_num, all_seqs = [], af2_preds=["AB", "B"],
                         mpnn_temp="0.1", mpnn_version="s_48_020", mpnn_chains=None, bidir=False):
    print('**creating new seqs with MPNN')
    print('AF2:', af2_preds, 'MPNN CHAINS:', mpnn_chains)
    pool = startseqs.copy()
    pdb_dirs = []
    #create a directory within running dir to run protein mpnn
    output_folder = run_dir + "MPNN_" + str(iter_num) + "/"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    for dsobj, j in zip(startseqs, range(len(startseqs))):
        key_seq = dsobj.get_sequence_string()
        print(len(key_seq))
        if mpnn_chains is not None:
            pdb = None
            for i, pred in zip(range(len(af2_preds)), af2_preds):
                if pred in mpnn_chains:
                    if not pdb:
                        pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                    else:
                        chains, residues, resindices = get_coordinates_pdb(pdb)
                        new_pdb = scored_seqs[key_seq]["data"][0]["pdb"][i]
                        chains_new, residues_new, resindices_new = get_coordinates_pdb(new_pdb)
                        chains_mod = []
                        chain_ind = 0
                        while len(chains_mod)<len(chains_new):
                            if chain_names[chain_ind] not in chains:
                                    chains_mod.append(chain_names[chain_ind])
                            chain_ind+=1
                        print(chains, chains_new, chains_mod, '*****')
                        for cm, cn in zip(chains_mod, chains_new):
                            new_pdb = change_chainid_pdb(new_pdb, old_chain=cn, new_chain=cm)
                        
                        pdb = append_pdbs(pdb, new_pdb)
                            
        else:
            pdb = scored_seqs[key_seq]["data"][0]["pdb"][0]
            
        #pdb = scored_seqs[key_seq]["data"][0]["pdb"][0]
        pdb_dir = output_folder + "seq_" + str(j) + "/"
        if not os.path.isdir(pdb_dir):
            os.makedirs(pdb_dir)
        with open(output_folder + "seq_" + str(j) + "/seq_" + str(j) + ".pdb", "w") as pdbf:
            pdbf.write(str(pdb))
        pdb_dirs.append(pdb_dir)

    k=0
    inf = 0
    while len(pool) < num_seqs:
        print('JSON DATA FOR MPNN:', json.dumps(startseqs[k].jsondata))
        results = mutate_by_protein_mpnn(pdb_dirs[k], startseqs[k], mpnn_temp, mpnn_version=mpnn_version, bidir=bidir)
        dsobj = startseqs[k]
        for result in results:
            inf +=1
            print(pool, num_seqs)
            seq = result[-1][-1].strip().split("/")
            newseq_sequence = "".join(seq)
            newseq_sequence_check = ",".join(seq)
            newseqobj = DesignSeqMSD(seq=newseq_sequence, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric, jdata=dsobj.jsondata)
            print(newseq_sequence_check, all_seqs)
            if newseq_sequence_check not in all_seqs:
                pool.append(newseqobj)
            if inf>=50:
                print("too many mpnn runs without generating a new sequence, using random mutation")
                newseq = dsobj.mutate()
                if ",".join([newseq.jsondata["sequence"][chain] for chain in newseq.jsondata["sequence"]]) not in all_seqs:
                    pool.append(newseq)
        k+=1
        if k>=len(startseqs):
            k=0

    return pool


def create_new_seqs_mpnn_old(startseqs, scored_seqs, num_seqs, run_dir, iter_num, all_seqs = [], mpnn_temp="0.1", mpnn_version="s_48_020"):
    pool = startseqs.copy()
    exampleds = None
    pdb_dirs = []
    
    #create a directory within running dir to run protein mpnn
    output_folder = run_dir + "MPNN_" + str(iter_num) + "/"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    for dsobj, j in zip(startseqs, range(len(startseqs))):
        key_seq = dsobj.get_sequence_string()
        #print(key_seq)
        pdb = scored_seqs[key_seq]["pdb"][0]
        pdb_dir = output_folder + "seq_" + str(j) + "/"
        if not os.path.isdir(pdb_dir):
            os.makedirs(pdb_dir)
        with open(output_folder + "seq_" + str(j) + "/seq_" + str(j) + ".pdb", "w") as pdbf:
            pdbf.write(str(pdb))
        pdb_dirs.append(pdb_dir)

    k=0
    while len(pool) < num_seqs:
        results = mutate_by_protein_mpnn(pdb_dirs[k], startseqs[k], mpnn_temp, mpnn_version=mpnn_version)
        dsobj = startseqs[k]
        for result in results:
            print(pool, num_seqs)
            seq = result[-1][-1].strip().split("/")
            newseq_sequence = "".join(seq)
            newseq_sequence_check = ",".join(seq)
            newseqobj = DesignSeq(seq=newseq_sequence, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
            if newseq_sequence_check not in all_seqs:
                pool.append(newseqobj)
        k+=1
        if k>=len(startseqs):
            k=0

    return pool

def write_outputs(seqs, scores, opfilename):
    with open(opfilename, "w") as opf:
        for seq, score in zip(seqs, scores):
            opf.write(str(seq) + "\t" + str(score) + "\n")

    print("done writing ", opfilename)


if __name__ == "__main__":
    rand_seqs = generate_random_seqs(20, [59])
    #print(rand_seqs)
    #path = "/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/tests/pd1_threehelix_test/"
    #dsobj = DesignSeq(jsonfile=path+"resfile.json")
    #print(read_starting_seqs(path+"starting_seqs.txt", dsobj))
