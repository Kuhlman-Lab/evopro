#import standard packages
import multiprocessing as mp
import collections
import random
import string
import time
import sys, os
import re
from DesignSeq import DesignSeq
from typing import Sequence, Union

#import custom packages
from folddesign.utils.distributor import Distributor

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2
from functools import partial
import numpy as np
all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def af2_init(proc_id: int, arg_file: str, lengths: Sequence[Union[str, Sequence[str]]], fitness_fxn):
    print('initialization of process', proc_id)

    os.environ['TF_FORCE_UNITED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '2.0'
    os.environ['TF_XLA_FLAGS'] = '--tf_xla_cpu_global_jit'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(proc_id)

    import jax
    from features import (getRawInputs, getChainFeatures, getInputFeatures)
    from setup import (getAF2Parser, QueryManager, getOutputDir)
    from model import (getModelNames, getModelRunner, predictStructure, getRandomSeeds)

    parser = getAF2Parser()
    args = parser.parse_args([f'@{arg_file}'])    

    output_dir = getOutputDir(out_dir=args.output_dir)

    print(f'lengths: {lengths}')
    sequences = []

    for chain in lengths:
        sequences.append("".join(random.choices(all_aas, k=chain)))

    sequences = [sequences]

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
        output_dir=output_dir)

    print('finished generating raw inputs')
    print(raw_inputs)

    model_names = getModelNames(
        first_n_seqs=len(queries[0][1]),
        last_n_seqs=len(queries[-1][1]),
        use_ptm=args.use_ptm, num_models=args.num_models)

    query_features = []
    for query in queries:
        sequences = query[1]
    
        features_for_chain = getChainFeatures(
            sequences=sequences,
            raw_inputs=raw_inputs,
            use_templates=args.use_templates)

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
                model_runner=model_runner,
                feature_dict=input_features,
                run_multimer=run_multimer)
            print(f'Model {model_name} took {time.time()-t} sec on GPU {proc_id}.')

        models.append((model_name, model_runner))

    af2_partial = partial(af2, arg_file=arg_file, proc_id=proc_id, fitness_fxn=fitness_fxn, compiled_runners=models)

    return af2_partial

def af2_init_2models(proc_id: int, arg_file: str, lengths: Sequence[Union[str, Sequence[str]]], fitness_fxn):
    print('initialization of process', proc_id)

    os.environ['TF_FORCE_UNITED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '2.0'
    os.environ['TF_XLA_FLAGS'] = '--tf_xla_cpu_global_jit'
    os.environ['CUDA_VISIBLE_DEVICES'] = str(proc_id)

    import jax
    from features import (getRawInputs, getChainFeatures, getInputFeatures)
    from setup import (getAF2Parser, QueryManager, getOutputDir)
    from model import (getModelNames, getModelRunner, predictStructure, getRandomSeeds)
    from utils.query_utils import generate_random_sequences

    parser = getAF2Parser()
    args = parser.parse_args([f'@{arg_file}'])

    output_dir = getOutputDir(out_dir=args.output_dir)

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
        output_dir=output_dir)

    print('finished generating raw inputs')
    print(raw_inputs)

    model_names = getModelNames(
        first_n_seqs=len(queries[0][1]),
        last_n_seqs=len(queries[-1][1]),
        use_ptm=args.use_ptm, num_models=args.num_models)
    print('model names', model_names)

    query_features = []
    for query in queries:
        sequences = query[1]

        features_for_chain = getChainFeatures(
            sequences=sequences,
            raw_inputs=raw_inputs,
            use_templates=args.use_templates)

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
                model_runner=model_runner,
                feature_dict=input_features,
                run_multimer=run_multimer)
            print(f'Model {model_name} took {time.time()-t} sec on GPU {proc_id}.')

        models.append((model_name, model_runner))

    af2_partial = partial(af2, arg_file=arg_file, proc_id=proc_id, fitness_fxn=fitness_fxn, compiled_runners=models)

    return af2_partial

def generate_random_seqs(num_seqs, lengths):
    oplist = []
    chainnames = list(string.ascii_uppercase)
    for i in range(num_seqs):
        sequences = []
        for chain in lengths:
            sequences.append("".join(random.choices(all_aas, k=chain)))
        newseq = {chainnames[i]:seq for i,seq in zip(range(len(sequences)), sequences)}
        oplist.append(DesignSeq(seq=newseq))

    return oplist

#NEEDS REWRITE
def read_starting_seqs(seqfile, dsobj, mut_only=True):

    if mut_only:
        newseqs = []
        seqs = []
        with open(seqfile,"r") as seqf:
            for l in seqf:
                seqs.append(l.strip())
        for seq in seqs:
            newseqforobj = {chain: list(dsobj.seq[chain]) for chain in dsobj.seq}
            for mut_id, aa_sub in zip(dsobj.mutable, seq):
                chain = re.split('(\d+)', mut_id)[0]
                aa = int(re.split('(\d+)', mut_id)[1])
                newseqforobj[chain][aa-1] = aa_sub
            newseq_joined = {chain: "".join(newseqforobj[chain]) for chain in newseqforobj}
            newseqobj = DesignSeq(seq=newseq_joined, mutable=dsobj.mutable)
            newseqs.append(newseqobj)    

    else:#THIS PART DOES NOT WORK YET
        with open(jsonfile, "r") as inpf:
            data = json.load(inpf)

    return newseqs

def create_new_seqs(startseqs, num_seqs, crossoverpercent = 0.2, symmetry=None):
    """takes starting sequences and creates a pool of size num_seqs by mutation and crossover"""
    pool = startseqs.copy()

    while len(pool)<num_seqs*(1-crossoverpercent):
        obj = random.choice(startseqs)
        newseq = obj.mutate(symmetry=symmetry)
        if newseq not in pool:
            pool.append(newseq)
    #print(pool)
    while len(pool)<num_seqs:
        if len(startseqs) > 1:
            oldseqs = random.sample(startseqs, 2)
        else:
            oldseqs = random.sample(pool, 2)
        newseq = oldseqs[0].crossover(oldseqs[1])
        if newseq not in pool:
            pool.append(newseq)

    return pool

def mutate_by_protein_mpnn(pdb_dir, dsobj):
    sys.path.append('/proj/kuhl_lab/proteinmpnn/run/')
    from run_protein_mpnn import run_protein_mpnn_func
    results = run_protein_mpnn_func(pdb_dir, dsobj.jsondata)

    return results

def create_new_seqs_mpnn(startseqs, scored_seqs, num_seqs, run_dir, iter_num):
    pool = startseqs.copy()
    exampleds = None
    # create a directory within running dir to run protein mpnn
    output_folder = run_dir + "MPNN_" + str(iter_num) + "/"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    for dsobj, j in zip(startseqs, range(len(startseqs))):
        key_seq = " ".join([dsobj.seq[chain] for chain in dsobj.seq])
        pdb = scored_seqs[key_seq][1][1][0]
        exampleds = dsobj
        with open(output_folder + "seq_" + str(j) + ".pdb", "w") as pdbf:
            pdbf.write(str(pdb))

    print(exampleds)
    results = mutate_by_protein_mpnn(output_folder, startseqs[0])
    for result in results:
        seq = result[-1][-1].strip().split("/")
        seq_dict={}
        for chain, chain_seq in zip(exampleds.seq, seq):
            seq_dict[chain] = chain_seq
        newseqobj = DesignSeq(seq=seq_dict, mutable=dsobj.mutable)
        pool.append(newseqobj)

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
