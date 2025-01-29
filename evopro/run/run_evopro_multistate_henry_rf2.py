import importlib
import math
import sys, os
import json
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.utils.plot_scores import plot_scores_general_dev
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import read_starting_seqs, create_new_seqs, mutate_by_protein_mpnn
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.plots import get_chain_lengths_rf2, plot_pae, plot_plddt
from evopro.utils.utils import compressed_pickle
from evopro.utils.pdb_parser import get_coordinates_pdb, change_chainid_pdb, append_pdbs
#from evopro.run.run_evopro_multistate_henry import DesignSeqBidir

sys.path.append('/proj/kuhl_lab/RosettaFold2/RoseTTAFold2/run/')
from run_rf2 import rf2_init

import string
chain_names = list(string.ascii_uppercase)
all_aas = ["A", "C",  "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

class DesignSeqBidir(DesignSeq):
    """DesignSeq object adapted to follow bidirectional symmetry constraints"""
    
    def mutate(self, mut_percent = 0.125, num_mut_choice = [-1, 0, 1], var=0, var_weights = [0.8, 0.1, 0.1]):
        all_aas = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

        #calculating number of mutants from mut_percent and maybe adding or subtracting one mutant for stochasticity
        #num_mut = round(len(str(self))*mut_percent) + random.choice([-1, 0, 1])
        num_mut=round(mut_percent*len(self.mutable.keys()))+random.choice(num_mut_choice)
        if num_mut<1:
            num_mut=1

        new_mut = copy.deepcopy(self.mutable) 

        mut_ids = random.sample(list(new_mut.keys()), len(new_mut.keys()))
        i = 0
        num_mut_curr = 0
        mutated = []
        while num_mut_curr < num_mut:
            mut_id = mut_ids[i]
            weights = new_mut[mut_id][0]["weights"]
            chain = re.split('(\d+)', mut_id)[0]
            aa = int(re.split('(\d+)', mut_id)[1])

            method = "sub"
            if var > 0:
                #print(len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])), len(self.sequence.keys()), var)
                if len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) >= len(self.sequence.keys()) + var:
                    method = random.choices(["sub", "del"], [var_weights[0], var_weights[2]])[0]

                elif len("".join([self.jsondata["sequence"][chain] for chain in self.jsondata["sequence"]])) <= len(self.sequence.keys()) - var:
                    method = random.choices(["sub", "insert"], [var_weights[0], var_weights[1]])[0]

                else:
                    method = random.choices(["sub", "insert", "del"], var_weights)[0]

            #print("mutating by", method, str(var), str(var_weights))
            old_res = new_mut[mut_id][-1]
            new_res = copy.deepcopy(old_res)
            if method == "sub":
                new_aa = random.choices(all_aas, weights)[0]
                if new_res["resid"][2] < 0:
                    print("Trying to mutate by substitution at a deletion. Mutating by insertion instead.")
                    new_res["resid"][2] = new_res["resid"][2] + 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if substitution, replace dict with sub
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)
            elif method == "insert":
                new_aa = random.choices(all_aas, weights)[0]
                new_res["resid"][2] = new_res["resid"][2] + 1
                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                #if insertion, append to list of dicts
                new_mut[mut_id].append(new_res)
                mutated.append(mut_id)
            elif method == "del":
                new_aa = ""
                new_res = copy.deepcopy(old_res)

                if new_res["resid"][2] < 0:
                    print("Trying to mutate by deletion at a deletion. Mutating by insertion instead.")
                    new_aa = random.choices(all_aas, weights)[0]
                    new_res["resid"][2] = new_res["resid"][2] + 1
                else:
                    new_res["resid"][2] = new_res["resid"][2] - 1

                new_res["resid"][3] = new_aa
                num_mut_curr += 1
                # if deletion, replace dict with no-aa dict
                new_mut[mut_id][-1] = new_res
                mutated.append(mut_id)

            #print("after mutation", new_mut[mut_id])
            i+=1

        newseqobj = DesignSeqBidir(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric)
        newseqobj._update_symmetric_positions(mutated)
        #newseqobj._check_symmetry()
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        newseqobj._check_symmetry()
        return newseqobj

    def crossover(self, otherDS, ncross = 1):
        chains = []
        for resid in self.mutable:
            chain = re.split('(\d+)', resid)[0]
            if chain not in chains:
                chains.append(chain)

        crossover_chain = random.choice(chains)
        mutated = []
        new_mut = {}
        
        mut_seq = [x for x in self.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        other_mut_seq = [x for x in otherDS.mutable if re.split('(\d+)', x)[0] == crossover_chain]
        new_mut = copy.deepcopy(self.mutable)

        s=0
        points = random.sample(range(len(mut_seq)), ncross)
        for i, mut_id, other_mut_id in zip(range(len(mut_seq)), mut_seq, other_mut_seq):
            options = [self.mutable[mut_id], otherDS.mutable[other_mut_id]]
            if i in points:
                if s==0:
                    s=1
                else:
                    s=0
            new_mut[mut_id] = options[s]
            if s==1:
                mutated.append(mut_id)

        newseqobj = DesignSeqBidir(sequence=self.sequence, mutable=new_mut, symmetric=self.symmetric)
        newseqobj._update_symmetric_positions(mutated)
        newseqobj._update_sequence()
        newseqobj._create_jsondata()
        newseqobj._check_symmetry()
        return newseqobj

    def _get_bidirectional_matches(self, res1):
        """Retrieve all possible bidirectional matches based on resid"""
        res_aa = res1["resid"][3]
        res_aa_idx = ALPHABET.index(res_aa)
        matches = BIDIR_MATRIX[res_aa_idx]
        aa_matches = []

        for n, aa in enumerate(ALPHABET[:-1]):
            # adjust to add weights or allowed/disallowed residues
            valid_match = matches[n] * res1["weights"][n]
            if valid_match:
                aa_matches.append(aa)

        return aa_matches
    
    def _update_symmetric_positions(self, mutated):
        """Instead of traditional symmetry, this uses bidirectional symmetry rules.
        So 'symmetric' residues must be bidirectional codon/anticodon matches."""
        for mut_check in mutated:
            if mut_check in self.symmetric:
                for mut_id in self.symmetric[mut_check]:
                    chain = re.split('(\d+)', mut_id)[0]
                    index = int(re.split('(\d+)', mut_id)[1])
                    if len(self.mutable[mut_check]) != len(self.mutable[mut_id]):
                        self.mutable[mut_id] = copy.deepcopy(self.mutable[mut_check])
                        for pos in self.mutable[mut_id]:
                            pos["resid"][0] = chain
                            pos["resid"][1] = index

                    for res1, res2 in zip(self.mutable[mut_check], self.mutable[mut_id]):
                        # check if bidirectional match, if not, randomly choose a match
                        pairs = self._get_bidirectional_matches(res1)
                        if res2["resid"][3] not in pairs:
                            # update residue by randomly choosing a bidir pair match
                            res2["resid"][3] = random.choice(pairs)
                        if res1["resid"][2] != res2["resid"][2]:
                            # insertion/deletion code - leave as-is
                            res2["resid"][2] = res1["resid"][2]


def create_new_seqs_mpnn_henry_bidir(startseqs, scored_seqs, num_seqs, run_dir, iter_num, all_seqs = [], af2_preds=["AB", "B"],
                         mpnn_temp="0.1", mpnn_version="s_48_020", mpnn_chains=None, bias_AA_dict=None, bias_by_res_dict=None, 
                         bidir=False):
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
        results = mutate_by_protein_mpnn(pdb_dirs[k], startseqs[k], mpnn_temp, mpnn_version=mpnn_version, bias_AA_dict=bias_AA_dict, bias_by_res_dict=bias_by_res_dict, bidir=bidir)
        dsobj = startseqs[k]
        for result in results:
            inf +=1
            #print(pool, num_seqs)
            seq = result[-1][-1].strip().split("/")
            newseq_sequence = "".join(seq)
            newseq_sequence_check = ",".join(seq)
            newseqobj = DesignSeqBidir(seq=newseq_sequence, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
            # newseqobj = DesignSeq(seq=newseq_sequence, sequence=dsobj.sequence, mutable=dsobj.mutable, symmetric=dsobj.symmetric)
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


def run_genetic_alg_multistate_henry(run_dir, af2_flags_file, score_func, startingseqs, poolsizes = [],
                               num_iter = 50, n_workers=1, mut_percents=None, contacts=None,
                               mpnn_temps="0.1", mpnn_version="s_48_020", skip_mpnn=[], mpnn_iters=None, 
                               mpnn_chains=None, bias_AA_dict=None, bias_by_res_dict=None,
                               repeat_af2=True, af2_preds=[], crossover_percent=0.2, vary_length=0,
                               write_pdbs=False, plot=[], conf_plot=False, write_compressed_data=True, 
                               path_to_starting=None, bidir=False):

    num_af2=0
    if type(mpnn_temps) == str:
        mpnn_temps = ["0.1"] * num_iter
    
    lengths = []
    num_preds = len(af2_preds)
    print(af2_preds, num_preds)
    for chains in af2_preds:
        c = list(chains)
        lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])][0] for chain in c])

    print("Repeating RF2", repeat_af2)

    print("Initializing distributor")
    dist = Distributor(n_workers, rf2_init, af2_flags_file, lengths)

    scored_seqs = {}
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
            
    if not mpnn_chains:
        mpnn_chains = [af2_preds[0]]

    mpnn_counter = 0

    #start genetic algorithm iteration
    while curr_iter <= num_iter:
        all_seqs = []
        all_scores = []
        pool = newpool
        print("Current pool", pool)

        if curr_iter == 1:
            print("\nIteration 1: Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs(pool, 
                                   poolsizes[curr_iter-1], 
                                   crossover_percent=0, 
                                   all_seqs = list(scored_seqs.keys()), 
                                   vary_length=vary_length)

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters and curr_iter not in skip_mpnn:
            print("Iteration " + str(curr_iter) + ": refilling with ProteinMPNN at temperature " + str(mpnn_temps[mpnn_counter]) + ".")
            print('AF2 PREDS:', af2_preds, '\nMPNN CHAINS:', mpnn_chains)
            pool = create_new_seqs_mpnn_henry_bidir(pool, 
                                        scored_seqs, 
                                        poolsizes[curr_iter-1], 
                                        run_dir, 
                                        curr_iter, 
                                        all_seqs = list(scored_seqs.keys()), 
                                        af2_preds=af2_preds,
                                        mpnn_temp=mpnn_temps[mpnn_counter], 
                                        mpnn_version=mpnn_version, 
                                        mpnn_chains=mpnn_chains,
                                        bias_AA_dict=bias_AA_dict, 
                                        bias_by_res_dict=bias_by_res_dict, 
                                        bidir=bidir)
            mpnn_counter+=1
        
        #otherwise refilling pool with just mutations and crossovers
        else:
            print("\nIteration " + str(curr_iter) + ": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs(pool, 
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
        for results, dsobj in zip(results_packed, scoring_pool):
            if path_to_starting is not None:   
                scores.append(score_func(results, path_to_starting, dsobj))
            else:
                scores.append(score_func(results, dsobj, contacts=contacts))
        
        #adding sequences and scores into the dictionary
        for score, results, dsobj in zip(scores, results_packed, scoring_pool):
            key_seq = dsobj.get_sequence_string()
            #print(seqs)
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
                scored_seqs[key_seq]["data"].append({"score": score_all, "pdb": pdbs, "result": results})
                #don't need to sort when using all 5 for average score
                #scored_seqs[key_seq]["data"].sort(key=lambda x: x["score"][0][0])
                #print(scored_seqs[key_seq]["data"])
                #print(overall_scores)
                sum_score = [0 for _ in overall_scores]
                #print("Before", sum_score)
                for elem in scored_seqs[key_seq]["data"]:
                    #print("Element", elem, elem["score"])
                    #print(sum_score, elem["score"][0], elem["score"][1])
                    sum_score = [x+y for x,y in zip(sum_score, elem["score"][0])]
                avg_score = [x/len(scored_seqs[key_seq]["data"]) for x in sum_score]
                #print("After", sum_score, avg_score)
                scored_seqs[key_seq]["average"] = avg_score
            else:
                scored_seqs[key_seq] = {"data": [{"score": score_all, "pdb": pdbs, "result": results}]}
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
            else:
                sorted_scored_pool.append((key_seq, scored_seqs[key_seq]["average"]))
            
            if write_pdbs:
                print("Writing pdbs...")
                pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
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
        result = scored_seqs[key_seq]["data"][0]["result"][0]
        results = [scored_seqs[key_seq]["data"][x]["result"] for x in range(len(scored_seqs[key_seq]["data"]))]
        
        if write_compressed_data:
            for r, k in zip(results, range(len(results))):
                compressed_pickle(os.path.join(output_dir, "seq_" + str(j) + "_result_" + str(k)), r)

        if conf_plot:
            # Plot confidence.
            Ls = get_chain_lengths_rf2(result)

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
    
    print("Number of RosettaFold2 predictions: ", num_af2)
    dist.spin_down()
        
    plotting = plot_scores_general_dev(plot, seqs_per_iteration, output_dir)
    print("plots created at", str(plotting))


def getEvoProParserHenry() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run EvoPro. Adds bidirectional flag"""

    parser = FileArgumentParser(description='EvoPro runner script that can take '
                                'command-line arguments or read from an '
                                'argument flag file.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--input_dir',
                        default='.',
                        type=str,
                        help='Path to directory that contains input files. Default is ./')

    parser.add_argument('--num_iter',
                        default='50',
                        type=int,
                        help='Number of iterations of genetic algorithm. Default is 50.')

    parser.add_argument('--pool_size',
                        default='20',
                        type=int,
                        help='Size of "genetic pool", or the number of sequences evaluated per '
                        'iteration. Default is 20.')

    parser.add_argument('--pool_size_variable',
                        action='store_true',
                        help='Specify a file pool_sizes.txt with pool sizes for every iteration (or until pool size stops changing). '
                        'Defaults to False and uses constant pool size.')

    parser.add_argument('--num_gpus',
                        default='1',
                        type=int,
                        help='Number of gpus available. Default is 1.')

    parser.add_argument('--score_file',
                        default='',
                        type=str,
                        help='Path and file name of python script containing the score'
                        ' function used to evaluate fitness of the alphafold predictions. Required.')

    parser.add_argument('--score_func',
                        default='',
                        type=str,
                        help='Name of the score function in the score file'
                        ' used to evaluate fitness of the alphafold predictions. Required.')
    
    parser.add_argument('--score_func_2',
                        default=None,
                        type=str,
                        help='Name of the second score function in the score file'
                        ' used to evaluate fitness of the alphafold predictions after score_func_2 iterations has been reached. Optional.')
    
    parser.add_argument('--score_func_2_iteration',
                        default=30,
                        type=int,
                        help='Number of iterations after which scoring function is switched to score_func_2. Default is 30.')

    parser.add_argument('--rmsd_func',
                        default=None,
                        type=str,
                        help='Name of the rmsd function in the score file'
                        ' used to evaluate fitness of the alphafold predictions. Optional, requires stabilize_binder=True.')

    parser.add_argument('--rmsd_to_starting',
                        default=None,
                        type=str,
                        help='Name of the rmsd function in the score file and the path/filename to pdb for RMSD.'
                        ' used to evaluate rmsd to the starting scaffold using a U-shaped potential. ')
    
    parser.add_argument('--path_to_starting',
                        default=None,
                        type=str,
                        help='path/filename to pdb to pass to scoring function for RMSD to starting.')

    parser.add_argument('--define_contact_area',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be targeted for contacts. Default is None.')
    
    parser.add_argument('--bonus_contacts',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be given a bonus for making contacts, followed by the'
                        'distance cutoff. Default is None and 4A.')
    
    parser.add_argument('--penalize_contacts',
                        default=None,
                        type=str,
                        help='File defining residues on target interface to be given a penalty for making contacts, followed by the'
                        'distance cutoff. Default is None and 4A.')
    
    parser.add_argument('--no_repeat_af2',
                         action='store_true',
                         help='Use this flag to specify if you want AF2 to be multiple times on the same sequence, and the score averaged.'
                         'This means that all sequences will be rescored every iteration (like FoldDesign protocol) until each sequence has'
                         'been scored 5 times. Default is False.')

    parser.add_argument('--dont_write_compressed_data',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--write_pdbs',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--mpnn_iters',
                        default=None,
                        type=str,
                        help='Iteration numbers at which MPNN is used to refill the pool. Defaults to mpnn_freq')

    parser.add_argument('--mpnn_freq',
                        default='10',
                        type=int,
                        help='Protein MPNN is used to refill the pool once every _ iterations.'
                        'Default is 10.')

    parser.add_argument('--skip_mpnn',
                        default=None,
                        type=str,
                        help='Skip MPNN refilling in these iterations. Default is None.')
    
    parser.add_argument('--mpnn_temp',
                        default='0.1',
                        type=str,
                        help='Protein MPNN is used to refill the pool at this sampling temperature.'
                        'Default is 0.1.')
    
    parser.add_argument('--mpnn_temp_variable',
                        action='store_true',
                        help='Specify a file mpnn_temps.txt with temperatures for every call to MPNN (or until temp stops changing). '
                        'Defaults to False and uses constant MPNN temp.')
    
    parser.add_argument('--mpnn_version',
                        default="s_48_020",
                        type=str,
                        help='Model version used to run MPNN. Default is s_48_020 (soluble).')
    
    parser.add_argument('--mpnn_bias_AA',
                        default=None,
                        type=str,
                        help='Path to json file containing bias dictionary for MPNN. Default is None.')
    
    parser.add_argument('--mpnn_bias_by_res',
                        default=None,
                        type=str,
                        help='Path to json file containing per residue bias dictionary for MPNN. Default is None.')
    
    parser.add_argument('--mpnn_chains',
                        default=None,
                        type=str,
                        help='Chains concatenated into a single pdb for MPNN. Default is None and'
                        'it will use the first af2 prediction. Example: AB,B')

    parser.add_argument('--plot_confidences',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--plot_scores_avg',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--plot_scores_median',
                         action='store_true',
                         help='Default is False.')
    
    parser.add_argument('--plot_scores_top',
                         action='store_true',
                         help='Default is False.')

    parser.add_argument('--crossover_percent',
                        default='0.2',
                        type=float,
                        help='Fraction of pool refilled by crossover.'
                        'Default is 0.2.')

    parser.add_argument('--vary_length',
                        default='0',
                        type=int,
                        help='How much the length is allowed to vary. Default is 0.')
    parser.add_argument('--substitution_insertion_deletion_weights',
                        default=None,
                        type=str,
                        help='Specify probability of substitutions, insertions, and deletions (in that order) during mutation. Default is 0.8,0.1,0.1')

    parser.add_argument('--mutation_percents',
                        default='0.125',
                        type=str,
                        help='Number of mutations made as a percentage of sequence length.'
                        'Default is 0.125 for every iteration. If more than one value is provided, number of iterations will be split evenly and assigned.')
    
    parser.add_argument('--force_single_mutation_only',
                        action='store_true',
                        help='Default is False.')
    
    parser.add_argument('--af2_preds',
                        default="AB",
                        type=str,
                        help='Chain ID permutations to run through individual AF2 runs, separated by commas. Only used for multistate design. Default is None.')
    
    parser.add_argument('--af2_preds_extra',
                        default=None,
                        type=str,
                        help='Chain ID permutations to run through individual AF2 runs, separated by commas. Default is None.')

    parser.add_argument('--bidirectional', 
                        action='store_true', 
                        help='Enable bidirectional coding constraint in MPNN steps. Default is False.')
    return parser


if __name__ == "__main__":
    parser = getEvoProParserHenry()
    args = parser.parse_args(sys.argv[1:])

    input_dir = args.input_dir
    if not input_dir.endswith("/"):
        input_dir = input_dir + "/"

    if input_dir not in ['', None]:
            onlyfiles = [f for f in sorted(os.listdir(input_dir)) if os.path.isfile(
                         os.path.join(input_dir, f))]

            flagsfile=None
            resfile=None
            starting_seqs=None
            for filename in onlyfiles:
                if "flag" in filename and "rf2" in filename:
                    flagsfile = filename
                    break

            for filename in onlyfiles:
                if "residue" in filename and "spec" in filename:
                    resfile = filename
                    break

            if flagsfile is None:
                raise ValueError("Flags file for Rosettafold runs not provided.")
            if resfile is None:
                raise ValueError("Please provide a residue specifications file.")

            #print(input_dir + resfile)
            # dsobj1 = DesignSeq(jsonfile=input_dir + resfile)
            dsobj1 = DesignSeqBidir(jsonfile=input_dir + resfile)

            for filename in onlyfiles:
                if "starting" in filename and "seqs" in filename:
                    if "mut" in filename:
                        starting_seqs = read_starting_seqs(input_dir + filename, dsobj1, mut_only=True)
                        starting_seqs.append(dsobj1)
                    else:
                        starting_seqs = read_starting_seqs(input_dir + filename, dsobj1)
                        starting_seqs.append(dsobj1)
                    print('Starting sequences successfully parsed from %s' % filename)
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
            
            mpnnfile=None
            if args.mpnn_temp_variable:
                for filename in onlyfiles:
                    if "mpnn" in filename and "temp" in filename:
                        mpnnfile = filename
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

    print("Writing compressed data:", not args.dont_write_compressed_data)
    print("Writing pdb files every iteration:", args.write_pdbs)
    print("Varying protein length:", args.vary_length>0)
    print("Repeating RF2:", not args.no_repeat_af2)
    
    af2_preds = args.af2_preds.strip().split(",")
    if args.mpnn_chains:
        mpnn_chains = args.mpnn_chains.strip().split(",")
    else:
        mpnn_chains = None
    
    mpnn_temps = []
    if mpnnfile:
        with open(mpnnfile, "r") as mpf:
            for lin in mpf:
                mpnn_temps.append(lin.strip())
    else:
        for i in range(len(mpnn_iters)):
            mpnn_temps.append(args.mpnn_temp)
    
    while len(mpnn_temps) < len(mpnn_iters):
        mpnn_temps.append(mpnn_temps[-1])
        
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
                
    print(bias_by_res_dict)
    
    if args.bidirectional:
        print('Bidirectional coding enabled for MPNN steps!')
                                   
    run_genetic_alg_multistate_henry(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, mut_percents=mut_percents, contacts=contacts, 
        mpnn_temps=mpnn_temps, mpnn_version=args.mpnn_version, skip_mpnn=mpnn_skips, mpnn_iters=mpnn_iters, 
        mpnn_chains=mpnn_chains, bias_AA_dict=bias_AA_dict, bias_by_res_dict=bias_by_res_dict, 
        repeat_af2=not args.no_repeat_af2, af2_preds = af2_preds, crossover_percent=args.crossover_percent, vary_length=args.vary_length, 
        write_pdbs=args.write_pdbs, plot=plot_style, conf_plot=args.plot_confidences, write_compressed_data=not args.dont_write_compressed_data, 
        path_to_starting=path_to_starting, bidir=args.bidirectional)
