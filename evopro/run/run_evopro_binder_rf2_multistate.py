import importlib
import math
import sys, os
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.genetic_alg.DesignSeq import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.run.generate_json import parse_mutres_input
from evopro.genetic_alg.geneticalg_helpers import read_starting_seqs, create_new_seqs, create_new_seqs_mpnn
from evopro.user_inputs.inputs import getEvoProParser
from evopro.utils.plots import get_chain_lengths_rf2, plot_pae, plot_plddt
from evopro.utils.utils import compressed_pickle
from evopro.utils.plot_scores import plot_scores_general_dev

sys.path.append('/proj/kuhl_lab/RosettaFold2/RoseTTAFold2/run/')
from run_rf2 import rf2_init

def run_genetic_alg_multistate(run_dir, af2_flags_file, score_func, startingseqs, af2_flags_file_2 = None, poolsizes = [], 
                               num_iter = 50, n_workers=1, mut_percents=None, contacts=None, distance_cutoffs=None,
                               score_func_2=None, score_func_2_iter=30, mpnn_chains=None,
                               mpnn_temps=None, mpnn_version="s_48_020", skip_mpnn=[], mpnn_iters=None, 
                               repeat_af2=True, af2_preds=[], crossover_percent=0.2, vary_length=0, 
                               write_pdbs=False, plot=[], conf_plot=False, write_compressed_data=True):

    num_af2=0
    
    print("Repeating RF2", repeat_af2)
    
    lengths = []
    num_preds = len(af2_preds)
    print(af2_preds, num_preds)
    for chains in af2_preds:
        c = list(chains)
        lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])][0] for chain in c])

    print(lengths)

    print("Initializing distributor")
    dist = Distributor(n_workers, rf2_init, af2_flags_file, lengths)
    
    if af2_flags_file_2:
        print("Initializing distributor 2")
        dist2 = Distributor(n_workers, rf2_init, af2_flags_file_2, lengths)

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
        pdb_folder = output_dir + "pdbs_per_iter/"
        if not os.path.isdir(pdb_folder):
            os.makedirs(pdb_folder)

    mpnn_counter = 0
    #start genetic algorithm iteration
    while curr_iter <= num_iter:

        pool = newpool

        if curr_iter == 1:
            print("Iteration 1: Creating new sequences.")
            #do not use crossover to create the initial pool
            pool = create_new_seqs(pool, poolsizes[curr_iter-1], crossover_percent=0, all_seqs = list(scored_seqs.keys()), vary_length=vary_length)

        #using protein mpnn to refill pool when specified
        elif curr_iter in mpnn_iters and curr_iter not in skip_mpnn:
            print("Iteration " + str(curr_iter) + ": refilling with ProteinMPNN at temperature " + str(mpnn_temps[mpnn_counter]) + ".")
            pool = create_new_seqs_mpnn(pool, 
                                        scored_seqs, 
                                        poolsizes[curr_iter-1],
                                        run_dir, 
                                        curr_iter, 
                                        all_seqs = list(scored_seqs.keys()), 
                                        af2_preds=af2_preds,
                                        mpnn_temp=mpnn_temps[mpnn_counter], 
                                        mpnn_version=mpnn_version, 
                                        mpnn_chains=mpnn_chains)
            mpnn_counter+=1
        #otherwise refilling pool with just mutations and crossovers
        else:
            print("Iteration " + str(curr_iter) + ": refilling with mutation and " + str(crossover_percent*100) + "% crossover.")
            pool = create_new_seqs(pool, poolsizes[curr_iter-1], crossover_percent=crossover_percent, mut_percent=mut_percents[curr_iter-1], all_seqs = list(scored_seqs.keys()), vary_length=vary_length)
        
        if repeat_af2:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in repeat_af2_seqs]
        else:
            scoring_pool = [p for p in pool if p.get_sequence_string() not in scored_seqs]
        
        if af2_flags_file_2:
            print("Detected 2nd AF2 flags file. Running 2nd AF2 predictions.")
            work_list_1 = []
            work_list_2 = []
            results_packed = []
            for p in scoring_pool:
                work_list_1.append([[p.jsondata["sequence"][chain] for chain in af2_preds[0]]])
                work_list_2.append([[p.jsondata["sequence"][chain] for chain in af2_preds[1]]])
                if len(af2_preds) > 2:
                    for c in af2_preds[2:]:
                        work_list_2.append([[p.jsondata["sequence"][chain] for chain in af2_preds[2]]])
            num_af2 += len(work_list_1)
            num_af2 += len(work_list_2)
            print("work list 1", work_list_1)
            print("work list 2", work_list_2)
            results1 = dist.churn(work_list_1)
            results2 = dist2.churn(work_list_2)
            results2_packed = [results2[i:i+(len(af2_preds)-1)] for i in range(0, len(results2), len(af2_preds)-1)]

            for result1, result2, in zip(results1, results2_packed):
                while type(result1) == list:
                    result1 = result1[0]
                results_seq = [result1]
                for res in result2:
                    while type(res) == list:
                        res = res[0]
                    results_seq.append(res)
                results_packed.append(results_seq)

        else:
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

            #print("These should be equal", len(work_list_all), len(results_all), num_preds*len(results_packed))
        
        scores = []
        if not contacts:
            contacts=(None, None, None)
        if not distance_cutoffs:
            distance_cutoffs=(4, 4, 8)
            
        for results, dsobj in zip(results_packed, scoring_pool):
            scores.append(score_func(results, dsobj, contacts=contacts, distance_cutoffs=distance_cutoffs))
        
        #adding sequences and scores into the dictionary
        for score, results, dsobj in zip(scores, results_packed, scoring_pool):
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

        curr_iter+=1

    with open(output_dir+"seqs_and_scores.log", "w") as logf:
        for tuplist, i in zip(seqs_per_iteration, range(len(seqs_per_iteration))):
            logf.write(">iteration"+str(i+1)+"\n")
            for val in tuplist:
                logf.write(str(val[0])+"\t"+str(val[1])+"\n")
        
    for key_seq, j in zip(newpool_seqs, range(len(newpool_seqs))):
        pdbs = scored_seqs[key_seq]["data"][0]["pdb"]
        result = scored_seqs[key_seq]["data"][0]["result"]
        results = [scored_seqs[key_seq]["data"][x]["result"] for x in range(len(scored_seqs[key_seq]["data"]))]

        for res, k in zip(results, range(len(results))):
            if write_compressed_data:
                compressed_pickle(os.path.join(output_dir, "seq_" + str(j) + "_result_" + str(k)), res)
        
        if conf_plot:
            for r, c in zip(result, af2_preds):
                while type(r) == list:
                    r = r[0]

                # Plot confidence.
                Ls = get_chain_lengths_rf2(r)
                print("Ls", Ls)

                # Plot pAEs.
                ptm, iptm = None, None
                if 'ptm' in r:
                    ptm = r['ptm']
                if 'iptm' in r:
                    iptm = r['iptm']

                pae_fig = plot_pae(r['pae_output'][0], Ls, ptm=ptm, iptm=iptm)
                pae_fig.savefig(output_dir + "seq_" + str(j) + "_final_chain_"+str(c)+"_pae.png")

                # Plot pLDDTs.
                plddt_fig = plot_plddt(r['plddt'], Ls)
                plddt_fig.savefig(output_dir + "seq_" + str(j) + "_final_chain_"+str(c)+"_plddt.png")
        
        for pdb, chains in zip(pdbs, af2_preds):
            with open(output_dir + "seq_" + str(j) + "_final_model_1_chain"+str(chains)+".pdb", "w") as pdbf:
                        pdbf.write(str(pdb))

    print("Number of RosettaFold2 predictions: ", num_af2)
    dist.spin_down()
    dist2.spin_down()
    
    plotting = plot_scores_general_dev(plot, seqs_per_iteration, output_dir)
    print("plots created at", str(plotting))

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
            flagsfile2=None
            resfile=None
            starting_seqs=None
            for filename in onlyfiles:
                if "flag" in filename and "rf2" in filename:
                    if "_2" in filename:
                        flagsfile2 = filename
                    else:
                        flagsfile = filename
                    if flagsfile2 and flagsfile:
                        break

            for filename in onlyfiles:
                if "residue" in filename and "spec" in filename:
                    resfile = filename
                    break

            if flagsfile is None:
                raise ValueError("Flags file for RosettaFold2 runs not provided.")
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
        if args.score_func_2:
            scorefunc2 = getattr(mod, args.score_func_2)
        else:
            scorefunc2 = None
    except:
        raise ValueError("Invalid score function")

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
        
    if args.mpnn_chains:
        mpnn_chains = args.mpnn_chains.strip().split(",")
    else:
        mpnn_chains = None

    print("Writing compressed data", (not args.dont_write_compressed_data))
    
    if args.af2_preds:
        af2_preds = args.af2_preds.strip().split(",")
                                   
    run_genetic_alg_multistate(input_dir, input_dir + flagsfile, scorefunc, starting_seqs, af2_flags_file_2=flagsfile2, poolsizes=pool_sizes, 
        num_iter = args.num_iter, n_workers=args.num_gpus, mut_percents=mut_percents, contacts=contacts, distance_cutoffs=distance_cutoffs,
        score_func_2=scorefunc2, score_func_2_iter=args.score_func_2_iteration, mpnn_chains=mpnn_chains, 
        mpnn_temps=mpnn_temps, mpnn_version=args.mpnn_version, skip_mpnn=mpnn_skips, mpnn_iters=mpnn_iters, 
        repeat_af2=not args.no_repeat_af2, af2_preds = af2_preds, crossover_percent=args.crossover_percent, vary_length=args.vary_length, 
        write_pdbs=args.write_pdbs, plot=plot_style, conf_plot=args.plot_confidences, write_compressed_data=not args.dont_write_compressed_data)