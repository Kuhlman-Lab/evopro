import sys, os, time
import argparse
import omegaconf

sys.path.append("/proj/kuhl_lab/evopro2/")
from utils.inputs import parse_arguments
from utils.parsing import parse_input_sequences
from objects.sequence import DesignSeq
from objects.pool import Pool
from utils.distributor_utils import initialize_distributor, get_lengths
from scoring.scoring import score_overall
from utils.utils import compressed_pickle
from utils.plotting import plot_scores, plot_pae, plot_plddt

def main(conf) -> None:
    conf = parse_arguments(conf)
    run_genetic_alg(conf)

def run_genetic_alg(conf):
    all_start = time.time()
    #creating a design sequence object from the residue specifications file
    dsobj1 = DesignSeq(conf)
    
    #creating a pool by point mutating the design sequence object
    pool = Pool(conf, dsobj1)

    num_struct_preds=0
    
    lengths = []
    num_preds = len(conf.structure_prediction.structure_pred_chains)

    lengths = get_lengths(dsobj1, conf)

    print("Compiling structure prediction model for lengths:", lengths)

    print("Initializing distributor")
    dist = initialize_distributor(conf, lengths, pre_func=False)

    #creating an outputs directory
    output_dir = conf.flags.run_dir + "outputs/"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    #optional flags
    if conf.flags.write_all_pdbs:
        pdb_folder = output_dir + "pdbs_per_iter/"
        if not os.path.isdir(pdb_folder):
            os.makedirs(pdb_folder)
    
    log_file = output_dir+"scores.csv"
    with open(log_file, "w") as logf:
        logf.write("sequence chains,overall score,individual score_terms\n")
            
    if not conf.sequence_prediction.sequence_pred_chains:
        conf.sequence_prediction.sequence_pred_chains = [conf.structure_prediction.structure_pred_chains[0]]
    
    curr_iter = 1
    scores_per_iter = []
    #start genetic algorithm iteration
    while curr_iter <= conf.flags.num_iter:
        
        print("Starting iteration", curr_iter)
        start = time.time()
        
        pool.size = conf.flags.pool_sizes[curr_iter-1]
    
        #removing the worst scoring sequences from the pool
        if curr_iter>1:
            pool.purge()

            if curr_iter in conf.sequence_prediction.mpnn_iters:
                pool.refill(conf, seq_pred = "mpnn")
            else:
                pool.refill(conf, seq_pred = "random")
        
        print("\nCurrent pool at iter", str(curr_iter) + ":", pool)

        #extracting sequences from pool for structure prediction
        work_list_all = []
        for d in pool.pool:
            if not d.scored:
                formatted_inputs = parse_input_sequences(d, conf)
                work_list_all.extend(formatted_inputs)

        print("work list", work_list_all)
        num_struct_preds += len(work_list_all)

        #running the structure prediction model
        results = dist.churn(work_list_all)
        print("done churning through structure prediction model")

        results_all = []
        for result in results:
            while type(result) == list:
                result = result[0]
            results_all.append(result)
        
        #separating results by design sequence object
        results_packed = [results_all[i:i+num_preds] for i in range(0, len(results_all), num_preds)]
        work_list_packed = [work_list_all[i:i+num_preds] for i in range(0, len(work_list_all), num_preds)]
        
        scores = []
        for results, sequences, dsobj in zip(results_packed, work_list_packed, pool.pool):
            if not dsobj.scored:
                score, pdbs, results_parsed = score_overall(results, sequences, dsobj, conf)
                dsobj.set_score(conf, score, pdbs, results_parsed)
                scores.append(score)

        pool.sort_by_score()
        pool.log_scores(log_file, curr_iter)

        if conf.flags.write_all_pdbs:
            for i, d in enumerate(pool.pool):
                curr_folder = pdb_folder + "iter_" + str(curr_iter) + "/"
                os.makedirs(curr_folder, exist_ok=True)
                for pdb, chains in zip(pdbs, conf.structure_prediction.structure_pred_chains):
                    with open(curr_folder + "seq_" + str(i) + "_chain"+str(chains)+".pdb", "w") as pdbf:
                                    pdbf.write(str(pdb))

        scores_per_iter.append(scores)
        print("Iteration " + str(curr_iter) + " time:", time.time()-start, "seconds.")
        curr_iter+=1

    print("Finished iterative optimization.")
    print("Number of structure predictions made:", num_struct_preds)
    
    dist.spin_down()
    
    #writing outputs
    pool.purge()
    for i, d in enumerate(pool.pool):
        scores = d.get_scores()
        for j, elem in enumerate(scores):
            if not conf.flags.dont_write_compressed_data:
                compressed_pickle(output_dir + "seq_" + str(i) + "_result_" + str(j), elem[2])
            for pdb, result, chain in zip(elem[1], elem[2], conf.structure_prediction.structure_pred_chains):
                with open(output_dir + "seq_" + str(i) + "_pred_" + str(j) +  "_chain"+str(chain)+".pdb", "w") as pdbf:
                    pdbf.write(str(pdb))
                if not conf.flags.no_plots_output:
                    
                    # Plot pAEs.
                    ptm, iptm = None, None
                    if 'ptm' in result:
                        ptm = result['ptm']
                    if 'iptm' in result:
                        iptm = result['iptm']

                    Ls = get_lengths(dsobj1, conf)
                    # print(Ls)
                    if 'pae' in result:
                        pae_fig = plot_pae(result['pae'], Ls=Ls, ptm=ptm, iptm=iptm)
                        pae_fig.savefig(output_dir + "seq_" + str(i) + "_pred_" + str(j) +  "_chain" + str(chain) + "_pae.png")

                    # Plot pLDDTs. 
                    # TODO: residue plddt for af3
                    if 'plddt' in result:
                        plddt_fig = plot_plddt(result['plddt'], Ls=Ls)
                        plddt_fig.savefig(output_dir + "seq_" + str(i) + "_pred_" + str(j) +  "_chain" + str(chain) + "_plddt.png")
            
    print("Final sequences and structures written to", output_dir)
    
    # TODO: make plots of scores pretty
    if not conf.flags.no_plots_output:
        for score_type in conf.flags.plot_scores:
            plotting = plot_scores(scores_per_iter, output_dir + "scores_"+score_type+".png", stat_type=score_type)
        #print("plots created at", str(plotting))
    
    print("Finished EvoPro2 run. Took", time.time()-all_start, "seconds.")

if __name__ == "__main__":
    # main()
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    argparser.add_argument(
        "--config_file",
        type=str,
        default="./evopro_basic.yaml",
        help="path to yaml config file to load options",
    )
    
    args = argparser.parse_args()
    
    # Load config
    conf = omegaconf.OmegaConf.load(args.config_file)
    
    main(conf)