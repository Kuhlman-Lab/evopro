import sys, os
import importlib

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle
from evopro.score_funcs.score_funcs import score_plddt_confidence_overall

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def run_af2_dist(args, score_func=None):
    
    # Load sequences
    seqs_list = []
    with open(args.sequence_file, "r") as sf:
        for line in sf:
            seq = [x for x in line.strip("\n").split(",") if x]
            seqs_list.append([seq])
    
    print(seqs_list)
    if args.max_chains_length:
        lengths = [[int(x) for x in args.max_chains_length.split(",")]]
    else:
        lengths = [[len(x) for x in seqs_list[0].split(",")]]
        
    if not os.path.isfile(args.af2_flags_file):
        raise ValueError("Invalid path to af2 flags file. Please provide a valid af2 flags file.")
    
    if not score_func:
        score_func = score_plddt_confidence_overall
    
    print("Compiling AF2 models for lengths:", lengths)
    dist = Distributor(args.n_workers, af2_init, args.af2_flags_file, lengths)
    
    results = dist.churn(seqs_list)
    
    print("done churning")
    dist.spin_down()
    
    seqs_and_scores = {}
    data = {}
    scores = []
    for seq, result in zip(seqs_list, results):
        s = ",".join(seq[0])
        score = score_func(result[0])
        scores.append(score)
        
        seqs_and_scores[s] = (score[0])
        data[s] = (score[1], score[-1], score[-2])
        
        #print(s, score[0])
    
    sorted_seqs_and_scores = sorted(seqs_and_scores.items(), key=lambda x:x[1])
    
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    #print(sorted_seqs_and_scores)
    for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
        pdb = data[seq[0]][-1]
        r = data[seq[0]][-2]
        with open(os.path.join(args.output_dir,"seq_" + str(i) + "_model_1.pdb"), "w") as pdbf:
            pdbf.write(str(pdb))
        compressed_pickle(os.path.join(args.output_dir, "seq_" + str(i) + "_result"), r)
        #print(seq[0], data[seq[0]][0])
    
    with open(os.path.join(args.output_dir, "seqs_and_scores.csv"), "w") as opf:
        for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
            opf.write(str(seq[0]) + "\t" + str(seqs_and_scores[seq[0]]) + "\t" + str(data[seq[0]][0]) + "\n")
            
def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--af2_flags_file',
                        default='./af2.flags',
                        type=str,
                        help='Path to and name of af2.flags file.')
    
    parser.add_argument('--sequence_file',
                        default="sequences.csv",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--output_dir',
                        default="outputs",
                        type=str,
                        help='Path to and name of csv file containing sequences.')
    
    parser.add_argument('--max_chains_length',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, overall max length of each chain separated by commas.')
    
    parser.add_argument('--custom_score',
                        default=None,
                        type=str,
                        help='Score function to use to generate scores for each prediction (optional).')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for AlphaFold2. Default is 1.')

    return parser
            
if __name__ == "__main__":
    parser = getFlagParser()
    args = parser.parse_args(sys.argv[1:])
    
    score_func = None
    if args.custom_score:
        file = args.custom_score.split(" ")[0]
        function = args.custom_score.split(" ")[1]
        
        try:
            scorefile = file.rsplit("/", 1)
            scorepath = scorefile[0]
            scorefilename = scorefile[1].split(".")[0]

            sys.path.append(scorepath)
            mod = importlib.import_module(scorefilename)
            score_func = getattr(mod, function)
        except:
            raise ValueError("Invalid score function. Please provide a valid python file and function name within that file, separated by a space.")
        
    run_af2_dist(args, score_func=score_func)
