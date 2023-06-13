import importlib
import subprocess
import sys, os
import shutil

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def run_protein_mpnn(pdb_dir, jsonstring, mpnn_temp, mpnn_version="s_48_020", bidir=False):
    sys.path.append('/proj/kuhl_lab/proteinmpnn/run/')
    from run_protein_mpnn import run_protein_mpnn_func
    results = run_protein_mpnn_func(pdb_dir, jsonstring, sampling_temp=mpnn_temp, model_name=mpnn_version, bidir=bidir)

    return results

#run and collect results from protein mpnn
def run_mpnn_on_diff_backbone(pdb_backbone, jsonfile, output_dir, mpnn_temp="0.1", mpnn_version="s_48_020", num_seqs=100):
    mpnn_seqs = []
    mpnn_seqs_af2 = []
    with open(jsonfile, 'r') as f:
        jsondata = f.read()
        
    while len(mpnn_seqs)<num_seqs:
        new_seq = run_protein_mpnn(pdb_backbone, jsondata, mpnn_temp, mpnn_version=mpnn_version, bidir=False)
        print(new_seq)
        seq = new_seq[0][-1][-1].strip().split("/")
        newseq_sequence = ",".join(seq)
        newseq_sequence = "," + newseq_sequence
        if newseq_sequence not in mpnn_seqs:
            
            mpnn_seqs_af2.append([seq])
            
            mpnn_seqs.append(newseq_sequence)

    print(mpnn_seqs)
    with open(os.path.join(output_dir, "sequences.csv"), 'w') as f:
        f.write("\n".join(mpnn_seqs))
        
    return mpnn_seqs_af2
    
def run_af2_on_mpnn_seqs(diff_backbone, mpnn_seqs_list, score_func, af2_flags_file, output_dir,  n_workers=1):
    # mpnn_seqs_list is a list of multi chain sequences

    num_af2=0

    lengths = []
    lengths.append([len(x) for x in mpnn_seqs_list[0][0]])
    #lengths.append([[x + vary_length for x in startingseqs[0].get_lengths([chain])][0] for chain in c])

    print("Compiling AF2 models for lengths:", lengths)

    print("Initializing distributor")
    dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)
        
    print("work list", mpnn_seqs_list)
    num_af2 += len(mpnn_seqs_list)

    results = dist.churn(mpnn_seqs_list)
    
    print("done churning")
    dist.spin_down()
    #print(results)
    
    seqs_and_scores = {}
    data = {}
    scores = []
    for seq, result in zip(mpnn_seqs_list, results):
        s = ",".join(seq[0])
        score = score_func(result[0], diff_backbone)
        scores.append(score)
        
        seqs_and_scores[s] = (score[0])
        data[s] = (score[1], score[-1], score[-2])
        
        #print(s, score[0])
    
    sorted_seqs_and_scores = sorted(seqs_and_scores.items(), key=lambda x:x[1])

    print(sorted_seqs_and_scores)
    for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
        pdb = data[seq[0]][-1]
        r = data[seq[0]][-2]
        with open(os.path.join(output_dir,"seq_" + str(i) + "_model_1.pdb"), "w") as pdbf:
            pdbf.write(str(pdb))
        compressed_pickle(os.path.join(output_dir, "seq_" + str(i) + "_result"), r)
        print(seq[0], data[seq[0]][0])
    
    with open(os.path.join(output_dir, "seqs_and_scores.csv"), "w") as opf:
        for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
            opf.write(str(seq[0]) + "\t" + str(seqs_and_scores[seq[0]]) + "\t" + str(data[seq[0]][0]) + "\n")

    print("Number of AlphaFold2 predictions: ", num_af2)
    
def run_af2_on_mpnn_seqs_withmonomers(diff_backbone, mpnn_seqs_list, score_func, af2_flags_file, output_dir,  n_workers=1):
    # mpnn_seqs_list is a list of multi chain sequences

    num_af2=0

    lengths = []
    #print(mpnn_seqs_list[0][0])
    lengths.append([len(x) for x in mpnn_seqs_list[0][0]])
    lengths = lengths + [[len(x)] for x in mpnn_seqs_list[0][0]]
    num_preds = len(lengths)

    print("Number of AF2 predictions per sequence:", num_preds)
    print("Compiling AF2 models for lengths:", lengths)
        
    work_list = []
    for seq in mpnn_seqs_list:
        #print(seq[0])
        work_list.append(seq)
        work_list = work_list + [[[x]] for x in seq[0]]
        
    print("work list", work_list)
    num_af2 += len(work_list)

    print("Initializing distributor")
    dist = Distributor(n_workers, af2_init, af2_flags_file, lengths)
    results = dist.churn(work_list)
    
    print("done churning")
    dist.spin_down()

    results_parsed = []
    for result in results:
        while type(result) == list:
            result = result[0]
        results_parsed.append(result)
    
    seqs_packed = [work_list[i:i+num_preds] for i in range(0, len(work_list), num_preds)]
    results_packed = [results_parsed[i:i+num_preds] for i in range(0, len(results_parsed), num_preds)]
    
    seqs_and_scores = {}
    data = {}
    scores = []
    for seqs, results in zip(seqs_packed, results_packed):
        s = ",".join(seqs[0][0])
        score = score_func(results, diff_backbone)
        scores.append(score)
        #print(score)
        seqs_and_scores[s] = (score[0])
        data[s] = (score[1], score[-1], score[-2])
            
    sorted_seqs_and_scores = sorted(seqs_and_scores.items(), key=lambda x:x[1])

    names = ["complex", "monomerA", "monomerB", "monomerC", "monomerD"]

    for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
        pdbs = data[seq[0]][-1]
        rs = data[seq[0]][-2]
        for pdb, name in zip(pdbs, names):
            with open(os.path.join(output_dir, "seq_" + str(i) + "_" + name + ".pdb"), "w") as pdbf:
                pdbf.write(str(pdb))
        for r in rs:
            compressed_pickle(os.path.join(output_dir, "seq_" + str(i) + "_" + name + "_result"), r)
    
    with open(os.path.join(output_dir, "seqs_and_scores.csv"), "w") as opf:
        for seq,i in zip(sorted_seqs_and_scores, range(len(sorted_seqs_and_scores))):
            opf.write(str(seq[0]) + "\t" + str(seqs_and_scores[seq[0]]) + "\t" + str(data[seq[0]][0]) + "\n")

    print("Number of AlphaFold2 predictions: ", num_af2)


def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--input_dir',
                        default='./',
                        type=str,
                        help='Path to and name of directory with input files (pdb_backbone_name.pdb, af2.flags, residue_specs.json) to extract chains and sequences.')
    
    parser.add_argument('--pdb_backbone',
                        default=None,
                        type=str,
                        help='Name of PDB file (if only one) to generate sequences using MPNN and check RMSD against.')
    
    parser.add_argument('--pdb_dir',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, name of PDB directory to iterate over and generate sequences using MPNN and check RMSD against.')
    
    parser.add_argument('--custom_score',
                        default="/proj/kuhl_lab/evopro/evopro/score_funcs/diff_score_binder.py score_seq_diff",
                        type=str,
                        help='Name of PDB file to generate sequences using MPNN and check RMSD against.')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for AlphaFold2. Default is 1.')
    
    parser.add_argument('--num_seqs_mpnn',
                        default=5,
                        type=int,
                        help='Number of sequences to generate using MPNN for provided PDB backbone. Default is 5.')
    
    parser.add_argument('--mpnn_temp',
                        default="0.1",
                        type=str,
                        help='MPNN temperature to use for sequence generation. Default is 0.1.')
    
    parser.add_argument('--mpnn_version',
                        default="s_48_020",
                        type=str,
                        help='MPNN version to use for sequence generation. Default is s_48_020.')
    
    parser.add_argument('--af2_preds_monomers',
                         action='store_true',
                         help='Default is False.')

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
    
    if args.pdb_backbone:
        input_dir = args.input_dir
        output_dir = os.path.join(input_dir, "outputs/")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        pdb_backbone = os.path.join(input_dir, args.pdb_backbone)
        jsonfile = os.path.join(input_dir, "residue_specs.json")
        af2_flags_file = os.path.join(input_dir, "af2.flags")
        mpnn_seqs_list = run_mpnn_on_diff_backbone(input_dir, jsonfile, output_dir, mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, num_seqs=args.num_seqs_mpnn)
        if args.af2_preds_monomers:
            run_af2_on_mpnn_seqs_withmonomers(pdb_backbone, mpnn_seqs_list, score_func, af2_flags_file, output_dir, n_workers=args.n_workers)
        else:
            run_af2_on_mpnn_seqs(pdb_backbone, mpnn_seqs_list, score_func, af2_flags_file, output_dir, n_workers=args.n_workers)
    
    elif args.pdb_dir:
        input_dir = args.input_dir
        if input_dir == "./":
            input_dir = os.getcwd()
        pdb_dir = os.path.join(input_dir, args.pdb_dir)
        output_dir = os.path.join(input_dir, "outputs/")
        af2_flags_file = os.path.join(input_dir, "af2.flags")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]
        for pdb in pdbs:
            
            name = pdb.split("/")[-1].split(".")[0]
            mpnn_dir = os.path.join(output_dir, name)
            if not os.path.isdir(mpnn_dir):
                os.makedirs(mpnn_dir)
            
            pdb_backbone = os.path.join(mpnn_dir, "design.pdb")
            shutil.copy(pdb, pdb_backbone)
            try:
                shutil.copy(os.path.join(input_dir, "json.flags"), mpnn_dir)
            except:
                raise ValueError("No json.flags file found in input directory. Please provide a json.flags file in the input directory that is generalizable to all input PDBs.")

            #input_dir_return = os.getcwd()
            
            os.chdir(mpnn_dir)
            subprocess.run(["python", "/proj/kuhl_lab/evopro/evopro/run/generate_json.py", "@json.flags"])
            jsonfile = os.path.join(mpnn_dir, "residue_specs.json")
            mpnn_seqs_list = run_mpnn_on_diff_backbone(mpnn_dir, jsonfile, mpnn_dir, mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, num_seqs=args.num_seqs_mpnn)
            if args.af2_preds_monomers:
                run_af2_on_mpnn_seqs_withmonomers(pdb_backbone, mpnn_seqs_list, score_func, af2_flags_file, mpnn_dir, n_workers=args.n_workers)
            else:
                run_af2_on_mpnn_seqs(pdb_backbone, mpnn_seqs_list, score_func, af2_flags_file, mpnn_dir, n_workers=args.n_workers)
            
            print("Finished running AlphaFold2 on MPNN sequences in", os.getcwd())
            os.chdir(input_dir)
            