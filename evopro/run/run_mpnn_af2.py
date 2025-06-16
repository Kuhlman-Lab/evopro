import sys
import os
import json
import importlib
import shutil
import subprocess
import omegaconf

sys.path.append("/proj/kuhl_lab/evopro_public/evopro/")
from evopro.utils.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle
from evopro.utils.parsing import parse_results_af2

af2_path = "/proj/kuhl_lab/alphafold/run/"
mpnn_path = "/proj/kuhl_lab/LigandMPNN/"

def change_dir(path):
    os.chdir(path)
    
def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--pdb_dir',
                        default=None,
                        type=str,
                        help='Path to PDB directory to iterate over.')
    
    parser.add_argument('--max_chains_length',
                        default=None,
                        type=str,
                        help='Max length of each chain across all PDB inputs, separated by commas.')
    
    parser.add_argument('--custom_score',
                        default="/proj/kuhl_lab/evopro2/scoring/standard_score_funcs.py score_seq_diff_backbone_af2",
                        type=str,
                        help='Name of python file and function to calculate scores on each prediction.')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for distributed AlphaFold3. Default is 1.')

    parser.add_argument('--mpnn_yaml',
                        default="mpnn_basic.yaml",
                        type=str,
                        help='Name of mpnn yaml file to set ligandmpnn options.')
    
    parser.add_argument('--af2_flags',
                        default="af2.flags",
                        type=str,
                        help='Name of af2 flags file.')

    parser.add_argument('--mpnn_only',
                        default=None,
                        action='store_true',
                        help='Only run MPNN on input PDBs. Do not run AF2.')
    
    #TODO: not working yet
    parser.add_argument('--prediction_chains',
                        default=None,
                        type=str,
                        help='If multistate, which chains to include per prediction separated by commas. eg. AB,A')
    
    return parser

def run_ligmpnn(pdb_dir, mpnn_conf, jsonfile, design_run=True):
    sys.path.append(mpnn_path)
    from run_mpnn import main as run_mpnn
    conf = omegaconf.OmegaConf.load(mpnn_conf)

    jsondata = None
    with open(jsonfile, 'r') as f:
        jsondata = json.load(f)
    
    pdb_paths = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]
    
    new_seqs = run_mpnn(conf, design_run=design_run, json_data=jsondata, pdb_paths=pdb_paths)

    return new_seqs

def score_and_write_results(results, out_dir, score_func, backbone_file):
    with open(backbone_file, 'r') as f:
        backbone = f.read()
    scores = []
    for result, i in zip(results, range(len(results))):
        parsed_result, pdbs = parse_results_af2(result, None)
        pdb = pdbs[0]
        with open(os.path.join(out_dir, "seq"+str(i)+".pdb"), 'w') as f:
            f.write(pdb)
        compressed_pickle(os.path.join(out_dir, "seq" + str(i)), parsed_result)
        while type(parsed_result) == list:
            parsed_result = parsed_result[0]
        score = score_func(parsed_result, backbone)
        scores.append(score)
    
    return scores

def change_dir(path):
    os.chdir(path)

def main():
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
    
    curr_dir = os.getcwd()
    pdb_dir = os.path.join(curr_dir, args.pdb_dir)
    
    mpnn_yaml = os.path.join(curr_dir, args.mpnn_yaml)

    output_dir = os.path.join(curr_dir, "outputs/")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    else:
        old_output_dir = os.path.join(curr_dir, "old_outputs/")
        shutil.move(output_dir, old_output_dir)

    templates_dir = None
    if os.path.isdir(os.path.join(curr_dir, "templates/")):
        templates_dir = os.path.join(curr_dir, "templates/")

    pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]
    
    lengths = []
    l = args.max_chains_length.split(",")
    l = [int(x) for x in l]
    lengths.append(l)

    sys.path.append(af2_path)
    from run_af2 import af2_init

    if not args.mpnn_only:
        print("Compiling AF2 models for lengths:", lengths)
        dist = Distributor(args.n_workers, af2_init, args.af2_flags, lengths)
        all_scores = []
    
    for pdb in pdbs:
        print("Running MPNN on", pdb)
        name = pdb.split("/")[-1].split(".")[0]
        mpnn_dir = os.path.join(output_dir, name)
        if not os.path.isdir(mpnn_dir):
            os.makedirs(mpnn_dir)

        pdb_backbone = os.path.join(mpnn_dir, "design.pdb")
        shutil.copy(pdb, pdb_backbone)
        
        try:
            shutil.copy(os.path.join(curr_dir, "json.flags"), mpnn_dir)
        except:
            raise ValueError("No json.flags file found in input directory. Please provide a json.flags file in the input directory that is generalizable to all input PDBs.")

        if templates_dir:
            shutil.copytree(templates_dir, os.path.join(mpnn_dir, "templates/"))
            
        os.chdir(mpnn_dir)
        
        subprocess.run(["python", "/proj/kuhl_lab/evopro2/run/generate_json.py", "@json.flags"])
        jsonfile = os.path.join(mpnn_dir, "residue_specs.json")
        if args.mpnn_only:
            run_ligmpnn(mpnn_dir, mpnn_yaml, jsonfile, design_run=False)
        else:
            mpnn_seqs_list = run_ligmpnn(mpnn_dir, mpnn_yaml, jsonfile, design_run=True)
                        
            work_list = []
            for seq in mpnn_seqs_list:
                #TODO
                if args.prediction_chains:
                    pass
                else:
                    af2_seq = []
                    for s in seq.split(","):
                        af2_seq.append(s)
                        
                    work_list.append([af2_seq])

            print("work list", work_list)
            results = dist.churn(work_list)
            scores = score_and_write_results(results, mpnn_dir, score_func, pdb_backbone)
            with open(os.path.join(mpnn_dir, "scores.csv"), 'w') as f:
                f.write("sequence chains,overall score,individual score_terms\n")
                for i in range(len(scores)):
                    f.write(mpnn_seqs_list[i] + "," + str(scores[i]) + "\n")
                    all_scores.append((mpnn_seqs_list[i], scores[i], i, mpnn_dir))
        
        os.chdir(curr_dir)
    
    if not args.mpnn_only:
        dist.spin_down()
        all_scores.sort(key=lambda x: x[1][0])
        with open(os.path.join(curr_dir, "all_scores.csv"), 'w') as f:
            f.write("sequence,overall score,sequenceID,path\n")
            for score in all_scores:
                f.write(score[0] + "," + str(score[1]) + "," + str(score[2]) + "," + score[3] + "\n")
        
if __name__ == "__main__":
    main()