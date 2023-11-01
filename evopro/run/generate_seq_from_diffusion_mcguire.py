import importlib
import subprocess
import sys, os
import shutil
from functools import partial

sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle

sys.path.append('/proj/kuhl_lab/alphafold/run')
from run_af2 import af2_init

def run_protein_mpnn(pdb_dir, jsonstring, mpnn_temp, mpnn_version="v_48_020", bidir=False):
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
        #print(new_seq)
        seq = new_seq[0][-1][-1].strip().split("/")
        newseq_sequence = ",".join(seq)
        newseq_sequence = "," + newseq_sequence
        if newseq_sequence not in mpnn_seqs:
            
            mpnn_seqs_af2.append([seq])
            
            mpnn_seqs.append(newseq_sequence)

    #print(mpnn_seqs)
    os.mkdir("af2_inputs")
    with open(os.path.join(output_dir, "sequences.csv"), 'w') as f:
        f.write("\n".join(mpnn_seqs))

    shutil.move('sequences.csv','af2_inputs/')

    return mpnn_seqs_af2


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
    
    #BELOW FEATURE STILL HAS BUGS
    parser.add_argument('--pdb_dir',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, name of PDB directory to iterate over and generate sequences using MPNN and check RMSD against.')
    parser.add_argument('--max_chains_length',
                        default=None,
                        type=str,
                        help='If more than one PDB input backbone, max length of each chain separated by commas.')
    
    
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

def change_dir(path):
    os.chdir(path)

if __name__ == "__main__":
    
    #print("length of command line args", len(sys.argv), sys.argv)
    parser = getFlagParser()
    #print(parser)
    args = parser.parse_args(sys.argv[1:])
    print(args)
    
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
        else:
            shutil.rmtree(output_dir)
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
        
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]
        
        
        for pdb in pdbs:
            print("Running MPNN on", pdb)
            name = pdb.split("/")[-1].split(".")[0]
            mpnn_dir = os.path.join(output_dir, name)
            if not os.path.isdir(mpnn_dir):
                os.makedirs(mpnn_dir)
            
                pdb_backbone = os.path.join(mpnn_dir, "DiffBB.pdb")
                shutil.copy(pdb, pdb_backbone)
                try:
                    shutil.copy(os.path.join(input_dir, "json.flags"), mpnn_dir)
                except:
                    raise ValueError("No json.flags file found in input directory. Please provide a json.flags file in the input directory that is generalizable to all input PDBs.")

            #input_dir_return = os.getcwd()
            
                os.chdir(mpnn_dir)
                print("generating json")
                subprocess.run(["python", "/proj/kuhl_lab/evopro/evopro/run/generate_json_dev.py", "@json.flags"])
                jsonfile = os.path.join(mpnn_dir, "residue_specs.json")
                mpnn_seqs_list = run_mpnn_on_diff_backbone(mpnn_dir, jsonfile, mpnn_dir, mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version, num_seqs=args.num_seqs_mpnn)
            
                print(f"mpnn_seqs_list: {mpnn_seqs_list}")

                mod_mpnn_seqs_list = []
                for element in mpnn_seqs_list:
                    element.append(partial(change_dir, mpnn_dir))
            
                with open("pippacksequences.fasta", 'w') as f:
                    f.write(f'{mpnn_seqs_list[0][0][0]}/{mpnn_seqs_list[0][0][1]}\n')

                os.chdir(input_dir)
        
        print("done churning")
            
