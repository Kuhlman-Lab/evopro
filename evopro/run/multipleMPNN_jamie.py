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

def run_mpnn_on_diff_backbone(pdb_backbone, jsonfile, output_dir, mpnn_temp="0.1", mpnn_version="s_48_020"):
    with open(jsonfile, 'r') as f:
        jsondata = f.read()

    new_seq = run_protein_mpnn(pdb_backbone, jsondata, mpnn_temp, mpnn_version=mpnn_version, bidir=False)
    print(f"new_seq: {new_seq}")
    seq = new_seq[0][-1][-1].strip().split("/")
    print(f"seq: {seq}")
        
    af2seq = "," + seq[0] + seq[1]
    print(f"af2seq: {af2seq}")
        
    pippackseq = seq[0] + "/" + seq[1]
    print(f"pippackseq: {pippackseq}")

    os.mkdir("af2_inputs")
    with open(os.path.join(output_dir, "af2_inputs", "sequences.csv"), 'w') as f:
        f.write(af2seq)

    with open(os.path.join(output_dir,"pippacksequences.fasta"),"w") as f:
        f.write(pippackseq)


def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')
    parser.add_argument('--target_name',default='',help='name of the diffusion target')
    parser.add_argument('--hotspot_res',default='',help='one-letter code of hotspots (eg YFL)')
    parser.add_argument('--other',default='',help='other identifying info (date, #diffruns,etc)')
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

    if args.pdb_dir:
        input_dir = args.input_dir
        if input_dir == "./":
            input_dir = os.getcwd()
        pdb_dir = os.path.join(input_dir, args.pdb_dir)
        output_dir = os.path.join(input_dir, "outputs/")

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

    if args.num_seqs_mpnn:
        num_seqs_mpnn = args.num_seqs_mpnn

        pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]

        for pdb in pdbs:
            print("Running MPNN on", pdb)
            name = pdb.split("/")[-1].split(".")[0]
            mpnn_dir = os.path.join(output_dir, name)
            if not os.path.isdir(mpnn_dir):
                for i in range(0,num_seqs_mpnn):
                    mpnn_run_dir = os.path.join(f"{mpnn_dir}",f"MPNN{i}")
                    if not os.path.isdir(mpnn_run_dir):
                        os.makedirs(mpnn_run_dir)
                        pdb_backbone = os.path.join(mpnn_run_dir, "DiffBB.pdb")
                        shutil.copy(pdb, pdb_backbone)
                        shutil.copy(os.path.join(input_dir, "json.flags"), mpnn_run_dir)
            
                        os.chdir(mpnn_run_dir)
                        print("generating json")
                        print(os.getcwd())
                        subprocess.run(["python", "/proj/kuhl_lab/evopro/evopro/run/generate_json_jamie.py", "@json.flags"])
                        jsonfile = os.path.join(mpnn_run_dir, "residue_specs.json")
                        run_mpnn_on_diff_backbone(mpnn_run_dir, jsonfile, mpnn_run_dir, mpnn_temp=args.mpnn_temp, mpnn_version=args.mpnn_version)

                        os.chdir("../")

                os.chdir(input_dir)

        print("done churning")
