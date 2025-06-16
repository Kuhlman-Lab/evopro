import sys
import os
import json
import random
import importlib
import shutil
import subprocess
import omegaconf

sys.path.append("/proj/kuhl_lab/evopro_public/evopro/")
from evopro.utils.inputs import FileArgumentParser
from evopro.utils.parsing import parse_rna
from evopro.objects.sequence import DesignSeq
from evopro.utils.distributor import Distributor
from evopro.utils.utils import compressed_pickle
from evopro.utils.parsing import parse_results_af3
from evopro.utils.parsing_utils import get_coordinates_pdb_extended

af3_path = "/proj/kuhl_lab/alphafold3/run/"
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
    
    parser.add_argument('--max_tokens_length',
                        default=1000,
                        type=int,
                        help='Max total number of tokens in all chains for AF3 compilation.')
    
    parser.add_argument('--custom_score',
                        default="/proj/kuhl_lab/evopro2/scoring/standard_score_funcs.py score_seq_diff_backbone_plddt_only_af3",
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
    
    parser.add_argument('--af3_flags',
                        default="af3.flags",
                        type=str,
                        help='Name of af3 flags file.')
    
    
    parser.add_argument('--custom_msa_chains',
                        default=None,
                        type=str,
                        help='Name of chains for which custom MSA is being provided.')
    
    parser.add_argument('--custom_msa_paths',
                        default=None,
                        type=str,
                        help='Name of file paths containing custom MSA for each chain specified above.')

    parser.add_argument('--custom_template_chains',
                        default=None,
                        type=str,
                        help='Name of chains for which custom template is being provided.')
    
    parser.add_argument('--custom_template_paths',
                        default=None,
                        type=str,
                        help='Name of file paths containing custom template for each chain specified above.')
    
    parser.add_argument('--mpnn_only',
                        default=None,
                        action='store_true',
                        help='Only run MPNN on input PDBs. Do not run AF3.')
    
    return parser

def get_chain_length(chain, residues):
    count = 0
    for residue in residues:
        if residue[0] == chain:
            count += 1

    return count

def generate_dummy_seqfile(input_seqfile, jsonfile, pdbfile, outdir):
    chains, residues, resindices = get_coordinates_pdb_extended(pdbfile, fil=True)
    mutable_chain_lengths = {}

    with open(jsonfile, 'r') as f:
        jsonflags = f.readlines()

    mut_res = [x for x in jsonflags if "mut_res" in x][0].split(" ")[1].strip().split(",")
    for res in mut_res:
        chain = res[0]
        if chain not in mutable_chain_lengths:
            mutable_chain_lengths[chain] = get_chain_length(chain, residues)
    
    with open(input_seqfile, 'r') as f:
        lines = f.readlines()

    new_seqfile = ""

    for lin in lines:
        if lin:
            l = lin.strip().split(":")
            if l[0] in mutable_chain_lengths and l[1] == "protein":
                new_seqfile += l[0] + ":" + l[1] + ":"
                i = mutable_chain_lengths[l[0]]
                for j in range(i):
                    new_seqfile += "G"
                new_seqfile += "\n"
            else:
                new_seqfile += lin

    with open(os.path.join(outdir, "seqfile.txt"), 'w') as f:
        f.write(new_seqfile)

    # print(new_seqfile)
    return new_seqfile

def run_ligmpnn(pdb_dir, mpnn_conf, jsonfile, design_run=True):
    sys.path.append(mpnn_path)
    from run_mpnn import main as run_mpnn
    conf = omegaconf.OmegaConf.load(mpnn_conf)

    jsondata = None
    with open(jsonfile, 'r') as f:
        jsondata = json.load(f)
    
    pdb_paths = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]
    
    new_seqs = run_mpnn(conf, design_run=design_run, json_data=jsondata, pdb_paths=pdb_paths)
    # print("145: mpnn sequences")
    # print(new_seqs)
    return new_seqs

def get_formatted_input_af3(d, custom_msa_dict={}, custom_template_dict={}):
    data = []
    
    output_dict = {}
    output_dict['name'] = "design"
    seed = random.randint(0, 1000000000)
    output_dict['modelSeeds'] = [seed]

    c = list(d.chains)
    for chain in c:
        t = d.get_chain_type(chain)
        seq = d.get_chain_sequence(chain)
        mods = d.chains[chain].modifications
        modifications = []
        for mod in mods:
            if mod["chain"] == chain:
                modifications.append({"ptmType": mod["type"], "ptmPosition": int(mod["resid"])})
        if "sequences" not in output_dict:  
            output_dict["sequences"] = []
        
        if t == "ligand":
            if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["smiles"]:
                output_dict["sequences"][t]["id"].append(chain)
            else:
                seq_dict = {t:{'smiles': seq, "id": [chain]}}
                if chain in custom_msa_dict:
                    seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                if chain in custom_template_dict:
                    seq_dict[t]["templates"] = custom_template_dict[chain]

                output_dict["sequences"].append(seq_dict)
        
        elif t=="rna":
            if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                output_dict["sequences"][t]["id"].append(chain)
            else:
                seq_dict = {t:{'sequence': parse_rna(seq), "id": [chain]}}
                if chain in custom_msa_dict:
                    seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                if chain in custom_template_dict:
                    seq_dict[t]["templates"] = custom_template_dict[chain]

                output_dict["sequences"].append(seq_dict)
        elif t=="dna":
            if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                output_dict["sequences"][t]["id"].append(chain)
            else:
                seq_dict = {t:{'sequence': seq.upper(), "id": [chain]}}
                if chain in custom_msa_dict:
                    seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                if chain in custom_template_dict:
                    seq_dict[t]["templates"] = custom_template_dict[chain]

                output_dict["sequences"].append(seq_dict)
        
        else:
            if t in output_dict["sequences"] and seq == output_dict["sequences"][t]["sequence"]:
                output_dict["sequences"][t]["id"].append(chain)

            else:
                if len(modifications) > 0:
                    seq_dict = {t:{'sequence': seq, "modifications": modifications, "id": [chain]}}
                else:
                    seq_dict = {t:{'sequence': seq, "id": [chain]}}
                if chain in custom_msa_dict:
                    seq_dict[t]["unpairedMsa"] = custom_msa_dict[chain]
                if chain in custom_template_dict:
                    seq_dict[t]["templates"] = custom_template_dict[chain]

                output_dict["sequences"].append(seq_dict)
    
    data.append(json.dumps(output_dict))
    
    with open("data.json", "w") as f:
        f.write(json.dumps(output_dict, indent=4))
    return data

def score_and_write_results(results, out_dir, score_func, backbone_file):
    with open(backbone_file, 'r') as f:
        backbone = f.read()
    scores = []
    for result, i in zip(results, range(len(results))):
        parsed_result, pdbs = parse_results_af3(result, None)
        pdb = pdbs[0]
        with open(os.path.join(out_dir, "seq"+str(i)+".pdb"), 'w') as f:
            f.write(pdb)
        compressed_pickle(os.path.join(out_dir, "seq" + str(i)), parsed_result)
        while type(parsed_result) == list:
            parsed_result = parsed_result[0]
        score = score_func(parsed_result, backbone)
        scores.append(score)
    
    return scores

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
    
    #parse custom msas and templates
    custom_msa_dict = {}
    if args.custom_msa_chains and args.custom_msa_paths:
        chains = args.custom_msa_chains.split(",")
        paths = args.custom_msa_paths.split(",")
        for chain, path in zip(chains, paths):
            custom_msa_dict[chain] = path
    
    custom_template_dict = {}
    if args.custom_template_chains and args.custom_template_paths:
        chains = args.custom_template_chains.split(",")
        paths = args.custom_template_paths.split(",")
        for chain, path in zip(chains, paths):
            custom_template_dict[chain] = path
    
    curr_dir = os.getcwd()
    pdb_dir = os.path.join(curr_dir, args.pdb_dir)
    
    mpnn_yaml = os.path.join(curr_dir, args.mpnn_yaml)

    output_dir = os.path.join(curr_dir, "outputs/")
    if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

    templates_dir = None
    if os.path.isdir(os.path.join(curr_dir, "templates/")):
        templates_dir = os.path.join(curr_dir, "templates/")

    pdbs = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if os.path.isfile(os.path.join(pdb_dir, f)) and f.endswith(".pdb")]

    lengths = [int(args.max_tokens_length)]

    sys.path.append(af3_path)
    from evopro_utils import init_af3

    if not args.mpnn_only:
        dist = Distributor(args.n_workers, init_af3, args.af3_flags, lengths)
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
        generate_dummy_seqfile(os.path.join(curr_dir, "seqfile.txt"), os.path.join(curr_dir, "json.flags"), pdb, mpnn_dir)
        
        subprocess.run(["python", "/proj/kuhl_lab/evopro2/run/generate_json.py", "@json.flags"])
        jsonfile = os.path.join(mpnn_dir, "residue_specs.json")
        if args.mpnn_only:
            run_ligmpnn(mpnn_dir, mpnn_yaml, jsonfile, design_run=False)
        else:
            mpnn_seqs_list = run_ligmpnn(mpnn_dir, mpnn_yaml, jsonfile, design_run=True)
            
            d = DesignSeq(jsonfile=jsonfile)
            
            work_list = []
            for seq in mpnn_seqs_list:
                # print("334: inputs for DesignSeq")
                # print(d)
                # print(seq)
                d_str = str(d)
                # print(f"338: dstr {d_str}") 
                
                # here seq will contain a residue "X" for ptms
                # try replacing "X" with corresponding residue in d object
                replace_string = ""
                for i, val in enumerate(seq):
                    # replace X with corresponding char in other string
                    if val == "X":
                        replace_string += d_str[i]
                    else:
                        replace_string += val
                # print(f"348: replace: {replace_string}")    
                seq = replace_string

                    
                new_d = DesignSeq(template_desseq=d, sequence=seq)
                # reformat mpnn sequences for af3
                print("337: reformatting for af3")
                print(new_d)
                work_list.extend(get_formatted_input_af3(new_d, custom_msa_dict=custom_msa_dict, custom_template_dict=custom_template_dict))

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