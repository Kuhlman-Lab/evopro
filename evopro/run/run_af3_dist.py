import sys
import os
import json
import random
import importlib

sys.path.append("/proj/kuhl_lab/evopro2/")
from utils.inputs import FileArgumentParser
from utils.parsing import parse_rna,parse_results_af3
from utils.distributor import Distributor
from utils.utils import compressed_pickle
from utils.parsing_utils import constituents_of_modified_fasta
from objects.chemical import ptms_dict

    
def getFlagParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script."""

    parser = FileArgumentParser(description='Parser that can take flag options for script.',
                                fromfile_prefix_chars='@')

    parser.add_argument('--sequence_csv',
                        default=None,
                        type=str,
                        help='Path and filename to CSV file that contains sequences for AF3.')
    
    parser.add_argument('--max_tokens_length',
                        default=1000,
                        type=int,
                        help='Max total number of tokens in all chains for AF3 compilation.')
    
    parser.add_argument('--custom_score',
                        default="/proj/kuhl_lab/evopro2/scoring/mpnn_af3_scoring.py score_seq_diff_backbone_plddt_only_af3",
                        type=str,
                        help='Name of python file and function to calculate scores on each prediction.')
    
    parser.add_argument('--custom_score_args',
                        default=None,
                        type=str,
                        help='String of arguments to be passed to the score function.')

    parser.add_argument('--n_workers',
                        default=1,
                        type=int,
                        help='Number of GPUs available for distributed AlphaFold3. Default is 1.')
    
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
    
    parser.add_argument('--af3_path',
                        default="/proj/kuhl_lab/alphafold3/run/",
                        type=str,
                        help='Path to AF3 installation.')
    
    return parser

def parse_protein_chains(csv_file_path):
    """
    Parse a CSV file containing protein chain information.
    
    Parameters:
    csv_file_path (str): Path to the CSV file.
    
    Returns:
    list: A list of dictionaries, each containing:
        - 'name': Protein name
        - 'chains': Dictionary of chain IDs mapped to their amino acid sequences
    """
    protein_list = {}
    
    with open(csv_file_path, 'r') as file:
        f = file.readlines()
        for line in f:
            if not line:  # Skip empty rows
                continue
            
            l = line.strip("\n").split(",")
                
            protein_name = l[0]
            chains = []
            
            # Parse chains starting from index 1 (after the name)
            for i in range(1, len(l)):
                chains.append(l[i])
            
            protein_list[protein_name] = chains
    
    return protein_list

def change_dir(path):
    os.chdir(path)

def get_chain_length(chain, residues):
    count = 0
    for residue in residues:
        if residue[0] == chain:
            count += 1

    return count

def get_chain_type(seq):
    dna_chars = ["a", "t", "g", "c"]
    rna_chars = ["b", "u", "h", "d"]
    protein_chars = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    if seq[0] in dna_chars and seq[1] in dna_chars and seq[2] in dna_chars and seq[-1] in dna_chars:
        return "dna"
    elif seq[0] in rna_chars and seq[1] in rna_chars and seq[2] in rna_chars and seq[-1] in rna_chars:
        return "rna"
    elif seq[0] in protein_chars and seq[1] in protein_chars and seq[2] in protein_chars and seq[-1] in protein_chars:
        return "protein"
    else:
        return "ligand"
    

def get_ptms(constituent_list):
    mods = []
    sequence = ""
    for aa, i in zip(constituent_list, range(len(constituent_list))):
        if len(aa) == 1:
            sequence += aa
        else:
            mods.append((aa, i+1))
            sequence = sequence + ptms_dict[aa]
    return mods, sequence
    

def get_formatted_input_af3(name, seqs, custom_msa_dict={}, custom_template_dict={}):
    data = []
    
    output_dict = {}
    output_dict['name'] = name
    seed = random.randint(0, 1000000000)
    output_dict['modelSeeds'] = [seed]

    num_chains = len(seqs)
    for sequence in seqs:
        s = sequence.strip(", \t").split(":")
        chain = s[0]
        t = s[1]
        seq = s[2]
        consts = constituents_of_modified_fasta(seq, t)
        
        modifications = []
        #parse sequence for mods
        mods, seq = get_ptms(consts)
        for mod in mods:
           modifications.append({"ptmType": mod[0], "ptmPosition": int(mod[1])})
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
    
    # with open("data.json", "w") as f:
    #     f.write(json.dumps(output_dict, indent=4))
    return data

def score_and_write_results(results, out_dir, score_func, name, args):

    scores = []
    for result, i in zip(results, range(len(results))):
        parsed_result, pdbs = parse_results_af3([result], None)
        pdb = pdbs[0]
        with open(os.path.join(out_dir, name+".pdb"), 'w') as f:
            f.write(pdb)
        compressed_pickle(os.path.join(out_dir, name), parsed_result)
        while type(parsed_result) == list:
            parsed_result = parsed_result[0]
        score = score_func(parsed_result, reslist=args.custom_score_args)
        scores.append((name, score))
    
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

    sequence_list = parse_protein_chains(args.sequence_csv)
        
    output_dir = os.path.join(curr_dir, "outputs/")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    lengths = [int(args.max_tokens_length)]

    sys.path.append(args.af3_path)
    from evopro_utils import init_af3

    dist = Distributor(args.n_workers, init_af3, args.af3_flags, lengths)
    all_scores = []
    
    work_list = []
    
    for name in sequence_list:
        work_list.extend(get_formatted_input_af3(name, sequence_list[name], custom_msa_dict=custom_msa_dict, custom_template_dict=custom_template_dict))

    print("work list", work_list)
    results = dist.churn(work_list)
    
    for name, result in zip(sequence_list, results):
        scores = score_and_write_results(result, output_dir, score_func, name, args)
        all_scores.extend(scores)

    dist.spin_down()
    all_scores.sort(key=lambda x: x[1][0])
    with open(os.path.join(output_dir, "all_scores.csv"), 'w') as f:
        f.write("sequence,score\n")
        for score in all_scores:
            f.write(score[0] + "," + str(score[1]) + "\n")
    
if __name__ == "__main__":
    main()