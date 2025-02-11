import sys

#SET PATH TO YOUR EVOPRO INSTALLATION HERE
#sys.path.append("/proj/kuhl_lab/evopro/")
sys.path.append("/nas/longleaf/home/amritan/Desktop/kuhlmanlab/evopro_temp/evopro/")
from evopro.utils.inputs import FileArgumentParser
from evopro.utils.aa_utils import three_to_one, one_to_three
from evopro.utils.pdb_parser import get_coordinates_pdb_old
import json
import re
import sys

def getBiasParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run generate_json"""

    parser = FileArgumentParser(description='Script that can take a sequence and amino acid bias numbers'
                                ' and convert to json file format for input to EvoPro', 
                                fromfile_prefix_chars='@')

    parser.add_argument('--pdb',
                        default=None,
                        type=str,
                        help='Path to and name of PDB file to extract chains and sequences.')
    parser.add_argument('--sequence_file',
                        default=None,
                        type=str,
                        help='Path to and name of text file to extract chains and sequences. Only provide if there is no PDB file.')
    parser.add_argument('--bias_res',
                        default='',
                        type=str,
                        help='PDB chain and residue numbers to mutate, separated by commas')
    parser.add_argument('--custom_bias',
                        #        A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   X
                        default='0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0',
                        type=str,
                        help='Custom setting for residue biases. Default is even across all 20 amino acids.')
    parser.add_argument('--default_bias',
                        default='0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0',
                        type=str,
                        help='Default setting for residue biases. Default is even across all 20 amino acids.')
    parser.add_argument('--output',
                        default='bias.json',
                        type=str,
                        help='path and name of output json file')

    return parser

def parse_seqfile(filename):
    ids = {}
    chains = {}
    with open(filename, "r") as f:
        i=0
        for lin in f:
            l = lin.strip().split(":")
            chains[l[0]] = l[1]
    
    for chain in chains:
        reslist = list(chains[chain])
        for aa, i in zip(reslist, range(len(reslist))):
            resid = chain + str(i+1)
            extended_resid = chain + "_" + one_to_three(aa) + "_" + str(i+1)
            ids[resid] = (extended_resid, chain, i+1)

    return ids, chains

def parse_pdbfile(filename):
    chains, residues, _ = get_coordinates_pdb_old(filename, fil=True)

    pdbids = {}
    chain_seqs = {}

    for chain in chains:
        chain_seqs[chain]=[]

    res_index_chain = 1
    chain_test = chains[0]
    for residue in residues:
        num_id = int(residue.split("_")[-1])
        chain = residue.split("_")[0]
        if chain != chain_test:
            res_index_chain = 1
            chain_test = chain
        pdbid = chain+str(num_id)
        pdbids[pdbid] = (residue, chain, res_index_chain)

        aa = three_to_one(residue.split("_")[1])
        chain_seqs[chain].append(aa)
        res_index_chain += 1

    for chain in chains:
        chain_seqs[chain] = "".join([x for x in chain_seqs[chain] if x is not None])
    
    return pdbids, chain_seqs

def generate_bias_json(chain_seqs, bias_res, custom_bias, default_bias, opf):

    print(chain_seqs)
    chains_bias = {}
    for chain in chain_seqs:
        chain_bias = []
        for i in range(len(chain_seqs[chain])):
            res_id = chain + str(i+1)
            if res_id in bias_res:
                chain_bias.append(custom_bias)
            else:
                chain_bias.append(default_bias)
        chains_bias[chain] = chain_bias
                

    jsonobj = json.dumps(chains_bias)

    with open(opf, "w") as outfile:
        outfile.write(jsonobj)

def parse_mutres_input(mutresstring, pdbids):
    mutres_temp = mutresstring.strip().split(",")
    mutres_temp = [x.strip() for x in mutres_temp if x]
    mutres = []
    for elem in mutres_temp:
        if "-" not in elem:
            mutres.append(elem)
        else:
            start, finish = elem.split("-")
            chain = re.split('(\d+)', start)[0]
            s = int(re.split('(\d+)', start)[1]) 
            f = int(re.split('(\d+)', finish)[1]) 
            for i in range(s, f+1):
                mutres.append(chain + str(i))
    return mutres
    

if __name__=="__main__":
    parser = getBiasParser()
    args = parser.parse_args(sys.argv[1:])
    
    if args.sequence_file:
        pdbids, chain_seqs = parse_seqfile(args.sequence_file)
        
    else:
        pdbids, chain_seqs = parse_pdbfile(args.pdb)
    
    bias_res = parse_mutres_input(args.bias_res, pdbids)
    
    custom_bias = [float(x) for x in args.custom_bias.split(" ")]
    default_bias = [float(x) for x in args.default_bias.split(" ")]
    
    print(bias_res, custom_bias, default_bias)

    generate_bias_json(chain_seqs, bias_res, custom_bias, default_bias, args.output)
