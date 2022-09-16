import json
import re
import sys
from inputs import FileArgumentParser

def generate_json(pdbfile, mut_res, opf, default):
    from folddesign.utils.pdb_parser import get_coordinates_pdb
    from folddesign.utils.aa_utils import three_to_one

    chains, residues, resindices = get_coordinates_pdb(pdbfile, fil=True)

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

    mutable = []
    for resind in mut_res:
        if resind in pdbids:
            mutable.append({"chain":pdbids[resind][1], "resid": pdbids[resind][2], "WTAA": three_to_one(pdbids[resind][0].split("_")[1]), "MutTo": default})

    dictionary = {"sequence" : chain_seqs, "designable": mutable}
    jsonobj = json.dumps(dictionary, indent = 4)

    with open(opf, "w") as outfile:
        outfile.write(jsonobj)
    
def getPDBParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run generate_json"""

    parser = FileArgumentParser(description='Script that can take a PDB and PDB residue numbers'
                                ' and convert to json file format for input to FDD', 
                                fromfile_prefix_chars='@')

    parser.add_argument('--pdb',
                        default='',
                        type=str,
                        help='Path to and name of PDB file')
    parser.add_argument('--mut_res',
                        default='',
                        type=str,
                        help='PDB chain and residue numbers to mutate, separated by commas')
    parser.add_argument('--default_mutres_setting',
                        default='all',
                        type=str,
                        help='Default setting for residues to mutate. Individual ones can be changed manually. Default is all')
    parser.add_argument('--output',
                        default='',
                        type=str,
                        help='path and name of output json file')
    return parser

def parse_mutres_input(mutresstring):
    mutres_temp = mutresstring.strip().split(",")
    mutres_temp = [x for x in mutres_temp if x]
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
    parser = getPDBParser()
    args = parser.parse_args(sys.argv[1:])
    mutres = parse_mutres_input(args.mut_res)
    generate_json(args.pdb, mutres, args.output, args.default_mutres_setting)

