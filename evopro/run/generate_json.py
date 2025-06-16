import argparse
import pprint
import re
import sys
import math
from typing import Optional, Sequence

from evopro.objects.chemical import to1letter, to3letter, alphabet, ptms, ptms_dict
from evopro.utils.parsing_utils import get_coordinates_pdb_extended, constituents_of_modified_fasta


class FileArgumentParser(argparse.ArgumentParser):
    """Overwrites default ArgumentParser to better handle flag files."""

    def convert_arg_line_to_args(self, arg_line: str) -> Optional[Sequence[str]]:
        """ Read from files where each line contains a flag and its value, e.g.
        '--flag value'. Also safely ignores comments denoted with '#' and
        empty lines.
        """

        # Remove any comments from the line
        arg_line = arg_line.split('#')[0]

        # Escape if the line is empty
        if not arg_line:
            return None

        # Separate flag and values
        split_line = arg_line.strip().split(' ')

        # If there is actually a value, return the flag value pair
        if len(split_line) > 1:
            return [split_line[0], ' '.join(split_line[1:])]
        # Return just flag if there is no value
        else:
            return split_line

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

def parse_seqfile(filename):
    ids = {}
    chains = {}
    with open(filename, "r") as f:
        i=0
        for lin in f:
            # loop over lines of seqfile.txt
            l = lin.strip().split(":")
            if len(l[2]) > 1:
                # chain[A] = [type, sequence]
                chains[l[0]] = [l[1],"".join(l[2:])]
            else:
                chains[l[0]] = [l[1],l[2]]
    
    for chain in chains:
        # parse ligand chains
        if chains[chain][0] == "ligand":
            reslist = chains[chain][1]
            resid = chain + "1"
            extended_resid = chain + "_LIG_1"
            ids[resid] = (extended_resid, chain, 1)
        # parse non-ligand chains
        else:
            reslist = constituents_of_modified_fasta(chains[chain][-1], chains[chain][0])
            # print("70: parsing non-ligand chains")
            # print(chains[chain][-1])
            # print(chains[chain][0])
            # print(reslist)
            for aa, i in zip(reslist, range(len(reslist))):
                resid = chain + str(i+1)
                if len(aa) == 1:
                    extended_resid = chain + "_" + to3letter[aa] + "_" + str(i+1)
                else:
                    extended_resid = chain + "_" + aa + "_" + str(i+1)
                ids[resid] = (extended_resid, chain, i+1)

    return ids, chains

def parse_pdbfile(filename):
    chains, residues, _ = get_coordinates_pdb_extended(filename, fil=True)

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

        try:
            aa = to1letter[residue.split("_")[1]]
        except:
            aa = "X"
        
        chain_seqs[chain].append(aa)
        res_index_chain += 1

    for chain in chains:
        chain_seqs[chain] = "".join([x for x in chain_seqs[chain] if x is not None])

    for chain_id in chain_seqs:
        sequence = chain_seqs[chain_id]
        type = "protein"
        if "a" in sequence or "t" in sequence or "c" in sequence or "g" in sequence:
            type = "dna"
        elif "b" in sequence or "u" in sequence or "d" in sequence or "h" in sequence:
            type = "rna"
        elif "X" in sequence:
            type = "ligand"
        chain_seqs[chain_id] = [type, sequence]
        
    return pdbids, chain_seqs

def generate_json(pdbids, chain_seqs, mut_res, opf, default, symmetric_res):
    if "all" in default:
        if default == "all":
            include = alphabet
        else:
            omit = list(default.split("-")[1])

            
        include = [x for x in alphabet if x not in omit]
        default = "".join(include)
    
    mods = []
    for pdbid in pdbids:
        pdbid = pdbids[pdbid][0]
        chain, res, num = pdbid.split("_")
        if res not in to1letter and res not in ["LIG", "dna"]:
            # add ptm modifications
            if res in ptms:
                mods.append({"chain":chain, "resid":num, "type":res})
            else:
                raise ValueError("Invalid residue type: " + res + ". Please use a valid residue type or PTM. You can add more PTM types to ptms list in objects/chemical.py.")
            
    mutable = []
    if "~" in "".join(mut_res):
        try:
            fixed_res = []
            fixed_seqs = [seq.strip("~") for seq in mut_res]
            #print(fixed_seqs)
            for c in chain_seqs:
                for fixed in fixed_seqs:
                    
                    if fixed in chain_seqs[c][-1]:
                        indices = list(find_all(chain_seqs[c][-1], fixed))

                        for i in indices:
                            
                            for j in range(i, len(fixed)+i):
                                fixed_res.append(pdbids[str(c) + str(j+1)])
            
            for pdbid in pdbids:
                if pdbids[pdbid] not in fixed_res:
                    mutable.append({"chain":pdbids[pdbid][1], "resid": pdbids[pdbid][2], "WTAA": to1letter[pdbids[pdbid][0].split("_")[1]], "MutTo": default})
        except:
            raise ValueError("Invalid mut_res specification. Try modifying the input syntax.")
                
    else:
        for resind in mut_res:
            if "*" in resind:
                try:
                    chain = resind.split("*")[0]
                    for pdbid in pdbids:
                        if pdbid.startswith(chain):
                            mutable.append({"chain":pdbids[pdbid][1], "resid": pdbids[pdbid][2], "WTAA": to1letter[pdbids[pdbid][0].split("_")[1]], "MutTo": default})
                except:
                    raise ValueError("Invalid specification. Try using the asterisk after the chain ID.")
            #starts with residue number specified (inclusive)
            elif "<G" in resind:
                residue = resind.split("<")[0]
                chain = re.split('(\\d+)', residue)[0]
                num_id = int(re.split('(\\d+)', residue)[1])
                chain_seq = chain_seqs[chain][-1]
                for i in range(num_id-1, len(chain_seq)):
                    if chain_seq[i] == "G":
                        mutable.append({"chain":chain, "resid": i+1, "WTAA": "G", "MutTo": default})
                    else:
                        break
            
            elif "<" in resind:
                try:
                    residue = resind.split("<")[0]
                    chain = re.split('(\\d+)', residue)[0]
                    num_id = int(re.split('(\\d+)', residue)[1])
                    for pdbid in pdbids:
                        chain_compare = re.split('(\\d+)', pdbid)[0]
                        num_id_compare = int(re.split('(\\d+)', pdbid)[1])
                        if chain == chain_compare and num_id<=num_id_compare:
                            mutable.append({"chain":pdbids[pdbid][1], "resid": pdbids[pdbid][2], "WTAA": to1letter[pdbids[pdbid][0].split("_")[1]], "MutTo": default})
                except:
                    raise ValueError("Invalid mut_res specification. Try using the less than sign after the residue ID.")
                
            elif resind in pdbids:
                mutable.append({"chain": pdbids[resind][1], "resid": pdbids[resind][2], "WTAA": to1letter[pdbids[resind][0].split("_")[1]], "MutTo": default})

    symmetric = []
    for symmetry in symmetric_res:
        values = list(symmetry.values())
        for tied_pos in zip(*values):
            skip_tie = False
            sym_res = []
            for pos in tied_pos:
                res_id = pos[0] + str(pos[1])
                if pdbids[res_id][0] == 'X':
                    skip_tie = True
                    break
                else:
                    sym_res.append( pdbids[res_id][1] + str(pdbids[res_id][2]) )

            if not skip_tie:
                symmetric.append(sym_res)

    chain_seqs_mod = {}
    for chain in chain_seqs:
        if chain_seqs[chain][0] == "ligand":
            chain_seqs_mod[chain] = {"sequence":chain_seqs[chain][-1], "type":chain_seqs[chain][0]}
        else:
            seq_input = constituents_of_modified_fasta(chain_seqs[chain][-1], chain_seqs[chain][0])
            # loop over seq and convert ptms to their 1 letter code
            for i, val in enumerate(seq_input):
                if val in ptms:
                    # print(f"237: found ptm {val}, converting to {ptms[val]}")
                    seq_input[i] = ptms_dict[val]
            # print("240, printing seq input")
            # print(seq_input)
            chain_seqs_mod[chain] = {"sequence":seq_input, "type":chain_seqs[chain][0]}
            # chain_seqs_mod[chain] = {"sequence":constituents_of_modified_fasta(chain_seqs[chain][-1], chain_seqs[chain][0]), "type":chain_seqs[chain][0]}
    # here for modified residues, the sequence should be written using the original residue
    #   chain_seqs_mod should be dictionary of chains written in 1 letter code
    dictionary = {"chains" : chain_seqs_mod, "modifications" : mods, "designable": mutable, "symmetric": symmetric}
    
    # write json to file with human-friendly formatting
    jsonobj = pprint.pformat(dictionary, compact=True, sort_dicts=False).replace("'",'"')

    with open(opf, "w") as outfile:
        outfile.write(jsonobj)
    
def getPDBParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run generate_json"""

    parser = FileArgumentParser(description='Script that can take a PDB and PDB residue numbers'
                                ' and convert to json file format for input to EvoPro', 
                                fromfile_prefix_chars='@')

    parser.add_argument('--pdb_file',
                        default=None,
                        type=str,
                        help='Path to and name of PDB file to extract chains and sequences.')
    parser.add_argument('--sequence_file',
                        default=None,
                        type=str,
                        help='Path to and name of text file to extract chains and sequences. Only provide if there is no PDB file.')
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
    parser.add_argument('--symmetric_res',
                        default='',
                        type=str,
                        help='PDB chain and residue numbers to force symmetry separated by a colon')
    return parser

def parse_mutres_input(mutresstring):
    mutres_temp = mutresstring.strip().split(",")
    mutres_temp = [x.strip() for x in mutres_temp if x]
    mutres = []
    for elem in mutres_temp:
        if "-" not in elem:
            mutres.append(elem)
        else:
            start, finish = elem.split("-")
            chain = re.split('(\\d+)', start)[0]
            s = int(re.split('(\\d+)', start)[1]) 
            f = int(re.split('(\\d+)', finish)[1]) 
            for i in range(s, f+1):
                mutres.append(chain + str(i))
    return mutres

def _check_res_validity(res_item):
    split_item = [item for item in re.split('(\\d+)', res_item) if item]

    if len(split_item) != 2:
        raise ValueError(f'Unable to parse residue: {res_item}.')
    return (split_item[0], int(split_item[1]))

def _check_range_validity(range_item):
    
    split_range = range_item.split('-')
    if len(split_range) != 2:
        raise ValueError(f'Unable to parse residue range: {range_item}')

    start_item, finish_item = split_range[0], split_range[1]

    s_chain, s_idx = _check_res_validity(start_item)
    f_chain, f_idx = _check_res_validity(finish_item)
    if s_chain != f_chain:
        raise ValueError(f'Residue ranges cannot span multiple chains: {range_item}')
    if s_idx >= f_idx:
        raise ValueError(f'Residue range starting index must be smaller than the ending index: '
                             f'{range_item}')

    res_range = []
    for i in range(s_idx, f_idx + 1):
        res_range.append( (s_chain, i) )

    return res_range

def _check_range_validity_asterisk(range_item, pdbids):
    res_range = []
    chain = range_item.split("*")[0]
    
    for pdbid in pdbids:
        chain_id = re.split('(\\d+)', pdbid)[0]
        if chain == chain_id:
            res_range.append( (chain_id, pdbids[pdbid][2]) )
    
    if len(res_range) == 0:
        raise ValueError(f'Unable to parse residue range: {range_item}')

    return res_range

def _check_symmetry_validity(symmetric_item, pdbids):
    split_item = symmetric_item.split(':')

    symmetry_dict = {}
    for subitem in split_item:
        if '-' in subitem:
            res_range = _check_range_validity(subitem)
            symmetry_dict[subitem] = res_range
        elif '*' in subitem:
            res_range = _check_range_validity_asterisk(subitem, pdbids)
            symmetry_dict[subitem] = res_range
        else:
            res_ch, res_idx = _check_res_validity(subitem)
            symmetry_dict[subitem] = [(res_ch, res_idx)]

    item_lens = [len(symmetry_dict[key]) for key in symmetry_dict]
    if math.floor(sum([l == item_lens[0] for l in item_lens])/len(item_lens)) != 1:
        raise ValueError(f'Tied residues and residue ranges must be of the same '
                             f'size for forcing symmetry: {symmetric_item}')

    return symmetry_dict

def parse_symmetric_res(symmetric_str, pdbids):

    symmetric_str = [s for s in symmetric_str.strip().split(",") if s]
    #print(symmetric_str)

    symmetric_res = []
    for item in symmetric_str:
        if ":" not in item:
            raise ValueError(f'No colon detected in symmetric res: {item}.')

        symmetry_dict = _check_symmetry_validity(item, pdbids)
        symmetric_res.append(symmetry_dict)
    #print(symmetric_res)

    return symmetric_res
    

if __name__=="__main__":
    parser = getPDBParser()
    args = parser.parse_args(sys.argv[1:])
    
    if args.sequence_file:
        pdbids, chain_seqs = parse_seqfile(args.sequence_file)
        
    else:
        pdbids, chain_seqs = parse_pdbfile(args.pdb_file)
        
    symres = parse_symmetric_res(args.symmetric_res, pdbids)
    mutres = parse_mutres_input(args.mut_res)

    generate_json(pdbids, chain_seqs, mutres, args.output, args.default_mutres_setting, symres)