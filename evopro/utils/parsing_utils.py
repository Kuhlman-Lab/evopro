import re
import string
from string import ascii_letters
from rdkit import Chem

def constituents_of_modified_fasta(x: str, chain_type: str):
    """
    Accepts amino acid and RNA/DNA inputs: 'agtc', 'AGT(ASP)TG', etc. Does not accept SMILES strings.
    Returns constituents, e.g, [A, G, T, ASP, T, G] or None if string is incorrect.
    Everything in returned list is single character, except for blocks specified in brackets.
    """
    if chain_type == "protein":
        x = x.strip().upper()
    # it is a bit strange that digits are here, but [NH2] was in one protein
    allowed_chars = ascii_letters + "()" + string.digits
    if not all(letter in allowed_chars for letter in x):
        return None

    current_modified: str | None = None

    constituents = []
    for letter in x:
        if letter == "(":
            if current_modified is not None:
                return None  # double open bracket
            current_modified = ""
        elif letter == ")":
            if current_modified is None:
                return None  # closed without opening
            if len(current_modified) <= 1:
                return None  # empty modification: () or single (K)
            constituents.append(current_modified)
            current_modified = None
        else:
            if current_modified is not None:
                current_modified += letter
            else:
                if letter not in ascii_letters:
                    return None  # strange single-letter residue
                constituents.append(letter)
    if current_modified is not None:
        return None  # did not close bracket
    return constituents

def count_nonH_atoms(smiles):
    """
    Count the number of non-hydrogen atoms in a molecule given its SMILES string.
    
    Parameters:
    smiles (str): SMILES string representing the molecule
    
    Returns:
    int: Number of non-hydrogen atoms
    None: If SMILES string is invalid
    
    Example:
    >>> count_nonH_atoms("CC(=O)OH")  # Acetic acid
    4  # (2 carbons, 2 oxygens)
    """
    try:
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Get the number of atoms (RDKit by default excludes implicit hydrogens)
        return mol.GetNumAtoms()
        
    except Exception as e:
        print(f"Error processing SMILES string: {e}")
        return None
    
def get_tokens(resstring, dsobj):
    if not resstring:
        return []
    resstring = str(resstring)
    res_temp = resstring.strip().split(",")
    res_temp = [x.strip() for x in res_temp if x]
    res = []
    for elem in res_temp:
        if "-" in elem:
            start, finish = elem.split("-")
            chain = re.split('(\d+)', start)[0]
            s = int(re.split('(\d+)', start)[1]) 
            f = int(re.split('(\d+)', finish)[1]) 
            for i in range(s, f+1):
                res.append(chain + str(i))
        
        elif "*" in elem:
            chain = elem.split("*")[0]
            #modify this to get all residues in chain, not tokens
            for i in range(1, dsobj.get_lengths(chain)+1):
                res.append(chain + str(i))
        
        elif "~" in elem:
            motif = elem.split("~")[0]
            #TODO: find motif and set all other residues mutable
            raise NotImplementedError("Motif-based design not yet implemented")
            
        elif "<G" in elem:
            residue = elem.split("<G")[0]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
                
        elif "<" in elem:
            residue = elem.split("<")[0]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
            chain_length = dsobj.get_lengths(chain)
            for i in range(num_id, chain_length+1):
                res.append(chain + str(i))
        
        elif ">" in elem:
            residue = elem.split(">")[1]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
            chain_length = dsobj.get_lengths(chain)
            for i in range(1, num_id+1):
                res.append(chain + str(i))
        
        else:
            res.append(elem)
    return res

def get_residues(resstring, dsobj):
    if not resstring:
        return []
    resstring = str(resstring)
    res_temp = resstring.strip().split(",")
    res_temp = [x.strip() for x in res_temp if x]
    res = []
    for elem in res_temp:
        if "-" in elem:
            start, finish = elem.split("-")
            chain = re.split('(\d+)', start)[0]
            s = int(re.split('(\d+)', start)[1]) 
            f = int(re.split('(\d+)', finish)[1]) 
            for i in range(s, f+1):
                res.append(chain + str(i))
        
        elif "*" in elem:
            chain = elem.split("*")[0]
            #modify this to get all residues in chain, not tokens
            for i in range(1, dsobj.get_residue_lengths(chain)+1):
                res.append(chain + str(i))
        
        elif "~" in elem:
            motif = elem.split("~")[0]
            #TODO: find motif and set all other residues mutable
            raise NotImplementedError("Motif-based design not yet implemented")
            
        elif "<G" in elem:
            residue = elem.split("<G")[0]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
                
        elif "<" in elem:
            residue = elem.split("<")[0]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
            chain_length = dsobj.get_residue_lengths(chain)
            for i in range(num_id, chain_length+1):
                res.append(chain + str(i))
        
        elif ">" in elem:
            residue = elem.split(">")[1]
            chain = re.split('(\d+)', residue)[0]
            num_id = int(re.split('(\d+)', residue)[1])
            chain_length = dsobj.get_residue_lengths(chain)
            for i in range(1, num_id+1):
                res.append(chain + str(i))
        
        else:
            res.append(elem)
    return res

def get_coordinates_pdb_tokens(pdb, fil=False):
    lines = []
    chains = []
    residues = {}
    residueindices = {}
    
    if fil:
        with open(pdb, "r") as f:
            pdb_split = f.read().split("\n")
    else:
        pdb_split = pdb.split("\n")
        
    i = 0
    j = 0
    pdb_split = [x for x in pdb_split if x]
    hetatom_chains = []
    for lin in pdb_split:
        x = lin[30:38].strip(' ')
        y = lin[38:46].strip(' ')
        z = lin[46:54].strip(' ')
        l = lin.strip().split()
        
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            if 'HETATM' in l[0]:
                # For HETATM records, create a unique residue number for each atom
                # but maintain the original format of chain+residue_number
                if l[4] not in hetatom_chains:
                    hetatom_chains.append(l[4])
                    j=0
                resid = l[4] + str(int(l[5]) + j)  # Using atom number to make unique
                j+=1
            else:
                # For regular ATOM records, keep original residue grouping
                resid = l[4] + l[5]
                
            atominfo = (l[1], l[2], (x, y, z))
            
            if l[4] not in chains:
                chains.append(l[4])
                
            if resid not in residues:
                residues[resid] = [atominfo]
            else:
                residues[resid].append(atominfo)
            
            # Assign an index to the residue if it doesn't have one
            if resid not in residueindices:
                residueindices[resid] = i
                i = i+1 
    # print(residues)

    return chains, residues, residueindices

def get_coordinates_pdb(pdb, fil=False):
    lines = []
    chains = []
    residues = {}
    residueindices = {}
    if fil:
        with open(pdb,"r") as f:
            pdb_split = f.read().split("\n")

    else:
        pdb_split = pdb.split("\n")
        
    i=0
    pdb_split = [x for x in pdb_split if x]
    for lin in pdb_split:
        x = lin[30:38].strip(' ')
        y = lin[38:46].strip(' ')
        z = lin[46:54].strip(' ')
        l = lin.strip().split()
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            resid = l[4]+l[5]
            atominfo = (l[1], l[2], (x, y, z))
            if l[4] not in chains:
                chains.append(l[4])
            if resid not in residues:
                residues[resid] = [atominfo]
            else:
                residues[resid].append(atominfo)
            if resid not in residueindices:
                residueindices[resid] = i
                i = i+1 

    return chains, residues, residueindices

#resid in format chain_aatype_resind
def get_coordinates_pdb_extended(pdb, fil = False):
    chains = []
    residues = {}
    residueindices = {}
    if fil:
        with open(pdb,"r") as f:
            pdb_split = f.read().split("\n")
            
    else:
        pdb_split = pdb.split("\n")
        
    i=0
    pdb_split = [x for x in pdb_split if x]
    for lin in pdb_split:
        x = lin[30:38].strip(' ')
        y = lin[38:46].strip(' ')
        z = lin[46:54].strip(' ')
        l = lin.strip().split()
        if l and ('ATOM' in l[0] or 'HETATM' in l[0]):
            resid = lin[21] + '_' + lin[17:20].strip() + "_" + lin[22:26].strip()

            atominfo = (l[1], l[2], (x, y, z))
            if lin[21] not in chains:
                chains.append(lin[21])
            if resid not in residues:
                residues[resid] = [atominfo]
            else:
                residues[resid].append(atominfo)
            if resid not in residueindices:
                residueindices[resid] = i
                i = i+1

    return chains, residues, residueindices

def get_ca_coordinates_pdb(pdb_str):
    residues = {}
    coords = []

    pdb_split = pdb_str.split("\n")
    i=0
    pdb_split = [x for x in pdb_split if x]
    for lin in pdb_split:
        x = lin[30:38].strip(' ')
        y = lin[38:46].strip(' ')
        z = lin[46:54].strip(' ')
        l = lin.strip().split()
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            if l[2] == "CA" or l[2] == "C1'":

                coord = [float(x), float(y), float(z)]
                coords.append(coord)

    return coords

def get_token_coordinates_pdb(pdb_str):
    """Gets coordinates of atoms according to token in the pdb file. 
    For proteins, gets CA atoms. For nucleic acids, gets C1' atoms. For ligands, gets all atoms.

    Args:
        pdb_str (string): pdb file as a string

    Returns:
        list: list of coordinates of atoms
    """
    coords = []

    pdb_split = pdb_str.split("\n")
    pdb_split = [x for x in pdb_split if x]
    
    for lin in pdb_split:
        x = lin[30:38].strip(' ')
        y = lin[38:46].strip(' ')
        z = lin[46:54].strip(' ')
        l = lin.strip().split()
        
        # Check if line starts with ATOM or HETATM
        if 'ATOM' in l[0] or 'HETATM' in l[0]:
            # For ATOM records, only include CA atoms
            if l[0] == 'ATOM' and (l[2] == "CA" or l[2] == "C1'"):
                coord = [float(x), float(y), float(z)]
                coords.append(coord)
            # For HETATM records, include all atoms
            elif l[0] == 'HETATM':
                coord = [float(x), float(y), float(z)]
                coords.append(coord)

    return coords

def parse_weights(weights_string, chains):
    weights_string = str(weights_string).strip()
    weights = [float(x.strip()) for x in weights_string.split(",")]
    if len(weights) != len(chains):
        print("warning: number of weights does not match number of scores. using first weight for all scores")
        print("weights[0]", weights[0])
        weights = [weights[0]]*len(chains)
    return weights