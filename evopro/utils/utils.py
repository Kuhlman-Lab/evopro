""" 
Utility functions, includes saving and loading data and others.
"""

# Standard imports.
import bz2
import pickle
import _pickle as cPickle
import hashlib
from typing import Any, Dict
from Bio.PDB import MMCIFParser, PDBIO
from io import StringIO

def cif_to_pdb(mmcif_str):
    parser = MMCIFParser()
    cif_fh = StringIO(mmcif_str) 
    structure = parser.get_structure("structure", cif_fh)
    
    # Truncate residue names to 3 letters
    for model in structure:
        for chain in model:
            for residue in chain:
                if len(residue.resname) > 3:
                    residue.resname = residue.resname[:3]
    
    io = PDBIO()
    io.set_structure(structure)
    output = StringIO()
    io.save(output)
    pdb_string = output.getvalue()
    
    return pdb_string

def full_pickle(title: str, data: Any) -> None:
    """
    Saves the 'data' with the 'title' and adds the extension .pkl.
    """
    pikd = open(title + '.pkl', 'wb')
    pickle.dump(data, pikd)
    pikd.close()


def loosen(file: str) -> Any:
    """
    Loads and returns a pickled object.
    """
    pikd = open(file, 'rb')
    data = pickle.load(pikd)
    pikd.close()
    return data


def compressed_pickle(title: str, data: Any) -> None:
    """
    Pickle a file and then compress it into a file with extension .pbz2.
    """

    with bz2.BZ2File(title + '.pbz2', 'w') as f:
        cPickle.dump(data, f)


def decompress_pickle(file: str) -> Any:
    """
    Load any compressed pickle file.
    """
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data


def get_hash(x: str) -> str:
    """
    Looks up and returns a hash of given string.
    """
    return hashlib.sha1(x.encode()).hexdigest()


def print_timing(timing: Dict[str, float]) -> None:
    """
    Prints timing results (stored in dict) in a prettier format.
    """
    for k, v in timing.items():
        print(f'{k} took {v:.2f} sec.')
        
def merge_related_lists(list_of_pairs):
    """
    Takes a list of two-element lists and merges lists that share common elements.
    
    Args:
        list_of_pairs: List of lists, where each inner list contains two elements
        
    Returns:
        List of lists, where related elements are grouped together
    """
    # Initialize result list to store merged groups
    merged_groups = []
    
    # Convert pairs to sets for easier comparison and merging
    sets = [set(pair) for pair in list_of_pairs]
    
    while sets:
        current = sets.pop(0)  # Take the first set
        merged = False
        
        # Check against existing merged groups
        for group in merged_groups:
            if current & group:  # If there's any overlap
                group.update(current)  # Merge current into existing group
                merged = True
                break
                
        if not merged:
            merged_groups.append(current)
            
        # Check if any remaining sets can be merged with existing groups
        i = 0
        while i < len(sets):
            merged_with_existing = False
            for group in merged_groups:
                if sets[i] & group:  # If there's any overlap
                    group.update(sets[i])
                    sets.pop(i)
                    merged_with_existing = True
                    break
            if not merged_with_existing:
                i += 1
    
    # Convert sets back to sorted lists
    return [sorted(list(group)) for group in merged_groups]