""" 
Utility functions, includes saving and loading data and others.
"""

# Standard imports.
import bz2
import pickle
import _pickle as cPickle
import hashlib
from typing import Any, Dict


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
    
