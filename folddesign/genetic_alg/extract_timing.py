import bz2
import os
import pickle
import _pickle as cPickle
from typing import Any, Dict

def decompress_pickle(file: str) -> Any:
    """
    Load any compressed pickle file.
    """
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data

def print_timing(timing: Dict[str, float]) -> None:
    """
    Prints timing results (stored in dict) in a prettier format.
    """
    for k, v in timing.items():
        print(f'{k} took {v:.2f} sec.')

if __name__=="__main__":
    import glob, os

    #times = []
    for filename in glob.iglob('/pine/scr/a/m/amritan/kuhlmanlab/fdd/fdd/genetic_alg/**', recursive=False):
        if os.path.isdir(filename): # filter dirs
            time_file = filename+"/timing.pbz2"
            if os.path.isfile(time_file):
                p = decompress_pickle(time_file)
                print_timing(p)
                #times.append(float(p["predict_model_1_ptm_0"]))
    #print("\nAvg =", sum(times)/len(times))
