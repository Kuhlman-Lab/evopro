from pyrosetta import *
import os
import numpy as np
init()

from pyrosetta.rosetta.core.scoring import get_score_function

from pyrosetta.rosetta.core.pose import total_energy_from_pose

from Bio import SeqIO

def get_filenames(dirname):
    names = []
    for filename in os.listdir(dirname):
        f = os.path.join(dirname, filename)
        # checking if it is a file
        if os.path.isfile(f):
            #checking if the file is a pdb file
            if f.endswith('.pdb'):
                without_packed = f.split("_packed_")[0]
                names.append(without_packed)

    return names

def get_energy(path):
    pose = pose_from_file(path)
    score_fn = get_score_function(True)
    #score_fn.show(pose) #run if troubleshooting
    energy = total_energy_from_pose(pose)
    return energy


def get_seq_number(name, dirname):
    file_without_dir = name.split("/")[1]
    avgs = {}
    files = [filename for filename in os.listdir(dirname) if filename.startswith(file_without_dir)]
    if len(files) == 1:
        after_packed = files[0].split("_packed_")[1]
        num = after_packed[0]
        return num

    for file in files:
        path = dirname + file
        after_packed = file.split("_packed_")[1]
        num = after_packed[0]

        if num not in avgs:
            energy = get_energy(path)
            avgs[num] = [1, energy]
        else:
            lst = avgs.get(num)
            lst[0] += 1
            energy = get_energy(path)
            lst[1] += energy

    for key in avgs:
        vals = avgs[key]
        total = vals[1]
        counter = vals[0]
        avgs[key] = total/counter

    lowest_energy = min(avgs.values())
    best_sequence_number = [key for key in avgs if avgs[key] == lowest_energy]

    return best_sequence_number[0]


def get_best_pdb_energy(name, dirname):
    file_without_dir = name.split("/")[1]
    avgs = {}
    files = [filename for filename in os.listdir(dirname) if filename.startswith(file_without_dir)]
    if len(files) == 1:
        after_packed = files[0].split("_packed_")[1]
        num_minus_pdb = after_packed.split(".")[0]
        return num_minus_pdb

    for file in files:
        path = dirname + file
        after_packed = file.split("_packed_")[1]
        num_minus_pdb = after_packed.split(".")[0]

        if num_minus_pdb not in avgs:
            energy = get_energy(path)
            avgs[num_minus_pdb] = [1, energy]
        else:
            lst = avgs.get(num_minus_pdb)
            lst[0] += 1
            energy = get_energy(path)
            lst[1] += energy

    lowest_energy = min(avgs.values())
    best_pdb = [key for key in avgs if avgs[key] == lowest_energy]

    return best_pdb[0]

def choose_seq(unique_names, dirname):
    for name in unique_names:
        print(name)
        best_sequence_number = get_seq_number(name, dirname)
        best_pdb = get_best_pdb_energy(name, dirname)
        seq_path = "seqs/" + name.split("/")[1] + ".fa"
        pdb_path = "check_output/" + name.split("/")[1] + "_packed_" + best_pdb + ".pdb"
        new_seq_path = "best_seqs/" + name.split("/")[1] + ".fa"
        new_pdb_path = "best_output/" + name.split("/")[1] + ".pdb"

        os.system(f'cp {pdb_path} {new_pdb_path}') 

        seqs = list(SeqIO.parse(seq_path,"fasta"))
        #sequence = seqs[int(best_sequence_number)].seq
        sequence = seqs[int(best_sequence_number)]

        with open(new_seq_path, "w") as output_handle:
            SeqIO.write(sequence, output_handle, "fasta")

# # Sample implementation
# dirname = "check_output/"
# names = get_filenames(dirname)
# unique_names = np.unique(names)
# choose_seq(unique_names, dirname)