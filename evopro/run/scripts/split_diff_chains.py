import os
import argparse

from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


import warnings
warnings.filterwarnings("ignore") 

def get_seq_from_pdb(pdbfilename):
    seq = ""
    with open(pdbfilename, "r") as f:
        for record in SeqIO.parse(f, 'pdb-atom'):
            seq += "," + record.seq
    
    return seq

def split_chains_diffusion_output(pdb_dir = './'):

    pdbs = [x for x in os.listdir(pdb_dir) if x.endswith(".pdb")]
    
    for pdb in pdbs:
        full_pdb = os.path.join(pdb_dir, pdb)
        seq = get_seq_from_pdb(os.path.join(pdb_dir, pdb)).strip(",")
        name = pdb[:-4]
        
        index = seq.rfind("GGG") + 3
        print(seq, index, seq[index:])
        
        structure = PDBParser().get_structure("", full_pdb)
        res_to_change = []

        for model in structure:
            
            for chains in model:
                for residues in chains:                    
                    if residues.get_id()[1] > index:
                        res_to_change.append(residues)

        for model in structure:
            for chain in model:
                [chain.detach_child(res.get_id()) for res in res_to_change]
        
        my_chain = Chain("B")

        model.add(my_chain)

        for res in res_to_change:
            my_chain.add(res)
            
        io = PDBIO()
        io.set_structure(model)
        savename = os.path.join(pdb_dir, "design{}_reformatted.pdb".format(name))

        io.save(savename,  write_end = True, preserve_atom_numbering = True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, default='./pdbs/')
    args = parser.parse_args()
    
    split_chains_diffusion_output(args.pdb_dir)
