import argparse
import numpy as np
import sys
sys.path.append('/proj/kuhl_lab/pdb_utils')
from renumber_pdb import from_pdb_string, to_pdb
import os
import glob
import re
from Bio.PDB import PDBParser, PDBIO

def gap_or_ungap_chains(pdb_string, mode="gap", gap_size=100, chains=""):
    protein = from_pdb_string(pdb_string)

    if not chains:
        chains = "".join(list(protein.chain_id_mapping.keys()))

    if mode == "gap":
        n_chains = max(protein.chain_index)
        for i in range(n_chains + 1):
            n_res = np.sum(protein.chain_index == i)
            protein.residue_index[protein.chain_index == i] = np.arange(1, n_res + 1)

        max_res_num = max(protein.residue_index[protein.chain_index == protein.chain_id_mapping[chains[0]]])
        offset = max_res_num + gap_size - 1
        for chain in chains[1:]:
            protein.residue_index[protein.chain_index == protein.chain_id_mapping[chain]] += offset
            protein.chain_index[protein.chain_index == protein.chain_id_mapping[chain]] = protein.chain_id_mapping[chains[0]]

            max_res_num = max(protein.residue_index[protein.chain_index == protein.chain_id_mapping[chains[0]]])
            offset = max_res_num + gap_size - 1
            
        pdb_str = to_pdb(protein)
        #output_filename = pdb_path[:-4] + "_gapped.pdb"
##############
    else:
        assert len(chains) == 1, "Can only ungap 1 chain at a time. To do multiple chains run this multiple times."

        merged_residue_index = protein.residue_index[protein.chain_index == protein.chain_id_mapping[chains]]
        index_diff = abs(merged_residue_index[:-1] - merged_residue_index[1:])
        gaps = np.where(index_diff == gap_size)[0]

        split_residue_index = []
        for i in range(len(gaps)):
            if i == 0:
                split_residue_index.append(merged_residue_index[:gaps[i] + 1])
            else:
                split_residue_index.append(merged_residue_index[gaps[i - 1] + 1:gaps[i] + 1])
        split_residue_index.append(merged_residue_index[gaps[-1] + 1:])

        new_residue_index = [np.arange(1, len(split) + 1) for split in split_residue_index]
        new_chain_index = []
        counter = 1
        for i in range(len(gaps)):
            if i == 0:
                split_residue_index.append(merged_residue_index[:gaps[i] + 1])
            else:
                split_residue_index.append(merged_residue_index[gaps[i - 1] + 1:gaps[i] + 1])
        split_residue_index.append(merged_residue_index[gaps[-1] + 1:])

        new_residue_index = [np.arange(1, len(split) + 1) for split in split_residue_index]
        new_chain_index = []
        counter = 1
        for i in range(len(split_residue_index)):
            if i == 0:
                new_chain_index.append(protein.chain_id_mapping[chains] * np.ones(split_residue_index[i].shape))
            else:
                new_chain_index.append((max(protein.chain_index) + counter) * np.ones(split_residue_index[i].shape))
                counter += 1

        protein.residue_index[protein.chain_index == protein.chain_id_mapping[chains]] = np.concatenate(new_residue_index)
        protein.chain_index[protein.chain_index == protein.chain_id_mapping[chains]] = np.concatenate(new_chain_index)

        pdb_str = to_pdb(protein)
        #output_filename = pdb_path[:-4] + "_gapped.pdb"

    #with open(output_filename, 'w') as f:
    #    f.write(pdb_str)
    #print(f"Output saved to {output_filename}")
    return pdb_str


def keep_first_model_only(input_pdb_path, output_pdb_path):
    parser = PDBParser()
    structure = parser.get_structure("structure", input_pdb_path)

    first_model = next(structure.get_models())  # Get the first model

    # Initialize PDBIO and set the structure to the first model
    io = PDBIO()
    io.set_structure(first_model)

    # Save the first model to a new file
    io.save(output_pdb_path)

if __name__ == "__main__":
    folder = '/work/users/z/y/zyyao/118.nomultimer_multiple_pippack_test02/outputs/diff_out1_0/02_templates_multimer'
    for pdb in glob.glob(os.path.join(folder, 'design_*.pdb')):
        seq_number_match = re.search(r'design_(\d+).pdb', os.path.basename(pdb))
        if seq_number_match:
            seq_number = seq_number_match.group(1)
            new_filename = f"ds{seq_number}.pdb"  # Filename without path
            new_filepath = os.path.join(folder, new_filename)  # Corrected to include the path
            output_filename = gap_or_ungap_chains(pdb, mode="gap")  # Assumes this returns a path to the new file
            keep_first_model_only(output_filename, new_filepath)  # Use new_filepath here
            os.remove(pdb)  # Optionally, remove the original file
            os.remove(output_filename)
        else:
            print(f"No match for {os.path.basename(pdb)}. Skipping file.")