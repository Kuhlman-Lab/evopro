import os
import sys
import numpy as np
import pandas as pd
import glob
import csv
import argparse
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.score_funcs.score_funcs import score_plddt_confidence
from evopro.score_funcs.score_funcs_efficient import score_contacts_pae_weighted_efficient
from evopro.utils.pdb_parser import get_coordinates_pdb

def score_pae_interaction(pae, pdb, chain1 = "A", chain2 = "B"):

    _, residues, _ = get_coordinates_pdb(pdb)    
    reslist2 = [x for x in residues if x.startswith(chain2)]   
    target_length = len(reslist2)
    pae_interaction1 = np.mean( pae[:target_length,target_length:] )
    pae_interaction2 = np.mean( pae[target_length:,:target_length] )
    
    return None, ( pae_interaction1 + pae_interaction2 ) / 2


def get_average_plddt_for_each_chain(plddt, pdb, chain1 = "A", chain2 = "B"):
    _, residues, resindices = get_coordinates_pdb(pdb)    

    reslist1 = [x for x in residues.keys() if x.startswith(chain1)]
    reslist2 = [x for x in residues.keys() if x.startswith(chain2)]
    #score full complex structure
    plddt_full = score_plddt_confidence(plddt, residues, resindices, rf2=True)
    #then, score chain A and B
    plddt_chainA = score_plddt_confidence(plddt, reslist1, resindices, rf2=True)
    plddt_chainB = score_plddt_confidence(plddt, reslist2, resindices, rf2=True)
    return plddt_full, plddt_chainA, plddt_chainB

if  __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent_dir", type=str, default='outputs/')
    args = parser.parse_args()
    
    parent_dir = args.parent_dir
    data_list = []
    info_output_list = []
    # Loop through all relevant folders and files
    folders = glob.glob(parent_dir + 'query*')
    
    sequence_file_path = os.path.join(parent_dir, 'sequences.csv')
    
    try:
        sequence_data = pd.read_csv(sequence_file_path, header=None)
    except FileNotFoundError:
        print(f"sequences.csv not found in {parent_dir}")

    for folder in folders:
        print(folder)
        for filepath in glob.glob(os.path.join(folder, '*.npz')):
            # Extract pLDDT
            print (filepath)
            result = np.load(filepath)
            PLDDT = result['lddt']
            PAE = result['pae']
            #average_plddt = np.average(result['plddt'][0:Ls[0]])
            path_parts = folder.split('/')
            last_part = path_parts[-1]
            seq_number = last_part.split('query_')[1]
            # Find the matching line in sequences.csv, find matching f"seq_{seq_number}_model_1.pdb"
            sequence_line = None
            if seq_number:
                seq_number_int = int(seq_number)  # Convert to integer for comparison
                if seq_number_int < len(sequence_data):
                    sequence_line = sequence_data.iloc[seq_number_int].to_list()
                    print(sequence_line)
                    binder_sequence = None
                    target_sequence = None
                    if sequence_line and len(sequence_line) == 3:
                        binder_sequence = sequence_line[1]  # Second element as binder sequence
                        target_sequence = sequence_line[2]  # Third element as target sequence
                    elif sequence_line and len(sequence_line) == 7:
                        target_sequence = sequence_line[1]
                        binder_sequence = sequence_line[4]
                    
                    else:
                        print(f"sequence.csv do not have correct 3 element as [nan, 'binder_sequence','target_sequence']  {seq_number} in {sequence_file_path}")
                else:
                    print(f"No matching sequence for seq_number {seq_number} in {sequence_file_path}")

            else:
                # Handle error case if pattern doesn't match
                continue
            
            # Extract RMSD and PAE, PLDDT
            pdb_filename = "S_00_pred.pdb"
            pdb_filepath = os.path.join(folder, pdb_filename)

            if os.path.exists(pdb_filepath):
                with open(pdb_filepath, 'r') as pdbf:
                    pdb_model = pdbf.read()

                _, residues, resindices = get_coordinates_pdb(pdb_model)
                #reslist1 = [x for x in residues.keys() if x.startswith("A")]
                #reslist2 = [x for x in residues.keys() if x.startswith("B")]
                reslist1 = [x for x in residues.keys() if x.startswith("A") or x.startswith("B") or x.startswith("C")]
                reslist2 = [x for x in residues.keys() if x.startswith("D")]
                pairs, interchain_PAE_new = score_pae_interaction(PAE, pdb_model, chain1 = "A", chain2 = "D")
                pairs, interchain_PAE_old = score_contacts_pae_weighted_efficient(result, pdb_model, reslist1, reslist2, rf2=True)
                average_plddt, binder_plddt, target_plddt = get_average_plddt_for_each_chain(result, pdb_model, chain1 = "A", chain2 = "D")
                average_plddt, binder_plddt, target_plddt = average_plddt*100, binder_plddt*100, target_plddt*100

            # Append to data
            data_list.append([binder_sequence, average_plddt, binder_plddt, interchain_PAE_old, interchain_PAE_new, folder, pdb_filepath, target_sequence ])

    # Sort & Write to CSV
    data_dicts = []
    for item in data_list:
        data_dict = {
            'binder_sequence': item[0],
            'average_plddt': item[1],
            'binder_plddt': item[2],
            'interchain_PAE_old': item[3],
            'interchain_PAE_new': item[4],
            'folder': item[5],
            'pdb_filepath': item[6],
            'target_sequence': item[7]
        }
        data_dicts.append(data_dict)

    data_dicts.sort(key=lambda x: (-x['average_plddt'], x['interchain_PAE_new']))  # Sort by average_plddt and interchain_PAE 

    with open('full_list_newPAE.csv', 'w', newline='') as csvfile1:
        csv_writer = csv.writer(csvfile1)
        csv_writer.writerow(["binder_sequence", "average_plddt", "binder_plddt", "interchain_PAE_old", "interchain_PAE_new", "folder", "pdb_filepath", "target_sequence"])
        for item in data_dicts:  
            row = [item['binder_sequence'], item['average_plddt'], item['binder_plddt'], item['interchain_PAE_old'], item['interchain_PAE_new'], item['folder'], item['pdb_filepath'], item['target_sequence']]
            csv_writer.writerow(row)