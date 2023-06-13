1. creating json file from pdb:

specify pdb file and residues to mutate in pdb_to_json_flags.txt
"python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @pdb_to_json.txt" from directory

2. make sure residue_specs.json file was created correctly, create starting_seqs.txt file (optional)

3. running evopro:

specify input folder and score function in evopro.flags (See /proj/kuhl_lab/evopro/evopro/user_inputs/inputs.py for all flags options)

"sbatch run_ep.sh" from input directory

4. analyze:

log files for scores, pdb files for output sequences/structures written to outputs directory
