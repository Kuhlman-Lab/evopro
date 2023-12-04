1. creating json file from pdb:

specify pdb file and residues to mutate in json.flags
"python /path/to/evopro/evopro/run/generate_json.py @json.flags" from directory
TO GENERATE BIAS JSON:
"python /path/to/evopro/evopro/run/generate_mpnn_bias_dict.py @bias.flags" from directory

2. make sure residue_specs.json file was created correctly, create starting_seqs.txt file (optional)

3. running evopro:

specify input folder and score function in evopro.flags. customize other flags.

"sbatch run_ep.sh" from input directory

4. analyze:

log files and plots for scores, pdb files for output sequences/structures written to outputs directory
