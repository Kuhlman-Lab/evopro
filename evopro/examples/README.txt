1. creating json file from pdb or sequences:

specify pdb file/seqfile and residues to mutate in json.flags
"python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @json.flags" from directory

2. make sure residue_specs.json file was created correctly, create starting_seqs.txt file (optional)

3. running evopro:

specify input folder, score function, and other options in evopro.flags

"sbatch run_ep.sh" from input directory

4. analyze:

log files for scores, pdb files for output sequences/structures written to input directory
