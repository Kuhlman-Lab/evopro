EvoPro is a genetic algorithm-based protein binder optimization pipeline, used for "in silico evolution". To run this pipeline, you will also need to clone our AlphaFold2 and ProteinMPNN repositories (link here).
After cloning all directories, they need to be linked. The path to AF2 and ProteinMPNN directories should be set at the top of the script evopro/evopro/run_geneticalg_gpus.py where indicated.
Additionally, you will need to set up the af2_mpnn conda environment.

1. Creating residue_specs.json file from pdb or sequence file:

specify pdb/sequence file and residues to mutate in json.flags
"python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @json.flags" from directory

2. Make sure residue_specs.json file was created correctly/make changes as needed.

3. Add af2.flags file to the running directory and set options for MSA and template generation.

4. Running evopro:
Specify input folder and score function in evopro.flags (See evopro/evopro/user_inputs/inputs.py for all flags options and descriptions)

Commands needed to run EvoPro can be found in test_ep.sh in the examples directories. The main one is:
"python /path/evopro/evopro/run/run_geneticalg_gpus.py @evopro.flags", which should be run on a system with GPUs available and conda environment activated.

5. Analyze results:
log files for scores, pdb files for output sequences/structures written to outputs directory
