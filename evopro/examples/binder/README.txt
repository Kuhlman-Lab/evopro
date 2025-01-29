EvoPro for binder design
1. Place target protein sequence and scaffold binder sequence in seqfile.txt as separate chains.

2. Specify which regions of the scaffold are allowed to mutate (and to what amino acids) in json.flags (Also symmetry, if needed). Then, generate the "residue_specs.json" file by running the command below:

python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @json.flags

Manually check the "residue_specs.json" file that was generated and modify if needed.

3. Select region on target protein of interest for binding and specify in evopro.flags under --define_contact_area flag. Check and modify other flags if needed, including the score file/function which dictates how the designs are evaluated. Flag options can be found starting on line 39 of /proj/kuhl_lab/evopro/evopro/user_inputs/inputs.py.

4. Provide af2.flags for AlphaFold2 runs (or RF2, based on Evopro version). Here, we recommend running AF2 on your target protein and using the MSA generated as a custom MSA. To do this, use "--msa_mode single_sequence" and "--custom_msa_path /path/to/MSA/here/". You can also supply custom templates here using "--use_templates" and "--custom_template_path /path/to/templates/" if the MSA alone is not sufficient to fold the target.

IMPORTANT: if you do not use single sequence mode for the MSA the AF2 runs in every iteration will attempt to query an external server to calculate the MSA, greatly increasing the computational time and subjecting the run to fail due to errors in server querying.

5. Once all the input files have been assembled, run the script run_evopro.sh. Do a test run first to figure out how many iterations it takes for the scores to plateau.
