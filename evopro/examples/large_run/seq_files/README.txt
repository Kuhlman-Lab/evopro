1. generating directories for large runs:

specify "representative directory" with files to copy to each directory, number of replicates for each pair of sequences, and other options in largerun.flags.
"python /proj/kuhl_lab/evopro/evopro/run/create_dirs_largerun.py @largerun.flags" from directory

2. running evopro:

"sbatch submit_array.sh" from main directory to start a job array of all directories

