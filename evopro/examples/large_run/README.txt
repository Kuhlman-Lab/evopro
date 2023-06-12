1. generating directories for large runs:
The script reads the .txt files and generates each individual run directory through combinations of the sequences in these files. So for each chain in chainA.txt, a directory is created with each chain in chainB.txt, and so on.
specify "representative directory" with files to copy to each directory, number of replicates for each pair of sequences, and other options in largerun.flags.
"python /proj/kuhl_lab/evopro/evopro/run/create_dirs_largerun.py @largerun.flags" from directory to generate the directories.

2. running evopro:

"sbatch submit_array.sh" from main directory to start a job array of all directories

