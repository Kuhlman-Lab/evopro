1. generating directories for large runs:

specify "representative directory" with files to copy to each directory, number of replicates for each pair of sequences, and other options in largerun.flags.
"python /path/to/evopro/evopro/run/create_dirs_largerun.py @largerun.flags" from directory

2. running evopro:

"python /path/to/evopro/evopro/run/run_largerun.py" from main directory to start jobs in all directories

