import os
import sys
import shutil
import itertools
import string
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser

alphabet = list(string.ascii_uppercase)
array_job_file = "/proj/kuhl_lab/evopro/evopro/data/submit_array.sh"

def getFlagsParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script for large scale EvoPro runs"""

    parser = FileArgumentParser(description='get arguments to run script for large scale EvoPro runs',
                                fromfile_prefix_chars='@')

    #provide text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc.
    parser.add_argument('--sequence_files',
                        default='chainA.txt chainB.txt',
                        type=str,
                        help='Path to and name of sequence files in order of chains, separated by spaces.')
    
    #provide af2.flags, evopro.flags, json.flags, run_evopro.sh, templates folder(optional), other pdbs(optional) in this directory
    parser.add_argument('--representative_directory',
                        default='./',
                        type=str,
                        help='Path to directory for files that are needed to run EvoPro, including af2.flags, evopro.flags etc.')
    
    parser.add_argument('--num_replicates',
                        default=3,
                        type=int,
                        help='number of independent trajectories to set up per sequence combination')
    return parser

def parse_sequences(filename):
    sequences = []
    with open(filename, "r") as fil:
        for lin in fil:
            l = lin.strip()
            sequences.append(l)

    return sequences

if __name__=="__main__":
    parser = getFlagsParser()
    args = parser.parse_args(sys.argv[1:])

    main_dir = os.getcwd()

    ip_files = args.sequence_files.strip().split(" ")
    seqs = []
    for fil in ip_files:
        seqs.append(parse_sequences(fil))
    
    seq_permutations = list(itertools.product(*seqs))

    num_pairs = len(seq_permutations)
    num_jobs =  num_pairs * args.num_replicates

    for perm, i in zip(seq_permutations, range(len(seq_permutations))):
        #create a directory for each combination of sequences
        dir_name = "pair" + str(i+1)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        else:
            shutil.rmtree(dir_name)
            os.makedirs(dir_name)

        #within this, make a directory for each replicate
        for j in range(args.num_replicates):
            run_name = "run" + str(j+1)
            os.mkdir(os.path.join(dir_name, run_name))

            #write a seqfile to use to generate the json file
            with open(os.path.join(dir_name, run_name, "seqfile.txt"), "w") as seqf:
                for elem, chain in zip(perm, alphabet[:len(perm)]):
                    seqf.write(chain + ":" + elem + "\n")
            
            #and copy all other needed files and directories from representative directory into each replicate
            onlyfiles = [x for x in os.listdir(args.representative_directory) if os.path.isfile(os.path.join(args.representative_directory, x))]
            onlydirs = [x for x in os.listdir(args.representative_directory) if os.path.isdir(os.path.join(args.representative_directory, x))]
            #print(onlyfiles, onlydirs)
            for f in onlyfiles:
                shutil.copy(os.path.join(args.representative_directory, f), os.path.join(dir_name, run_name))
            for d in onlydirs:
                shutil.copytree(os.path.join(args.representative_directory, d), os.path.join(dir_name, run_name, d))

            #generate json for each dir
            os.chdir(os.path.join(main_dir, dir_name, run_name))
            os.system("python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @json.flags")
            os.chdir(main_dir)
            
    #generate bash script to create array job of all dirs that is run manually
    l = []
    with open(array_job_file, "r") as jobf:
        for lin in jobf:
            if "--array" in lin:
                l.append("#SBATCH --array=1-" + str(num_jobs) +"\n")
            elif lin.startswith("$N_PAIRS"):
                l.append("$N_PAIRS=" + str(num_pairs) + "\n")
            elif lin.startswith("$N_REPS"):
                l.append("$N_REPS=" + str(args.num_replicates) + "\n")
            else:
                l.append(lin)
    
    print("writing job file at:", main_dir + "/submit_array.sh")
    with open(main_dir + "/submit_array.sh", "w") as jobf:
        for lin in l:
            jobf.write(lin)

    
