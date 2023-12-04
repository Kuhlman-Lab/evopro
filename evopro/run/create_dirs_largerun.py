<<<<<<< HEAD
=======
import sys
#SET PATH TO YOUR EVOPRO INSTALLATION HERE
#sys.path.append("/proj/kuhl_lab/evopro/")
sys.path.append("/nas/longleaf/home/amritan/Desktop/kuhlmanlab/evopro_temp/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser
>>>>>>> backup-stable
import os
import sys
import shutil
import itertools
import string
<<<<<<< HEAD
sys.path.append("/proj/kuhl_lab/evopro/")
from evopro.user_inputs.inputs import FileArgumentParser

alphabet = list(string.ascii_uppercase)
array_job_file = "/proj/kuhl_lab/evopro/evopro/data/submit_array.sh"
=======

alphabet = list(string.ascii_uppercase)
>>>>>>> backup-stable

def getFlagsParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script for large scale EvoPro runs"""

    parser = FileArgumentParser(description='get arguments to run script for large scale EvoPro runs',
                                fromfile_prefix_chars='@')

<<<<<<< HEAD
    #provide text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc.
    parser.add_argument('--sequence_files',
                        default='chainA.txt chainB.txt',
                        type=str,
                        help='Path to and name of sequence files in order of chains, separated by spaces.')
=======
    parser.add_argument('--sequence_file',
                        default='chains.txt',
                        type=str,
                        help='Path to and name of sequence file. Each line in the sequence file will get a set of directories. Default is chains.txt')
    
    #provide text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc.
    parser.add_argument('--sequence_files',
                        default=None,
                        type=str,
                        help='Path to and name of sequence files in order of chains, separated by spaces. Provide here text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc')
>>>>>>> backup-stable
    
    #provide af2.flags, evopro.flags, json.flags, run_evopro.sh, templates folder(optional), other pdbs(optional) in this directory
    parser.add_argument('--representative_directory',
                        default='./',
                        type=str,
                        help='Path to directory for files that are needed to run EvoPro, including af2.flags, evopro.flags etc.')
    
    parser.add_argument('--num_replicates',
                        default=3,
                        type=int,
                        help='number of independent trajectories to set up per sequence combination')
<<<<<<< HEAD
=======
    
    parser.add_argument('--custom_pdb_dir',
                        default=None,
                        type=str,
                        help='directory of unique pdb files to use in each directory, numbered to match number of pairs')
    parser.add_argument('--custom_pdb_prefix',
                        default="scaffold",
                        type=str,
                        help='prefix of unique pdb files to use in each directory, numbered to match number of pairs')
>>>>>>> backup-stable
    return parser

def parse_sequences(filename):
    sequences = []
    with open(filename, "r") as fil:
        for lin in fil:
            l = lin.strip()
            sequences.append(l)

    return sequences

<<<<<<< HEAD
=======
def parse_multimer_sequences(filename):
    sequences = []
    with open(filename, "r") as fil:
        for lin in fil:
            l = lin.strip().split(",")
            sequences.append(l)

    return sequences

>>>>>>> backup-stable
if __name__=="__main__":
    parser = getFlagsParser()
    args = parser.parse_args(sys.argv[1:])

    main_dir = os.getcwd()
<<<<<<< HEAD

    ip_files = args.sequence_files.strip().split(" ")
    seqs = []
    for fil in ip_files:
        seqs.append(parse_sequences(fil))
    
    seq_permutations = list(itertools.product(*seqs))

    num_pairs = len(seq_permutations)
    num_jobs =  num_pairs * args.num_replicates

    for perm, i in zip(seq_permutations, range(len(seq_permutations))):
=======
    
    if args.sequence_files:

        ip_files = args.sequence_files.strip().split(" ")
        ch_seqs = []
        seqs = []
        for fil in ip_files:
            ch_seqs.append(parse_sequences(fil))
        
        seq_permutations = list(itertools.product(*ch_seqs))
        for perm in seq_permutations:
            seq_split = []
            for elem in perm:
                seq_split.extend(elem.strip().split(","))
            seqs.append(seq_split)
        
    else:
        ip_files = args.sequence_file.strip().split(" ")
        seqs = []
        for fil in ip_files:
            seqs.append(parse_multimer_sequences(fil))
        
        seqs = seqs[0]

    num_pairs = len(seqs)
    num_jobs =  num_pairs * args.num_replicates
    
    custom_pdbs = False
    if args.custom_pdb_dir:
        custom_pdbs = True
        custom_pdbs_list = [os.path.join(main_dir, args.custom_pdb_dir, args.custom_pdb_prefix + "_" + str(i+1) + ".pdb") for i in range(num_pairs)]
        print(custom_pdbs_list)

    for perm, i in zip(seqs, range(num_pairs)):
>>>>>>> backup-stable
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
<<<<<<< HEAD

=======
            
>>>>>>> backup-stable
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
<<<<<<< HEAD
=======
            if custom_pdbs:
                shutil.copy(custom_pdbs_list[i], os.path.join(dir_name, run_name, args.custom_pdb_prefix + ".pdb"))
>>>>>>> backup-stable

            #generate json for each dir
            os.chdir(os.path.join(main_dir, dir_name, run_name))
            os.system("python /proj/kuhl_lab/evopro/evopro/run/generate_json.py @json.flags")
            os.chdir(main_dir)
<<<<<<< HEAD
            
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

    
=======

    print("Use script at run_largerun.py to run jobs.")

>>>>>>> backup-stable
