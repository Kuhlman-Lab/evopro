import os
import sys
import shutil
import itertools
import string
from evopro.utils.inputs import FileArgumentParser

alphabet = list(string.ascii_uppercase)

def getFlagsParser() -> FileArgumentParser:
    """Gets an FileArgumentParser with necessary arguments to run script for large scale EvoPro runs"""

    parser = FileArgumentParser(description='get arguments to run script for large scale EvoPro runs',
                                fromfile_prefix_chars='@')

    parser.add_argument('--sequence_file',
                        default='chains.txt',
                        type=str,
                        help='Path to and name of sequence file. Each line in the sequence file will get a set of directories. Default is chains.txt')
    
    #provide text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc.
    parser.add_argument('--sequence_files',
                        default=None,
                        type=str,
                        help='Path to and name of sequence files in order of chains, separated by spaces. Provide here text files with options for each chain eg.--sequence_files chainA.txt chainB.txt etc')
    
    #provide af2.flags, evopro.flags, json.flags, run_evopro.sh, templates folder(optional), other pdbs(optional) in this directory
    parser.add_argument('--representative_directory',
                        default='./',
                        type=str,
                        help='Path to directory for files that are needed to run EvoPro, including af2.flags, evopro.yaml, mpnn.yaml etc.')
    
    parser.add_argument('--num_replicates',
                        default=3,
                        type=int,
                        help='number of independent trajectories to set up per sequence combination')
    
    parser.add_argument('--custom_pdb_dir',
                        default=None,
                        type=str,
                        help='directory of unique pdb files to use in each directory, numbered to match number of pairs')
    
    parser.add_argument('--custom_pdb_prefix',
                        default="scaffold",
                        type=str,
                        help='prefix of unique pdb files to use in each directory, numbered to match number of pairs')
    parser.add_argument('--evopro_path',
                        default="/proj/kuhl_lab/evopro2/",
                        type=str,
                        help='path to the evopro direcotory containing the run/generate_json.py script')
    return parser

def parse_sequences(filename):
    sequences = []
    with open(filename, "r") as fil:
        for lin in fil:
            l = lin.strip()
            sequences.append(l)

    return sequences

def parse_multimer_sequences(filename):
    sequences = []
    with open(filename, "r") as fil:
        for lin in fil:
            l = lin.strip().split(",")
            sequences.append(l)

    return sequences

if __name__=="__main__":
    parser = getFlagsParser()
    args = parser.parse_args(sys.argv[1:])

    main_dir = os.getcwd()
    
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
            if custom_pdbs:
                shutil.copy(custom_pdbs_list[i], os.path.join(dir_name, run_name, args.custom_pdb_prefix + ".pdb"))

            #generate json for each dir
            os.chdir(os.path.join(main_dir, dir_name, run_name))
            os.system("python " + str(args.evopro_path) + "run/generate_json.py @json.flags")
            os.chdir(main_dir)