"""Script that consolidates multiple directories of PDB files into one.
The PDB files are renamed using prefixes of the directory names they came from."""

import os
import argparse
import shutil

def copy_file(oldfile, newfile):
    with open(oldfile, "r") as f:
        lines = f.readlines()

    with open(newfile, "w") as f:
        for line in lines:
            f.write(line)

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent_dir", type=str, default='./')
    parser.add_argument("--new_directory_name", type=str, default='all_pdbs')
    parser.add_argument("--prefix", type=str, default='pair')
    parser.add_argument("--file_prefix", type=str, default='outputs/')
    parser.add_argument("--copy_pbz", action='store_true', help='Default is False.')
    args = parser.parse_args()
    
    parent_dir = args.parent_dir
    newdirname = args.new_directory_name
    prefix = args.prefix
    file_prefix = args.file_prefix
    
    onlydirs = [os.path.join(parent_dir, x) for x in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, x)) and x.startswith(prefix)]

    newdir = os.path.join(parent_dir, newdirname)
    if not os.path.exists(newdir):
        os.mkdir(newdir)
        print("Created directory " + newdir)
    
    for dirname in onlydirs:
        subdirs = [os.path.join(dirname, x) for x in os.listdir(dirname) if os.path.isdir(os.path.join(dirname, x))]
        #print(subdirs)
        if len(subdirs) == 0:
            subdirs = [dirname]
            print("No subdirectories found in " + dirname + ". Using " + dirname + " as the subdirectory.")
        for subdir in subdirs:
            print(subdir)
            files = []
            try:
                files = [os.path.join(subdir, file_prefix, x) for x in os.listdir(os.path.join(subdir, file_prefix)) if x.endswith(".pdb")]
                if args.copy_pbz:
                    files+=[os.path.join(subdir, file_prefix, x) for x in os.listdir(os.path.join(subdir, file_prefix)) if x.endswith("result_0.pbz2")]
                print(files)
                for oldfile in files:
                    #os.path.basename(os.path.normpath(dirname))
                    #print(os.path.basename(os.path.normpath(dirname)))
                    #print(os.path.basename(os.path.normpath(subdir)))
                    print(oldfile)
                    print(os.path.basename(oldfile))
                    newfilename = os.path.basename(os.path.normpath(dirname)) + "_" + os.path.basename(os.path.normpath(subdir)) + "_" + os.path.basename(oldfile)
                    print("here", newfilename)
                    newfile = os.path.join(newdir, newfilename)
                    print(oldfile, newfile)
                    shutil.copyfile(oldfile, newfile)
                    #copy_file(oldfile, newfile)
            except:
                print("No pdb or pbz files found in " + subdir)