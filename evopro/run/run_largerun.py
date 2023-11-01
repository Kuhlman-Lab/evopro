import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--main_dir", type=str, default=None)
    parser.add_argument("--dir_prefix", type=str, default='pair')
    parser.add_argument("--num_dirs", type=str, default='1-10')
    parser.add_argument("--subdir_prefix", type=str, default='run')
    parser.add_argument("--num_subdirs", type=int, default=5)
    parser.add_argument("--run_file_name", type=str, default="run_evopro_rf2.sh")

    args = parser.parse_args()
    
    if not args.main_dir:
        main_dir = os.getcwd()
    
    numdirs_temp = args.num_dirs.strip().split(",")
    numdirs_temp = [x.strip() for x in numdirs_temp if x]
    numdirs = []
    for elem in numdirs_temp:
        if "-" not in elem:
            numdirs.append(elem)
        else:
            start, finish = elem.split("-")
            s = int(start)
            f = int(finish)
            for i in range(s, f+1):
                numdirs.append(i)
                
    for numdir in numdirs:
        for i in range(args.num_subdirs):
            os.chdir(os.path.join(main_dir, args.dir_prefix + str(numdir), args.subdir_prefix + str(i+1)))
            print(args.dir_prefix + str(numdir), args.subdir_prefix + str(i+1))
            os.system("sbatch " + args.run_file_name)