#!/bin/bash

#SBATCH -p volta-gpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=16g
#SBATCH -t 08-00:00:00
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:2
#SBATCH --mail-type=end
#SBATCH --mail-user=amritan@ad.unc.edu

source ~/.bashrc
conda activate af2
module add cuda
python run_geneticalg_gpus.py @fdd_flags_testing.txt
