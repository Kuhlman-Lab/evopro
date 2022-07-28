#!/bin/bash

#SBATCH -p volta-gpu
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=16g
#SBATCH -t 02-00:00:00
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:2
#SBATCH --mail-type=end
#SBATCH --mail-user=amritan@ad.unc.edu

source ~/.bashrc
conda activate af2
module add cuda
python /pine/scr/a/m/amritan/kuhlmanlab/folddesign/folddesign/genetic_alg/run_geneticalg_gpus.py @fdd_flags.txt
