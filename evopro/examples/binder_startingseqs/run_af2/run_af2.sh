#!/bin/bash

#SBATCH -J af2_aid5
#SBATCH -p volta-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH -t 09:00:00
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1

source ~/.bashrc
module load gcc
module load cuda
conda activate af2
python /proj/kuhl_lab/alphafold/run/run_af2.py @af2.flags
