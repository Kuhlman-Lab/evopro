#!/bin/bash

#SBATCH -p kuhlab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH -t 10:00:00
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=amritan@ad.unc.edu

source ~/.bashrc
conda activate af2_mpnn
module load gcc
module add cuda
python /proj/kuhl_lab/evopro/evopro/run/run_geneticalg_gpus.py @evopro.flags
