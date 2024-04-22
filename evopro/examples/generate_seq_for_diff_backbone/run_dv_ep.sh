#!/bin/bash

#SBATCH -J 147
#SBATCH -p volta-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 5-
#SBATCH --mem=32g
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=gprida@unc.edu

source ~/.bashrc
conda activate af2_mpnn
module load gcc
module load cuda
python /proj/kuhl_lab/evopro/evopro/run/run_geneticalg_gpus.py @evopro.flags
