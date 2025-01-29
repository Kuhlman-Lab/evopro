#!/bin/bash

#SBATCH -J binder-example
#SBATCH -p a100-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 6-00:00:00
#SBATCH --mem=80g
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=amritan@ad.unc.edu

source ~/.bashrc
conda activate af2_mpnn
module load gcc
module load cuda
python /proj/kuhl_lab/evopro/evopro/run/run_evopro_binder_dev.py @evopro.flags
