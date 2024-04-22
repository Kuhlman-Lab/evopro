#!/bin/bash

#SBATCH -p kuhlab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=40g
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=amritan@ad.unc.edu

source ~/.bashrc
conda activate rf2_mpnn
module load gcc
module load cuda
python /proj/kuhl_lab/evopro/evopro/run/run_evopro_binder_rf2.py @evopro.flags
