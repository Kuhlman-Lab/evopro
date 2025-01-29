#!/bin/bash

#SBATCH -p volta-gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH -t 20:00:00
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1

source ~/.bashrc
conda activate af2_mpnn
module load gcc
module load cuda

python /proj/kuhl_lab/evopro/evopro/run/generate_json_msd_henry.py @json.flags
cp ./pdbs/msd/msd.pdb .
python /proj/kuhl_lab/evopro/evopro/run/run_evopro_multistate_henry.py @evopro.flags
