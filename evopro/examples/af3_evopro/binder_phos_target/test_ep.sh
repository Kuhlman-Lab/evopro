#!/bin/bash
  
#SBATCH -p kuhlab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=80g
#SBATCH -t 04:00:00
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1

source ~/.bashrc
conda activate evopro2_af3
module load gcc
module load cuda
#python /proj/kuhl_lab/evopro2/run/generate_json.py @json.flags
python /proj/kuhl_lab/evopro2/run/run_evopro.py --config_file evopro_basic.yaml
