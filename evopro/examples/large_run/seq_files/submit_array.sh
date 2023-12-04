#!/bin/bash
# submit_array.sh

#SBATCH -p kuhlab
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 6-00:00:00
#SBATCH --mem=32g
#SBATCH --constraint=rhel8
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=onyen@unc.edu
#SBATCH --array=1-15

N_PAIRS=3
N_REPS=5
cd $SLURM_SUBMIT_DIR/pair(($SLURM_ARRAY_TASK_ID%$N_PAIRS)+1)/run(($SLURM_ARRAY_TASK_ID%$N_REPS)+1)
./run_evopro.sh
