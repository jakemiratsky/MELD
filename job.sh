#!/bin/bash 
#SBATCH -p gpu
#SBATCH -q wildfire 
#SBATCH --gres=gpu:4
#SBATCH -n 4
#SBATCH -o meld.log

# Other parameters for slurm may be added 
# For MELD, the number of tasks much match the number of allocated GPUs!

module purge
module load anaconda/py3
source activate meld 

python setup_MELD.py

mpirun launch_remd 
