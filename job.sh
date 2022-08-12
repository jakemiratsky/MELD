#!/bin/bash 
#SBATCH -N 16
#SBATCH -p asinghargpu1
#SBATCH -q wildfire 
#SBATCH --gres=gpu:1
#SBATCH -n 16
#SBATCH -o meld.log

# Other parameters for slurm may be added 
# For MELD, the number of tasks much match the number of allocated GPUs!

module purge
module load anaconda/py3
source activate meld 

python setup_MELD.py

mpirun launch_remd 
