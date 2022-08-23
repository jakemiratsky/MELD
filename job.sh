#!/bin/bash 
#SBATCH -N 16
#SBATCH -p asinghargpu1
#SBATCH -q asinghargpu1
#SBATCH --gres=gpu:1
#SBATCH -n 16 #This will only assign one core per task... may need to use c... try -c 8 next run
#SBATCH -o meld.log

# Other parameters for slurm may be added 
# For MELD, the number of tasks much match the number of allocated GPUs!

module purge
module load anaconda/py3
source activate meld 

python setup_MELD.py

mpirun launch_remd 
