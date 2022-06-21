#!/bin/bash/
#SBATCH -p serial 
#SBATCH -q normal 
#SBATCH -t 0-01:00:00
module load ambertools/15
pdb4amber -i [file], -d > [outputfile] 
