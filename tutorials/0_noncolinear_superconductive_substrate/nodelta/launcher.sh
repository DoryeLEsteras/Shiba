#!/bin/bash -l
#SBATCH --output=task.out.%j
#SBATCH --error=task.err.%j
#SBATCH -D ./

# Job Name
#SBATCH --job-name=tutncnodelta
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --partition=general
#SBATCH --ntasks-per-node=72
#SBATCH --time=0-07:29:00

#Libraries and modules
module load intel/21.5.0
module load  impi/2021.5

wstm.py -input input
