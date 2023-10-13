#!/bin/bash -l
#SBATCH --output=task.out.%j
#SBATCH --error=task.err.%j
#SBATCH -D ./

# Job Name
#SBATCH --job-name=tut3nodelta
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --partition=general
#SBATCH --ntasks-per-node=72
#SBATCH --time=0-06:59:00

#Libraries and modules
module load intel/21.5.0
module load  impi/2021.5
module load mkl/2022.2
module load qe/6.8
module load wannier90/3.1.0
 module load anaconda/3/2023.03

wstm.py -input input
