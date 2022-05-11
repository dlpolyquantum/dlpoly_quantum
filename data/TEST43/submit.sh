#!/bin/bash -I
#SBATCH --job-name=test
#SBATCH --output=myoutput%j.out
#SBATCH --error=jobError.out
#SBATCH --time=30-00:00:00
#SBATCH -p shakib
#SBATCH -N 2
#SBATCH --ntasks-per-node 24
#SBATCH --mem=20G
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=dkl27@njit.edu

module purge
module load intel impi
module list

echo JOB STARTED AT:
 date

 mpirun -n 48 dlpoly.x

echo JOB FINISHED AT:
 date
echo ALL DONE!!!

