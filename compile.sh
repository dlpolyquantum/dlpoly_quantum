#!/bin/bash
#SBATCH -N 1
##SBATCH -p shakib
#SBATCH --ntasks-per-node 24
#SBATCH --mem=60G
#SBATCH --mail-user=momeni@njit.edu
#SBATCH --mail-type=ALL

START_TIME=`date`
echo The calculation started at: ${START_TIME}...

module purge
module load intel impi
module list

cd $SLURM_SUBMIT_DIR/source
make clean

make dlpoly

