#!/bin/bash 
#SBATCH --nodes=30
#SBATCH --ntasks=1200
#SBATCH --time=20:00:00
#SBATCH --job-name bubble
#SBATCH --output=mpi_ex_%j.txt

cd $SLURM_SUBMIT_DIR

module load gcc/7.3.0
module load python/2.7.14
module load openmpi/3.1.1
mpirun ./bucket2D bucket2D.xml

