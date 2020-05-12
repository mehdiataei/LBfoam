#!/bin/bash 
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --time=20:00:00
#SBATCH --job-name bubble3D
#SBATCH --output=mpi_ex_%j.txt

cd $SLURM_SUBMIT_DIR

module load gcc/7.3.0
module load python/3.6.5
module load openmpi/3.1.1
mpirun ./bucket3D bucket3D.xml

