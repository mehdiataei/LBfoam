#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=20:00:00
#SBATCH --job-name bubble3D
#SBATCH --output=mpi_ex_%j.txt

cd $SLURM_SUBMIT_DIR

module load gcc/7.3.0
module load python/2.7.14
module load openmpi/3.1.0rc3
mpirun ./bucket3D bucket3D.xml

