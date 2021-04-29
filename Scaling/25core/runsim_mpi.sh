#!/bin/bash
#SBATCH --job-name=mpi-25cores-sim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=03:00:00

module load StdEnv/2020
module load gcc/9.3.0
module load mpi4py
module load scipy-stack/2020b
module load python/3.8.2

echo "Running MPI job"
mpiexec -n 25 python mpi_blocking_many.py
echo "Done."
