#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

module load StdEnv/2020
module load python/3.8.2

echo "Running MPI job"
python sim_algs_parallel_many.py
echo "Done."
