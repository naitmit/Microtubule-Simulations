#!/bin/bash
# cd ./rank0/
# rm -r *png
# echo "Clearing ./rank0/ pictures"
# cd ../rank1/
# rm -r *png
# echo "Clearing ./rank1/ pictures"
# cd ../
# echo "Running script"
mpiexec -n 4 python mpi_blocking_many.py
# mpiexec -n 2 python mpietc.py
