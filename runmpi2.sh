#!/bin/bash
# cd ./rank0/
# rm -r *png
# echo "Clearing ./rank0/ pictures"
# cd ../rank1/
# rm -r *png
# echo "Clearing ./rank1/ pictures"
# cd ../
# echo "Running script"
mpiexec -n 2 python mpi_blocking_2.py
# mpiexec -n 2 python mpietc.py
