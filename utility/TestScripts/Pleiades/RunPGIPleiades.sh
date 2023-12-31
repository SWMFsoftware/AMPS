#!/bin/csh

cd PGI/AMPS
module purge
module load comp-pgi/19.10 
module load mpi-sgi/mpt
module load mpi-hpe/mpt.2.17r13
module load tecplot/2017r2
make $1 MPIRUN="mpiexec -n 8 dplace -s1 -c 8-15" >>& test_amps.log

cd ../..
echo Done > AmpsTestPGIComplete 
