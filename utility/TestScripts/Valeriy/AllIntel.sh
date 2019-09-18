#!/bin/csh

set WorkDir = $HOME
source $WorkDir/module/intel 

echo -n "Compiling Intel....."

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS; 
make test_compile >>& test_amps.log

echo " done."

echo -n "Executing tests Intel....."
make TESTMPIRUN4="mpirun -np 4" MPIRUN="export DYLD_LIBRARY_PATH=/Users/ccmc/boost/lib:/opt/intel/compilers_and_libraries/mac/lib;mpirun -np 4" TESTMPIRUN1="export DYLD_LIBRARY_PATH=/Users/ccmc/boost/lib:/opt/intel/compilers_and_libraries/mac/lib;mpirun -np 1" test_run >>& test_amps.log
echo " done."
