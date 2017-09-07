#!/bin/csh

set WorkDir = $HOME

source $WorkDir/Tmp_AMPS_test/AMPS/utility/TestScripts/CompilerSetup/set_mpi_intel.valeriy
echo -n "Compiling Intel....."

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS; 
make test_compile >>& test_amps.log

echo " done."

echo -n "Executing tests Intel....."
make TESTMPIRUN4="mpirun -np 4" MPIRUN="mpirun -np 4" TESTMPIRUN1="export DYLD_LIBRARY_PATH=/Users/ccmc/boost/lib:/opt/intel/compilers_and_libraries/mac/lib;mpirun -np 1" test_run >>& test_amps.log
echo " done."
