#!/bin/csh

set WorkDir = $HOME
source $WorkDir/module/intel 

echo -n "Executing tests Intel....."
make TESTMPIRUN4="mpirun -np 4"  MPIRUN="mpirun -np 8" TESTMPIRUN1="mpirun -np 1" test_run >>& test_amps.log

echo " done."

cd $WorkDir/Tmp_AMPS_test
echo Done > AmpsTestIntelComplete
