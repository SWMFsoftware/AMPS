#!/bin/csh
# This script checks out the latest version of the AMPS,
# executes the AMPS tests
# This script is meant to be used by the AMPS developers.

# The script can be executed by simply typing
#
# ./run_test_amps.sh
#
# To run the AMPS tests perodically, you can use the crontab facility.
# Type 'crontab -e' to add a new entry to your crontab.
# Here is an example entry for nightly runs at 12:30 am:
#
# 30 0 * * * $HOME/bin/run_test_amps.sh

# NOTE:
# the body of the script is divided into blocks of the format
#  #>BlockName ##############
#  #command_1               #
#  # ...                    #
#  #command_last           <#
# certain blocks will be uncommented at the test installation,

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

#Set the working directory
set WorkDir = /glade2/scratch2/vtenishe  

#Go to your home directory
cd $WorkDir

#Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

#Remove the previous test directory if necessary
rm -rf AMPS */AMPS

#Checkout the latest code version
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS 
cd AMPS
cvs co -D "`date +%m/%d/%Y` 23:20" AMPS_data 
cd ..

#Undata model data files 
rsync -r amps@tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST/ $WorkDir/AMPS_DATA_TEST

#Create separate folders for different compilers
rm -rf GNU
mkdir -p GNU;   cp -r AMPS GNU/;

rm -rf Intel
mkdir -p Intel; cp -r AMPS Intel/;

rm -rf PGI
mkdir -p PGI;   cp -r AMPS PGI/;

#Install AMPS
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                         
./Config.pl -install -compiler=gfortran,gcc_mpicc -f-link-option=-lmpi_cxx   >& test_amps.log    

cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                       
./Config.pl -f-link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >& test_amps.log

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                         
./Config.pl -f-link-option=-lc++ -install -compiler=pgf90,pgccmpicxx      >& test_amps.log    

#Compile the tests and waite completion of the compilation procedure
cd $WorkDir/Tmp_AMPS_test
rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete

AMPS/utility/TestScripts/Yellowstone/CompileGNUYellowstone.sh & 
AMPS/utility/TestScripts/Yellowstone/CompilePGIYellowstone.sh &
AMPS/utility/TestScripts/Yellowstone/CompileIntelYellowstone.sh &

while ((! -f AmpsCompilingIntelComplete) || (! -f AmpsCompilingGNUComplete) || (! -f AmpsCompilingPGIComplete))
  sleep 60
end

rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete

# Execute the tests 
rm -f AmpsTestComplete

foreach job (test_amps.yellowstone.all.*job)                
  #execute the next part of the tests
  bsub -K -J "AMPS_test[1-3]" < $job
  wait

  while ((! -f AmpsTestGNUComplete) || (! -f AmpsTestIntelComplete) || (! -f AmpsTestPGIComplete))
    sleep 60
  end

  rm -rf AmpsTestGNUComplete
  rm -rf AmpsTestIntelComplete
  rm -rf AmpsTestPGIComplete
end







