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

# init command to load modules      
#>Pleiades ############################
#source /usr/share/modules/init/csh  <#

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

# set the working directory
set WorkDir = $HOME  

#>Pleiades ############################
#set WorkDir = /nobackup/`whoami`    <#

#>Yellowstone ###########
#set WorkDir =         <#

#>Stampede #############
 set WorkDir = $SCRATCH

# Go to your home directory
cd $WorkDir

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

# Remove the previous test directory if necessary
rm -rf */AMPS

# Checkout the latest code version
set CheckoutTime = `date`

# Checkout the latest code version
if (-e AMPS) then
  cd AMPS_Legacy; git pull
  cd ../BATL; git pull
  cd ../AMPS; git pull
  cd SWMF_data; git pull

  cd $WorkDir/Tmp_AMPS_test 
else
  gitclone AMPS
  gitclone AMPS_Legacy
  gitclone BATL

  cd AMPS
  gitclone SWMF_data

  cd ../
endif 

# create separate folders for different compilers
#>GNUAll ############################
mkdir -p GNU;   cp -r AMPS GNU/;  
cp -r BATL GNU/AMPS/

#>IntelAll ##########################
mkdir -p Intel; cp -r AMPS Intel/;
cp -r BATL Intel/AMPS/
#>PGIAll ############################
#mkdir -p PGI;   cp -r AMPS PGI/;  <#

# install AMPS
#>GNUAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                        #
echo AMPS was checked out on $CheckoutTime > test_amps.log
./Config.pl -install -compiler=gfortran,gcc_mpicc    >>& test_amps.log    
utility/TestScripts/BuildTest.pl -test-run-time=30

#>IntelAll #################################################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
echo AMPS was checked out on $CheckoutTime > test_amps.log
./Config.pl -install -compiler=ifortmpif90,iccmpicxx >>& test_amps.log    
utility/TestScripts/BuildTest.pl -test-run-time=30

#>Valeriy ##################################################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
#./Config.pl -link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >>& test_amps.log<#
#>PGIAll ###################################################################
#cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                        #
#./Config.pl -install -compiler=pgf90,pgccmpicxx      >& test_amps.log    <#

# Update data files for test at supercomputers
#>Pleiades>Yellowstone>Stampede #########################
#rsync -r amps@tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST/ $WorkDir/AMPS_DATA_TEST


cd $WorkDir

if (-e AMPS_DATA_TEST) then
  cd AMPS_DATA_TEST
  git pull
else
  git clone amps@tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST
endif

# copy job files to the AMPS directory on supercomputers
# #>Pleiades ###############################################
# #cp AMPS/utility/TestScripts/test_amps.pleiades.*.job . <#
# #>Stampede ###############################################
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS 
rm -f utility/TestScripts/test_amps.stampede.all.overtime*.job
utility/TestScripts/BuildTest.pl -test-run-time=30  

rm -rf ../../test_amps.stampede.*.job
cp utility/TestScripts/test_amps.stampede.*.job ../.. 

# Compile AMPS tests
cd $WorkDir/Tmp_AMPS_test
rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete

AMPS/utility/TestScripts/Stampede/CompileGNUPStampede.sh &
AMPS/utility/TestScripts/Stampede/CompileIntelStampede.sh &

#waite untill all compilation is complete
while ((! -f AmpsCompilingIntelComplete) || (! -f AmpsCompilingGNUComplete) )
  sleep 60
end

rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete

echo Compiling of AMPS is completed

# Run test
#>Stampede ####################################
foreach job (test_amps.stampede.all*.job) # 
  rm -rf AmpsTestComplete
  @ iTest=0

  sbatch $job

# while ($iTest<2) 
#   sbatch $job

#   echo $job is submitted 
#   set JobCompletedFlag='false'

#   while ($JobCompletedFlag == 'false')
#     sleep 60

#     if (`squeue -u vtenishe | grep AMPS_tes` == "") then
#       set JobCompletedFlag='true'
#     endif
#   end

#   echo $job is completed
#   sleep 120

#   if (-e AmpsTestDone) then 
#     break
#   endif 

#   @ iTest++ 
# end
end

echo The test routine is completed

