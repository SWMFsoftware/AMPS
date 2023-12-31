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
source /usr/share/modules/init/csh  

# source shell run commands (ONLY csh AND tcsh SHELLS ARE USED FOR NOW)
# to set CVSROOT or CVS_RSH variables
# note: it is better to have these variables set in the beginning of rc file
source $HOME/.cshrc

# set the working directory
set WorkDir = $HOME  

#>Pleiades ############################
set WorkDir = /nobackup/`whoami`    

#>Yellowstone ###########
#set WorkDir =         <#

#>Stampede #############
# set WorkDir = $WORK <#

# Go to your home directory
cd $WorkDir

# Create a temporary directory for the tests
mkdir -p Tmp_AMPS_test
cd Tmp_AMPS_test

# Remove previous job scripts
rm -f *.job

# Remove the previous test directory if necessary
rm -rf */AMPS

# Checkout the latest code version
set CheckoutTime = `date`

if (-e AMPS) then
  cd AMPS_Legacy; git pull
  cd ../BATL; git pull
  cd ../AMPS; git pull
  cd SWMF_data; git pull

  cd ../../
else
  gitclone AMPS_Legacy
  gitclone AMPS
  gitclone BATL

  cd AMPS
  gitclone SWMF_data

  cd ../
endif


# Update data files for test at supercomputers
#>Pleiades>Yellowstone>Stampede #########################

#rsync -P -a -c -r amps@tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST/ $WorkDir/AMPS_DATA_TEST 

cd $WorkDir 

if (-e AMPS_DATA_TEST) then   
  cd AMPS_DATA_TEST
  git pull
else 
  git clone tower-left.engin.umich.edu:/Volumes/Data01/AMPS_DATA_TEST
endif

cd $WorkDir/Tmp_AMPS_test

#cp /home3/vtenishe/Table /nobackupp8/vtenishe/Tmp_AMPS_test/AMPS/MakefileTest

# create separate folders for different compilers
#>GNUAll ############################
rm -rf GNU
mkdir -p GNU;   cp -r AMPS GNU/;  
cp -r BATL GNU/AMPS/
#>IntelAll ##########################
rm -rf Intel
mkdir -p Intel; cp -r AMPS Intel/;
cp -r BATL Intel/AMPS/
#>PGIAll ############################
rm -rf PGI
mkdir -p PGI;   cp -r AMPS PGI/;  
cp -r BATL PGI/AMPS/

# install AMPS
#>GNUAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                                        #
echo AMPS was checked out on $CheckoutTime > test_amps.log
./Config.pl -install -compiler=gfortran,gcc_mpicc    >>& test_amps.log    
./utility/TestScripts/BuildTest.pl -test-run-time=60

#>IntelAll #################################################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                      #
echo AMPS was checked out on $CheckoutTime > test_amps.log
./Config.pl -install -compiler=ifortmpif90,iccmpicxx >>& test_amps.log    
./utility/TestScripts/BuildTest.pl -test-run-time=60

#>PGIAll ###################################################################
cd $WorkDir/Tmp_AMPS_test/PGI/AMPS                                        #
echo AMPS was checked out on $CheckoutTime > test_amps.log
./Config.pl -install -compiler=pgf90,pgccmpicxx -cpp-compiler=pgc++ -link-option=-L/nasa/sgi/mpt/2.15r20/lib,-lmpi++,-lmpi       >>& test_amps.log    
./utility/TestScripts/BuildTest.pl -test-run-time=60

# copy job files to the AMPS directory on supercomputers
#>Pleiades ###############################################
cp utility/TestScripts/test_amps.pleiades.*.job ../..
#>Stampede ###############################################
#cp AMPS/utility/TestScripts/test_amps.stampede.*.job . <#


# Exeptions
echo -n "Set Exeptions....."

#>Valeriy ##################################################################
#cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                                                                      #
#./Config.pl -link-option=-lc++ -install -compiler=ifort,iccmpicxx -link-option=-cxxlib >>& test_amps.log <# 

#>Pleiades ##############################################
cd $WorkDir/Tmp_AMPS_test/Intel/AMPS                   #
./Config.pl -cpp-link-option=-lmpi,-lmpi++,-lpthread 
./Config.pl -f-link-option=-lmpi,-lmpi++                     #
./Config.pl -cpplib-rm=-lmpi_cxx
./Config.pl -noopenmp
./Config.pl -cpp-compiler=icpc

cd $WorkDir/Tmp_AMPS_test/GNU/AMPS                     #
./Config.pl -cpplib-rm=-lmpi_cxx
./Config.pl -f-link-option=-lmpi++                  

cd $WorkDir/Tmp_AMPS_test/PGI/AMPS
./Config.pl -f-link-option=-pgf90libs,-pgc++libs

cd $WorkDir/Tmp_AMPS_test
rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete
rm -f AmpsTestComplete

echo " done."

# compile AMPS tests
$HOME/bin/CompileGNUPleiades.sh &
$HOME/bin/CompileIntelPleiades.sh &
$HOME/bin/CompilePGIPleiades.sh &

# Run test
# Super computers
#>Pleiades ####################################
cd $WorkDir/Tmp_AMPS_test

#waite untill all compilation is complete
while ((! -f AmpsCompilingIntelComplete) || (! -f AmpsCompilingGNUComplete) || (! -f AmpsCompilingPGIComplete)) 
  sleep 60
end

rm -f AmpsCompilingIntelComplete
rm -f AmpsCompilingGNUComplete
rm -f AmpsCompilingPGIComplete

echo Compiling of AMPS is completed

#########################################
#Execute the tests 
source $HOME/.cshrc
#$HOME/bin/RunAllPleiades.sh
#perl $HOME/bin/RunAllPleiades.pl

# Run test
foreach job (test_amps.pleiades.all*.job) #
   qsub $job
end


