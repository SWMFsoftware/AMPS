#!/usr/bin/perl

use strict;
use warnings;

my @TestTable=("test_batl-reader.diff/test_batl-reader_check", 
"test_binary-reader.diff/test_binary-reader_check", 
"test_CCMC.diff/test_CCMC_check", 
"test_CCMC-Individual_Trajectories.diff/test_CCMC-Individual_Trajectories_check", 
"test_CG[?].diff/test_CG_check", 
"test_CG-BATSRUS-dust-coupling.diff/test_CG-BATSRUS-dust-coupling_check", 
"test_CG-dust-test\\[?\\].diff/test_CG-dust-test_check", 
"test_COPS.diff/test_COPS_check", 
"test_lightwave-current-on.diff/test_lightwave-current-on_check", 
"test_Sphere-drag-coefficient\\[?\\].diff/test_Sphere-drag-coefficient_check", 
"test_fast-wave.diff/test_fast-wave_check", 
"test_fast-wave-openmp.diff/test_fast-wave-openmp_check", 
"test_lightwave-current-off.diff/test_lightwave-current-off_check", 
"test_lightwave-current-off-openmp.diff/test_lightwave-current-off-openmp_check");  

my ($fname,$target,$pair,$size,$cmd);

foreach $pair (@TestTable) {
  ($fname,$target) = split /\//, $pair, 2;

  my @files=glob($fname);

  foreach $fname (@files) {
    $size = -s "$fname"; 
   
    if ($size!=0) {
      print "Found nont-zero diff file ($fname) - recheck it\n";
       
      if (-e "test_amps.check") {
        $cmd="make ".$target." >> test_amps.check";
      }
      else {
        $cmd="make ".$target." > test_amps.check";
      }

      print "$cmd\n";
     
      system($cmd); 
    }
  }
}

