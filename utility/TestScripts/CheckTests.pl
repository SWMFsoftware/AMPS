#!/usr/bin/perl

use strict;
use warnings;

my @TestTable=("test_batl-reader.diff/test_batl-reader_check", 
"test_binary-reader.diff/test_binary-reader_check", 
"test_CCMC.diff/test_CCMC_check", 
"test_CCMC-Individual_Trajectories.diff/test_CCMC-Individual_Trajectories_check", 
"test_CG\\[?\\].diff/test_CG_check", 
"test_CG-BATSRUS-dust-coupling.diff/test_CG-BATSRUS-dust-coupling_check", 
"test_CG-dust-test\\[?\\].diff/test_CG-dust-test_check", 
"test_COPS.diff/test_COPS_check", 
"test_lightwave-current-on.diff/test_lightwave-current-on_check", 
"test_Sphere-drag-coefficient\\[?\\].diff/test_Sphere-drag-coefficient_check", 
"test_fast-wave.diff/test_fast-wave_check", 
"test_fast-wave-openmp.diff/test_fast-wave-openmp_check", 
"test_lightwave-current-off.diff/test_lightwave-current-off_check", 
"test_lightwave-current-off-openmp.diff/test_lightwave-current-off-openmp_check",
"test_photon-test.diff/test_photon-test_check",
"test_X-Wing.diff/test_X-Wing_check",
"test_Moon.diff/test_Moon_check",
"test_Moon-restart.diff/test_Moon-restart_check",
"test_Individual_Trajectories--Boris\\[?\\].diff/test_Individual_Trajectories--Boris_check",
"test_Individual_Trajectories--Boris-relativistic\\[?\\].diff/test_Individual_Trajectories--Boris-relativistic_check",
"test_Individual_Trajectories--Markidis2010\\[?\\].diff/test_Individual_Trajectories--Markidis2010_check",
"test_Individual_Trajectories--guiding-center\\[?\\].diff/test_Individual_Trajectories--guiding-center_check",
"test_Drag-Coefficient--Specular-Reflection\\[?\\].diff/test_Drag-Coefficient--Specular-Reflection_check",
"test_Mars-Hot-Carbon.diff/test_Mars-Hot-Carbon_check", 
"test_Earth-Impulse-Source.diff/test_Earth-Impulse-Source_check",  
"test_InterpolateAMR-AMPS-CenterBased-test.diff/test_InterpolateAMR-AMPS-test_check", 
"test_InterpolateAMR-AMPS-CornerBased-test.diff/test_InterpolateAMR-AMPS-test_check", 
"test_InterpolateAMR-AMPS-test/test_InterpolateAMR-AMPS-test_check",
"test_Periodic-BC-NonUniform.diff/test_Periodic-BC-NonUniform_check",
"test_Mars-ions.diff/test_Mars-ions_check"); 

my ($fname,$target,$pair,$size,$cmd);

foreach $pair (@TestTable) {
  ($fname,$target) = split /\//, $pair, 2;

  my @files=glob($fname);

  foreach $fname (@files) {
    $size = -s "$fname"; 
   
    if (defined $size) {
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
}

