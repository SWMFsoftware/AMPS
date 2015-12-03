#!/usr/bin/perl
#$Id$
#reads the "Exosphere" section of the input file and modify the exosphere model of the PIC code

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);

use ampsConfigLib;


#my $command= $];
#print "Perl version : ".$command;



#arguments: 
#$ARGV[0] -> the name of the intput file 
#$ARGV[1] -> the working directory


#print "Process the exospehre model input:\n";


my $InputFileName=$ARGV[0]; #"cg.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName=$InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/; #"cg.input.Assembled.Species";
my $WorkingSourceDirectory="srcTemp";  #$ARGV[1];

$ampsConfigLib::WorkingSourceDirectory=$WorkingSourceDirectory;

my $line;
my $LineOriginal;

my $InputFileLineNumber=0;
my $FileName;
my $InputLine;
my $InputComment;
my $s0;
my $s1;
my $s2;
my $spec;


#get the list of the species                                                   
my $TotalSpeciesNumber;
my @SpeciesList;

open (SPECIES,"<$SpeciesFileName") || die "Cannot open file $SpeciesFileName\\
n";
$TotalSpeciesNumber=<SPECIES>;
@SpeciesList=<SPECIES>;
close (SPECIES);

chomp($TotalSpeciesNumber);
foreach (@SpeciesList) {
    chomp($_);
}


my @BjornSourceRate=(0)x$TotalSpeciesNumber;
my @UniformSourceRate=(0)x$TotalSpeciesNumber;
my @JetSourceRate=(0)x$TotalSpeciesNumber;

#add the model header to the pic.h
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
print PIC_H "#include \"Comet.h\"\n";
close (PIC_H);

open (InputFile,"<","$InputFileName") || die "Cannot find file \"$InputFileName\"\n";

while ($line=<InputFile>) {
  ($InputFileLineNumber,$FileName)=split(' ',$line);
  $line=<InputFile>;
  
  $LineOriginal=$line;
  chomp($line);
  
  
  ($InputLine,$InputComment)=split('!',$line,2);
  $InputLine=uc($InputLine);
  chomp($InputLine);
 
  $InputLine=~s/[=():]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  if ($InputLine eq "SODIUMSTICKINGPROBABILITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
     
    if ($InputLine eq "CONST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_","_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_","main/Comet.cpp");
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_","$InputLine","main/Comet.cpp");
    } 
    elsif ($InputLine eq "YAKSHINSKIY2005SS") {
      ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_","_EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_","main/Comet.cpp");
    }
    else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  elsif ($InputLine eq "SODIUMREEMISSIONFRACTION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    
    ampsConfigLib::RedefineMacro("_EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_","$InputLine","main/Comet.cpp");
  }

  #change the starting locations at which tracking of the partucle trajecoties begis
  elsif ($InputLine eq "TRACINGSURFACERADIUS") { #TracingSurfaceRadius
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
     ampsConfigLib::ChangeValueOfVariable("const double TracingSurfaceRadius",$InputLine,"main/Comet.h"); 
  }
  
  #the range of the dust grain speed (used for evaluation of the local time step)
  elsif ($InputLine eq "DUSTVELOCITYLIMIT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "MIN") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity",$InputLine,"main/Comet.cpp");   
    }
    elsif ($InputLine eq "MAX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity",$InputLine,"main/Comet.cpp");   
    }     
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }     
  }
  
  #redefine a value of a macro 
  elsif ($InputLine eq "DEFINE") {
    my ($macro,$value,$s0,$s1);
    
    ($InputLine,$InputComment)=split('!',$line,2);
    ($s0,$macro,$value,$s1)=split(' ',$InputLine,4);
    
    $s0=$macro;    
    $s0=~s/[()=]/ /g;
    ($s0,$s1)=split(' ',$s0,2);

    if (!defined $value) {
      $value=" ";
    }

    ampsConfigLib::AddLine2File("\n#undef $s0\n#define $macro $value\n","main/Comet.dfn");    
  }
  
  #define the mesh sigrature that will be used in the simuation (for generation of the particular mesh settings and loadging particular model data files)
  elsif ($InputLine eq "MESHSIGNATURE") {
    $line=~s/[=():]/ /g;
    ($InputLine,$line)=split(' ',$line,2);
    ($InputLine,$line)=split(' ',$line,2);
    
    ampsConfigLib::ChangeValueOfVariable("char Comet::Mesh::sign\\[_MAX_STRING_LENGTH_PIC_\\]","\"".$InputLine."\"","main/Comet.cpp");
  }
  
  #mean density of the dust grains (DustMeanDensity)
  elsif ($InputLine eq "DUSTMEANDENSITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2); 

    ampsConfigLib::ChangeValueOfVariable("double ElectricallyChargedDust::MeanDustDensity",$InputLine,"models/dust/Dust.cpp");
  } 

  #set the drag coefficient (GrainDragCoefficient)
  elsif ($InputLine eq "DUSTDRAGCOEFFICIENT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

    ampsConfigLib::ChangeValueOfVariable("double GrainDragCoefficient",$InputLine,"main/Comet.h");
  } 

  
  #forces that will be accounted during the simulation
  elsif ($InputLine eq "FORCES") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
    if ($InputLine eq "GRAVITY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__GRAVITY_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__GRAVITY_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }   
    elsif ($InputLine eq "FRAMEROTATION") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__FRAME_ROTATION_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__FRAME_ROTATION_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "RADIATIONPRESSURE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__RADIATION_PRESSURE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__RADIATION_PRESSURE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "LORENTZFORCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__LORENTZ_FORCE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__LORENTZ_FORCE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }       
    elsif ($InputLine eq "DRAGFORCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "DRAGFORCETANGENTIALCOMPONENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE__TANGENTIAL_COMPONENT_","_PIC_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_CG_DUST_FORCE_MODE__DRAG_FORCE__TANGENTIAL_COMPONENT_","_PIC_MODE_OFF_","main/Comet.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
  #dust charging processes that will be modeled
  elsif ($InputLine eq "DUSTCHARGING") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "ELECTRONCOLLECTIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ELECTRON_COLLECTION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ELECTRON_COLLECTION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "IONCOLLECTIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ION_COLLECTION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__ION_COLLECTION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "PHOTOELECTRONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__PHOTO_ELECTRON_EMISSION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    elsif ($InputLine eq "SECONDARYELECTRONEMISSIONCURRENT") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if  ($InputLine eq "ON") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_","_PIC_MODE_ON_","models/dust/Dust.dfn");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING__SECONDARY_ELECTRON_EMISSION__MODE_","_PIC_MODE_OFF_","models/dust/Dust.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }
    
    elsif ($InputLine eq "MODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);

      if  ($InputLine eq "OFF")  {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__OFF_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__OFF_","pic/picGlobal.dfn");
      }
      elsif  ($InputLine eq "TIMEDEPENDENT") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__TIME_DEPENDENT_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "EQUILIBRIUM") {
        ampsConfigLib::RedefineMacro("_DUST__CHARGING_MODE_","_DUST__CHARGING_MODE__EQUILIBRIUM_","models/dust/Dust.dfn");
        ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_","pic/picGlobal.dfn");
      } 
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
    }    
    
    else {
      die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }    
  
  #gravity mode  
  elsif ($InputLine eq "GRAVITY3D") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_3DGRAVITY__MODE_","_3DGRAVITY__MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_3DGRAVITY__MODE_","_3DGRAVITY__MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "SAVEOUTPUTBINARY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_SAVING_BINARY_OUTPUT_MODE_","_SAVING_BINARY_OUTPUT_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_SAVING_BINARY_OUTPUT_MODE_","_SAVING_BINARY_OUTPUT_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "READNEUTRALSFROMBINARY") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_READ_NEUTRALS_FROM_BINARY_MODE_","_READ_NEUTRALS_FROM_BINARY_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_READ_NEUTRALS_FROM_BINARY_MODE_","_READ_NEUTRALS_FROM_BINARY_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "NUMBEROFNEUTRALSTOREAD") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("int Comet::CometData::nNeutrals","$InputLine","main/CometData.cpp");
  }
  elsif ($InputLine eq "RADIATIVECOOLINGMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "CROVISIER") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIATIVECOOLING__MODE_","_PIC_MODEL__RADIATIVECOOLING__MODE__CROVISIER_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIATIVECOOLING__MODE_","_PIC_MODEL__RADIATIVECOOLING__MODE__OFF_","pic/picGlobal.dfn");
      }
  }
  elsif ($InputLine eq "RADIALVELOCITYMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIAL_VELOCITY_MODE_","_PIC_MODEL__RADIAL_VELOCITY_MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__RADIAL_VELOCITY_MODE_","_PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_","pic/picGlobal.dfn");
      }
  }
  elsif ($InputLine eq "SAMPLEBACKFLUXMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_SAMPLE_BACKFLUX_MODE_","_SAMPLE_BACKFLUX_MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_SAMPLE_BACKFLUX_MODE_","_SAMPLE_BACKFLUX_MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "COMPUTEMAXIMUMLIFTABLESIZEMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE_","_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE__ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE_","_COMPUTE_MAXIMUM_LIFTABLE_SIZE_MODE__OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "HELIOCENTRICDISTANCE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double HeliocentricDistance","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double HeliocentricDistance","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "SUBSOLARPOINTAZIMUTH") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointAzimuth","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointAzimuth","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "SUBSOLARPOINTZENITH") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointZenith","$InputLine","main/Comet.cpp");
      ampsConfigLib::ChangeValueOfVariable("double subSolarPointZenith","$InputLine","main/main.cpp");
  }
  elsif ($InputLine eq "NDIST") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("static int ndist","$InputLine","main/Comet.h");
  }
  elsif ($InputLine eq "COMETTEMPERATUREMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "BJORN") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__BJORN_","main/Comet.dfn");
      }
      elsif ($InputLine eq "CONSTANT") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__CONSTANT_","main/Comet.dfn");
      }
      elsif ($InputLine eq "ANALYTICAL") {
	  ampsConfigLib::RedefineMacro("_COMET_TEMPERATURE_MODE_","_COMET_TEMPERATURE_MODE__ANALYTICAL_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "BJORNPRODUCTIONRATEUSERDEFINED") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_OFF_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "ANALYTICAL") {
	  ampsConfigLib::RedefineMacro("_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_","_BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ANALYTICAL_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "MODELSOURCEDFMS") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_MODEL_SOURCE_DFMS_","_MODEL_SOURCE_DFMS_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_MODEL_SOURCE_DFMS_","_MODEL_SOURCE_DFMS_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "MULTISPECIESANALYTICALMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_MULTISPECIES_ANALYTICAL_MODE_","_MULTISPECIES_ANALYTICAL_MODE_ON_","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_MULTISPECIES_ANALYTICAL_MODE_","_MULTISPECIES_ANALYTICAL_MODE_OFF_","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "TRACKINGSURFACEELEMENTMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_TRACKING_SURFACE_ELEMENT_MODE_","_TRACKING_SURFACE_ELEMENT_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_TRACKING_SURFACE_ELEMENT_MODE_","_TRACKING_SURFACE_ELEMENT_MODE_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "DFMSRATIOMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_DFMS_RATIO_MODE_","_DFMS_RATIO_MODE_ON_ ","main/Comet.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_DFMS_RATIO_MODE_","_DFMS_RATIO_MODE_OFF_ ","main/Comet.dfn");
      }
  }
  elsif ($InputLine eq "BJORNPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $BjornSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Bjorn_SourceRate\\[\\]",\@BjornSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "UNIFORMPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $UniformSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Uniform_SourceRate\\[\\]",\@UniformSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "JETPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $spec=$InputLine;
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      $JetSourceRate[getSpeciesNumber($spec)]=$InputLine;
      ampsConfigLib::ChangeValueOfArray("static double Jet_SourceRate\\[\\]",\@JetSourceRate,"main/Comet.h");
  }
  elsif ($InputLine eq "DUSTMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      if ($InputLine eq "ON") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__MODE_","_PIC_MODEL__DUST__MODE__ON_","pic/picGlobal.dfn");
      }
      elsif ($InputLine eq "OFF") {
	  ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__MODE_","_PIC_MODEL__DUST__MODE__OFF_","pic/picGlobal.dfn");
          ampsConfigLib::RedefineMacro("_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_","_PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__OFF_","pic/picGlobal.dfn");
      }
      else {
	  die "Unknown option\n";
      }
  }
  elsif ($InputLine eq "DUSTRMIN") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeMin","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "DUSTRMAX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeMax","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "NDUSTRADIUSGROUPS") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("int DustSampleIntervals","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "DUSTTOTALMASSPRODUCTIONRATE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustTotalMassProductionRate","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "POWERLAWINDEX") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("double DustSizeDistribution","$InputLine","main/Comet.cpp");
  }
  elsif ($InputLine eq "#ENDBLOCK") {
      last;
  }
   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
}


#=============================== Determine the species number  =============================                                                               
sub getSpeciesNumber {
    my $res=-1;
    my $species=$_[0];

    for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
	if ($species eq $SpeciesList[$i]) {
	    $res=$i;
	    last;
	}
    }

    if ($res eq -1) {
	die "Cannot fild species $_\n";
    }

    return $res;
}


=comment
#=============================== Change a definition of a macro in a source code  =============================
sub RedefineMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$WorkingSourceDirectory/$File") || die "Cannot open file $WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro /) {
#      print "\#define $Macro $Value\n";
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}
=cut
