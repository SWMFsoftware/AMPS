#!/usr/bin/perl
#$Id$
#reads the "Exosphere" section of the input file and modify the exosphere model of the PIC code

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use Cwd qw(cwd);
use FindBin qw($Bin);

use lib cwd,$Bin,"$Bin/input";
use ampsConfigLib;
use CometReactionNetwork;
use CometPhysicsDatabase;


#my $command= $];
#print "Perl version : ".$command;



#arguments: 
#$ARGV[0] -> the name of the intput file 
#$ARGV[1] -> the working directory


#print "Process the exospehre model input:\n";


my $InputFileName=$ARGV[0]; #"cg.input.Assembled.Block";   #$ARGV[0];
my $SpeciesFileName=$InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/; #"cg.input.Assembled.Species";
my $WorkingSourceDirectory=$ARGV[1];

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
# State below belongs to the new generic comet block.  The parser retains the
# older Hartley 2 keywords later in this file for DIM3 compatibility, while
# these flags make every dimension-sensitive or physics-sensitive quantity
# explicit before generated source is accepted.
my $Dimension=3;
my $NucleusRadiusIsSet=0;
my $DomainRadialExtentIsSet=0;
my $DomainOuterRadiusIsSet=0;
my $RadialCellNumberIsSet=0;
my $SurfaceResolutionIsSet=0;
my $HeliocentricDistanceAUIsSet=0;
my $ChemistryModeIsSet=0;
my $ConfiguredGeometry="";
my $PhysicsDatabaseName="NONE";
my $Rubin2011OHGroundRate;
my $Rubin2011MethanolH2COBranch;
my $RadialSimulationModeIsSet=0;
my $RadialCourantNumberIsSet=0;
my $RadialInnerBoundaryTemperatureIsSet=0;
my $RadialMaxIterationsIsSet=0;
my $RadialSamplingLengthIsSet=0;
my $RadialConvergenceToleranceIsSet=0;
my $RadialConvergenceWindowsIsSet=0;


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


# Legacy arrays remain for the historical triangulated Hartley 2 source
# routines.  The Model arrays are their geometry-neutral replacements and are
# indexed by the same compile-time SpeciesList order.
my @BjornSourceRate=(0)x$TotalSpeciesNumber;
my @UniformSourceRate=(0)x$TotalSpeciesNumber;
my @JetSourceRate=(0)x$TotalSpeciesNumber;
my @SurfaceProductionRate=(0)x$TotalSpeciesNumber;
my @SurfaceDistributionMode=(0)x$TotalSpeciesNumber;
my @SurfaceTemperature=(-1.0)x$TotalSpeciesNumber;
my @SurfaceBulkSpeed=(-1.0)x$TotalSpeciesNumber;
my @ExtendedProductionRate=(0)x$TotalSpeciesNumber;
my @ExtendedDistributionMode=(0)x$TotalSpeciesNumber;
my @RadialCharacteristicSpeed=(-1.0)x$TotalSpeciesNumber;
my $ReactionNetwork=CometReactionNetwork->new(species=>\@SpeciesList);
my $PhysicsDatabase=CometPhysicsDatabase->new(species=>\@SpeciesList);

#add the model header to the pic.h
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
print PIC_H "#include \"Comet.h\"\n";
print PIC_H "#include \"CometChemistry.h\"\n";
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
 
  $InputLine=~s/[=():,]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  
  # These records define the generic comet case.  Each value is written into
  # the generated working source tree with ampsConfigLib, matching the normal
  # AMPS compile-time input model.  Species-bearing records resolve names here
  # so runtime callbacks use integer indices only.
  if ($InputLine eq "DIMENSION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    $Dimension=$InputLine;

    if ($Dimension eq "1") {
      ampsConfigLib::RedefineMacro("DIM","1","pic/picGlobal.dfn");
      ampsConfigLib::RedefineMacro("_AMR_SYMMETRY_MODE_","_AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_","meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_AMR__CUT_CELL__MODE_","_AMR__CUT_CELL__MODE__OFF_","meshAMR/meshAMRdef.h");
      ampsConfigLib::RedefineMacro("_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_","_USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_OFF_","meshAMR/meshAMRdef.h");
      ampsConfigLib::AddLineAfterMarker2File("#include \"Exosphere_DefaultSurfaceDataStructure.h\"\n","//MARKER--add-definitions-after-this-line","meshAMR/meshAMR_UserDefinitions.h");

      # The second-order mover still contains an unfinished spherical branch.
      # Select the implemented first-order radial mover for a 1D build.
      ampsConfigLib::RedefineMacroFunction("_PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_","(ptr,LocalTimeStep,node) PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory(ptr,LocalTimeStep,node);","pic/picGlobal.dfn");
      ampsConfigLib::RedefineMacroFunction("_PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_","(ptr,LocalTimeStep,node,FirstBoundaryFlag) PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(ptr,LocalTimeStep,node,FirstBoundaryFlag);","pic/picGlobal.dfn");
      ampsConfigLib::RedefineMacroFunction("_PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_FACE_","(ptr,LocalTimeStep,node,FirstBoundaryFlag,InjectionFace) PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(ptr,LocalTimeStep,node,FirstBoundaryFlag);","pic/picGlobal.dfn");

      ampsConfigLib::SubstituteCodeLine('^MAINLIBOBJ =',"MAINLIBOBJ = RadialMesh.o RadialMeshProfile.o RadialSimulation.o CometModel.o CometChemistry.o CometPhysics.o","main/makefile");
      ampsConfigLib::SubstituteCodeLine('^MAINOBJ =',"MAINOBJ = RadialMain.o","main/makefile");
      ampsConfigLib::SubstituteCodeLine('#include "Comet.h"',"// The radial mesh build does not use the 3D Comet model header.","pic/pic.h");
      ampsConfigLib::ChangeValueOfVariable("int Comet::Model::GeometryMode","Comet::Model::Geometry::Spherical","main/CometModel.cpp");
      $ConfiguredGeometry="SPHERICAL";
    }
    elsif ($Dimension eq "3") {
      # The source defaults select the existing three-dimensional model.
      ampsConfigLib::ChangeValueOfVariable("int Comet::Model::GeometryMode","Comet::Model::Geometry::Triangulated","main/CometModel.cpp");
      $ConfiguredGeometry="TRIANGULATED";
    }
    else {
      die "The generic comet application currently executes Dimension=1 or Dimension=3; the axisymmetric Dimension=2 backend is designed but not implemented, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  elsif ($InputLine eq "GEOMETRY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;

    my %mode=(
      "SPHERICAL"=>"Comet::Model::Geometry::Spherical",
      "AXISYMMETRIC"=>"Comet::Model::Geometry::Axisymmetric",
      "TRIANGULATED"=>"Comet::Model::Geometry::Triangulated");
    die "Unknown comet geometry '$InputLine', line=$InputFileLineNumber ($InputFileName)\n"
        unless exists $mode{$InputLine};
    ampsConfigLib::ChangeValueOfVariable("int Comet::Model::GeometryMode",$mode{$InputLine},"main/CometModel.cpp");
    $ConfiguredGeometry=$InputLine;
  }
  elsif ($InputLine eq "HELIOCENTRICDISTANCEAU") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "HeliocentricDistanceAU must be a positive number, line=$InputFileLineNumber ($InputFileName)\n"
        unless isPositiveNumber($InputLine);
    ampsConfigLib::ChangeValueOfVariable("double Comet::Model::HeliocentricDistanceAU",$InputLine,"main/CometModel.cpp");
    ampsConfigLib::ChangeValueOfVariable("double Comet::Chemistry::HeliocentricDistanceAU",$InputLine,"main/CometChemistry.cpp");
    $HeliocentricDistanceAUIsSet=1;
  }
  elsif ($InputLine eq "SURFACEMESHFILE") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*SurfaceMeshFile\s*[=:]\s*//i;
    my ($file,$unit)=split(/[\s,]+/,$value,2);
    $file=trim($file);
    $unit="KM" unless defined $unit && $unit ne '';
    $unit=uc(trim($unit));
    my %scale=("M"=>"1.0","KM"=>"1.0e3");
    die "SurfaceMeshFile unit must be M or KM, line=$InputFileLineNumber ($InputFileName)\n"
        unless exists $scale{$unit};
    ampsConfigLib::ChangeValueOfVariable("char Comet::Model::SurfaceMeshFileName\\[_MAX_STRING_LENGTH_PIC_\\]","\"$file\"","main/CometModel.cpp");
    ampsConfigLib::ChangeValueOfVariable("double Comet::Model::SurfaceMeshLengthUnit",$scale{$unit},"main/CometModel.cpp");
  }
  elsif ($InputLine eq "SURFACESOURCE") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*SurfaceSource\s*[=:]\s*//i;
    my @field=map {trim($_)} split(/\|/,$value,-1);
    die "SurfaceSource expects SPECIES | DISTRIBUTION | RATE_S-1 | TEMPERATURE_K | BULK_SPEED_M_S, line=$InputFileLineNumber ($InputFileName)\n"
        unless @field==5;
    my ($name,$distribution,$rate,$temperature,$speed)=@field;
    my $species=getSpeciesNumber(uc($name));
    my %mode=("DISABLED"=>0,"UNIFORM"=>1,"COSINE"=>2,"TABULATED"=>3,
        "LEGACYBJORN"=>4,"LEGACYJET"=>5);
    $distribution=uc($distribution);
    die "Unknown surface distribution '$distribution', line=$InputFileLineNumber ($InputFileName)\n"
        unless exists $mode{$distribution};
    die "Invalid surface production rate '$rate', line=$InputFileLineNumber ($InputFileName)\n"
        unless isNonNegativeNumber($rate);
    $temperature=parseOptionalNumber($temperature,"surface temperature");
    $speed=parseOptionalNumber($speed,"surface bulk speed");
    $SurfaceProductionRate[$species]=$rate;
    $SurfaceDistributionMode[$species]=$mode{$distribution};
    $SurfaceTemperature[$species]=$temperature;
    $SurfaceBulkSpeed[$species]=$speed;
  }
  elsif ($InputLine eq "EXTENDEDSOURCE") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*ExtendedSource\s*[=:]\s*//i;
    my @field=map {trim($_)} split(/\|/,$value,-1);
    die "ExtendedSource expects SPECIES | DISTRIBUTION | TOTAL_RATE_S-1, line=$InputFileLineNumber ($InputFileName)\n"
        unless @field==3;
    my ($name,$distribution,$rate)=@field;
    my $species=getSpeciesNumber(uc($name));
    my %mode=("DISABLED"=>0,"RADIALPROFILE"=>1,"GRAINSUBLIMATION"=>2,
        "USERDEFINED"=>3);
    $distribution=uc($distribution);
    die "Unknown extended distribution '$distribution', line=$InputFileLineNumber ($InputFileName)\n"
        unless exists $mode{$distribution};
    die "Invalid extended production rate '$rate', line=$InputFileLineNumber ($InputFileName)\n"
        unless isNonNegativeNumber($rate);
    $ExtendedProductionRate[$species]=$rate;
    $ExtendedDistributionMode[$species]=$mode{$distribution};
  }
  elsif ($InputLine eq "CHEMISTRYMODE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    die "ChemistryMode must be ON or OFF, line=$InputFileLineNumber ($InputFileName)\n"
        unless $InputLine eq "ON" || $InputLine eq "OFF";
    ampsConfigLib::ChangeValueOfVariable("bool Comet::Chemistry::Enabled",
        $InputLine eq "ON" ? "true" : "false","main/CometChemistry.cpp");
    if ($InputLine eq "ON") {
      ampsConfigLib::RedefineMacro("_PIC_PHOTOLYTIC_REACTIONS_MODE_",
          "_PIC_PHOTOLYTIC_REACTIONS_MODE_ON_","pic/picGlobal.dfn");
      ampsConfigLib::RedefineMacroFunction(
          "_PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_",
          "(ptr,FirstParticleCell,node) Comet::Chemistry::ProcessParticle(ptr,FirstParticleCell,node);",
          "pic/picGlobal.dfn");
    }
    else {
      ampsConfigLib::RedefineMacro("_PIC_PHOTOLYTIC_REACTIONS_MODE_",
          "_PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_","pic/picGlobal.dfn");
    }
    $ChemistryModeIsSet=1;
  }
  elsif ($InputLine eq "CHEMISTRYREACTION") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*ChemistryReaction\s*[=:]\s*//i;
    $ReactionNetwork->add_reaction($value,
        "$FileName:$InputFileLineNumber");
  }
  elsif ($InputLine eq "PHYSICSDATABASE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    die "PhysicsDatabase must be NONE or RUBIN2011, line=$InputFileLineNumber ($InputFileName)\n"
        unless $InputLine eq "NONE" || $InputLine eq "RUBIN2011";
    $PhysicsDatabaseName=$InputLine;
  }
  elsif ($InputLine eq "NEUTRALCOLLISION") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*NeutralCollision\s*[=:]\s*//i;
    $PhysicsDatabase->add_collision($value,
        "$FileName:$InputFileLineNumber");
  }
  elsif ($InputLine eq "RUBIN2011OHGROUNDRATE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "Rubin2011OHGroundRate must be numeric, line=$InputFileLineNumber ($InputFileName)\n"
        unless isNonNegativeNumber($InputLine);
    $Rubin2011OHGroundRate=$InputLine+0.0;
  }
  elsif ($InputLine eq "RUBIN2011METHANOLH2COBRANCH") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    die "Rubin2011MethanolH2COBranch must be ON or OFF, line=$InputFileLineNumber ($InputFileName)\n"
        unless $InputLine eq "ON" || $InputLine eq "OFF";
    $Rubin2011MethanolH2COBranch=$InputLine;
  }
  elsif ($InputLine eq "RADIALSIMULATIONMODE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
    die "RadialSimulationMode must be ON or OFF, line=$InputFileLineNumber ($InputFileName)\n"
        unless $InputLine eq "ON" || $InputLine eq "OFF";
    ampsConfigLib::ChangeValueOfVariable(
        "bool Comet::RadialSimulation::Enabled",
        $InputLine eq "ON" ? "true" : "false","main/RadialSimulation.cpp");
    if ($InputLine eq "ON") {
      ampsConfigLib::RedefineMacro(
          "_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_",
          "_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_",
          "pic/picGlobal.dfn");
    }
    $RadialSimulationModeIsSet=1;
  }
  elsif ($InputLine eq "RADIALCOURANTNUMBER") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialCourantNumber must be in (0,1], line=$InputFileLineNumber ($InputFileName)\n"
        unless isPositiveNumber($InputLine) && $InputLine<=1.0;
    ampsConfigLib::ChangeValueOfVariable(
        "double Comet::RadialSimulation::CourantNumber",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialCourantNumberIsSet=1;
  }
  elsif ($InputLine eq "RADIALCHARACTERISTICSPEED") {
    my $value=$LineOriginal;
    $value=~s/!.*$//;
    $value=~s/^\s*RadialCharacteristicSpeed\s*[=:]\s*//i;
    my @field=map {trim($_)} split(/\|/,$value,-1);
    die "RadialCharacteristicSpeed expects SPECIES | SPEED_M_S, line=$InputFileLineNumber ($InputFileName)\n"
        unless @field==2;
    my $species=getSpeciesNumber(uc($field[0]));
    die "invalid radial characteristic speed '$field[1]', line=$InputFileLineNumber ($InputFileName)\n"
        unless isPositiveNumber($field[1]);
    $RadialCharacteristicSpeed[$species]=$field[1]+0.0;
  }
  elsif ($InputLine eq "RADIALINNERBOUNDARYTEMPERATURE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialInnerBoundaryTemperature must be positive, line=$InputFileLineNumber ($InputFileName)\n"
        unless isPositiveNumber($InputLine);
    ampsConfigLib::ChangeValueOfVariable(
        "double Comet::RadialSimulation::InnerBoundaryTemperature",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialInnerBoundaryTemperatureIsSet=1;
  }
  elsif ($InputLine eq "RADIALMAXITERATIONS") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialMaxIterations must be a positive integer, line=$InputFileLineNumber ($InputFileName)\n"
        unless defined $InputLine && $InputLine=~/^\d+$/ && $InputLine>0;
    ampsConfigLib::ChangeValueOfVariable(
        "long int Comet::RadialSimulation::MaxIterations",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialMaxIterationsIsSet=1;
  }
  elsif ($InputLine eq "RADIALSAMPLINGLENGTH") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialSamplingLength must be a positive integer, line=$InputFileLineNumber ($InputFileName)\n"
        unless defined $InputLine && $InputLine=~/^\d+$/ && $InputLine>0;
    ampsConfigLib::ChangeValueOfVariable(
        "long int Comet::RadialSimulation::SamplingLength",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialSamplingLengthIsSet=1;
  }
  elsif ($InputLine eq "RADIALCONVERGENCETOLERANCE") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialConvergenceTolerance must be in (0,1], line=$InputFileLineNumber ($InputFileName)\n"
        unless isPositiveNumber($InputLine) && $InputLine<=1.0;
    ampsConfigLib::ChangeValueOfVariable(
        "double Comet::RadialSimulation::ConvergenceTolerance",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialConvergenceToleranceIsSet=1;
  }
  elsif ($InputLine eq "RADIALCONVERGENCEWINDOWS") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    die "RadialConvergenceWindows must be a positive integer, line=$InputFileLineNumber ($InputFileName)\n"
        unless defined $InputLine && $InputLine=~/^\d+$/ && $InputLine>0;
    ampsConfigLib::ChangeValueOfVariable(
        "int Comet::RadialSimulation::RequiredConvergedWindows",$InputLine,
        "main/RadialSimulation.cpp");
    $RadialConvergenceWindowsIsSet=1;
  }
  elsif ($InputLine eq "NUCLEUSRADIUS") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("double Comet::RadialModel::NucleusRadius",$InputLine,"main/RadialMesh.cpp");
    $NucleusRadiusIsSet=1;
  }
  elsif ($InputLine eq "DOMAINRADIALEXTENT") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("double Comet::RadialModel::DomainRadialExtent",$InputLine,"main/RadialMesh.cpp");
    $DomainRadialExtentIsSet=1;
  }
  elsif ($InputLine eq "DOMAINOUTERRADIUS") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("double Comet::RadialModel::DomainOuterRadius",$InputLine,"main/RadialMesh.cpp");
    $DomainOuterRadiusIsSet=1;
  }
  elsif ($InputLine eq "RADIALCELLNUMBER") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("int Comet::RadialModel::RequestedCellNumber",$InputLine,"main/RadialMesh.cpp");
    $RadialCellNumberIsSet=1;
  }
  elsif ($InputLine eq "NUCLEUSSURFACERESOLUTION") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    ampsConfigLib::ChangeValueOfVariable("double Comet::RadialModel::SurfaceResolution",$InputLine,"main/RadialMesh.cpp");
    $SurfaceResolutionIsSet=1;
  }
  elsif ($InputLine eq "SODIUMSTICKINGPROBABILITY") {
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


  #type of the gas velocity distribution at the nucleus 
  elsif ($InputLine eq "INITIALGASVELOCITYDISTRIBUTION") { 
    ($InputLine,$InputComment)=split(' ',$InputComment,2);

     if ($InputLine eq "UNIFORM") { 
       ampsConfigLib::RedefineMacro("_COMET_GAS_INJECTION_VELOCITY_DIRECTION_MODE_","_COMET_GAS_INJECTION_VELOCITY_DIRECTION_MODE__UNIFORM_","main/Comet.dfn");
     }
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
  
  #dust initial velocity model
  elsif ($InputLine eq "DUSTINITIALVELOCITY") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    if ($InputLine eq "MODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      
      if ($InputLine eq "CONSTANT") {
        ampsConfigLib::ChangeValueOfVariable("int Comet::DustInitialVelocity::VelocityModelMode","Comet::DustInitialVelocity::Mode::ConstantVelocity","main/Comet.cpp");
      }
      elsif ($InputLine eq "ROTATIONBODY") {
        ampsConfigLib::ChangeValueOfVariable("int Comet::DustInitialVelocity::VelocityModelMode","Comet::DustInitialVelocity::Mode::RotationBody","main/Comet.cpp");        
      }
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }
    }
    elsif ($InputLine eq "INJECTIONCONSTANTSPEED") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      ampsConfigLib::ChangeValueOfVariable("double Comet::DustInitialVelocity::InjectionConstantVelocity","$InputLine","main/Comet.cpp");   
    } 
    elsif ($InputLine eq "ROTATIONPERIOD") {
      my @all;
      
      $line=$LineOriginal;
      chomp($line);
      $line=~s/[=:!]/ /g;
      
      @all=split(' ',$line);
      ampsConfigLib::ChangeValueOfVariable("double Comet::DustInitialVelocity::RotationPeriod","$all[2]","main/Comet.cpp");    
    } 
    elsif ($InputLine eq "ROTATIONAXIS") {
      my @Axis;
      my @t;
      
      $InputComment=~s/,/ /g;
      @t=split(' ',$InputComment);
      
      #normalize the rotation axis vector
      my $l=sqrt($t[0]*$t[0]+$t[1]*$t[1]+$t[2]*$t[2]);
      
      push(@Axis,$t[0]/$l);
      push(@Axis,$t[1]/$l);
      push(@Axis,$t[2]/$l);
      
      ampsConfigLib::ChangeValueOfArray("double Comet::DustInitialVelocity::RotationAxis\\[3\\]",\@Axis,"main/Comet.cpp");
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
  elsif ($InputLine eq "NUMERICALLOSSRATEINCREASE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      ampsConfigLib::ChangeValueOfVariable("const double NumericalLossRateIncrease","$InputLine","main/Comet.h");
      ampsConfigLib::ChangeValueOfVariable("const double NumericalLossRateIncrease","$InputLine","models/exosphere/Exosphere_Chemistry.cpp");
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
  
  
  #dust surface injection corection mode
  elsif ($InputLine eq "DUSTINJECTIONCORRECTIONMODE") {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
#     print "$InputLine\n";
      
      if (($InputLine eq "INJECTIONRATE") || ($InputLine eq "CORRECTIONFLUXRELATIVEH2O")) {
        if  ($InputLine eq "INJECTIONRATE") {
          ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__INJECTION_RATE_ ","main/Comet.dfn");
        }
        elsif ($InputLine eq "CORRECTIONFLUXRELATIVEH2O") {
          ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__RELATIVE_TO_H2O_ ","main/Comet.dfn");
        }
        
        #get the file name 
        my ($s1,$s2,$s3,$s4,$s5);
 
        $LineOriginal=~s/[=():,]/ /g;
        ($s1,$s2,$s3,$s4,$s5)=split(' ',$LineOriginal,5);
        chomp($LineOriginal);
        
        ampsConfigLib::ChangeValueOfVariable("char DustInjectionRateCorrectionFactorFileName\\[\\]","\"".$s4."\"","main/Comet.cpp");
      }
      elsif ($InputLine eq "OFF") {
        ampsConfigLib::RedefineMacro("_DUST_INJECTION_RATE_CORRECTION_MODE_","_DUST_INJECTION_RATE_CORRECTION_MODE__OFF_ ","main/Comet.dfn");
      }
      else {
        die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
      }    
  }   
  elsif ($InputLine eq "#ENDBLOCK") {
      last;
  }
   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
}

die "HeliocentricDistanceAU is required by the generic comet model\n"
    unless $HeliocentricDistanceAUIsSet;
die "Dimension=1 requires Geometry=Spherical\n"
    if $Dimension eq "1" && $ConfiguredGeometry ne "SPHERICAL";
die "Dimension=3 requires Geometry=Triangulated\n"
    if $Dimension eq "3" && $ConfiguredGeometry ne "TRIANGULATED";

# Named databases are expanded only after the complete application block has
# been read.  This makes the input keywords order independent while retaining
# strict validation of every paper-specific choice.  Explicit
# ChemistryReaction records may be used together with a named database; the
# reaction parser rejects duplicate IDs.
if ($PhysicsDatabaseName eq "RUBIN2011") {
  $PhysicsDatabase->load_rubin2011(
      network=>$ReactionNetwork,
      oh_ground_rate=>$Rubin2011OHGroundRate,
      methanol_h2co_branch=>$Rubin2011MethanolH2COBranch);
}

die "ChemistryMode=ON requires at least one ChemistryReaction\n"
    if $ChemistryModeIsSet && $ReactionNetwork->reaction_count()==0 &&
       chemistryIsEnabled("$WorkingSourceDirectory/main/CometChemistry.cpp");

for (my $species=0;$species<$TotalSpeciesNumber;$species++) {
  die "surface production for $SpeciesList[$species] requires a distribution\n"
      if $SurfaceProductionRate[$species]>0.0 &&
         $SurfaceDistributionMode[$species]==0;
  die "extended production for $SpeciesList[$species] requires a distribution\n"
      if $ExtendedProductionRate[$species]>0.0 &&
         $ExtendedDistributionMode[$species]==0;
}

ampsConfigLib::ChangeValueOfArray("double Comet::Model::SurfaceProductionRate\\[PIC::nTotalSpecies\\]",\@SurfaceProductionRate,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("int Comet::Model::SurfaceDistributionMode\\[PIC::nTotalSpecies\\]",\@SurfaceDistributionMode,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("double Comet::Model::SurfaceTemperature\\[PIC::nTotalSpecies\\]",\@SurfaceTemperature,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("double Comet::Model::SurfaceBulkSpeed\\[PIC::nTotalSpecies\\]",\@SurfaceBulkSpeed,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("double Comet::Model::ExtendedProductionRate\\[PIC::nTotalSpecies\\]",\@ExtendedProductionRate,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("int Comet::Model::ExtendedDistributionMode\\[PIC::nTotalSpecies\\]",\@ExtendedDistributionMode,"main/CometModel.cpp");
ampsConfigLib::ChangeValueOfArray("double Comet::RadialSimulation::CharacteristicSpeed\\[PIC::nTotalSpecies\\]",\@RadialCharacteristicSpeed,"main/RadialSimulation.cpp");
$ReactionNetwork->write_header(
    "$WorkingSourceDirectory/main/CometChemistry.generated.h");
$PhysicsDatabase->write_collision_header(
    "$WorkingSourceDirectory/main/CometPhysics.generated.h");

if ($Dimension eq "1") {
  die "Dimension=1 requires NucleusRadius\n" unless $NucleusRadiusIsSet;
  die "Dimension=1 requires RadialCellNumber\n" unless $RadialCellNumberIsSet;
  die "Dimension=1 requires NucleusSurfaceResolution\n" unless $SurfaceResolutionIsSet;
  die "Dimension=1 requires exactly one of DomainRadialExtent and DomainOuterRadius\n"
      unless ($DomainRadialExtentIsSet+$DomainOuterRadiusIsSet)==1;

  # The mesh-only mode intentionally needs no numerical evolution controls.
  # Once a particle run is selected, every timestep and termination input must
  # be explicit; no paper-independent numerical defaults are inserted here.
  if (radialSimulationIsEnabled(
      "$WorkingSourceDirectory/main/RadialSimulation.cpp")) {
    die "RadialSimulationMode=ON requires RadialCourantNumber\n"
        unless $RadialCourantNumberIsSet;
    die "RadialSimulationMode=ON requires RadialInnerBoundaryTemperature\n"
        unless $RadialInnerBoundaryTemperatureIsSet;
    die "RadialSimulationMode=ON requires RadialMaxIterations\n"
        unless $RadialMaxIterationsIsSet;
    die "RadialSimulationMode=ON requires RadialSamplingLength\n"
        unless $RadialSamplingLengthIsSet;
    die "RadialSimulationMode=ON requires RadialConvergenceTolerance\n"
        unless $RadialConvergenceToleranceIsSet;
    die "RadialSimulationMode=ON requires RadialConvergenceWindows\n"
        unless $RadialConvergenceWindowsIsSet;
    for (my $species=0;$species<$TotalSpeciesNumber;$species++) {
      die "RadialSimulationMode=ON requires RadialCharacteristicSpeed for $SpeciesList[$species]\n"
          unless $RadialCharacteristicSpeed[$species]>0.0;
    }
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

sub trim {
    my ($value)=@_;
    $value="" unless defined $value;
    $value=~s/^\s+//;
    $value=~s/\s+$//;
    return $value;
}

sub isNonNegativeNumber {
    my ($value)=@_;
    return 0 unless defined $value &&
        $value=~/^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?$/;
    return $value>=0.0;
}

sub isPositiveNumber {
    my ($value)=@_;
    return isNonNegativeNumber($value) && $value>0.0;
}

sub parseOptionalNumber {
    my ($value,$description)=@_;
    $value=trim($value);
    return -1.0 if uc($value) eq "UNSPECIFIED";
    die "Invalid $description '$value', line=$InputFileLineNumber ($InputFileName)\n"
        unless isNonNegativeNumber($value);
    return $value;
}

sub chemistryIsEnabled {
    my ($path)=@_;
    open(my $in,'<',$path) or die "Cannot open $path\n";
    while (my $entry=<$in>) {
      if ($entry=~/^bool\s+Comet::Chemistry::Enabled\s*=\s*true\s*;/) {
        close($in);
        return 1;
      }
    }
    close($in);
    return 0;
}

sub radialSimulationIsEnabled {
    my ($path)=@_;
    open(my $in,'<',$path) or die "Cannot open $path\n";
    while (my $entry=<$in>) {
      if ($entry=~/^bool\s+Comet::RadialSimulation::Enabled\s*=\s*true\s*;/) {
        close($in);
        return 1;
      }
    }
    close($in);
    return 0;
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
