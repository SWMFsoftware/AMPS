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


my $InputFileName=$ARGV[0]; #"protostellar-nebula.input.Assembled.Block"; #$ARGV[0];  # moon.input.Assembled.Block";
my $SpeciesFileName=$InputFileName; $SpeciesFileName =~ s/\.Block$/.Species/;
my $WorkingSourceDirectory=$ARGV[1];   #"srcTemp"; #$ARGV[1];   # srcTemp

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


#add the model header to the pic.h
open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
print PIC_H "#include \"ProtostellarNebula.h\"\n";
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
 
  $InputLine=~s/[=()]/ /g;
  ($InputLine,$InputComment)=split(' ',$InputLine,2);
  $InputLine=~s/ //g;
  

  if ($InputLine eq "SPHEREINSIDEDOMAIN") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    $InputLine=~s/ //g;
     
    if ($InputLine eq "ON") {
      ampsConfigLib::ChangeValueOfVariable("static const bool SphereInsideDomain","true","main/main_lib.cpp");
    } 
    elsif ($InputLine eq "OFF") {
      ampsConfigLib::ChangeValueOfVariable("static const bool SphereInsideDomain","false","main/main_lib.cpp");
    }
    else {
      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
    }
  }
  
#  elsif ($InputLine eq "EPD_FLUX") {  #the injection flux of the high energy ions
#    ($InputLine,$InputComment)=split('!',$line,2);
#    chomp($InputLine);
#    $InputLine=~s/[=();]/ /g;
    
#    ($s0,$s1,$s2)=split(' ',$InputLine,3);
#    ampsConfigLib::ChangeValueOfVariable("const static double EPD_Flux",$s1,"main/Europa.h");   
#  }

#  elsif ($InputLine eq "UNIMOLECULARREACTION") { 
#    ($InputLine,$InputComment)=split(' ',$InputComment,2);
#    $InputLine=~s/ //g;
#     
#    if ($InputLine eq "PHOTOIONIZATION") {
#      ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_\n#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_\n\n","main/UserDefinition.PIC.h");   
#      ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_\n#define _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(t0,t1,t2,t3,t4) Europa::ExospherePhotoionizationReactionProcessor(t0,t1,t2,t3,t4);\n","main/UserDefinition.PIC.h");
#      ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_\n#define _PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(t0,t1,t2,t3) Europa::ExospherePhotoionizationLifeTime(t0,t1,t2,t3);\n","main/UserDefinition.PIC.h");
#    }    
#    elsif ($InputLine eq "GENERICTRANSFORMATION") {
#      my $ReactionProcessor;
#      
#      ($InputLine,$InputComment)=split(' ',$InputComment,2);
#      if ($InputLine eq "FUNC") {    
#        $line=~s/[=()]/ /g;
#        ($ReactionProcessor,$line)=split(' ',$line,2);
#        ($ReactionProcessor,$line)=split(' ',$line,2);
#        ($ReactionProcessor,$line)=split(' ',$line,2);
#        ($ReactionProcessor,$line)=split(' ',$line,2);
#      }
#      else {
#        die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
#      }
      
#      ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_\n#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_\n\n","main/UserDefinition.PIC.h");
#      ampsConfigLib::AddLine2File("#undef _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_\n#define _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_\n\n","main/UserDefinition.PIC.h");
#      ampsConfigLib::AddLine2File("#undef _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_\n#define _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(t0,t1,t2,t3,t4,t5,t6,t7) $ReactionProcessor(t0,t1,t2,t3,t4,t5,t6,t7);\n","main/UserDefinition.PIC.h"); 
#    }
#    elsif ($InputLine eq "OFF") {
#      ampsConfigLib::AddLine2File("#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_\n#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_\n\n","main/UserDefinition.PIC.h");
#    }
#    else {
#      die "The option is not recognized, line=$InputFileLineNumber ($InputFileName)\n";
#    }
#  }  
  
  elsif ($InputLine eq "FORCES") {    
    $InputComment=~s/[,]/ /g;
    
    while (defined $InputComment) {
      ($InputLine,$InputComment)=split(' ',$InputComment,2);
      $InputLine=~s/ //g;
      
      
      if ($InputLine eq "OFF") {
        #do nothing
      }
      elsif ($InputLine eq "GRAVITY") {
        ampsConfigLib::RedefineMacro("_FORCE_GRAVITY_MODE_","_PIC_MODE_ON_","main/ProtostellarNebula.h");
      }
      elsif ($InputLine eq "LORENTZ") {
        ampsConfigLib::RedefineMacro("_FORCE_LORENTZ_MODE_","_PIC_MODE_ON_","main/ProtostellarNebula.h");
      }
      elsif ($InputLine eq "FRAMEROTATION") {
        ampsConfigLib::RedefineMacro("_FORCE_FRAMEROTATION_MODE_","_PIC_MODE_ON_","main/ProtostellarNebula.h");
      }
      else {
        die "The option '$InputLine' is not recognized";
      }
    }
  }

 elsif($InputLine eq "STELLARANGULARVELOCITY"){
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("static const double Omega",$s1,"main/ProtostellarNebula.h");        
 }
  
 
 
  elsif ($InputLine eq "#ENDBLOCK") {
    last;
  }   
  else {
    die "Option is unknown, line=$InputFileLineNumber ($InputFileName)\n";
  }
  
  
  
  #extract and redefine the ProtostellarNebula's model related variables from '.ampConfig.Settings' ('.ampConfig.Settings. is created by Config.pl)
  if (-e ".ampsConfig.Settings") {
    my @Settings;
  
    open (AMPSSETTINGS,".ampsConfig.Settings") || die "Cannot open file\n";
    @Settings=<AMPSSETTINGS>;
    close (AMPSSETTINGS);
  
    
    foreach (@Settings){
#      if (/^ICESLOCATION=(.*)$/i) {ampsConfigLib::ChangeValueOfVariable("const char IcesLocationPath\\[\\]","\"".$1."\"","main/main_lib.cpp"); next};
    }    
  }

  
}



