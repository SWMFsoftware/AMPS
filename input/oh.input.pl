#!/usr/bin/perl
#$Id$
#reads the corresponding section of the input file and modifies source code

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use Scalar::Util qw/looks_like_number/;

use ampsConfigLib;


#my $command= $];
#print "Perl version : ".$command;



#arguments: 
#$ARGV[0] -> the name of the intput file 
#$ARGV[1] -> the working directory


#print "Process the exospehre model input:\n";


my $InputFileName=$ARGV[0]; #"europa.input.Assembled.Block"; #$ARGV[0];  # moon.input.Assembled.Block";
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


#get the list of the speces 
my $TotalSpeciesNumber;
my @SpeciesList;

open (SPECIES,"<$SpeciesFileName") || die "Cannot open file $SpeciesFileName\n";
$TotalSpeciesNumber=<SPECIES>;
@SpeciesList=<SPECIES>;
close (SPECIES);
      
chomp($TotalSpeciesNumber);
foreach (@SpeciesList) {
  chomp($_);
}


#add the model header to the pic.h
#open (PIC_H,">>$WorkingSourceDirectory/pic/pic.h") || die "Cannot open $WorkingSourceDirectory/pic/pic.h";
#print PIC_H "#include \"Europa.h\"\n";
#close (PIC_H);

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
  

  if ($InputLine eq "INJECTIONVELOCITY") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[=();,]/ /g;
      
      my @v;
      
      ($InputLine,$InputComment)=split(' ',$InputLine,2);
      @v=split(' ',$InputComment);
      
      ampsConfigLib::ChangeValueOfArray("double OH::InjectionVelocity\\[3\\]",\@v,"main/OH.cpp");  
  }

  elsif ($InputLine eq "INJECTIONNDENSITY") { 
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("double OH::InjectionNDensity",$s1,"main/OH.cpp");   
  }
  elsif ($InputLine eq "INJECTIONTEMPERATURE") {
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("double OH::InjectionTemperature",$s1,"main/OH.cpp");   
  }

  elsif ($InputLine eq "USERGLOBALTIMESTEP") {
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("double OH::UserGlobalTimeStep",$s1,"main/OH.cpp");   
  }

  elsif ($InputLine eq "VELOCITYDISTRIBUTIONSAMPLING"){
      my $ReadSampleLocations=0;
      my @SampleLocationTable;
      my $t;
      
      $InputComment=~s/[(),]/ /g;
      
      while (defined $InputComment) {
	  ($InputLine,$InputComment)=split(' ',$InputComment,2);
	  
	  if ($ReadSampleLocations == 1) { #reading of the sample location coordinates is in progress
	      if (looks_like_number($InputLine)) {
		  push(@SampleLocationTable,$InputLine);
	      } 
	      else {
		  $ReadSampleLocations=0;
	      }         
	  }
	  
	  if ($ReadSampleLocations == 0) { #'$InputLine' is not a number
	      if ($InputLine eq "ON") {
              ampsConfigLib::ChangeValueOfVariable("const bool OH::Sampling::DistributionFunctionSample::Use","true","main/OH_sample_distribution_function.cpp");

          }
          elsif ($InputLine eq "OFF") {
              ampsConfigLib::ChangeValueOfVariable("const bool OH::Sampling::DistributionFunctionSample::Use","false","main/OH_sample_distribution_function.cpp");
          }
          elsif ($InputLine eq "VMIN") {
            ($InputLine,$InputComment)=split(' ',$InputComment,2);
            
            if (defined $InputLine) {
              ampsConfigLib::ChangeValueOfVariable("double OH::Sampling::DistributionFunctionSample::vMin",$InputLine,"main/OH_sample_distribution_function.cpp");
            }
            else {
              warn ("Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)");
              die "#1 Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          }
          elsif ($InputLine eq "VMAX") {
            ($InputLine,$InputComment)=split(' ',$InputComment,2);
            
            if (defined $InputLine) {
              ampsConfigLib::ChangeValueOfVariable("double OH::Sampling::DistributionFunctionSample::vMax",$InputLine,"main/OH_sample_distribution_function.cpp");
            }
            else {
              warn ("Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)");
              die "#2 Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          }
          elsif ($InputLine eq "NSAMPLEINTERVALS") {
            ($InputLine,$InputComment)=split(' ',$InputComment,2);
            
            if (defined $InputLine) {
              ampsConfigLib::ChangeValueOfVariable("long int OH::Sampling::DistributionFunctionSample::nSampledFunctionPoints",$InputLine,"main/OH_sample_distribution_function.cpp");
            }
            else {
              warn ("Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)");
              die "#3 Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
            }
          }
          elsif ($InputLine eq "X") {
            #start reading of the sample location coordinates
            $ReadSampleLocations=1;
          }
          else {
            warn ("Cannot recognize the option (line=$InputLine, nline=$InputFileLineNumber)");
            die "#4 $InputLine: Cannot recognize line $InputFileLineNumber ($line) in $InputFileName.Assembled\n";
          }
        }
        
      }
      
      #insert the number of the sample points into the code
      my ($i,$j);
      my @Table;
      
      $j=0;
      
      for ($i=0;$i<(scalar @SampleLocationTable);$i++) {
        if ($j==0) {
          $t="{";
        }
                
        $t=$t."$SampleLocationTable[$i]";
        $j++;
        
        if ($j==3) {          
          $t=$t."}";
          $j=0;
         
          push(@Table,"$t");
        }
        else {
          $t=$t.",";
        }
      }
      
      ampsConfigLib::ChangeValueOfArray("double OH::Sampling::DistributionFunctionSample::SamplingLocations\\[\\]\\[3\\]",\@Table,"main/OH_sample_distribution_function.cpp");
      ampsConfigLib::ChangeValueOfVariable("int OH::Sampling::DistributionFunctionSample::nSampleLocations",(scalar @Table),"main/OH_sample_distribution_function.cpp");      

  }

  elsif ($InputLine eq "DOMAINXMAX") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[=();,]/ /g;
      
      my @x;
      
      ($InputLine,$InputComment)=split(' ',$InputLine,2);
      @x=split(' ',$InputComment);
      
      ampsConfigLib::ChangeValueOfArray("double OH::DomainXMax\\[3\\]",\@x,"main/OH.cpp");  
  }

  elsif ($InputLine eq "DOMAINXMIN") {
      ($InputLine,$InputComment)=split('!',$line,2);
      chomp($InputLine);
      $InputLine=~s/[=();,]/ /g;
      
      my @x;
      
      ($InputLine,$InputComment)=split(' ',$InputLine,2);
      @x=split(' ',$InputComment);
      
      ampsConfigLib::ChangeValueOfArray("double OH::DomainXMin\\[3\\]",\@x,"main/OH.cpp");  
  }

  elsif ($InputLine eq "DOMAINDXMIN") {
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("double OH::DomainDXMin",$s1,"main/OH.cpp");   
  }


  elsif ($InputLine eq "DOMAINDXMAX") {
    ($InputLine,$InputComment)=split('!',$line,2);
    chomp($InputLine);
    $InputLine=~s/[=();]/ /g;
    
    ($s0,$s1,$s2)=split(' ',$InputLine,3);
    ampsConfigLib::ChangeValueOfVariable("double OH::DomainDXMax",$s1,"main/OH.cpp");   
  }
  
  #the number of the ENAs source regions for sampling of the number density of particles produced in each individual source region
  elsif ($InputLine eq "SOURCEREGIONNUMBER") {
    ($InputLine,$InputComment)=split(' ',$InputComment,2);
    
    ampsConfigLib::ChangeValueOfVariable("int OH::Sampling::OriginLocation::nSampledOriginLocations",$InputLine,"main/OH.cpp")
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
  my $spec=$_[0];
  
  for (my $i=0;$i<$TotalSpeciesNumber;$i++) {
    if ($spec eq $SpeciesList[$i]) {
      $res=$i;
      last;
    }
  } 
  
  if ($res eq -1) {
    die "Cannot find species $spec\n";
  }
  
  return $res;
}

