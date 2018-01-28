#!/usr/bin/perl
#$Id$
# The scrpit builds Makefile.test specific for the current machine
use strict;

# the hostname of the current machine
my $hostname = `hostname -s`;
chomp($hostname);

# check hostname for supercomputers
my $hostname_full = `hostname -f`;

my $Stampede = 'stampede';
if ($hostname_full =~ m/login2/) {
  $hostname = $Stampede;
}

my $Pleiades = 'pleiades';
if($hostname_full =~ m/pfe(.*)\.nas\.nasa\.gov/) {
  $hostname = $Pleiades;
}

my $Yellowstone = 'yellowstone';
if($hostname_full =~ m/yslogin(.*)/){
  $hostname = $Yellowstone;
}

if($hostname_full =~ m/srbwks2014-0079.engin.umich.edu/) {
  $hostname = "valeriy";
}

#path to the Makefile.test source
my $path="MakefileTest";

#table of the nightly tests
my $fTable="$path/Table";

#content of the table
my @Table=read_content($fTable);

my @Table_NoSpaces;

foreach (@Table) {
  if ($_ ne "\n") {
    push(@Table_NoSpaces,$_);
  }
}

@Table=@Table_NoSpaces;

#base for the Makefile.test
my $fBase="$path/makefile.test.base";
my @Base=read_content($fBase);

#individual tests for the Makefile.test
my $fApp="$path/makefile.test.app";
my @App=read_content($fApp);

my $fFinal="Makefile.test";
my @Final;
my @FinalApps;

#parameters of tests
my $Name;
my $Keys;
my $CustomRefSolutionPaths;
my $Outs;
my ($Refs,$ExeptionCodes);
my $Time;
my $TimeTotal = 0;
my $TimeLimit = 45;
my $iOvertime = 0;
#non-generic targets of tests 
my @Tars=('','','','');
#additional parameter to take into account location of input file
my $PathName;

#hostName only print once
my $HostNameFlag=0;

#determine whether time limit is set in the command line 
foreach (@ARGV) {
  if (/^-test-run-time=(.*)$/i) {
    $TimeLimit=$1;

    print "New Test run time limit is set of $TimeLimit \n";
  }
} 

#process table with test description
while (@Table) {
  my $refTable;
  my $refTars;

  ($refTable,$Name,$Time,$Keys,$CustomRefSolutionPaths,$Outs,$refTars,$Refs,$ExeptionCodes) = get_next_test(@Table);

  @Table   = @$refTable;
  @Tars = @$refTars;
  $PathName=$Name;
  if ($Name=~ m/(.*)\/(.*)$/) {$Name=$2;}
  $Outs="test_$Name" unless($Outs);
  next unless($Name);
  $TimeTotal+=$Time;
  &check_overtime;
   

  $Keys=~s/\s+$//; #remove spaces from the end of the line
  $CustomRefSolutionPaths=~s/\s+$//;
    
  for (my $i = 0; $i<=@Base-1; $i++) {
    #general part with all tests
    if ($Base[$i]=~m/<APP.*?>/) {
      
      if ($Base[$i]=~m/MAKE/) {
        $Final[$i]=$Final[$i]."\nifneq (\$(TEST<APP>-EXEPTIONCODE),SKIP)\n".$Base[$i]."endif";
      }
      else {
        $Final[$i]=$Final[$i].$Base[$i];
      }
            
      $Final[$i]=~s/<APP>/$Name/g;
      $Final[$i]=~s/<APPPATH>/$PathName/g;
      $Final[$i]=~s/<APPKEYS>/$Keys/g;
      $Final[$i]=~s/<APPCUSTOMREFSOLUTIONPATHS>/$CustomRefSolutionPaths/g;
      $Final[$i]=~s/<APPOUTS>/$Outs/g;
      $Final[$i]=~s/<APPREF>/$Refs/g;
      $Final[$i]=~s/<APPEXEPTIONCODE>/$ExeptionCodes/g;
    }
    else {
      if ($Base[$i]=~m/<HOST.*?>/ ) {
        if ($HostNameFlag==0) {
          #hostname only prints once
          $Final[$i]=$Final[$i].$Base[$i];
          $Final[$i]=~s/<HOSTNAME>/$hostname/g;
          $HostNameFlag=1;
        }
      }
      else {
        $Final[$i]=$Base[$i];
      }
    }
  }

  #application specific blocks
  my @lines = @App;
  #index of target blocks(1-Compile,2-Rundir,3-Run,4-Check)
  my $iTar = '';

  for (my $i = 0; $i<=@lines-1; $i++) {
    #substitute generic names with test's parameters
    $lines[$i]=~s/<APP>/$Name/g;
    $lines[$i]=~s/<APPPATH>/$PathName/g;
    $lines[$i]=~s/<APPKEYS>/$Keys/g;
    $lines[$i]=~s/<APPCUSTOMREFSOLUTIONPATHS>/$CustomRefSolutionPaths/g;
    $lines[$i]=~s/<APPOUTS>/$Outs/g;

    #check if inside a target block
    $iTar='' unless($lines[$i]=~m/^\t(.*)/);
    if ($lines[$i]=~m/^test\_$Name\_compile:/&& $Tars[0]) {$iTar=1;next;}
    if ($lines[$i]=~m/^test\_$Name\_rundir:/ && $Tars[1]) {$iTar=2;next;}
    if ($lines[$i]=~m/^test\_$Name\_run:/    && $Tars[2]) {$iTar=3;next;}
    if ($lines[$i]=~m/^test\_$Name\_check:/  && $Tars[3]) {$iTar=4;next;}

    if ($iTar) {
      #substitute the whole target
      $lines[$i]=$Tars[$iTar-1];
      $Tars[$iTar-1]='';
    }
  }

  push(@FinalApps,@lines);
}

#build the final Makefile
push(@Final,"\n\n");
push(@Final,@FinalApps);

#write it
write_content($fFinal,"include");
&write_content($fFinal,@Final);

#return 1;

#------------------------------------------------------------------------------
sub read_content {
  #read the content of a file
  open(my $fh,'<',$_[0]);
  my @lines = <$fh>;
  close($fh);
  @lines;
}

#------------------------------------------------------------------------------
sub write_content {
  # write the content of a file
  open(my $fh,'>',shift(@_));
  print $fh @_;
  close($fh);
}

#------------------------------------------------------------------------------
sub get_next_test{
  #returns test description extracted from content of fTable
  my $IsTestBlock='';

  #number of lines to remove
  my $nRemove=0;

  #test description
  my $Name=''; my $Keys=''; my $CustomRefSolutionPaths=''; my $Outs=''; my $Time=0;
  my $Refs='';
  my $ExeptionCodes='';

  #target block flag
  my $iTar='';

  #non-generic target descriptions
  my @Tars=('','','','');

  #error flag
  my $ErrorRead='';

  foreach my $line (@_) {
    #remove spaces at the end of $line
    $line =~ s/\s*$//g;

    if ($line =~ m/<#$/) {
      $IsTestBlock='0';

      #time for test MUST be set and be positive
      $ErrorRead='1' unless ($Time > 0);
      last;
    }

    if ($IsTestBlock) {
      #extract test description: Name, Key, Outs
      if ($line =~ m/Name=(.*)/) {
        unless ($Name) {$Name=$1;}
        else {
          #there can be only ONE definition for NAME
          $ErrorRead='1';
        }
      }
      elsif ($line =~ m/Time=([0-9]+)/) {$Time+=$1;}
      elsif ($line =~ m/Keys=(.*)/) {$Keys="$Keys$1,," if($1);}
      elsif ($line =~ m/CustomRefSolutionPaths=(.*)/) {$CustomRefSolutionPaths="$CustomRefSolutionPaths$1,," if($1);}
      elsif ($line =~ m/Outs=(.*)/) {$Outs="$Outs$1,," if($1);}
      elsif ($line =~ m/Ref=(.*)/) {$Refs="$Refs$1,," if($1);}  
      elsif ($line =~ m/ExeptionCode=(.*)/) {$ExeptionCodes="$ExeptionCodes$1,," if($1);}

      #non-generic target description header
      elsif ($line =~ m/Compile=(.*)/) {
        $iTar = 1; $Tars[0].=',,' if $Tars[0];
      }
      elsif ($line =~ m/Rundir=(.*)/ ) {
        $iTar = 2; $Tars[1].=',,' if $Tars[1];
      }
      elsif ($line =~ m/Run=(.*)/    ) {
        $iTar = 3; $Tars[2].=',,' if $Tars[2];
      }
      elsif ($line =~ m/Check=(.*)/  ) {
        $iTar = 4; $Tars[3].=',,' if $Tars[3];
      }

      #extract target description
      elsif ($line =~ m/^>>>(.*)/){$Tars[$iTar-1].="\t$1\n"if($1);}
      elsif ($line =~ m/^<<<(.*)/){$Tars[$iTar-1].="$1\n"if($1);}
      else {$ErrorRead='1';}
    }

    if ($line =~ m/^#>/) {
      $IsTestBlock='1';
    }

    $nRemove++;
  }

  #remove processed lines
  splice(@_,0,$nRemove+1);

  unless ($ErrorRead) {
    my @Keys    = split(/,,/,$Keys)   if($Keys);   $Keys   ='';
    my @CustomRefSolutionPaths = split(/,,/,$CustomRefSolutionPaths)   if($CustomRefSolutionPaths);   $CustomRefSolutionPaths   ='';
    my @Outs    = split(/,,/,$Outs)   if($Outs);   $Outs   ='';
    my @Compile = split(/,,/,$Tars[0])if($Tars[0]);$Tars[0]='';
    my @Rundir  = split(/,,/,$Tars[1])if($Tars[1]);$Tars[1]='';
    my @Run     = split(/,,/,$Tars[2])if($Tars[2]);$Tars[2]='';
    my @Check   = split(/,,/,$Tars[3])if($Tars[3]);$Tars[3]='';
    my @Refs    = split(/,,/,$Refs)   if($Refs);   $Refs   ='';
    my @ExeptionCodes    = split(/,,/,$ExeptionCodes)   if($ExeptionCodes);   $ExeptionCodes   ='';


    unless ($ErrorRead) {
      #process parameters for this machine
      $Name = process_option($Name);

      #process Keys
      foreach my $Key (@Keys) {
        $Key = process_option($Key);
        $Keys="$Keys$Key " if($Key);
      }

      #process CustomRefSolutionPaths
      foreach my $CustomRefSolutionPath (@CustomRefSolutionPaths) {
        $CustomRefSolutionPath = process_option($CustomRefSolutionPath);
        $CustomRefSolutionPaths="$CustomRefSolutionPaths$CustomRefSolutionPath " if($CustomRefSolutionPath);
      }
 
      #process Outs
      foreach my $Out (@Outs){
        $Out = process_option($Out);
        $Outs="$Outs$Out " if($Out);
      }

      #process Refs
      foreach my $Ref (@Refs) {
        $Ref = process_option($Ref);

        #the reference change is requested for this machine
        if ($Ref) {
          my ($nRef,$Compiler,$TestName); 

          $Ref=~s/[}{]/ /g;
          ($nRef,$Compiler)=split(' ',$Ref); 

          $TestName=$Name;
          if ($TestName=~ m/(.*)\/(.*)$/) {$TestName=$2;}
          $Refs=$Refs."\n\nifeq (\$(COMPILE.c),$Compiler)\nTEST".$TestName."-REF=[$nRef]\nendif\n";
        }
      }
      
      #process Elexptions
      foreach my $ExeptionCode (@ExeptionCodes) {
        $ExeptionCode  = process_option($ExeptionCode);

        #the individual exeption code is requested for this machine
        if ($ExeptionCode) {
          my ($Code,$Compiler,$TestName); 

          $ExeptionCode=~s/[}{]/ /g;
          ($Code,$Compiler)=split(' ',$ExeptionCode); 

          $TestName=$Name;
          if ($TestName=~ m/(.*)\/(.*)$/) {$TestName=$2;}
          $ExeptionCodes=$ExeptionCodes."\n\nifeq (\$(COMPILE.c),$Compiler)\nTEST".$TestName."-EXEPTIONCODE=$Code\nendif\n";
        }
      }

      #process targets
      foreach my $Tar (@Compile){
        $Tar    = process_option($Tar);
        $Tars[0]= "$Tar"if($Tar);
      }

      foreach my $Tar (@Rundir){
        $Tar    = process_option($Tar);
        $Tars[1]= "$Tar"if($Tar);
      }

      foreach my $Tar (@Run){
        $Tar    = process_option($Tar);
        $Tars[2]= "$Tar"if($Tar);
      }

      foreach my $Tar (@Check){
        $Tar    = process_option($Tar);
        $Tars[3]= "$Tar"if($Tar);
      }

    }
    (\@_,$Name,$Time,$Keys,$CustomRefSolutionPaths,$Outs,\@Tars,$Refs,$ExeptionCodes);
  }
  else {
    @Tars=('',0,'','','');
    (\@_,'','','',\@Tars,$Refs,$ExeptionCodes);
  }
}

#------------------------------------------------------------------------------
sub process_option{
  my $Machines;
  my $OptionName=$_[0];
    
  if ($OptionName =~ m/(.*)\t?@@\((.*)\)\n?/s) {
    # check if current machine is in the list of
    $OptionName = $1;
    $Machines   = $2;

    # first, count number of occurences
    my $nOccur = () = $Machines =~ /$hostname/g;

    if ($nOccur > 1) {return('');}
    if ($nOccur==0) {
      #check if all machines are removed
	    if($Machines =~ m/^0/){return('');}
    }
    
    if ($nOccur==1) {
      if ($Machines =~ m/\+$hostname/) {
      #do nothing
	    }
	    elsif ($Machines =~ m/-$hostname/) {return('');}
	    else {return('');}
    }
  }

  return ($OptionName);
}

#------------------------------------------------------------------------------
sub check_overtime {
  # for supercomputers time in development queue is limited to 2 hours;
  # check if test exceed the limit and create separate job-files in the case
  #..........................................................................
  # check if this is a supercomputer
  return unless(($hostname eq $Stampede) || ($hostname eq $Pleiades)|| ($hostname eq $Yellowstone));

  # check if went overtime
  return if ($TimeTotal < $TimeLimit);
    
  # number of this overtime case
  $iOvertime++;
  # reset time counter
  $TimeTotal -= $TimeLimit;

  # find line in final Makefile.test where to unsert a new target
  for (my $i = 0; $i<=@Final-1; $i++) {
    next unless ($Final[$i] =~ m/^test_run:/);
	
    #current line is the one with default target name
	  #append new target name at the next line
    $Final[$i+1] .= "\ntest_run_OVERTIME$iOvertime:\n";	
    last;
  }

  # create a new job file for overtimed tests for each compiler
  my $path = 'utility/TestScripts';
 
  my $original="$path/test_amps.$hostname.all.job";
  my $overtime="$path/test_amps.$hostname.all.overtime$iOvertime.job";
	
  if (-e $original) {
	  system("cp $original $overtime");
    my @lines = read_content("$overtime");
	    
    foreach my $line (@lines) {
      $line =~ s/test_run/test_run_OVERTIME$iOvertime/;
	  }
	    
	  &write_content($overtime, @lines);
	}
}
