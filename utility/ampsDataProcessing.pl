#!/usr/bin/perl
#$Id$

use strict;
use warnings;

#the structure of the array FileMask: [mask,rm,preplot,send]  
my @FileMask;

#the structure of the array FileList: [FileName,rm,preplot,send]
my @FileList;

#constants
use constant _TRUE_ => 0;
use constant _FALSE_ => 1;

#variables used for the parallel execution of the script 
my $nTotalThreads=1;
my $ThisThread=0;
my @childs;

my $host;
my $dir=".";
my $wait=_FALSE_;

#read the argument line
foreach (@ARGV) {
  
  if (/^-np=(.*)/i) {$nTotalThreads=$1; next};
  if (/^-host=(.*)/i) {$host=$1; next};
  if (/^-dir=(.*)/i) {$dir=$1; next};
  if (/^-wait/i) {$wait=_TRUE_; next;}
  
  if (/^-send=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_FALSE_, 'preplot'=>_FALSE_,'send'=>_TRUE_});next};
  if (/^-send-rm=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_FALSE_,'send' =>_TRUE_});next};
  
  if (/^-preplot=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_TRUE_,'send'=>_TRUE_});next};
  
  if (/^-rm=(.*)/i) {push(@FileMask, {'mask'=>$1, 'rm'=>_TRUE_, 'preplot'=>_FALSE_,'send'=>_FALSE_});next};
  
  if (/^-help/i) {
    print "Example:\n";
    print "./ampsDataProcessing.pl -np=1 -host=tower-left.engin.umich.edu:/Volumes/data/EUROPA-LOCAL-RUN--TEMP -preplot=\'*.sdf\' -dir=./\n"; 
    print "The argument line:\n";
    print "-help             -> print the list of the arguments\n";
    print "-wait   ->  the script will wait for the new files being generated untill it is killed by user\n";
    print "-np [number]      -> the number of threads\n\n";

    print "-send=[mask for file search, separated by ',']    -> files send to the remote host without preprocessing, the files are NOT removed\n";
    print "-send-rm=[mask for file search, separated by ','] -> files send to the remote host without preprocessing, the files are will be removed\n\n";
    
    print "-preplot=[mask for file search, separated by ',']    -> preplot the data files and send; after the files sent ther are removed\n";

    print "-rm=[mask for file search, separated by ',']   -> remove the files\n";
    
    print "-host   -> the name of the remote host and directory where the files will be copied\n";
    
    print "-dir   -> the directory where the output files will be created\n";
    exit;
  }
}



print "Starting post processing script\n";
print "Total threads:  $nTotalThreads\n";
 
do {
  if (-d $dir) {
    #create the file list
    for (my $i=0;$i<=$#FileMask;$i++) {
      my @mask;
      my @files;
      
      @mask=split(',',$FileMask[$i]{'mask'});
      
      foreach (@mask) {
        @files=glob "$dir"."/"."$_";
        
        foreach (@files) {
          push(@FileList,{'file'=>$_},'rm'=>$FileMask[$i]{'rm'}, 'preplot'=>$FileMask[$i]{'preplot'}, 'send'=>$FileMask[$i]{'send'});
        }      
      }     
    }
	
    #start slave processes	
    if (@FileList) {
      for (my $count=0;$count<$nTotalThreads;$count++) {
        my $pid = fork();
      
        if ($pid) {
          # parent
          #print "pid is $pid, parent $$\n";
          push(@childs, $pid);
        } elsif ($pid == 0) {
           # child
           ProcessDataFiles($count);
           exit 0;
        } else {
           die "couldnt fork: $!\n";
        } 
      }
     
      foreach (@childs) {
        my $tmp = waitpid($_, 0);
      }
    }
  }
  
   sleep(120);
}
while ($wait == _TRUE_);
 
print "Done.\n";

#=============================== Process Data Files =============================
sub ProcessDataFiles {
  #process the data files
  
  for (my $i=0;$i<=$#FileList;$i++) {
    if ($i%$nTotalThreads==$ThisThread) {
      #process the file
      my $fname=$dir."/".$FileList[$i]{'file'};
      my $rm=$dir."/".$FileList[$i]{'rm'};
      my $send=$dir."/".$FileList[$i]{'send'};
      
      #preplot the data file
      if (-e $fname) {
        `preplot $fname`;
        `rm -f $fname`;
        $fname=~s/.dat$/*plt/;
        
        $rm=_TRUE_;
        $send=_TRUE_;
      }
      
      #send the data file
      if ((-e $fname) && ($send == _TRUE_)) {
        `scp $fname $host`;
      }
      
      #remove the data file
      if ((-e $fname) && ($rm == _TRUE_)) {
        `rm -f $fname`;
      }      
    }
  }
}
 


