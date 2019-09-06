#!/usr/bin/perl
#remove a column from a file 

use strict;
use warnings;
use POSIX qw(strftime);
use List::Util qw(first);
use IO::File;
use File::Copy;


#arguments:
#$ARGV[0] -> the column that need to be removed
#$ARGV[1] -> the name of the data file

my $ColumnOffset=$ARGV[0];

open (fIn,"<$ARGV[1]");
open (fOut,">$ARGV[1].tmp");

my $line;

while ($line=<fIn>) {
  my @l = split(/ /, $line);

  splice(@l,$ColumnOffset,1);
  print fOut "@l";
}

close (fIn);
close (fOut);

move ("$ARGV[1].tmp","$ARGV[1]");
