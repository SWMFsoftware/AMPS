#!/usr/bin/perl


push @INC, ".";

use strict;

#parse the command line: -fname=the file name pattern -n=number of threads -out = the number of the index 'out'
my ($fname,$nthreads,$out);


foreach (@ARGV) {
  if (/^-fname=(.*)/i) {
    $fname=$1;
  }
 
  if (/^-n=(.*)$/i) {
    $nthreads=$1;
  }
  
  if (/^-out=(.*)$/i) {
    $out=$1;
  }

  if (/^-h/i) {
    print "ConOutput.pl -fname=[the file name pattern] -n=[number of threads] -out =[the number of the index \'out\']\n";
    exit;
  }
}


#read the number of cells and corners in the resulted output file
my @mesh_size_files;
my @mesh_data_files;
my @mesh_connectivity_files;
my @corners_number;
my @cells_number;
my $line;

opendir(DIR, "."); 
my @all_files=readdir(DIR);




#@mesh_size_files = grep(/$fname\.mesh-size\.out=$out\.tread=.*\.tmp$/,@all_files);
#@mesh_data_files = grep(/$fname\.data\.out=$out\.tread=.*\.tmp$/,@all_files);
#@mesh_connectivity_files= grep(/$fname\.connectivity\.out=$out\.tread=.*\.tmp$/,@all_files);


for (my $i=0;$i<$nthreads;$i++) {
  my $fname_full;

  $fname_full=$fname."\.mesh-size\.tread=".$i."\.tmp"; 
  push @mesh_size_files,$fname_full;

  $fname_full=$fname."\.data\.tread=".$i."\.tmp";
  push @mesh_data_files,$fname_full;

print "$fname_full\n";


 $fname_full=$fname."\.connectivity\.tread=".$i."\.tmp";

print "$fname_full\n";


  push @mesh_connectivity_files,$fname_full;
}


print @mesh_size_files; 

print @mesh_data_files;
print @mesh_connectivity_files;

print "!!!!!!!!\n";
#mesh.mesh-size.out=0.tread=0.tmp


my $nTotalCells=0;
my $nTotalCorners=0;

for (my $i=0;$i<scalar @mesh_size_files;$i++) {
  open (fIn,"<",$mesh_size_files[$i]); 

  $line=<fIn>;
  $line=~s/Corners=//;
  push @corners_number,$line;
  $nTotalCorners+=$line;

  $line=<fIn>;
  $line=~s/Cells=//g; 
  push @cells_number,$line;
  $nTotalCells+=$line;

print "line=$line\n";
print "nTotalCorners= $nTotalCorners\n";

  close fIn;
}

print @corners_number;
print @cells_number;
print "$nTotalCorners\n";
print "$nTotalCells\n";

#create the resulting output file and copy the header
open (fOut,">","$fname\.all\.dat");

open (fIn,"<","$fname.\header.tmp");
$line=<fIn>;
print fOut  "$line\n";
print fOut "ZONE N=$nTotalCorners, E=$nTotalCells, DATAPACKING=POINT, ZONETYPE=FEBRICK\n";
close fIn;



#copy the data content 
for (my $i=0;$i<scalar @mesh_data_files;$i++) {
  open (fIn,"<",$mesh_data_files[$i]);

  while ($line=<fIn>) {
    print fOut  "$line";
  }

  close fIn; 
}

#combine the connectivity list  
my $offset=0;
my ($s0,$s1,$s2,$s3,$s4,$s5,$s6,$s7);

for (my $i=0;$i<scalar @mesh_connectivity_files;$i++) {
  open (fIn,"<",$mesh_connectivity_files[$i]);

  while ($line=<fIn>) {
    ($s0,$s1,$s2,$s3,$s4,$s5,$s6,$s7)=split(' ',$line);

    $s0+=1+$offset;
    $s1+=1+$offset;
    $s2+=1+$offset;
    $s3+=1+$offset; 
    $s4+=1+$offset; 
    $s5+=1+$offset;
    $s6+=1+$offset;
    $s7+=1+$offset; 

    print fOut  "$s0 $s1 $s2 $s3 $s4 $s5 $s6 $s7\n";
  } 

# $offset+=$cells_number[$i];

$offset+=$corners_number[$i]; 

print "offset=$offset\n";
}



close fOut;
