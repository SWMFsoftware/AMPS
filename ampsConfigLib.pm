#$Id$

use strict;
use warnings;

use POSIX qw(strftime);
use List::Util qw(first);

package ampsConfigLib;
our @EXPORT=qw($WorkingSourceDirectory);

my $WorkingSourceDirectory;


sub printwd {
  print "$ampsConfigLib::WorkingSourceDirectory\n";
}

#=============================== Change a value of a variable in the code  =============================
sub ChangeValueOfVariable {
  my $Variable=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  
  foreach (@FileContent) {
    if ($_=~/($Variable)/) {
      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."=".$Value.";\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub ChangeValueOfVariableNoSemiColon {
  my $Variable=$_[0];
  my $Value=$_[1];
  my $File=$_[2];

  my @FileContent;

  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  @FileContent=<FILEIN>;
  close (FILEIN);

  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");


  foreach (@FileContent) {
    if ($_=~/($Variable)/) {
      my $t=$Variable;

      $t=~s/\\//g;
      $_=$t."=".$Value."\n";
    }

    print FILEOUT "$_";
  }

  close (FILEOUT);
}

#=============================== Change a definition of a macro in a source code  =============================
sub RedefineMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro /) {
      $_="\#define $Macro $Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub RedefineMacroFunction {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
  
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";
  
  foreach (@FileContent) {
    if ($_=~/\#define $Macro/) {
      $_="\#define $Macro$Value\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

sub AddMacro {
  my $Macro=$_[0];
  my $Value=$_[1];
  my $File=$_[2];
  
  InsertText2File("//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n#undef $Macro\n#define $Macro $Value\n\n",$File);
}
  
#=============================== Change a value of a array in the code  =============================
sub ChangeValueOfArray {
  my $Variable=$_[0];
  my @Value=@{$_[1]};
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($Variable)/) {
      my $t=$Variable;
      
      $t=~s/\\//g;
      $_=$t."={";
      
      for (my $i=0;$i<1+$#Value;$i++) {
        $_=$_."$Value[$i]";
        
        if ($i ne $#Value) {
          $_=$_.",";
        }
      }
      
      $_=$_."};\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT); 
}

#===============================  Get the first appearence of a string $_[0] in the array $_[1] ==
sub GetElementNumber {
  my $str=$_[0];
  my @List=@{$_[1]};
  
  my $nspec=0;
  my $found=0;
  
  foreach (@List) {
    if ($List[$nspec] eq $str) {
      $found=1; #the species has been found in the list of those used in the simulation
      last;
    }
    else {
      $nspec++;
    }
  }
  
  if ($found == 0) {
    $nspec=-1;
  }

  return $nspec;
}

#===============================  Substitute a line in a source file =============================
sub SubstituteCodeLine {
  my $oldLinstKey=$_[0];
  my $newLine=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    if ($_=~m/($oldLinstKey)/) {
      $_=$newLine."\n";
    }
    
    print FILEOUT "$_";
  }
  
  close (FILEOUT);
}

#===============================  Insert a line into a file after a marker ==========================
sub AddLineAfterMarker2File {
  my $NewLine=$_[0];
  my $Marker=$_[1];
  my $File=$_[2];
  
  my @FileContent;
  
  open (FILEIN,"<$ampsConfigLib::WorkingSourceDirectory/$File") || die "Cannot open file $ampsConfigLib::WorkingSourceDirectory/$File\n";  
  @FileContent=<FILEIN>;
  close (FILEIN);
 
  open (FILEOUT,">$ampsConfigLib::WorkingSourceDirectory/$File");
  
  foreach (@FileContent) {
    print FILEOUT "$_";
    
    if ($_=~m/($Marker)/) {
      print FILEOUT  "//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n";
      print FILEOUT "$NewLine\n";
    }
  }
  
  close (FILEOUT);
}  


#===============================  RECURSIVE SUBSTITUDE OF A STRING IN THE SOURCE DIRECTORY ==========
sub RecursiveSubstitute {
  my $init=$_[0];
  my $final=$_[1];
  my $dir=$_[2];
  
  my $fname;
  my @FileList;
  
  opendir(DIR,$dir) or die "Cannot open directory $dir\n";
  @FileList=readdir(DIR);
  closedir(DIR);
  
  foreach $fname (@FileList) {
    #check if '$fname' is a directory    
    if (($fname eq '.')||($fname eq '..')||($fname eq 'CVS')||($fname eq 'RCS')) {
      next;
    }
    
    if (-d "$dir/$fname") {
      #fname is a directory -> call the function recursevely
      RecursiveSubstitute($init,$final,"$dir/$fname");
    }
    else {
      #fname is a file -> make the substitution
      my @lines;
      
      open (FILE,"<$dir/$fname") || die "Cannot open file $dir/$fname\n";
      @lines=<FILE>;
      close (FILE);
      
      foreach (@lines) {
        $_=~s/$init/$final/g;
      }
      
      open (FILE,">$dir/$fname") || die "Cannot open file $dir/$fname\n";
      print FILE @lines;
      close (FILE);
    }
  }
}

sub AddLine2File {
  my $fname=$_[1];
  my $newline=$_[0];
  
  InsertText2File("//The following changes are added while reading input file. [".POSIX::strftime("%m/%d/%Y", localtime)."]\n$newline\n",$fname);
}

#=============================== Add text to a header of the build tree  =============================
# The settings of the application input are written into the build copies of the headers (picGlobal.dfn,
# Exosphere.dfn, meshAMRdef.h, ...).  Text appended after the #endif that closes a header's include guard is
# applied again by every repeated #include of the header (meshAMRcutcell.h and Exosphere.h include
# picGlobal.dfn after pic.h, for example).  In the middle of a translation unit it then undoes whatever was set
# in between, e.g. the nightly-test settings, and the code before and after that point, or different
# translation units, are compiled with different values of the same macro.  InsertText2File therefore puts
# the text, in the order of the calls,
#   - before the line that contains $MARKER:APPLICATION-SETTINGS$ (picGlobal.dfn: the nightly-test settings
#     that follow the marker take precedence over the application input), otherwise
#   - before the #endif that closes the file's include guard (#ifndef X / #define X ... #endif), otherwise
#   - at the end of the file (no include guard: as before; a missing file is created, as by the former '>>').
sub InsertText2File {
  my $Text=$_[0];
  my $File=$_[1];
  my $Path="$ampsConfigLib::WorkingSourceDirectory/$File";
  my @FileContent;

  if (! -e $Path) {
    open (FILEOUT,">>$Path") || die "Cannot open file $Path\n";
    print FILEOUT $Text;
    close (FILEOUT);
    return;
  }

  open (FILEIN,"<$Path") || die "Cannot open file $Path\n";
  @FileContent=<FILEIN>;
  close (FILEIN);

  $Text.="\n" unless ($Text=~/\n$/); #the text is followed by the guard's '#endif' line
  splice(@FileContent,InsertionLine(\@FileContent),0,$Text);

  open (FILEOUT,">$Path") || die "Cannot open file $Path\n";
  print FILEOUT @FileContent;
  close (FILEOUT);
}

#the index of the line before which InsertText2File places new text
sub InsertionLine {
  my @Lines=@{$_[0]};
  my @Directives; #[line index, directive, argument] of the preprocessor directives and code lines, without comments
  my $InComment=0;

  for (my $i=0;$i<=$#Lines;$i++) {
    return $i if ($Lines[$i]=~/\$MARKER:APPLICATION-SETTINGS\$/);
  }

  for (my $i=0;$i<=$#Lines;$i++) {
    my $l=$Lines[$i];

    if ($InComment) {
      next unless ($l=~s/^.*?\*\///);
      $InComment=0;
    }

    $l=~s/\/\*.*?\*\///g;
    $InComment=1 if ($l=~s/\/\*.*$//);
    $l=~s/\/\/.*$//;

    if ($l=~/^\s*#\s*(\w+)\s*(.*?)\s*$/) {push(@Directives,[$i,$1,$2]);}
    elsif ($l=~/\S/) {push(@Directives,[$i,"",""]);}
  }

  #the include guard: the last directive is #endif, its #ifndef X is followed by '#define X'
  return scalar(@Lines) if ((scalar(@Directives)==0)||($Directives[-1][1] ne "endif"));

  my $depth=0;

  for (my $n=$#Directives;$n>=0;$n--) {
    my $d=$Directives[$n][1];

    if ($d eq "endif") {$depth++;}
    elsif (($d eq "if")||($d eq "ifdef")||($d eq "ifndef")) {$depth--;}

    if ($depth==0) {
      my ($Guard)=split(' ',$Directives[$n][2]);

      if (($d eq "ifndef")&&(defined $Guard)&&($n<$#Directives)&&($Directives[$n+1][1] eq "define")&&($Directives[$n+1][2] eq $Guard)) {
        return $Directives[-1][0];
      }

      last;
    }
  }

  return scalar(@Lines);
}

1;
