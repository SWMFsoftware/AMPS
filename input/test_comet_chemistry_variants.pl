#!/usr/bin/env perl
# Tests of the chemistry-variant records of input/CometReactionNetwork.pm
# (srcHartley2/CHEMISTRY_VARIANTS.md, section 5).  The ChemistryReaction,
# ChemistryGroup and ChemistryRateTable records of input/i3.conf must be
# accepted and produce the expected generated tables; each listed defect must
# be rejected with its own message; and records without the optional fields
# must keep generating the reaction table of the build before variants
# existed (the first reactions and terms unchanged, new ones appended).
use strict;
use warnings;
use File::Basename qw(dirname);
use File::Temp qw(tempdir);
use lib dirname(__FILE__);
use CometReactionNetwork;

my $input_dir=dirname(__FILE__);
my @species=qw(H H2 O OH H2O CO CO2);
my ($pass,$fail)=(0,0);

sub check {
  my ($ok,$name)=@_;
  printf "%s %s\n",$ok ? 'PASS' : 'FAIL',$name;
  $ok ? $pass++ : $fail++;
}

sub build {
  my (@lines)=@_;
  my $network=CometReactionNetwork->new(species=>\@species);
  for my $line (@lines) {
    (my $value=$line)=~s/!.*$//;
    $value=~s/^\s+//;
    if ($value=~s/^ChemistryReaction\s*=\s*//i) {$network->add_reaction($value,'test');}
    elsif ($value=~s/^ChemistryGroup\s*=\s*//i) {$network->add_group($value,'test');}
    elsif ($value=~s/^ChemistryRateTable\s*=\s*//i) {$network->add_rate_table($value,'test');}
  }
  return $network;
}

sub header_section {
  my ($text,$start)=@_;
  my ($section)=$text=~/\Q$start\E[^\n]*\n(.*?)\n\};/s;
  return [split /\n/,$section];
}

open(my $in,'<',"$input_dir/i3.conf") or die "cannot read i3.conf: $!\n";
my @records=grep {/^\s*Chemistry(Reaction|Group|RateTable)\s*=/} <$in>;
close($in);

my $dir=tempdir(CLEANUP=>1);
my $network=build(@records);
$network->write_header("$dir/full.h");
open(my $full,'<',"$dir/full.h") or die;
my $full_text=do {local $/; <$full>};
close($full);
check($network->has_table_rates(),'i3.conf has table-rate reactions');
check($full_text=~/static const std::size_t GroupCount=1;/,'one group generated');
check($full_text=~/\{"OH_LOSS",0,2,GroupVariants0\}/,'OH_LOSS default STATE_RESOLVED');
check($full_text=~/\{"OH_TAU_SA1982","SA1982_OH_LIFETIME",83300,0.05/,'rate table record');
my $variants=header_section($full_text,'Comet::Chemistry::cCompiledReactionVariant ReactionVariants[]');
check(@$variants==19,'one variant row per reaction');
check($variants->[3]=~/\{0,0,-1,0\}/ && $variants->[12]=~/\{0,0,-1,0\}/,
    'state-resolved OH rows');
check($variants->[17]=~/\{0,1,0,1\}/ && $variants->[18]=~/\{0,1,0,2\}/,
    'SA1982_VH rows: dissociation and ionization roles');

# Without the variant records the same reactions must generate exactly the
# leading part of the new tables.
my @legacy=map {s/ \| GROUP=OH_LOSS:STATE_RESOLVED//r}
    grep {!/Chemistry(Group|RateTable)|SA1982/} @records;
my $old=build(@legacy);
$old->write_header("$dir/legacy.h");
open(my $legacy_in,'<',"$dir/legacy.h") or die;
my $legacy_text=do {local $/; <$legacy_in>};
close($legacy_in);
for my $start ('ChemicalSolver::cStoichiometricTerm Terms[]',
    'ChemicalSolver::cReaction Reactions[]') {
  my ($a,$b)=(header_section($legacy_text,$start),header_section($full_text,$start));
  check(join("\n",@$a) eq join("\n",@{$b}[0..$#$a]),"unchanged prefix of $start");
}

my @defects=(
  ['variant without reaction',sub {grep {!/SA1982 \|/ && !/OH_O_H_SA1982|OH_ION_SA1982/} @_},
      qr/has no reaction|not used by any reaction/],
  ['undeclared group',sub {map {s/GROUP=OH_LOSS:SA1982_VH/GROUP=OHX:SA1982_VH/r} @_},
      qr/undeclared ChemistryGroup 'OHX'/],
  ['unlisted variant',sub {map {s/OH_LOSS:SA1982_VH/OH_LOSS:FOO/r} @_},qr/variant 'FOO'/],
  ['default not a variant',sub {map {s/DEFAULT=STATE_RESOLVED/DEFAULT=BAR/r} @_},
      qr/is not one of its VARIANTS/],
  ['undeclared table',sub {map {s/TABLE:OH_TAU_SA1982:DISSOCIATION/TABLE:NOPE:DISSOCIATION/r} @_},
      qr/undeclared ChemistryRateTable 'NOPE'/],
  ['unknown role',sub {map {s/:DISSOCIATION/:FOO/r} @_},qr/role 'FOO' is not defined/],
  ['ionization with products',sub {map {s/OH_ION_SA1982 \| OH \| NONE/OH_ION_SA1982 | OH | O + H/r} @_},
      qr/IONIZATION table rate requires PRODUCTS=NONE/],
  ['table rate with RREF_AU=2',sub {map {/OH_O_H_SA1982/ ?
      s/\| 0 \| 0 \| 0 \| 1 \| 0.36/| 0 | 0 | 0 | 2 | 0.36/r : $_} @_},
      qr/requires LAW=HELIOCENTRIC and RREF_AU=1/],
  ['two parents in one group',sub {map {/OH_ION_SA1982/ ? s/\| OH \| NONE/| H2O | NONE/r : $_} @_},
      qr/must share one parent species/],
  ['unknown reaction option',sub {map {/OH_O_H_SA1982/ ? s/RATE=/RAT=/r : $_} @_},
      qr/unknown option 'RAT'/],
  ['unknown table model',sub {map {s/MODEL=SA1982_OH_LIFETIME/MODEL=XYZ/r} @_},qr/unknown MODEL=XYZ/],
  ['negative SCALE_S',sub {map {s/SCALE_S=0.833e5/SCALE_S=-1/r} @_},qr/SCALE_S .* must be positive/],
  ['unused table',sub {(@_,"ChemistryRateTable=SPARE | MODEL=SA1982_OH_LIFETIME | SCALE_S=1 | ION_FRACTION=0\n")},
      qr/'SPARE' .* not used/],
  ['duplicate group',sub {(@_,"ChemistryGroup=OH_LOSS | DEFAULT=A | VARIANTS=A,B\n")},
      qr/duplicate ChemistryGroup 'OH_LOSS'/],
  ['single-variant group',sub {map {s/VARIANTS=STATE_RESOLVED,SA1982_VH/VARIANTS=STATE_RESOLVED/r} @_},
      qr/at least two variants/],
);
for my $defect (@defects) {
  my ($name,$mutate,$pattern)=@$defect;
  my $error='';
  eval {my $bad=build($mutate->(@records)); $bad->write_header("$dir/bad.h"); 1} or $error=$@;
  check($error=~$pattern,"rejects $name".($error eq '' ? ' (accepted!)' : ''));
}

printf "%s test_comet_chemistry_variants passed=%d failed=%d\n",
    $fail==0 ? 'PASS' : 'FAIL',$pass,$fail;
exit($fail==0 ? 0 : 1);
