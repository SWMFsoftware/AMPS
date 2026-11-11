#!/usr/bin/perl

use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib $Bin;
use CometTotalReactionForcingPolicy;

my @species=qw(H2O OH H C2 C CO2);
my $reactions=[
  {id=>'H2O_OH_H',rate=>9.0,reactants=>[
      {species=>0,coefficient=>1}]},
  {id=>'H2O_ION',rate=>1.0,reactants=>[
      {species=>0,coefficient=>1}]},
  {id=>'C2_2C',rate=>1.0e-7,reactants=>[
      {species=>3,coefficient=>1}]},
  {id=>'C2_ION',rate=>9.0e-7,reactants=>[
      {species=>3,coefficient=>1}]},
  {id=>'CO2_DISABLED',rate=>0.0,reactants=>[
      {species=>5,coefficient=>1}]},
  {id=>'OH_BIMOLECULAR',rate=>1.0,reactants=>[
      {species=>1,coefficient=>1},{species=>2,coefficient=>1}]},
  {id=>'H_NONUNIT',rate=>1.0,reactants=>[
      {species=>2,coefficient=>2}]},
];

sub must_fail {
  my ($description,$callback,$pattern)=@_;
  my $error='';
  eval {$callback->();};
  $error=$@;
  die "$description was accepted\n" if $error eq '';
  die "$description produced unexpected error: $error\n"
      if defined($pattern) && $error!~$pattern;
}

# Omitted mode and explicit OFF are strict default-off identity policies.
for my $explicit_off (0,1) {
  my $policy=CometTotalReactionForcingPolicy->new(species=>\@species);
  $policy->set_mode('off','mode') if $explicit_off;
  my $resolved=$policy->validate_configuration_and_resolve(
      reactions=>$reactions);
  die "off total-reaction policy configured a species\n"
      if grep {$_!=0} @{$resolved->{configured}};
}

# Stage 8 is parent-species based: naming C2 enables its complete physical
# total rate (neutral and sink channels together) without changing either row.
my $policy=CometTotalReactionForcingPolicy->new(species=>\@species);
$policy->set_mode('ON','mode');
$policy->add_species_record(
    'ChemistryTotalReactionForcing(SPECIES:C2)=ON','C2 record');
my $resolved=$policy->validate_configuration_and_resolve(
    reactions=>$reactions);
die "C2 total-reaction policy flags are wrong\n"
    unless join(',',@{$resolved->{configured}}) eq '0,0,0,1,0,0';

# Input order is irrelevant, while generated order always follows species.
my $many=CometTotalReactionForcingPolicy->new(species=>\@species);
$many->add_species_record(
    'ChemistryTotalReactionForcing(SPECIES:C2):on','second');
$many->add_species_record(
    'ChemistryTotalReactionForcing(SPECIES:H2O)=ON','first');
$many->set_mode('on','mode');
my $many_resolved=$many->validate_configuration_and_resolve(
    reactions=>$reactions);
die "multiple total-reaction policy flags are wrong\n"
    unless join(',',@{$many_resolved->{configured}}) eq '1,0,0,1,0,0';

must_fail('records with omitted mode',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:C2)=ON','record');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/require ChemistryTotalReactionForcingMode=ON/);
must_fail('records with explicit off',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->set_mode('off','mode');
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:C2)=ON','record');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/require ChemistryTotalReactionForcingMode=ON/);
must_fail('on mode without records',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/requires at least one/);
must_fail('duplicate mode',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->set_mode('on','first'); $p->set_mode('off','second');
},qr/duplicate/);
must_fail('duplicate species record',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:C2)=ON','first');
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:C2)=ON','second');
},qr/duplicate/);
must_fail('case-mismatched species',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:c2)=ON','case');
},qr/unknown/);
must_fail('unknown species',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:UNKNOWN)=ON','unknown');
},qr/unknown/);
must_fail('per-species OFF record',sub {
  my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
  $p->add_species_record(
      'ChemistryTotalReactionForcing(SPECIES:C2)=OFF','off');
},qr/must be ON/);
for my $species_name (qw(CO2 OH H)) {
  must_fail("ineligible species $species_name",sub {
    my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
    $p->set_mode('on','mode');
    $p->add_species_record(
        "ChemistryTotalReactionForcing(SPECIES:$species_name)=ON",'row');
    $p->validate_configuration_and_resolve(reactions=>$reactions);
  },qr/no positive-rate unit-unimolecular reaction/);
}
for my $line (
    'ChemistryTotalReactionForcing=ON',
    'ChemistryTotalReactionForcing(REACTION:C2_2C)=ON',
    'ChemistryTotalReactionForcing(SPECIES:)=ON',
    'ChemistryTotalReactionForcing(SPECIES:C2)=MAYBE') {
  must_fail("malformed record '$line'",sub {
    my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
    $p->add_species_record($line,'syntax');
  });
}

# The descriptor adapter is a checked API boundary.  Cover duplicate IDs,
# invalid rates/species/coefficients, and a missing reactant array.
for my $case (
    [{id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>1}]},
     {id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>1}]}],
    [{id=>'X',rate=>-1,reactants=>[{species=>0,coefficient=>1}]}],
    [{id=>'X',rate=>1,reactants=>[{species=>99,coefficient=>1}]}],
    [{id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>0}]}],
    [{id=>'X',rate=>1}]) {
  must_fail('invalid reaction descriptor',sub {
    my $p=CometTotalReactionForcingPolicy->new(species=>\@species);
    $p->validate_configuration_and_resolve(reactions=>$case);
  });
}

# Generated output is self-contained and preserves species-table order.
my $directory=tempdir(CLEANUP=>1);
my $header="$directory/CometTotalReactionForcingPolicy.generated.h";
$many->write_header($header,$many_resolved);
open(my $stream,'<',$header) or die "cannot read generated header: $!\n";
my $text=do {local $/; <$stream>};
close($stream);
die "generated header is not self-contained\n"
    unless $text=~/#include <cstddef>/;
die "generated species count is wrong\n"
    unless $text=~/SpeciesCount=6;/;
die "generated forcing flags are wrong\n"
    unless $text=~/Configured\[\]=\{1,0,0,1,0,0\};/;

must_fail('wrong generated species count',sub {
  $many->write_header("$directory/bad-length.h",{configured=>[1]});
},qr/species count/);
must_fail('invalid generated flag',sub {
  $many->write_header("$directory/bad-flag.h",
      {configured=>[0,0,0,2,0,0]});
},qr/invalid generated/);

# Empty species tables still need a legal C++ sentinel while retaining the
# semantic SpeciesCount=0 contract.
my $empty=CometTotalReactionForcingPolicy->new(species=>[]);
my $empty_resolved=$empty->validate_configuration_and_resolve(reactions=>[]);
my $empty_header="$directory/empty.h";
$empty->write_header($empty_header,$empty_resolved);
open(my $empty_stream,'<',$empty_header) or die "cannot read empty header: $!\n";
my $empty_text=do {local $/; <$empty_stream>};
close($empty_stream);
die "empty generated header lacks a safe sentinel\n"
    unless $empty_text=~/SpeciesCount=0;/ &&
        $empty_text=~/Configured\[\]=\{0\};/;

print "Comet total-reaction-forcing policy test passed\n";
