#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use CometPhysicsDatabase;
use CometReactionNetwork;

# Use the exact 27 neutral species represented in Rubin et al. (2011).  This
# test exercises name-to-index resolution independently of the AMPS build
# preprocessor, which makes a missing or misspelled catalog species fail fast.
my @species=qw(H H2 C CH CH2 N CH3 NH CH4 NH2 O NH3 OH H2O C2 C2H
    C2H2 CN HCN C2H4 CO C2H5 HCO C2H6 H2CO CH3OH CO2);
my $network=CometReactionNetwork->new(species=>\@species);
my $physics=CometPhysicsDatabase->new(species=>\@species);

# The lower published OH endpoint is a deterministic parser-test value; the
# production Halley input must state its own selection explicitly.
$physics->load_rubin2011(
    network=>$network,
    oh_ground_rate=>3.0e-6,
    methanol_h2co_branch=>'ON');
die "Rubin2011 parent collision table should contain 66 unordered pairs\n"
    unless $physics->collision_count()==66;

my $chemistry='/tmp/CometChemistry.rubin2011.generated.test.h';
my $collisions='/tmp/CometPhysics.rubin2011.generated.test.h';
$network->write_header($chemistry);
$physics->write_collision_header($collisions);

open(my $cin,'<',$chemistry) or die "cannot open '$chemistry': $!\n";
my $chemistry_text=do {local $/; <$cin>};
close($cin);
die "Rubin2011 chemistry is missing the methanol sensitivity branch\n"
    unless $chemistry_text=~/CH3OH_H2CO_H2/;
die "Rubin2011 chemistry is missing a photoionization sink\n"
    unless $chemistry_text=~/H2O_ION/;
die "Rubin2011 chemistry is missing the selected OH endpoint\n"
    unless $chemistry_text=~/OH_O_H_GROUND/;

open(my $pin,'<',$collisions) or die "cannot open '$collisions': $!\n";
my $collision_text=do {local $/; <$pin>};
close($pin);
die "Rubin2011 collision header has the wrong database name\n"
    unless $collision_text=~/DatabaseName\[\]=\"RUBIN2011\"/;
die "Rubin2011 collision header is missing the Crifo temperature law\n"
    unless $collision_text=~/TemperaturePowerLaw/;
die "Rubin2011 collision header is missing the published 50--300 K range\n"
    unless $collision_text=~/1\.66.*?300.*?-0\.599999.*?50.*?300/s;

# The Halley production policy preserves all 66 Rubin rows, adds the 18
# non-overlapping Combi rows, and estimates only the remaining 294 pairs from
# CHEMKIN diameters.  The generated source strings make those tiers auditable.
my $hierarchy_network=CometReactionNetwork->new(species=>\@species);
my $hierarchy=CometPhysicsDatabase->new(species=>\@species);
$hierarchy->load_rubin2011(
    network=>$hierarchy_network,
    oh_ground_rate=>3.0e-6,
    methanol_h2co_branch=>'OFF');
$hierarchy->complete_with_combi1996_chemkin_lorentz_hs();
die "completed collision table should contain all 378 unordered pairs\n"
    unless $hierarchy->collision_count()==378;
my $hierarchy_header='/tmp/CometPhysics.hierarchy.generated.test.h';
$hierarchy->write_collision_header($hierarchy_header);
open(my $hierarchy_in,'<',$hierarchy_header)
    or die "cannot open '$hierarchy_header': $!\n";
my $hierarchy_text=do {local $/; <$hierarchy_in>};
close($hierarchy_in);
die "completed collision header has the wrong database name\n"
    unless $hierarchy_text=~/DatabaseName\[\]="RUBIN2011\+COMBI1996\+CHEMKIN_LORENTZ_HS"/;
my $combi_count=()=$hierarchy_text=~/Combi \(1996\), reproduced in Fougere dissertation Table 2\.1/g;
die "completed table should contain 18 non-overlapping Combi rows, found $combi_count\n"
    unless $combi_count==18;
my $estimate_count=()=$hierarchy_text=~/CHEMKIN Table 5\.1 diameter with Lorentz mixing; geometric hard-sphere estimate/g;
die "completed table should contain 294 CHEMKIN/Lorentz estimates, found $estimate_count\n"
    unless $estimate_count==294;
die "completion replaced the Rubin/Crifo H2O self-collision law\n"
    unless $hierarchy_text=~/TemperaturePowerLaw.*?Crifo \(1989\), as listed in Rubin et al\. \(2011\), Table 3/s;

# Check one row from each priority tier numerically. Species indices follow
# @species above: H=0, C=2, O=10, OH=12, H2O=13, and CO=20.
sub find_collision {
  my ($database,$species0,$species1)=@_;
  for my $entry (@{$database->{collisions}}) {
    return $entry if ($entry->{species0}==$species0 &&
                      $entry->{species1}==$species1) ||
                     ($entry->{species0}==$species1 &&
                      $entry->{species1}==$species0);
  }
  return undef;
}
sub relatively_equal {
  my ($a,$b)=@_;
  my $scale=abs($a)>abs($b) ? abs($a) : abs($b);
  return abs($a-$b)<=1.0e-12*$scale;
}

my $rubin_h2o_co=find_collision($hierarchy,13,20);
die "Rubin H2O-CO priority row was not preserved\n"
    unless defined $rubin_h2o_co &&
        $rubin_h2o_co->{source}=~/Rubin et al\./ &&
        relatively_equal($rubin_h2o_co->{sigma},3.2e-19);
my $combi_h_o=find_collision($hierarchy,0,10);
die "Combi H-O row is missing or incorrect\n"
    unless defined $combi_h_o && $combi_h_o->{source}=~/Combi \(1996\)/ &&
        relatively_equal($combi_h_o->{sigma},1.2e-19);
my $chemkin_c_c=find_collision($hierarchy,2,2);
my $expected_c_c=3.14159265358979323846*3.298*3.298*1.0e-20;
die "CHEMKIN/Lorentz C-C estimate is missing or incorrect\n"
    unless defined $chemkin_c_c &&
        $chemkin_c_c->{source}=~/CHEMKIN Table 5\.1/ &&
        relatively_equal($chemkin_c_c->{sigma},$expected_c_c);

# A second completion request is an input error rather than a way to append
# duplicate estimates.
my $completion_failed=0;
eval {$hierarchy->complete_with_combi1996_chemkin_lorentz_hs();};
$completion_failed=1 if $@;
die "collision completion accepted a duplicate request\n"
    unless $completion_failed;

# Table 4's OH interval is intentionally not defaulted by the database.
my $failed=0;
eval {
  my $n=CometReactionNetwork->new(species=>\@species);
  my $p=CometPhysicsDatabase->new(species=>\@species);
  $p->load_rubin2011(network=>$n,methanol_h2co_branch=>'ON');
};
$failed=1 if $@;
die "Rubin2011 accepted an unspecified OH interval value\n" unless $failed;

# Explicit collision records carry comet-specific or supplemental values in
# input rather than C++. Check both supported forms and their generated name.
my $custom=CometPhysicsDatabase->new(species=>\@species);
$custom->add_collision(
    'H2O | CO2 | CONSTANT | 3.4e-19 | parser test','unit test');
$custom->add_collision(
    'H2O | H2O | TEMPERATURE_POWER | 1.66e-19 | 300 | -0.6 | 50 | 300 | parser test',
    'unit test');
my $custom_header='/tmp/CometPhysics.custom.generated.test.h';
$custom->write_collision_header($custom_header);
open(my $custom_in,'<',$custom_header)
    or die "cannot open '$custom_header': $!\n";
my $custom_text=do {local $/; <$custom_in>};
close($custom_in);
die "custom collision header has the wrong database name\n"
    unless $custom_text=~/DatabaseName\[\]=\"CUSTOM\"/;
die "custom collision header has the wrong number of entries\n"
    unless $custom_text=~/CollisionCount=2/;

print "Comet physics database test passed\n";
