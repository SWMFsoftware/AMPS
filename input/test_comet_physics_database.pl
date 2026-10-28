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
