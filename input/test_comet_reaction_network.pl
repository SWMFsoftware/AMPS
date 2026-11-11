#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use CometReactionNetwork;

# A small synthetic network covers all input laws and the explicit NONE sink
# without tying the parser test to a particular comet or AMPS build tree.
my @species=('H2O','OH','H','CO2');
my $network=CometReactionNetwork->new(species=>\@species);
$network->add_reaction(
    'H2O_PHOTO | H2O | OH + H | HELIOCENTRIC | 1.03e-5 | 0 | 0 | 1 | 3.41',
    'unit test');
$network->add_reaction(
    'OH_LOSS | 2*OH | H2O + H | ARRHENIUS | 2.0e-12 | 0.5 | 100 | 1 | 0',
    'unit test');
$network->add_reaction(
    'H2O_IONIZATION | H2O | NONE | HELIOCENTRIC | 4.04e-7 | 0 | 0 | 1 | 0',
    'unit test');
die "unexpected reaction count\n" unless $network->reaction_count()==3;

my $output='/tmp/CometChemistry.generated.test.h';
$network->write_header($output);
open(my $in,'<',$output) or die "cannot open '$output': $!\n";
my $text=do {local $/; <$in>};
close($in);
die "generated table is missing the photolysis reaction\n"
    unless $text=~/H2O_PHOTO/;
die "generated table is missing the bimolecular reaction\n"
    unless $text=~/OH_LOSS/;
die "generated table is missing the zero-product loss reaction\n"
    unless $text=~/H2O_IONIZATION.*?,1,.*?,0,/s;

my $failed=0;
eval {
  $network->add_reaction(
      'BAD | UNKNOWN | H | CONSTANT | 1 | 0 | 0 | 1 | 0','unit test');
};
$failed=1 if $@;
die "unknown species was accepted\n" unless $failed;

print "Comet reaction-network parser test passed\n";
