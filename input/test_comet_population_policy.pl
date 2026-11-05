#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use CometPopulationPolicy;
use CometReactionNetwork;

my @species=qw(H2O OH H C2H2 C2 CO2);
my $network=CometReactionNetwork->new(species=>\@species);
$network->add_reaction(
    'H2O_OH_H | H2O | OH + H | CONSTANT | 1 | 0 | 0 | 1 | 0',
    'population-policy test');
$network->add_reaction(
    'C2H2_C2_H2 | C2H2 | C2 + H | CONSTANT | 1 | 0 | 0 | 1 | 0',
    'population-policy test');
$network->add_reaction(
    'ZERO_UNUSED | CO2 | C2 | CONSTANT | 0 | 0 | 0 | 1 | 0',
    'population-policy test');
my $reachable=$network->reachable_species_from([0,3]);
die "transitive chemistry reachability is incorrect\n"
    unless $reachable->[0] && $reachable->[1] && $reachable->[2] &&
           $reachable->[3] && $reachable->[4] && !$reachable->[5];

my $policy=CometPopulationPolicy->new(species=>\@species);
$policy->add_range_record(
    'RadialParticleNumberPerCellRange=800 | 3200','global test');
$policy->add_lock_record(
    'RadialParticleWeightLock=DIRECT_SOURCE','direct group test');
$policy->add_range_record(
    'RadialParticleNumberPerCellRange(CHEMICAL_ONLY)=50 | 200',
    'chemical group test');
$policy->add_range_record(
    'RadialParticleNumberPerCellRange(C2)=100 | 400',
    'explicit species test');
$policy->add_lock_record('RadialParticleWeightLock=CO2','explicit lock test');
my $resolved=$policy->resolve(
    direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);

# Explicit species wins over the derived group, and a species outside the
# reaction closure still accepts an explicit lock.  The direct group wins over
# the global range for both directly supplied parents.
die "direct-source lock did not override the global range\n"
    unless $resolved->{locked}->[0] && $resolved->{locked}->[3] &&
        $resolved->{source}->[0]==CometPopulationPolicy::POLICY_DIRECT_SOURCE;
die "chemical-only range was not resolved\n"
    unless !$resolved->{locked}->[1] &&
        $resolved->{minimum}->[1]==50 && $resolved->{maximum}->[1]==200 &&
        $resolved->{source}->[1]==CometPopulationPolicy::POLICY_CHEMICAL_ONLY;
die "explicit species did not override the chemical-only group\n"
    unless $resolved->{minimum}->[4]==100 &&
        $resolved->{maximum}->[4]==400 &&
        $resolved->{source}->[4]==CometPopulationPolicy::POLICY_SPECIES;
die "explicit lock was not resolved\n"
    unless $resolved->{locked}->[5] &&
        $resolved->{source}->[5]==CometPopulationPolicy::POLICY_SPECIES;

sub must_fail {
  my ($description,$callback)=@_;
  my $failed=0;
  eval {$callback->();};
  $failed=1 if $@;
  die "$description was accepted\n" unless $failed;
}

must_fail('duplicate selector',sub {
  my $p=CometPopulationPolicy->new(species=>\@species);
  $p->add_range_record('RadialParticleNumberPerCellRange=10 | 20','one');
  $p->add_range_record('RadialParticleNumberPerCellRange=20 | 40','two');
});
must_fail('range/lock contradiction',sub {
  my $p=CometPopulationPolicy->new(species=>\@species);
  $p->add_range_record(
      'RadialParticleNumberPerCellRange(CHEMICAL_ONLY)=10 | 20','range');
  $p->add_lock_record('RadialParticleWeightLock=CHEMICAL_ONLY','lock');
});
must_fail('unknown selector',sub {
  my $p=CometPopulationPolicy->new(species=>\@species);
  $p->add_range_record(
      'RadialParticleNumberPerCellRange(UNKNOWN)=10 | 20','unknown');
});
must_fail('invalid range',sub {
  my $p=CometPopulationPolicy->new(species=>\@species);
  $p->add_range_record('RadialParticleNumberPerCellRange=10 | 19','range');
});
must_fail('global lock',sub {
  my $p=CometPopulationPolicy->new(species=>\@species);
  $p->add_lock_record('RadialParticleWeightLock=GLOBAL','lock');
});

print "Comet population-policy parser test passed\n";
