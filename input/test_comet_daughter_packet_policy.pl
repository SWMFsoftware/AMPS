#!/usr/bin/perl

use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use POSIX qw(LONG_MAX);
use lib $Bin;
use CometDaughterPacketPolicy;

my @species=qw(H2O OH H C2H2 C2 CO2);
my @chemical_only=(0,1,1,0,1,0);
my $reactions=[
  {id=>'H2O_OH_H',rate=>1.0,products=>[
      {species=>1,coefficient=>1},{species=>2,coefficient=>1}]},
  {id=>'C2H2_C2_2H',rate=>1.0,products=>[
      {species=>4,coefficient=>1},{species=>2,coefficient=>2}]},
  {id=>'GLOBAL_ONLY',rate=>1.0,products=>[
      {species=>0,coefficient=>1}]},
  {id=>'GROUP_ONLY',rate=>1.0,products=>[
      {species=>1,coefficient=>1}]},
  {id=>'CO2_SINK',rate=>1.0,products=>[]},
  {id=>'DISABLED',rate=>0.0,products=>[
      {species=>4,coefficient=>1}]},
];

sub must_fail {
  my ($description,$callback)=@_;
  my $failed=0;
  eval {$callback->();};
  $failed=1 if $@;
  die "$description was accepted\n" unless $failed;
}

# Omission is the compatibility path: all reactions retain one packet and no
# mode-specific hard limit is required.
my $default=CometDaughterPacketPolicy->new(species=>\@species);
my $default_resolved=$default->validate_configuration_and_resolve(
    reactions=>$reactions,chemical_only=>\@chemical_only,
    transaction_limit=>24);
die "default packet policy did not retain one packet\n"
    unless join(',',@{$default_resolved->{packet_count}}) eq '1,1,1,1,1,1';

my $policy=CometDaughterPacketPolicy->new(species=>\@species);
$policy->set_mode('on','mode');
$policy->set_maximum(8,'maximum');
$policy->add_count_record('ChemistryDaughterPacketCount=2','global');
$policy->add_count_record(
    'ChemistryDaughterPacketCount(CHEMICAL_ONLY)=3','group');
$policy->add_count_record(
    'ChemistryDaughterPacketCount(SPECIES:H)=4','H');
$policy->add_count_record(
    'ChemistryDaughterPacketCount(SPECIES:C2)=8','C2');
$policy->add_count_record(
    'ChemistryDaughterPacketCount(REACTION:H2O_OH_H)=5','reaction');
my $resolved=$policy->validate_configuration_and_resolve(
    reactions=>$reactions,chemical_only=>\@chemical_only,
    transaction_limit=>24);

# Reaction ID wins for H2O.  C2H2 has explicit H=4 and C2=8 requests, so the
# approved complete-packet rule chooses their maximum.  Sinks are always one;
# the zero-rate reaction inherits the species policy but is never executed.
die "reaction override was not applied\n"
    unless $resolved->{packet_count}->[0]==5 &&
        $resolved->{policy_source}->[0]==
            CometDaughterPacketPolicy::POLICY_REACTION;
die "maximum coproduct species request was not applied\n"
    unless $resolved->{packet_count}->[1]==8 &&
        $resolved->{policy_source}->[1]==
            CometDaughterPacketPolicy::POLICY_SPECIES;
die "global fallback was not applied\n"
    unless $resolved->{packet_count}->[2]==2 &&
        $resolved->{policy_source}->[2]==
            CometDaughterPacketPolicy::POLICY_GLOBAL;
die "chemical-only fallback was not applied\n"
    unless $resolved->{packet_count}->[3]==3 &&
        $resolved->{policy_source}->[3]==
            CometDaughterPacketPolicy::POLICY_CHEMICAL_ONLY;
die "sink reaction did not remain one packet\n"
    unless $resolved->{packet_count}->[4]==1;
die "transaction counts do not include repeated products\n"
    unless join(',',@{$resolved->{transaction_count}}) eq '10,24,2,3,0,8';

# Resolution is independent of record order.
my $reordered=CometDaughterPacketPolicy->new(species=>\@species);
$reordered->add_count_record(
    'ChemistryDaughterPacketCount(REACTION:H2O_OH_H)=5','reaction');
$reordered->add_count_record(
    'ChemistryDaughterPacketCount(SPECIES:C2)=8','C2');
$reordered->add_count_record(
    'ChemistryDaughterPacketCount(SPECIES:H)=4','H');
$reordered->add_count_record(
    'ChemistryDaughterPacketCount(CHEMICAL_ONLY)=3','group');
$reordered->add_count_record('ChemistryDaughterPacketCount=2','global');
$reordered->set_maximum(8,'maximum');
$reordered->set_mode('ON','mode');
my $again=$reordered->validate_configuration_and_resolve(
    reactions=>$reactions,chemical_only=>\@chemical_only,
    transaction_limit=>24);
for my $field (qw(packet_count policy_source selector transaction_count)) {
  die "packet policy depends on input order ($field)\n"
      unless join('|',@{$again->{$field}}) eq
          join('|',@{$resolved->{$field}});
}

# Header emission preserves reaction order and has a zero-network sentinel.
my $directory=tempdir(CLEANUP=>1);
my $header="$directory/policy.generated.h";
$policy->write_header($header,$resolved);
open(my $in,'<',$header) or die "cannot read generated packet header: $!\n";
my $text=join('',<$in>);
close($in);
die "generated packet header has incorrect count/order\n"
    unless $text=~/#include <cstddef>/ && $text=~/ReactionCount=6;/ &&
        $text=~/PacketCount\[\]=\{5,8,2,3,1,8\};/ &&
        $text=~/PolicySource\[\]=\{4,3,1,2,0,3\};/;

my $empty=CometDaughterPacketPolicy->new(species=>\@species);
my $empty_resolved=$empty->validate_configuration_and_resolve(
    reactions=>[],chemical_only=>\@chemical_only,transaction_limit=>1);
my $empty_header="$directory/empty.generated.h";
$empty->write_header($empty_header,$empty_resolved);
open(my $empty_in,'<',$empty_header) or die "cannot read empty header: $!\n";
my $empty_text=join('',<$empty_in>);
close($empty_in);
die "empty generated packet header has no sentinel\n"
    unless $empty_text=~/ReactionCount=0;/ &&
        $empty_text=~/PacketCount\[\]=\{1\};/ &&
        $empty_text=~/PolicySource\[\]=\{0\};/;

must_fail('invalid mode',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('maybe','mode');
});
must_fail('duplicate mode',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','one'); $p->set_mode('off','two');
});
for my $invalid (qw(0 1 -2 2.5 Inf)) {
  must_fail("invalid maximum $invalid",sub {
    my $p=CometDaughterPacketPolicy->new(species=>\@species);
    $p->set_maximum($invalid,'maximum');
  });
}
for my $invalid ('0','-1','2.5','Inf') {
  must_fail("invalid count $invalid",sub {
    my $p=CometDaughterPacketPolicy->new(species=>\@species);
    $p->add_count_record("ChemistryDaughterPacketCount=$invalid",'count');
  });
}
must_fail('unknown species selector',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->add_count_record(
      'ChemistryDaughterPacketCount(SPECIES:BAD)=2','species');
});
for my $selector ('BAD:X','SPECIES:','REACTION:') {
  must_fail("malformed selector $selector",sub {
    my $p=CometDaughterPacketPolicy->new(species=>\@species);
    $p->add_count_record(
        "ChemistryDaughterPacketCount($selector)=2",'selector');
  });
}
for my $case (
    ['global','ChemistryDaughterPacketCount=2'],
    ['group','ChemistryDaughterPacketCount(CHEMICAL_ONLY)=2'],
    ['species','ChemistryDaughterPacketCount(SPECIES:C2)=2'],
    ['reaction','ChemistryDaughterPacketCount(REACTION:H2O_OH_H)=2']) {
  must_fail("duplicate $case->[0] selector",sub {
    my $p=CometDaughterPacketPolicy->new(species=>\@species);
    $p->add_count_record($case->[1],'one');
    $p->add_count_record($case->[1],'two');
  });
}
must_fail('unknown reaction selector',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(4,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(REACTION:BAD)=2','reaction');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('record with disabled mode',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->add_count_record('ChemistryDaughterPacketCount=2','global');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('enabled mode without maximum',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->add_count_record('ChemistryDaughterPacketCount=2','global');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('enabled mode without effective split',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record('ChemistryDaughterPacketCount=1','global');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('packet count above maximum',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(4,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(SPECIES:C2)=8','species');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('transaction above hard limit',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(SPECIES:C2)=8','species');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>23);
});
must_fail('transaction integer overflow',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(LONG_MAX,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(SPECIES:C2)='.LONG_MAX,'species');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>LONG_MAX);
});
must_fail('sink reaction split',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(REACTION:CO2_SINK)=2','sink');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('sink reaction explicit identity policy',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record('ChemistryDaughterPacketCount(GLOBAL)=2','global');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(REACTION:CO2_SINK)=1','sink');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('zero-rate reaction split',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(REACTION:DISABLED)=2','disabled');
  $p->validate_configuration_and_resolve(reactions=>$reactions,
      chemical_only=>\@chemical_only,transaction_limit=>24);
});
must_fail('zero-rate inherited rule is not an effective split',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->set_mode('on','mode'); $p->set_maximum(8,'maximum');
  $p->add_count_record(
      'ChemistryDaughterPacketCount(SPECIES:C2)=2','species');
  $p->validate_configuration_and_resolve(
      reactions=>[$reactions->[5]],chemical_only=>\@chemical_only,
      transaction_limit=>24);
});

# A reaction-level N=1 record intentionally suppresses a broader request.
my $suppression=CometDaughterPacketPolicy->new(species=>\@species);
$suppression->set_mode('on','mode');
$suppression->set_maximum(8,'maximum');
$suppression->add_count_record('ChemistryDaughterPacketCount(GLOBAL)=2',
    'global');
$suppression->add_count_record(
    'ChemistryDaughterPacketCount(REACTION:GLOBAL_ONLY)=1','suppression');
my $suppressed=$suppression->validate_configuration_and_resolve(
    reactions=>$reactions,chemical_only=>\@chemical_only,
    transaction_limit=>24);
die "reaction identity rule did not suppress the broader policy\n"
    unless $suppressed->{packet_count}->[2]==1 &&
        $suppressed->{policy_source}->[2]==
            CometDaughterPacketPolicy::POLICY_REACTION;

for my $invalid_reactions (
    [{id=>'ZERO_COEFFICIENT',rate=>1.0,
      products=>[{species=>0,coefficient=>0}]}],
    [{id=>'FRACTIONAL_COEFFICIENT',rate=>1.0,
      products=>[{species=>0,coefficient=>1.5}]}],
    [{id=>'BAD_SPECIES_INDEX',rate=>1.0,
      products=>[{species=>99,coefficient=>1}]}]) {
  must_fail('invalid reaction product descriptor',sub {
    my $p=CometDaughterPacketPolicy->new(species=>\@species);
    $p->validate_configuration_and_resolve(reactions=>$invalid_reactions,
        chemical_only=>\@chemical_only,transaction_limit=>24);
  });
}

must_fail('invalid generated header packet count',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->write_header("$directory/invalid-count.h",{
      packet_count=>[0],policy_source=>[0]});
});
must_fail('invalid generated header policy source',sub {
  my $p=CometDaughterPacketPolicy->new(species=>\@species);
  $p->write_header("$directory/invalid-source.h",{
      packet_count=>[1],policy_source=>[99]});
});

print "Comet daughter-packet policy test passed\n";
