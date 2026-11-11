#!/usr/bin/perl

use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib $Bin;
use CometChannelImportancePolicy;
use CometReactionNetwork;

my @species=qw(H2O OH H C2 C CO2);
my $reactions=[
  {id=>'H2O_OH_H',rate=>9.0,reactants=>[
      {species=>0,coefficient=>1}]},
  {id=>'H2O_ION',rate=>1.0,reactants=>[
      {species=>0,coefficient=>1}]},
  {id=>'C2_2C',rate=>1.0e-7,reactants=>[
      {species=>3,coefficient=>1}]},
  # A product-free sink is deliberately indistinguishable here: Stage 7
  # groups channels by the tracked parent and must include the sink in q*.
  {id=>'C2_ION',rate=>9.1e-7,reactants=>[
      {species=>3,coefficient=>1}]},
  {id=>'CO2_ONLY',rate=>2.0,reactants=>[
      {species=>5,coefficient=>1}]},
  {id=>'DISABLED',rate=>0.0,reactants=>[
      {species=>5,coefficient=>1}]},
  {id=>'BIMOLECULAR',rate=>1.0,reactants=>[
      {species=>1,coefficient=>1},{species=>2,coefficient=>1}]},
  {id=>'NONUNIT',rate=>1.0,reactants=>[
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

sub add_c2_family {
  my ($policy,$first,$second)=@_;
  $policy->add_probability_record(
      "ChemistryChannelBiasedProbability(REACTION:C2_2C)=$first",'C2 neutral');
  $policy->add_probability_record(
      "ChemistryChannelBiasedProbability(REACTION:C2_ION)=$second",'C2 sink');
}

# Omission and explicit OFF are strict identity paths.  They resolve every
# generated reaction to Configured=0/q*=0 without requiring any mode record.
for my $explicit_off (0,1) {
  my $policy=CometChannelImportancePolicy->new(species=>\@species);
  $policy->set_mode('off','mode') if $explicit_off;
  my $resolved=$policy->validate_configuration_and_resolve(
      reactions=>$reactions);
  die "off policy configured a reaction\n"
      unless join(',',@{$resolved->{configured}}) eq ('0,'x7).'0';
  die "off policy generated a nonzero probability\n"
      unless !grep {$_!=0.0} @{$resolved->{biased_probability}};
}

# A complete family may include a product-free sink.  Stable IDs are retained
# exactly and untouched H2O/CO2 families remain analog.
my $policy=CometChannelImportancePolicy->new(species=>\@species);
$policy->set_mode('ON','mode');
add_c2_family($policy,0.5,0.5);
my $resolved=$policy->validate_configuration_and_resolve(
    reactions=>$reactions);
die "C2 family flags are wrong\n"
    unless join(',',@{$resolved->{configured}}) eq '0,0,1,1,0,0,0,0';
die "C2 family probabilities are wrong\n"
    unless join(',',@{$resolved->{biased_probability}}) eq
        '0,0,0.5,0.5,0,0,0,0';

# Record order cannot change the reaction-indexed generated policy.
my $reordered=CometChannelImportancePolicy->new(species=>\@species);
$reordered->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:C2_ION)=0.5','sink');
$reordered->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5','neutral');
$reordered->set_mode('on','mode');
my $again=$reordered->validate_configuration_and_resolve(
    reactions=>$reactions);
die "record order changed the resolved policy\n"
    unless join(',',@{$again->{configured}}) eq
            join(',',@{$resolved->{configured}}) &&
        join(',',@{$again->{biased_probability}}) eq
            join(',',@{$resolved->{biased_probability}});

# A single positive channel has the normalized identity policy q*=1.  A
# disabled sibling is excluded from the complete positive-rate family.
my $single=CometChannelImportancePolicy->new(species=>\@species);
$single->set_mode('on','mode');
$single->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:CO2_ONLY)=1','single');
my $single_resolved=$single->validate_configuration_and_resolve(
    reactions=>$reactions);
die "single-channel identity policy failed\n"
    unless $single_resolved->{configured}->[4]==1 &&
        $single_resolved->{biased_probability}->[4]==1.0 &&
        $single_resolved->{configured}->[5]==0;

# A many-channel family exercises normalization beyond the two-channel case
# and accepts scientific notation and leading-dot input forms.
my $many_reactions=[
  {id=>'A',rate=>1.0,reactants=>[{species=>1,coefficient=>1}]},
  {id=>'B',rate=>2.0,reactants=>[{species=>1,coefficient=>1}]},
  {id=>'C',rate=>3.0,reactants=>[{species=>1,coefficient=>1}]},
];
my $many=CometChannelImportancePolicy->new(species=>\@species);
$many->set_mode('on','mode');
$many->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:A)=1e-12','tiny');
$many->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:B)=.25','middle');
$many->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:C)=0.749999999999','dominant');
my $many_resolved=$many->validate_configuration_and_resolve(
    reactions=>$many_reactions);
die "many-channel policy did not configure every sibling\n"
    unless join(',',@{$many_resolved->{configured}}) eq '1,1,1';

# Strict mode, selector, ID, family, rate, molecularity, probability, and
# normalization failures.  Each case uses a new object so no earlier failure
# can mask a later parser condition.
must_fail('records with omitted mode',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  add_c2_family($p,0.5,0.5);
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/require ChemistryChannelImportanceMode=ON/);
must_fail('records with explicit off',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('off','mode'); add_c2_family($p,0.5,0.5);
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/require ChemistryChannelImportanceMode=ON/);
must_fail('on mode without records',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/requires at least one/);
must_fail('duplicate mode',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','first'); $p->set_mode('off','second');
},qr/duplicate/);
must_fail('duplicate reaction record',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5','first');
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5','second');
},qr/duplicate/);
must_fail('case-mismatched reaction ID',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:c2_2c)=1','case');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/unknown/);
must_fail('unknown reaction ID',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:UNKNOWN)=1','unknown');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/unknown/);
must_fail('partial sibling family',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:C2_2C)=1','partial');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/missing reaction 'C2_ION'/);
must_fail('zero-rate target',sub {
  my $p=CometChannelImportancePolicy->new(species=>\@species);
  $p->set_mode('on','mode');
  $p->add_probability_record(
      'ChemistryChannelBiasedProbability(REACTION:DISABLED)=1','disabled');
  $p->validate_configuration_and_resolve(reactions=>$reactions);
},qr/zero-rate/);
for my $id (qw(BIMOLECULAR NONUNIT)) {
  must_fail("non-unit target $id",sub {
    my $p=CometChannelImportancePolicy->new(species=>\@species);
    $p->set_mode('on','mode');
    $p->add_probability_record(
        "ChemistryChannelBiasedProbability(REACTION:$id)=1",'nonunit');
    $p->validate_configuration_and_resolve(reactions=>$reactions);
  },qr/non-unit-unimolecular/);
}
for my $value ('0','-0.1','1.0001','NaN','Inf','') {
  must_fail("invalid probability '$value'",sub {
    my $p=CometChannelImportancePolicy->new(species=>\@species);
    $p->add_probability_record(
        "ChemistryChannelBiasedProbability(REACTION:C2_2C)=$value",'value');
  },qr/positive finite number/);
}
for my $pair ([0.4,0.5],[0.6,0.5]) {
  must_fail("non-normalized sum @$pair",sub {
    my $p=CometChannelImportancePolicy->new(species=>\@species);
    $p->set_mode('on','mode'); add_c2_family($p,@$pair);
    $p->validate_configuration_and_resolve(reactions=>$reactions);
  },qr/sum to/);
}
for my $line (
    'ChemistryChannelBiasedProbability=0.5',
    'ChemistryChannelBiasedProbability(SPECIES:C2)=0.5',
    'ChemistryChannelBiasedProbability(REACTION:)=0.5') {
  must_fail("malformed record '$line'",sub {
    my $p=CometChannelImportancePolicy->new(species=>\@species);
    $p->add_probability_record($line,'syntax');
  },qr/invalid/);
}

# Descriptor validation is a public boundary, not an assumption delegated to
# CometReactionNetwork.  Cover duplicate IDs, bad rate/species/coefficient,
# and a missing reactant array.
for my $case (
    [{id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>1}]},
     {id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>1}]}],
    [{id=>'X',rate=>-1,reactants=>[{species=>0,coefficient=>1}]}],
    [{id=>'X',rate=>1,reactants=>[{species=>99,coefficient=>1}]}],
    [{id=>'X',rate=>1,reactants=>[{species=>0,coefficient=>0}]}],
    [{id=>'X',rate=>1}]) {
  must_fail('invalid reaction descriptor',sub {
    my $p=CometChannelImportancePolicy->new(species=>\@species);
    $p->validate_configuration_and_resolve(reactions=>$case);
  });
}

# The generated header is self-contained and preserves reaction-table order.
my $directory=tempdir(CLEANUP=>1);
my $header="$directory/CometChannelImportancePolicy.generated.h";
$policy->write_header($header,$resolved);
open(my $stream,'<',$header) or die "cannot read generated header: $!\n";
my $text=do {local $/; <$stream>};
close($stream);
die "generated header is not self-contained\n"
    unless $text=~/#include <cstddef>/;
die "generated reaction count is wrong\n"
    unless $text=~/ReactionCount=8;/;
die "generated flags are wrong\n"
    unless $text=~/Configured\[\]=\{0,0,1,1,0,0,0,0\};/;
die "generated probabilities are wrong\n"
    unless $text=~/BiasedProbability\[\]=\{0,0,0\.5,0\.5,0,0,0,0\};/;

# Empty-network generation needs one valid C++ sentinel while retaining the
# semantic reaction count zero.
my $empty=CometChannelImportancePolicy->new(species=>\@species);
my $empty_resolved=$empty->validate_configuration_and_resolve(reactions=>[]);
my $empty_header="$directory/empty.h";
$empty->write_header($empty_header,$empty_resolved);
open(my $empty_stream,'<',$empty_header) or die "cannot read empty header: $!\n";
my $empty_text=do {local $/; <$empty_stream>};
close($empty_stream);
die "empty generated header lacks safe sentinels\n"
    unless $empty_text=~/ReactionCount=0;/ &&
        $empty_text=~/Configured\[\]=\{0\};/ &&
        $empty_text=~/BiasedProbability\[\]=\{0\};/;

# Header writer rejects malformed resolved data rather than emitting a C++
# table that only fails at application startup.
must_fail('mismatched generated arrays',sub {
  $policy->write_header("$directory/bad-length.h",{
      configured=>[1],biased_probability=>[]});
},qr/differ in length/);
must_fail('configured zero probability',sub {
  $policy->write_header("$directory/bad-value.h",{
      configured=>[1],biased_probability=>[0]});
},qr/invalid generated/);
must_fail('unconfigured nonzero probability',sub {
  $policy->write_header("$directory/bad-sentinel.h",{
      configured=>[0],biased_probability=>[0.5]});
},qr/invalid generated/);

# Exercise the production descriptor adapter itself, including the fact that
# it returns defensive copies rather than references into the reaction
# database.  Both a neutral channel and its NONE sink must remain siblings.
my $network=CometReactionNetwork->new(species=>\@species);
$network->add_reaction(
    'NET_C2_2C | C2 | 2*C | HELIOCENTRIC | 1.0e-7 | 0 | 0 | 1 | 3.6',
    'network neutral');
$network->add_reaction(
    'NET_C2_ION | C2 | NONE | HELIOCENTRIC | 9.0e-7 | 0 | 0 | 1 | 0',
    'network sink');
my $descriptors=$network->channel_importance_policy_reactions();
die "reaction-network importance descriptors are incomplete\n"
    unless @$descriptors==2 && $descriptors->[0]->{id} eq 'NET_C2_2C' &&
        $descriptors->[1]->{id} eq 'NET_C2_ION' &&
        $descriptors->[0]->{reactants}->[0]->{species}==3 &&
        $descriptors->[0]->{reactants}->[0]->{coefficient}==1.0;
my $network_policy=CometChannelImportancePolicy->new(species=>\@species);
$network_policy->set_mode('on','network mode');
$network_policy->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:NET_C2_2C)=0.5',
    'network neutral');
$network_policy->add_probability_record(
    'ChemistryChannelBiasedProbability(REACTION:NET_C2_ION)=0.5',
    'network sink');
my $network_resolved=$network_policy->validate_configuration_and_resolve(
    reactions=>$descriptors);
die "reaction-network descriptors did not resolve a complete family\n"
    unless join(',',@{$network_resolved->{configured}}) eq '1,1';
$descriptors->[0]->{id}='MUTATED';
$descriptors->[0]->{reactants}->[0]->{species}=0;
my $fresh=$network->channel_importance_policy_reactions();
die "reaction-network descriptor was not a defensive copy\n"
    unless $fresh->[0]->{id} eq 'NET_C2_2C' &&
        $fresh->[0]->{reactants}->[0]->{species}==3;

print "Comet channel-importance policy test passed\n";
