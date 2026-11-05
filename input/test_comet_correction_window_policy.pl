#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use CometCorrectionWindowPolicy;
use CometReactionNetwork;

my @species=qw(H2O OH H C2H2 C2 CO2);
my $network=CometReactionNetwork->new(species=>\@species);
$network->add_reaction(
    'H2O_OH_H | H2O | OH + H | CONSTANT | 1 | 0 | 0 | 1 | 0',
    'correction-window test');
$network->add_reaction(
    'C2H2_C2_H2 | C2H2 | C2 + H | CONSTANT | 1 | 0 | 0 | 1 | 0',
    'correction-window test');
my $reachable=$network->reachable_species_from([0,3]);

# Omission must preserve the pre-Stage-5 behavior.  The mode and hard split
# limit are parsed by this pure module so these semantics do not depend on a
# configured or compiled AMPS tree.
my $default=CometCorrectionWindowPolicy->new(species=>\@species);
die "correction windows are not disabled by default\n"
    if $default->is_enabled() || $default->has_maximum_splits();
die "empty default correction-window configuration was rejected\n"
    if defined $default->validate_configuration_and_resolve(
        direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);

my $policy=CometCorrectionWindowPolicy->new(species=>\@species);
$policy->set_mode('on','mode test');
$policy->set_maximum_splits('32','split-limit test');
$policy->add_window_record(
    'ParticleWeightCorrectionWindow=0.125 | 1 | 8','global test');
$policy->add_window_record(
    'ParticleWeightCorrectionWindow(DIRECT_SOURCE)=0.5 | 1 | 2',
    'direct group test');
$policy->add_window_record(
    'ParticleWeightCorrectionWindow(CHEMICAL_ONLY)=0.0625 | 0.25 | 4',
    'chemical group test');
$policy->add_window_record(
    'ParticleWeightCorrectionWindow(C2)=0.03125 | 0.125 | 2',
    'explicit species test');

die "global correction-window access is incorrect\n"
    unless $policy->has_global_window() && $policy->window_count()==4;
my ($global_lower,$global_target,$global_upper)=$policy->global_window();
die "global correction-window values are incorrect\n"
    unless $global_lower==0.125 && $global_target==1.0 &&
        $global_upper==8.0;

my $resolved=$policy->resolve(
    direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
die "correction-window mode or split limit was not retained\n"
    unless $policy->is_enabled() && $policy->maximum_splits()==32;
my $validated=$policy->validate_configuration_and_resolve(
    direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
die "validated correction-window resolution differs from direct resolution\n"
    unless join('|',@{$validated->{selector}}) eq
        join('|',@{$resolved->{selector}});

# Both direct parents use the derived group.  Chemical daughters use the
# chemical-only group except for C2, whose explicit record has precedence.
# CO2 is neither sourced nor reachable and therefore uses GLOBAL.
for my $species_index (0,3) {
  die "direct-source correction window was not resolved\n"
      unless $resolved->{lower}->[$species_index]==0.5 &&
          $resolved->{target}->[$species_index]==1.0 &&
          $resolved->{upper}->[$species_index]==2.0 &&
          $resolved->{source}->[$species_index]==
              CometCorrectionWindowPolicy::POLICY_DIRECT_SOURCE;
}
die "chemical-only correction window was not resolved\n"
    unless $resolved->{lower}->[1]==0.0625 &&
        $resolved->{target}->[1]==0.25 &&
        $resolved->{upper}->[1]==4.0 &&
        $resolved->{source}->[1]==
            CometCorrectionWindowPolicy::POLICY_CHEMICAL_ONLY;
die "explicit species did not override the chemical-only group\n"
    unless $resolved->{lower}->[4]==0.03125 &&
        $resolved->{target}->[4]==0.125 &&
        $resolved->{upper}->[4]==2.0 &&
        $resolved->{source}->[4]==
            CometCorrectionWindowPolicy::POLICY_SPECIES;
die "global correction window was not inherited\n"
    unless $resolved->{lower}->[5]==0.125 &&
        $resolved->{target}->[5]==1.0 &&
        $resolved->{upper}->[5]==8.0 &&
        $resolved->{source}->[5]==
            CometCorrectionWindowPolicy::POLICY_GLOBAL;

# Resolution must not depend on input record order.  This is important because
# group membership is derived only after the complete source/network input.
my $reordered=CometCorrectionWindowPolicy->new(species=>\@species);
$reordered->add_window_record(
    'ParticleWeightCorrectionWindow(C2)=0.03125 | 0.125 | 2','species');
$reordered->add_window_record(
    'ParticleWeightCorrectionWindow(CHEMICAL_ONLY)=0.0625 | 0.25 | 4',
    'chemical');
$reordered->add_window_record(
    'ParticleWeightCorrectionWindow=0.125 | 1 | 8','global');
$reordered->add_window_record(
    'ParticleWeightCorrectionWindow(DIRECT_SOURCE)=0.5 | 1 | 2','direct');
my $resolved_reordered=$reordered->resolve(
    direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
for my $field (qw(lower target upper source selector)) {
  die "correction-window resolution depends on record order ($field)\n"
      unless join('|',@{$resolved->{$field}}) eq
          join('|',@{$resolved_reordered->{$field}});
}

sub must_fail {
  my ($description,$callback)=@_;
  my $failed=0;
  eval {$callback->();};
  $failed=1 if $@;
  die "$description was accepted\n" unless $failed;
}

must_fail('duplicate selector',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->add_window_record(
      'ParticleWeightCorrectionWindow=0.25 | 1 | 4','one');
  $p->add_window_record(
      'ParticleWeightCorrectionWindow(GLOBAL)=0.5 | 1 | 2','two');
});
must_fail('invalid mode',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_mode('MAYBE','invalid mode');
});
must_fail('duplicate mode',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_mode('ON','one');
  $p->set_mode('OFF','two');
});
for my $invalid_limit (
    '0','1','2.5','-2','Inf','999999999999999999999999999999999999') {
  must_fail("invalid split limit '$invalid_limit'",sub {
    my $p=CometCorrectionWindowPolicy->new(species=>\@species);
    $p->set_maximum_splits($invalid_limit,'invalid limit');
  });
}
must_fail('duplicate split limit',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_maximum_splits(2,'one');
  $p->set_maximum_splits(4,'two');
});
must_fail('unknown selector',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->add_window_record(
      'ParticleWeightCorrectionWindow(UNKNOWN)=0.25 | 1 | 4','unknown');
});
for my $invalid (
    '0 | 1 | 2','-1 | 1 | 2','NaN | 1 | 2','Inf | 1 | 2',
    '0.25 | 0.25 | 2','0.25 | 2 | 2','0.25 | 1',
    '0.25 | 1 | 2 | 4','0.9 | 1 | 1.1') {
  must_fail("invalid correction window '$invalid'",sub {
    my $p=CometCorrectionWindowPolicy->new(species=>\@species);
    $p->add_window_record(
        "ParticleWeightCorrectionWindow=$invalid",'invalid');
  });
}
must_fail('missing global fallback',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->add_window_record(
      'ParticleWeightCorrectionWindow(C2)=0.25 | 1 | 4','species');
  $p->resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});
must_fail('enabled mode missing a window',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_mode('ON','mode');
  $p->set_maximum_splits(4,'limit');
  $p->validate_configuration_and_resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});
must_fail('enabled mode missing split limit',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_mode('ON','mode');
  $p->add_window_record(
      'ParticleWeightCorrectionWindow=0.25 | 1 | 4','global');
  $p->validate_configuration_and_resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});
must_fail('window record with disabled mode',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_mode('OFF','mode');
  $p->add_window_record(
      'ParticleWeightCorrectionWindow=0.25 | 1 | 4','global');
  $p->validate_configuration_and_resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});
must_fail('split limit with omitted mode',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->set_maximum_splits(4,'limit');
  $p->validate_configuration_and_resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});
must_fail('invalid direct-source flag vector',sub {
  $policy->resolve(
      direct_source=>[1],chemical_reachable=>$reachable);
});
must_fail('invalid reachability flag vector',sub {
  $policy->resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>[1]);
});
must_fail('direct-source window excluding correction one',sub {
  my $p=CometCorrectionWindowPolicy->new(species=>\@species);
  $p->add_window_record(
      'ParticleWeightCorrectionWindow=0.01 | 0.1 | 0.5','global');
  $p->resolve(
      direct_source=>[1,0,0,1,0,0],chemical_reachable=>$reachable);
});

print "Comet correction-window policy test passed\n";
