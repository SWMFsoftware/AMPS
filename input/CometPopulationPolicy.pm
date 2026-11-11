package CometPopulationPolicy;

use strict;
use warnings;

# Input-side resolver for the radial particle-population policy.  Parsing and
# resolution live in one small module so the application parser and focused
# tests use exactly the same selector grammar and precedence rules.
#
# A selector owns one policy at its precedence level: either a target range or
# a lock.  Consequently a selector cannot be both ranged and locked.  Policies
# are resolved only after the full source and reaction input has been read:
#
#   explicit species > derived DIRECT_SOURCE/CHEMICAL_ONLY group > GLOBAL.
#
# DIRECT_SOURCE and CHEMICAL_ONLY are deliberately derived from input data,
# not from a comet-specific species list.  Species outside both groups inherit
# the global range unless an explicit species policy is supplied.

use constant {
  POLICY_GLOBAL        => 0,
  POLICY_DIRECT_SOURCE => 1,
  POLICY_CHEMICAL_ONLY => 2,
  POLICY_SPECIES       => 3,
  CLASS_OTHER          => 0,
  CLASS_DIRECT_SOURCE  => 1,
  CLASS_CHEMICAL_ONLY  => 2,
};

sub new {
  my ($class,%args)=@_;
  die "CometPopulationPolicy requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my @species;
  my %species_index;
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=uc(_trim($args{species}->[$i]));
    die "empty species name in population-policy species list\n"
        if $name eq '';
    die "duplicate species '$name' in population-policy species list\n"
        if exists $species_index{$name};
    push @species,$name;
    $species_index{$name}=$i;
  }

  return bless {
    species=>\@species,
    species_index=>\%species_index,
    policies=>{},
  },$class;
}

sub add_range_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid RadialParticleNumberPerCellRange at $location\n"
      unless $record=~/^\s*RadialParticleNumberPerCellRange\s*(?:\(\s*([^)]+?)\s*\))?\s*[=:]\s*(.*?)\s*$/i;
  my ($selector,$value)=($1,$2);
  my @field=map {_trim($_)} split(/\|/,$value,-1);
  die "RadialParticleNumberPerCellRange expects MINIMUM | MAXIMUM at $location\n"
      unless @field==2;
  die "RadialParticleNumberPerCellRange values must be positive integers at $location\n"
      unless $field[0]=~/^\d+$/ && $field[0]>0 &&
             $field[1]=~/^\d+$/ && $field[1]>0;
  die "RadialParticleNumberPerCellRange maximum must be at least twice the minimum for power-of-two adaptation at $location\n"
      unless $field[1]>=2*$field[0];

  $selector=defined($selector) ? $self->_normalize_selector($selector) :
      'GLOBAL';
  $self->_add_policy($selector,{
      kind=>'range',minimum=>$field[0]+0,maximum=>$field[1]+0,
      location=>$location});
}

sub add_lock_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid RadialParticleWeightLock at $location\n"
      unless $record=~/^\s*RadialParticleWeightLock\s*[=:]\s*(.*?)\s*$/i;
  my $selector=$self->_normalize_selector($1);
  die "RadialParticleWeightLock requires a species, DIRECT_SOURCE, or CHEMICAL_ONLY selector at $location\n"
      if $selector eq 'GLOBAL';
  $self->_add_policy($selector,{kind=>'lock',location=>$location});
}

sub has_global_range {
  my ($self)=@_;
  return exists $self->{policies}->{GLOBAL} &&
      $self->{policies}->{GLOBAL}->{kind} eq 'range';
}

sub global_range {
  my ($self)=@_;
  die "global particle-number range is not configured\n"
      unless $self->has_global_range();
  my $entry=$self->{policies}->{GLOBAL};
  return ($entry->{minimum},$entry->{maximum});
}

sub policy_count {
  my ($self)=@_;
  return scalar keys %{$self->{policies}};
}

sub resolve {
  my ($self,%args)=@_;
  my $direct=$args{direct_source};
  my $reachable=$args{chemical_reachable};
  die "population-policy resolution requires direct-source flags\n"
      unless ref($direct) eq 'ARRAY' && @$direct==@{$self->{species}};
  die "population-policy resolution requires chemical-reachability flags\n"
      unless ref($reachable) eq 'ARRAY' && @$reachable==@{$self->{species}};
  die "population-policy resolution requires a global range\n"
      unless $self->has_global_range();

  my (@minimum,@maximum,@locked,@source,@class,@selector);
  for (my $species=0;$species<@{$self->{species}};$species++) {
    my $population_class=CLASS_OTHER;
    my $group;
    if ($direct->[$species]) {
      $population_class=CLASS_DIRECT_SOURCE;
      $group='DIRECT_SOURCE';
    }
    elsif ($reachable->[$species]) {
      $population_class=CLASS_CHEMICAL_ONLY;
      $group='CHEMICAL_ONLY';
    }

    my $selected='GLOBAL';
    $selected=$group if defined($group) &&
        exists $self->{policies}->{$group};
    my $species_selector=$self->{species}->[$species];
    $selected=$species_selector
        if exists $self->{policies}->{$species_selector};

    my $entry=$self->{policies}->{$selected};
    my $global=$self->{policies}->{GLOBAL};
    # Locked policies still carry the global range in generated tables.  The
    # range is inactive while locked, but retaining valid bounds keeps common
    # emergency diagnostics and startup reporting well defined.
    push @minimum,$entry->{kind} eq 'range' ?
        $entry->{minimum} : $global->{minimum};
    push @maximum,$entry->{kind} eq 'range' ?
        $entry->{maximum} : $global->{maximum};
    push @locked,$entry->{kind} eq 'lock' ? 1 : 0;
    push @source,$selected eq 'GLOBAL' ? POLICY_GLOBAL :
        $selected eq 'DIRECT_SOURCE' ? POLICY_DIRECT_SOURCE :
        $selected eq 'CHEMICAL_ONLY' ? POLICY_CHEMICAL_ONLY : POLICY_SPECIES;
    push @class,$population_class;
    push @selector,$selected;
  }

  return {
    minimum=>\@minimum,maximum=>\@maximum,locked=>\@locked,
    source=>\@source,class=>\@class,selector=>\@selector,
  };
}

sub _normalize_selector {
  my ($self,$raw)=@_;
  my $selector=uc(_trim($raw));
  die "empty population-policy selector\n" if $selector eq '';
  return $selector if $selector eq 'DIRECT_SOURCE' ||
      $selector eq 'CHEMICAL_ONLY' || $selector eq 'GLOBAL';
  die "unknown population-policy selector '$selector'\n"
      unless exists $self->{species_index}->{$selector};
  return $selector;
}

sub _add_policy {
  my ($self,$selector,$entry)=@_;
  if (exists $self->{policies}->{$selector}) {
    my $old=$self->{policies}->{$selector};
    my $old_location=$old->{location} ne '' ? " at $old->{location}" : '';
    die "duplicate population policy for selector '$selector'".
        "$old_location; a selector cannot have both a range and a lock\n";
  }
  $self->{policies}->{$selector}=$entry;
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
