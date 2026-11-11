package CometCorrectionWindowPolicy;

use strict;
use warnings;
use POSIX qw(isfinite LONG_MAX);

# Input-side resolver for the Stage-5 individual-correction windows.  This is
# intentionally separate from CometPopulationPolicy: population ranges govern
# a block's base W/dt carrier, whereas these windows govern the dimensionless
# correction carried by one particle.  Combining the records would make it
# too easy to apply a correction window while individual corrections are off.
#
# Selectors use the same deterministic precedence as the population policy:
#
#   explicit species > derived DIRECT_SOURCE/CHEMICAL_ONLY group > GLOBAL.
#
# Group membership is supplied only after all source and reaction records have
# been parsed.  The policy therefore contains no comet- or network-specific
# species classification and remains independent of input record order.

use constant {
  POLICY_GLOBAL        => 0,
  POLICY_DIRECT_SOURCE => 1,
  POLICY_CHEMICAL_ONLY => 2,
  POLICY_SPECIES       => 3,
};

sub new {
  my ($class,%args)=@_;
  die "CometCorrectionWindowPolicy requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my @species;
  my %species_index;
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=uc(_trim($args{species}->[$i]));
    die "empty species name in correction-window species list\n"
        if $name eq '';
    die "duplicate species '$name' in correction-window species list\n"
        if exists $species_index{$name};
    push @species,$name;
    $species_index{$name}=$i;
  }

  return bless {
    species=>\@species,
    species_index=>\%species_index,
    windows=>{},
    enabled=>0,
    mode_is_set=>0,
    maximum_splits=>undef,
  },$class;
}

sub set_mode {
  my ($self,$value,$location)=@_;
  $location='' unless defined $location;
  $value=uc(_trim($value));
  die "ParticleWeightCorrectionWindowMode must be ON or OFF at $location\n"
      unless $value eq 'ON' || $value eq 'OFF';
  die "duplicate ParticleWeightCorrectionWindowMode at $location\n"
      if $self->{mode_is_set};
  $self->{enabled}=($value eq 'ON') ? 1 : 0;
  $self->{mode_is_set}=1;
}

sub is_enabled {
  my ($self)=@_;
  # Omission is intentionally OFF so old inputs retain their exact behavior.
  return $self->{enabled};
}

sub set_maximum_splits {
  my ($self,$value,$location)=@_;
  $location='' unless defined $location;
  $value=_trim($value);
  die "ParticleWeightCorrectionWindowMaxSplits must be an integer of at least two at $location\n"
      unless $value=~/^\d+$/ && $value>=2 && $value<=LONG_MAX;
  die "duplicate ParticleWeightCorrectionWindowMaxSplits at $location\n"
      if defined $self->{maximum_splits};
  $self->{maximum_splits}=0+$value;
}

sub has_maximum_splits {
  my ($self)=@_;
  return defined $self->{maximum_splits};
}

sub maximum_splits {
  my ($self)=@_;
  die "particle-correction window split limit is not configured\n"
      unless $self->has_maximum_splits();
  return $self->{maximum_splits};
}

sub add_window_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid ParticleWeightCorrectionWindow at $location\n"
      unless $record=~/^\s*ParticleWeightCorrectionWindow\s*(?:\(\s*([^)]+?)\s*\))?\s*[=:]\s*(.*?)\s*$/i;
  my ($selector,$value)=($1,$2);
  my @field=map {_trim($_)} split(/\|/,$value,-1);
  die "ParticleWeightCorrectionWindow expects LOWER | TARGET | UPPER at $location\n"
      unless @field==3;

  my ($lower,$target,$upper)=map {
    _positive_finite_number($_,
        "ParticleWeightCorrectionWindow values",$location)
  } @field;
  die "ParticleWeightCorrectionWindow requires LOWER < TARGET < UPPER at $location\n"
      unless $lower<$target && $target<$upper;
  # Runtime splitting uses ceil(c/target) equal children.  Immediately above
  # the upper boundary this is the smallest split correction the controller
  # can create.  Reject a contradictory window here instead of discovering it
  # after a production run has begun.
  my $first_split_count=int($upper/$target)+1;
  die "ParticleWeightCorrectionWindow is incompatible with ceil(c/target) splitting at $location\n"
      if $lower>$upper/$first_split_count;

  $selector=defined($selector) ? $self->_normalize_selector($selector) :
      'GLOBAL';
  $self->_add_window($selector,{
      lower=>$lower,target=>$target,upper=>$upper,location=>$location});
}

sub has_global_window {
  my ($self)=@_;
  return exists $self->{windows}->{GLOBAL};
}

sub global_window {
  my ($self)=@_;
  die "global particle-correction window is not configured\n"
      unless $self->has_global_window();
  my $entry=$self->{windows}->{GLOBAL};
  return ($entry->{lower},$entry->{target},$entry->{upper});
}

sub window_count {
  my ($self)=@_;
  return scalar keys %{$self->{windows}};
}

sub resolve {
  my ($self,%args)=@_;
  my $direct=$args{direct_source};
  my $reachable=$args{chemical_reachable};
  die "correction-window resolution requires direct-source flags\n"
      unless ref($direct) eq 'ARRAY' && @$direct==@{$self->{species}};
  die "correction-window resolution requires chemical-reachability flags\n"
      unless ref($reachable) eq 'ARRAY' && @$reachable==@{$self->{species}};
  die "correction-window resolution requires a global window\n"
      unless $self->has_global_window();

  my (@lower,@target,@upper,@source,@selector);
  for (my $species=0;$species<@{$self->{species}};$species++) {
    my $group;
    if ($direct->[$species]) {
      $group='DIRECT_SOURCE';
    }
    elsif ($reachable->[$species]) {
      $group='CHEMICAL_ONLY';
    }

    my $selected='GLOBAL';
    $selected=$group if defined($group) &&
        exists $self->{windows}->{$group};
    my $species_selector=$self->{species}->[$species];
    $selected=$species_selector
        if exists $self->{windows}->{$species_selector};

    my $entry=$self->{windows}->{$selected};
    # Directly injected particles enter the model with correction one.  A
    # Stage-5 window must not immediately roulette or split that population;
    # changing direct-source statistics requires a separate, explicit model.
    die "correction window '$selected' for direct-source species " .
        "$self->{species}->[$species] must contain correction 1\n"
        if $direct->[$species] &&
           !($entry->{lower}<=1.0 && 1.0<=$entry->{upper});
    push @lower,$entry->{lower};
    push @target,$entry->{target};
    push @upper,$entry->{upper};
    push @source,$selected eq 'GLOBAL' ? POLICY_GLOBAL :
        $selected eq 'DIRECT_SOURCE' ? POLICY_DIRECT_SOURCE :
        $selected eq 'CHEMICAL_ONLY' ? POLICY_CHEMICAL_ONLY : POLICY_SPECIES;
    push @selector,$selected;
  }

  return {
    lower=>\@lower,target=>\@target,upper=>\@upper,
    source=>\@source,selector=>\@selector,
  };
}

sub validate_configuration_and_resolve {
  my ($self,%args)=@_;
  if (!$self->is_enabled()) {
    die "ParticleWeightCorrectionWindow records and split limit require ParticleWeightCorrectionWindowMode=ON\n"
        if $self->window_count()!=0 || $self->has_maximum_splits();
    return undef;
  }

  die "ParticleWeightCorrectionWindowMode=ON requires a global ParticleWeightCorrectionWindow\n"
      unless $self->has_global_window();
  die "ParticleWeightCorrectionWindowMode=ON requires ParticleWeightCorrectionWindowMaxSplits\n"
      unless $self->has_maximum_splits();
  return $self->resolve(%args);
}

sub _normalize_selector {
  my ($self,$raw)=@_;
  my $selector=uc(_trim($raw));
  die "empty correction-window selector\n" if $selector eq '';
  return $selector if $selector eq 'DIRECT_SOURCE' ||
      $selector eq 'CHEMICAL_ONLY' || $selector eq 'GLOBAL';
  die "unknown correction-window selector '$selector'\n"
      unless exists $self->{species_index}->{$selector};
  return $selector;
}

sub _add_window {
  my ($self,$selector,$entry)=@_;
  if (exists $self->{windows}->{$selector}) {
    my $old=$self->{windows}->{$selector};
    my $old_location=$old->{location} ne '' ? " at $old->{location}" : '';
    die "duplicate correction window for selector '$selector'".
        "$old_location\n";
  }
  $self->{windows}->{$selector}=$entry;
}

sub _positive_finite_number {
  my ($value,$description,$location)=@_;
  # Perl accepts spellings such as Inf and NaN numerically.  Restrict the
  # grammar before conversion so generated C++ receives a finite literal.
  die "$description must be positive finite numbers at $location\n"
      unless defined($value) &&
          $value=~/^[+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[Ee][+-]?\d+)?$/;
  my $number=0.0+$value;
  die "$description must be positive finite numbers at $location\n"
      unless $number>0.0 && isfinite($number);
  return $number;
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
