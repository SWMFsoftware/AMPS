package CometChannelImportancePolicy;

use strict;
use warnings;
use POSIX qw(isfinite);

# Resolve Stage-7 stochastic channel-importance requests into immutable target
# selection probabilities indexed by the generated reaction table.  Records
# use stable, case-sensitive reaction IDs; generated array positions are never
# an input interface.
#
# A policy is configured per complete unit-unimolecular parent family.  If one
# positive-rate channel for a parent is named, every positive-rate sibling --
# including product-free sinks -- must be named, every q* must be positive,
# and the family probabilities must sum to one.  This all-or-none rule prevents
# an omitted rare or loss channel from being silently assigned probability
# zero and acquiring an infinite likelihood ratio.

use constant NORMALIZATION_TOLERANCE => 1.0e-12;

sub new {
  my ($class,%args)=@_;
  die "CometChannelImportancePolicy requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my @species;
  my %seen;
  for my $raw (@{$args{species}}) {
    my $name=uc(_trim($raw));
    die "empty species in channel-importance policy\n" if $name eq '';
    die "duplicate species '$name' in channel-importance policy\n"
        if $seen{$name}++;
    push @species,$name;
  }

  return bless {
    species=>\@species,
    enabled=>0,
    mode_is_set=>0,
    records=>{},
    record_count=>0,
  },$class;
}

sub set_mode {
  my ($self,$value,$location)=@_;
  $location='' unless defined $location;
  $value=uc(_trim($value));
  die "ChemistryChannelImportanceMode must be ON or OFF at $location\n"
      unless $value eq 'ON' || $value eq 'OFF';
  die "duplicate ChemistryChannelImportanceMode at $location\n"
      if $self->{mode_is_set};
  $self->{enabled}=($value eq 'ON') ? 1 : 0;
  $self->{mode_is_set}=1;
}

sub is_enabled {
  my ($self)=@_;
  return $self->{enabled};
}

sub add_probability_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid ChemistryChannelBiasedProbability at $location\n"
      unless $record=~/^\s*ChemistryChannelBiasedProbability\s*\(\s*REACTION\s*:\s*([^)]+?)\s*\)\s*[=:]\s*(.*?)\s*$/i;

  # Preserve the reaction ID exactly.  The REACTION prefix is syntax and is
  # case-insensitive, but changing ID case would make input resolution depend
  # on a convention not present in the reaction database.
  my $id=_trim($1);
  die "empty channel-importance reaction ID at $location\n" if $id eq '';
  my $probability=_positive_probability($2,$location);
  die "duplicate channel-importance probability for reaction '$id' at $location\n"
      if exists $self->{records}->{$id};
  $self->{records}->{$id}={
    probability=>$probability,
    location=>$location,
  };
  $self->{record_count}++;
}

sub validate_configuration_and_resolve {
  my ($self,%args)=@_;
  my $reactions=$args{reactions};
  die "channel-importance resolution requires reaction descriptors\n"
      unless ref($reactions) eq 'ARRAY';

  my %reaction_index;
  my %family;
  for (my $i=0;$i<@$reactions;$i++) {
    my $reaction=$reactions->[$i];
    die "invalid channel-importance reaction descriptor at index $i\n"
        unless ref($reaction) eq 'HASH';
    my $id=$reaction->{id};
    die "invalid channel-importance reaction ID at index $i\n"
        unless defined($id) && $id ne '';
    die "duplicate reaction ID '$id' in channel-importance descriptors\n"
        if exists $reaction_index{$id};
    $reaction_index{$id}=$i;

    my $rate=$reaction->{rate};
    die "reaction '$id' has an invalid nominal rate\n"
        unless defined($rate) && _is_number($rate) &&
            isfinite(0.0+$rate) && $rate>=0.0;
    die "reaction '$id' has no reactant descriptor array\n"
        unless ref($reaction->{reactants}) eq 'ARRAY';
    for my $term (@{$reaction->{reactants}}) {
      die "reaction '$id' has an invalid reactant descriptor\n"
          unless ref($term) eq 'HASH';
      my $species=$term->{species};
      my $coefficient=$term->{coefficient};
      die "reaction '$id' has an invalid reactant species index\n"
          unless defined($species) && $species=~/^\d+$/ &&
              $species<@{$self->{species}};
      die "reaction '$id' has an invalid reactant coefficient\n"
          unless defined($coefficient) && _is_number($coefficient) &&
              isfinite(0.0+$coefficient) && $coefficient>0.0;
    }

    # Only the unit-unimolecular families supported by the particle executor
    # participate.  Positive general mass-action reactions remain valid in a
    # generic database but cannot be selected by this Stage-7 policy.
    if ($rate>0.0 && _is_unit_unimolecular($reaction)) {
      my $parent=$reaction->{reactants}->[0]->{species};
      push @{$family{$parent}},$i;
    }
  }

  if (!$self->is_enabled()) {
    die "channel-importance probability records require ChemistryChannelImportanceMode=ON\n"
        if $self->{record_count}!=0;
    return {
      configured=>[(0)x scalar(@$reactions)],
      biased_probability=>[(0.0)x scalar(@$reactions)],
    };
  }

  die "ChemistryChannelImportanceMode=ON requires at least one biased-probability record\n"
      if $self->{record_count}==0;

  # Validate every named target before family completion so an input error is
  # reported as unknown/disabled/non-unit rather than as an indirect missing
  # sibling error.
  my %touched_parent;
  for my $id (keys %{$self->{records}}) {
    die "unknown channel-importance reaction ID '$id'\n"
        unless exists $reaction_index{$id};
    my $reaction=$reactions->[$reaction_index{$id}];
    die "channel-importance policy targets zero-rate reaction '$id'\n"
        unless $reaction->{rate}>0.0;
    die "channel-importance policy targets non-unit-unimolecular reaction '$id'\n"
        unless _is_unit_unimolecular($reaction);
    $touched_parent{$reaction->{reactants}->[0]->{species}}=1;
  }

  my @configured=(0)x scalar(@$reactions);
  my @probability=(0.0)x scalar(@$reactions);
  for my $parent (sort {$a<=>$b} keys %touched_parent) {
    my $siblings=$family{$parent};
    die "channel-importance parent family $parent has no positive unit-unimolecular reactions\n"
        unless ref($siblings) eq 'ARRAY' && @$siblings;

    my $sum=0.0;
    for my $index (@$siblings) {
      my $id=$reactions->[$index]->{id};
      die "channel-importance family for parent '$self->{species}->[$parent]' is missing reaction '$id'\n"
          unless exists $self->{records}->{$id};
      my $value=$self->{records}->{$id}->{probability};
      $sum+=$value;
      $configured[$index]=1;
      $probability[$index]=$value;
    }
    die "channel-importance probabilities for parent '$self->{species}->[$parent]' sum to $sum instead of one\n"
        unless isfinite($sum) &&
            abs($sum-1.0)<=NORMALIZATION_TOLERANCE;
  }

  return {
    configured=>\@configured,
    biased_probability=>\@probability,
  };
}

sub write_header {
  my ($self,$path,$resolved)=@_;
  die "channel-importance header requires a resolved policy\n"
      unless ref($resolved) eq 'HASH';
  my $configured=$resolved->{configured};
  my $probability=$resolved->{biased_probability};
  die "channel-importance generated arrays differ in length\n"
      unless ref($configured) eq 'ARRAY' &&
          ref($probability) eq 'ARRAY' &&
          @$configured==@$probability;

  for (my $i=0;$i<@$configured;$i++) {
    die "invalid generated channel-importance configured flag at index $i\n"
        unless defined($configured->[$i]) &&
            ($configured->[$i]==0 || $configured->[$i]==1);
    my $value=$probability->[$i];
    die "invalid generated channel-importance probability at index $i\n"
        unless defined($value) && _is_number($value) &&
            isfinite(0.0+$value) &&
            (($configured->[$i] && $value>0.0 && $value<=1.0) ||
             (!$configured->[$i] && $value==0.0));
  }

  open(my $out,'>',$path) or die "cannot write '$path': $!\n";
  print $out "// Generated from the comet application input. Do not edit.\n";
  print $out "#ifndef _COMET_CHANNEL_IMPORTANCE_POLICY_GENERATED_H_\n";
  print $out "#define _COMET_CHANNEL_IMPORTANCE_POLICY_GENERATED_H_\n\n";
  print $out "#include <cstddef>\n\n";
  print $out "namespace Comet { namespace Chemistry { namespace ChannelImportanceGenerated {\n";
  print $out "static const std::size_t ReactionCount=".scalar(@$configured).";\n";
  my @generated_configured=@$configured ? @$configured : (0);
  my @generated_probability=@$probability ? @$probability : (0.0);
  print $out "static const unsigned char Configured[]={".
      join(',',@generated_configured)."};\n";
  print $out "static const double BiasedProbability[]={".
      join(',',map {_format_number($_)} @generated_probability)."};\n";
  print $out "}}} // namespace Comet::Chemistry::ChannelImportanceGenerated\n\n";
  print $out "#endif\n";
  close($out) or die "cannot close '$path': $!\n";
}

sub _is_unit_unimolecular {
  my ($reaction)=@_;
  return 0 unless ref($reaction->{reactants}) eq 'ARRAY' &&
      @{$reaction->{reactants}}==1;
  my $coefficient=$reaction->{reactants}->[0]->{coefficient};
  return defined($coefficient) && _is_number($coefficient) &&
      isfinite(0.0+$coefficient) && $coefficient==1.0;
}

sub _positive_probability {
  my ($value,$location)=@_;
  $value=_trim($value);
  die "ChemistryChannelBiasedProbability must be a positive finite number no greater than one at $location\n"
      unless _is_number($value);
  my $number=0.0+$value;
  die "ChemistryChannelBiasedProbability must be a positive finite number no greater than one at $location\n"
      unless isfinite($number) && $number>0.0 && $number<=1.0;
  return $number;
}

sub _is_number {
  my ($value)=@_;
  return defined($value) &&
      $value=~/^[+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[Ee][+-]?\d+)?$/;
}

sub _format_number {
  my ($value)=@_;
  return sprintf('%.17g',0.0+$value);
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
