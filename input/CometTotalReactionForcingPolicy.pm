package CometTotalReactionForcingPolicy;

use strict;
use warnings;
use POSIX qw(isfinite);

# Resolve the Stage-8 deterministic total-reaction-forcing policy into one
# immutable flag per generated species.  Stage 8 does not alter a physical
# rate and has no biased Bernoulli probability: an enabled parent is split
# into both physical branches, with carrier fractions P and 1-P.  Keeping the
# input policy boolean prevents a numerical tuning parameter from being
# mistaken for chemistry.
#
# Species symbols are stable, case-sensitive identifiers.  The surrounding
# parser may treat record keywords and ON/OFF case-insensitively, but it must
# never silently rewrite a species symbol.

sub new {
  my ($class,%args)=@_;
  die "CometTotalReactionForcingPolicy requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my (%species_index,@species);
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=_trim($args{species}->[$i]);
    die "empty species in total-reaction-forcing policy\n" if $name eq '';
    die "duplicate species '$name' in total-reaction-forcing policy\n"
        if exists $species_index{$name};
    push @species,$name;
    $species_index{$name}=$i;
  }

  return bless {
    species=>\@species,
    species_index=>\%species_index,
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
  die "ChemistryTotalReactionForcingMode must be ON or OFF at $location\n"
      unless $value eq 'ON' || $value eq 'OFF';
  die "duplicate ChemistryTotalReactionForcingMode at $location\n"
      if $self->{mode_is_set};
  $self->{enabled}=($value eq 'ON') ? 1 : 0;
  $self->{mode_is_set}=1;
}

sub is_enabled {
  my ($self)=@_;
  return $self->{enabled};
}

sub add_species_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid ChemistryTotalReactionForcing record at $location\n"
      unless $record=~/^\s*ChemistryTotalReactionForcing\s*\(\s*SPECIES\s*:\s*([^)]+?)\s*\)\s*[=:]\s*(.*?)\s*$/i;

  my $species=_trim($1);
  my $value=uc(_trim($2));
  die "ChemistryTotalReactionForcing species is empty at $location\n"
      if $species eq '';
  die "ChemistryTotalReactionForcing must be ON at $location\n"
      unless $value eq 'ON';
  die "unknown total-reaction-forcing species '$species' at $location\n"
      unless exists $self->{species_index}->{$species};
  die "duplicate total-reaction-forcing policy for species '$species' at $location\n"
      if exists $self->{records}->{$species};
  $self->{records}->{$species}={location=>$location};
  $self->{record_count}++;
}

sub validate_configuration_and_resolve {
  my ($self,%args)=@_;
  my $reactions=$args{reactions};
  die "total-reaction-forcing resolution requires reaction descriptors\n"
      unless ref($reactions) eq 'ARRAY';

  # Validate the descriptor API rather than trusting a particular database
  # implementation.  A future generic network can contain mass-action rows,
  # but the particle forcing executor may target only a positive-rate,
  # unit-unimolecular parent family.
  my %reaction_id;
  my @eligible=(0)x scalar(@{$self->{species}});
  for (my $i=0;$i<@$reactions;$i++) {
    my $reaction=$reactions->[$i];
    die "invalid total-reaction-forcing descriptor at index $i\n"
        unless ref($reaction) eq 'HASH';
    my $id=$reaction->{id};
    die "invalid total-reaction-forcing reaction ID at index $i\n"
        unless defined($id) && $id ne '';
    die "duplicate reaction ID '$id' in total-reaction-forcing descriptors\n"
        if exists $reaction_id{$id};
    $reaction_id{$id}=1;

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
    if ($rate>0.0 && _is_unit_unimolecular($reaction)) {
      my $parent=$reaction->{reactants}->[0]->{species};
      $eligible[$parent]=1;
    }
  }

  if (!$self->is_enabled()) {
    die "total-reaction-forcing records require ChemistryTotalReactionForcingMode=ON\n"
        if $self->{record_count}!=0;
    return {configured=>[(0)x scalar(@{$self->{species}})]};
  }
  die "ChemistryTotalReactionForcingMode=ON requires at least one species record\n"
      if $self->{record_count}==0;

  my @configured=(0)x scalar(@{$self->{species}});
  for my $species (keys %{$self->{records}}) {
    my $index=$self->{species_index}->{$species};
    die "total-reaction-forcing species '$species' has no positive-rate unit-unimolecular reaction\n"
        unless $eligible[$index];
    $configured[$index]=1;
  }
  return {configured=>\@configured};
}

sub write_header {
  my ($self,$path,$resolved)=@_;
  die "total-reaction-forcing header requires a resolved policy\n"
      unless ref($resolved) eq 'HASH' &&
          ref($resolved->{configured}) eq 'ARRAY';
  my $configured=$resolved->{configured};
  die "total-reaction-forcing generated species count is inconsistent\n"
      unless @$configured==@{$self->{species}};
  for (my $i=0;$i<@$configured;$i++) {
    die "invalid generated total-reaction-forcing flag at species $i\n"
        unless defined($configured->[$i]) &&
            ($configured->[$i]==0 || $configured->[$i]==1);
  }

  open(my $out,'>',$path) or die "cannot write '$path': $!\n";
  print $out "// Generated from the comet application input. Do not edit.\n";
  print $out "#ifndef _COMET_TOTAL_REACTION_FORCING_POLICY_GENERATED_H_\n";
  print $out "#define _COMET_TOTAL_REACTION_FORCING_POLICY_GENERATED_H_\n\n";
  print $out "#include <cstddef>\n\n";
  print $out "namespace Comet { namespace Chemistry { namespace TotalReactionForcingGenerated {\n";
  print $out "static const std::size_t SpeciesCount=".scalar(@$configured).";\n";
  my @sentinel=@$configured ? @$configured : (0);
  print $out "static const unsigned char Configured[]={".
      join(',',@sentinel)."};\n";
  print $out "}}} // namespace Comet::Chemistry::TotalReactionForcingGenerated\n\n";
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

sub _is_number {
  my ($value)=@_;
  return defined($value) &&
      $value=~/^[+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[Ee][+-]?\d+)?$/;
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
