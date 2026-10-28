package CometReactionNetwork;

use strict;
use warnings;

# Input-side builder for the geometry-independent C++ reaction tables.
# Species names are resolved here because AMPS fixes the species order while
# preprocessing an input file.  The emitted header then contains only numeric
# indices and immutable arrays, so runtime chemistry performs no text parsing.
#
# Reaction records have nine pipe-separated fields:
#   ID | REACTANTS | PRODUCTS | LAW | A | BETA | TA_K | RREF_AU | EXCESS_EV
# A side is a '+' separated list of NAME or COEFFICIENT*NAME terms.  The
# literal NONE is accepted only for products and represents loss from the
# tracked neutral system, for example photoionization.  Excess energy is
# converted from eV to joules at generation time.

sub new {
  my ($class,%args)=@_;
  die "CometReactionNetwork requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my %species_index;
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=uc(_trim($args{species}->[$i]));
    die "empty species name in SpeciesList\n" if $name eq '';
    die "duplicate species '$name' in SpeciesList\n"
        if exists $species_index{$name};
    $species_index{$name}=$i;
  }

  return bless {
    species_index=>\%species_index,
    reactions=>[],
    ids=>{},
  },$class;
}

sub add_reaction {
  my ($self,$text,$location)=@_;
  $location='' unless defined $location;
  $text=_trim($text);

  my @field=map {_trim($_)} split(/\|/,$text,-1);
  die "invalid ChemistryReaction at $location: expected " .
      "ID | REACTANTS | PRODUCTS | LAW | A | BETA | TA_K | RREF_AU | EXCESS_EV\n"
      unless @field==9;

  my ($id,$reactants,$products,$law,$a,$beta,$ta,$rref,$energy_ev)=@field;
  die "empty reaction ID at $location\n" if $id eq '';
  die "duplicate reaction ID '$id' at $location\n" if $self->{ids}->{$id};

  $law=uc($law);
  die "unknown rate law '$law' for reaction '$id' at $location\n"
      unless $law eq 'CONSTANT' || $law eq 'HELIOCENTRIC' ||
          $law eq 'ARRHENIUS';

  for my $entry ([$a,'A'],[$beta,'BETA'],[$ta,'TA_K'],
      [$rref,'RREF_AU'],[$energy_ev,'EXCESS_EV']) {
    die "invalid numeric $entry->[1] '$entry->[0]' for reaction '$id' at $location\n"
        unless _is_number($entry->[0]);
  }

  $a+=0.0,$beta+=0.0,$ta+=0.0,$rref+=0.0,$energy_ev+=0.0;
  die "negative A for reaction '$id' at $location\n" if $a<0.0;
  die "negative EXCESS_EV for reaction '$id' at $location\n"
      if $energy_ev<0.0;
  die "HELIOCENTRIC reaction '$id' requires RREF_AU>0 at $location\n"
      if $law eq 'HELIOCENTRIC' && $rref<=0.0;

  # Keep the parsed representation independent of Perl source ordering.  Flat
  # term offsets are assigned only in write_header after every input record and
  # optional named database have been loaded.
  my $reaction={
    id=>$id,
    reactants=>$self->_parse_side($reactants,$id,'reactant',$location),
    products=>$self->_parse_side($products,$id,'product',$location,1),
    law=>$law,
    a=>$a,
    beta=>$beta,
    ta=>$ta,
    rref=>$rref,
    energy_j=>$energy_ev*1.602176634e-19,
  };

  push @{$self->{reactions}},$reaction;
  $self->{ids}->{$id}=1;
}

sub reaction_count {
  my ($self)=@_;
  return scalar @{$self->{reactions}};
}

sub write_header {
  my ($self,$path)=@_;
  open(my $out,'>',$path) or die "cannot write '$path': $!\n";

  # Flatten both sides into one table.  Offsets and counts make the generated
  # C++ arrays compact while preserving arbitrary reaction molecularity.
  my @terms;
  my @reactions;
  for my $reaction (@{$self->{reactions}}) {
    my $reactant_offset=scalar @terms;
    push @terms,@{$reaction->{reactants}};
    my $product_offset=scalar @terms;
    push @terms,@{$reaction->{products}};
    push @reactions,{
      %$reaction,
      reactant_offset=>$reactant_offset,
      product_offset=>$product_offset,
    };
  }

  print $out "// Generated from the comet application input. Do not edit.\n";
  print $out "#ifndef _COMET_CHEMISTRY_GENERATED_H_\n";
  print $out "#define _COMET_CHEMISTRY_GENERATED_H_\n\n";
  print $out "namespace Comet { namespace Chemistry { namespace Generated {\n";
  print $out "static const std::size_t TermCount=".(scalar @terms).";\n";
  print $out "static const std::size_t ReactionCount=".(scalar @reactions).";\n";

  print $out "static const ChemicalSolver::cStoichiometricTerm Terms[] = {\n";
  if (@terms) {
    for my $term (@terms) {
      print $out "  {$term->{species},"._number($term->{coefficient})."},\n";
    }
  }
  else {
    print $out "  {0,1.0}\n";
  }
  print $out "};\n\n";

  print $out "static const ChemicalSolver::cReaction Reactions[] = {\n";
  if (@reactions) {
    for my $reaction (@reactions) {
      my $law="ChemicalSolver::RateLaw::$reaction->{law}";
      $law='ChemicalSolver::RateLaw::HeliocentricInverseSquare'
          if $reaction->{law} eq 'HELIOCENTRIC';
      $law='ChemicalSolver::RateLaw::Constant'
          if $reaction->{law} eq 'CONSTANT';
      $law='ChemicalSolver::RateLaw::Arrhenius'
          if $reaction->{law} eq 'ARRHENIUS';
      my $id=$reaction->{id};
      $id=~s/\\/\\\\/g;
      $id=~s/"/\\"/g;
      print $out "  {\"$id\",$reaction->{reactant_offset},".
          scalar(@{$reaction->{reactants}}).",$reaction->{product_offset},".
          scalar(@{$reaction->{products}}).",$law,".
          join(',',map {_number($_)} ($reaction->{a},$reaction->{beta},
              $reaction->{ta},$reaction->{rref},$reaction->{energy_j}))."},\n";
    }
  }
  else {
    print $out "  {\"none\",0,0,0,0,ChemicalSolver::RateLaw::Constant,".
        "0.0,0.0,0.0,1.0,0.0}\n";
  }
  print $out "};\n";
  print $out "}}} // namespace Comet::Chemistry::Generated\n\n";
  print $out "#endif\n";
  close($out) or die "cannot close '$path': $!\n";
}

sub _parse_side {
  my ($self,$text,$id,$kind,$location,$allow_empty)=@_;
  $text=_trim($text);
  # An empty product list is meaningful only when the input says NONE.  A
  # blank field is rejected because it is more likely a malformed record.
  return [] if $allow_empty && uc($text) eq 'NONE';
  die "reaction '$id' has no ${kind}s at $location\n" if $text eq '';

  my @result;
  for my $raw (split(/\+/,$text)) {
    my $term=_trim($raw);
    my ($coefficient,$name);
    if ($term=~/^([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)\s*\*\s*(\S+)$/) {
      ($coefficient,$name)=($1+0.0,$2);
    }
    elsif ($term=~/^(\S+)$/) {
      ($coefficient,$name)=(1.0,$1);
    }
    else {
      die "invalid $kind term '$term' in reaction '$id' at $location\n";
    }

    $name=uc($name);
    die "non-positive coefficient in $kind term '$term' of reaction '$id' at $location\n"
        unless $coefficient>0.0;
    die "unknown species '$name' in reaction '$id' at $location\n"
        unless exists $self->{species_index}->{$name};
    push @result,{
      species=>$self->{species_index}->{$name},
      coefficient=>$coefficient,
    };
  }

  return \@result;
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

sub _is_number {
  my ($value)=@_;
  return $value=~/^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?$/;
}

sub _number {
  my ($value)=@_;
  return sprintf('%.17g',$value+0.0);
}

1;
