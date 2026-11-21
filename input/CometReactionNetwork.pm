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
#
# Chemistry variants (srcHartley2/CHEMISTRY_VARIANTS.md).  Optional KEY=VALUE
# fields may follow the nine fields:
#   GROUP=<GROUP>:<VARIANT>        the reaction belongs to one variant of a
#                                  reaction group; it is active only when that
#                                  variant is selected at run time
#   RATE=TABLE:<TABLE>:<ROLE>      the runtime rate at 1 AU is computed from a
#                                  rate table (ChemistryRateTable) instead of a
#                                  chemistry.rate.<ID>_1au_s-1 input key
# Group and table records:
#   ChemistryGroup=<GROUP> | DEFAULT=<VARIANT> | VARIANTS=<V1>,<V2>,...
#   ChemistryRateTable=<TABLE> | MODEL=SA1982_OH_LIFETIME | SCALE_S=<s> | ION_FRACTION=<f>
# The checks in validate_variants() run before any header is written.  Records
# without the optional fields keep their exact previous meaning and order.

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
    groups=>[],        # in input order: {id,default,variants=>[...],location}
    group_index=>{},
    tables=>[],        # in input order: {id,model,scale_s,ion_fraction,location}
    table_index=>{},
  },$class;
}

# Rate-table models known to the runtime (CometRuntimeInput.cpp) and the roles
# a reaction may take in each.  SA1982_OH_LIFETIME: OH loss with the
# heliocentric-velocity-dependent lifetime of the legacy BLOCKDATA OHPHOT;
# DISSOCIATION is k_d = 1/tau_FUV + 1/(SCALE_S*OHTAUS(v_h)), IONIZATION is
# ION_FRACTION*k_d (srcHartley2/CHEMISTRY_VARIANTS.md, section 3).
my %TABLE_MODEL_ROLES=(
  SA1982_OH_LIFETIME=>{DISSOCIATION=>1,IONIZATION=>2},
);

sub add_group {
  my ($self,$text,$location)=@_;
  $location='' unless defined $location;
  my @field=map {_trim($_)} split(/\|/,_trim($text),-1);
  die "invalid ChemistryGroup at $location: expected " .
      "GROUP | DEFAULT=<VARIANT> | VARIANTS=<V1>,<V2>,...\n" unless @field==3;
  my $id=uc($field[0]);
  die "invalid ChemistryGroup name '$field[0]' at $location\n"
      unless $id=~/^[A-Z][A-Z0-9_]*$/;
  die "duplicate ChemistryGroup '$id' at $location\n"
      if exists $self->{group_index}->{$id};
  my %option=_parse_options([@field[1..2]],"ChemistryGroup '$id'",$location,
      {DEFAULT=>1,VARIANTS=>1});
  die "ChemistryGroup '$id' at $location requires DEFAULT and VARIANTS\n"
      unless defined($option{DEFAULT}) && defined($option{VARIANTS});
  my @variants=map {uc(_trim($_))} split(/,/,$option{VARIANTS},-1);
  my %seen;
  for my $variant (@variants) {
    die "invalid variant name '$variant' in ChemistryGroup '$id' at $location\n"
        unless $variant=~/^[A-Z][A-Z0-9_]*$/;
    die "duplicate variant '$variant' in ChemistryGroup '$id' at $location\n"
        if $seen{$variant}++;
  }
  die "ChemistryGroup '$id' at $location needs at least two variants\n"
      unless @variants>=2;
  my $default=uc($option{DEFAULT});
  die "DEFAULT=$default of ChemistryGroup '$id' at $location is not one of its VARIANTS\n"
      unless $seen{$default};
  push @{$self->{groups}},{id=>$id,default=>$default,variants=>\@variants,
      location=>$location};
  $self->{group_index}->{$id}=$#{$self->{groups}};
}

sub add_rate_table {
  my ($self,$text,$location)=@_;
  $location='' unless defined $location;
  my @field=map {_trim($_)} split(/\|/,_trim($text),-1);
  die "invalid ChemistryRateTable at $location: expected " .
      "TABLE | MODEL=<MODEL> | SCALE_S=<s> | ION_FRACTION=<f>\n" unless @field==4;
  my $id=uc($field[0]);
  die "invalid ChemistryRateTable name '$field[0]' at $location\n"
      unless $id=~/^[A-Z][A-Z0-9_]*$/;
  die "duplicate ChemistryRateTable '$id' at $location\n"
      if exists $self->{table_index}->{$id};
  my %option=_parse_options([@field[1..3]],"ChemistryRateTable '$id'",$location,
      {MODEL=>1,SCALE_S=>1,ION_FRACTION=>1});
  for my $key ('MODEL','SCALE_S','ION_FRACTION') {
    die "ChemistryRateTable '$id' at $location requires $key\n"
        unless defined $option{$key};
  }
  my $model=uc($option{MODEL});
  die "unknown MODEL=$model of ChemistryRateTable '$id' at $location; known: " .
      join(', ',sort keys %TABLE_MODEL_ROLES)."\n"
      unless exists $TABLE_MODEL_ROLES{$model};
  for my $key ('SCALE_S','ION_FRACTION') {
    die "invalid numeric $key '$option{$key}' of ChemistryRateTable '$id' at $location\n"
        unless _is_number($option{$key});
  }
  die "SCALE_S of ChemistryRateTable '$id' at $location must be positive\n"
      unless $option{SCALE_S}>0.0;
  die "ION_FRACTION of ChemistryRateTable '$id' at $location must be non-negative\n"
      unless $option{ION_FRACTION}>=0.0;
  push @{$self->{tables}},{id=>$id,model=>$model,scale_s=>$option{SCALE_S}+0.0,
      ion_fraction=>$option{ION_FRACTION}+0.0,location=>$location};
  $self->{table_index}->{$id}=$#{$self->{tables}};
}

# Check every group, variant, and rate-table reference after all records have
# been read (records may appear in any order).  Dies with the record location
# and the violated rule.
sub validate_variants {
  my ($self)=@_;
  my %members;     # group -> variant -> [reactions]
  my %table_used;
  for my $reaction (@{$self->{reactions}}) {
    my $where="reaction '$reaction->{id}' at $reaction->{location}";
    if (defined $reaction->{group}) {
      my $group=$reaction->{group};
      die "$where names undeclared ChemistryGroup '$group'\n"
          unless exists $self->{group_index}->{$group};
      my $record=$self->{groups}->[$self->{group_index}->{$group}];
      die "$where names variant '$reaction->{variant}', which ChemistryGroup '$group' does not list\n"
          unless grep {$_ eq $reaction->{variant}} @{$record->{variants}};
      # A variant selects mutually exclusive channels of one parent; a
      # unit-unimolecular reaction is the only kind the particle executor
      # treats as a parent channel.
      die "$where: grouped reactions must have exactly one reactant with coefficient 1\n"
          unless @{$reaction->{reactants}}==1 &&
              $reaction->{reactants}->[0]->{coefficient}==1.0;
      push @{$members{$group}->{$reaction->{variant}}},$reaction;
    }
    if (defined $reaction->{rate_table}) {
      my $table=$reaction->{rate_table};
      die "$where names undeclared ChemistryRateTable '$table'\n"
          unless exists $self->{table_index}->{$table};
      my $record=$self->{tables}->[$self->{table_index}->{$table}];
      my $roles=$TABLE_MODEL_ROLES{$record->{model}};
      die "$where: role '$reaction->{rate_role}' is not defined for MODEL=$record->{model}; known: " .
          join(', ',sort keys %$roles)."\n"
          unless exists $roles->{$reaction->{rate_role}};
      # The table gives the rate at 1 AU; the heliocentric law with RREF_AU=1
      # applies the (1 AU/r_h)^2 scaling.
      die "$where: a table rate requires LAW=HELIOCENTRIC and RREF_AU=1\n"
          unless $reaction->{law} eq 'HELIOCENTRIC' && $reaction->{rref}==1.0;
      die "$where: an IONIZATION table rate requires PRODUCTS=NONE\n"
          if $reaction->{rate_role} eq 'IONIZATION' && @{$reaction->{products}};
      die "$where: a DISSOCIATION table rate requires products\n"
          if $reaction->{rate_role} eq 'DISSOCIATION' && !@{$reaction->{products}};
      $table_used{$table}=1;
    }
  }
  for my $group (@{$self->{groups}}) {
    my %parent;
    for my $variant (@{$group->{variants}}) {
      die "variant '$variant' of ChemistryGroup '$group->{id}' at $group->{location} has no reaction\n"
          unless $members{$group->{id}}->{$variant};
      $parent{$_->{reactants}->[0]->{species}}=1
          for @{$members{$group->{id}}->{$variant}};
    }
    die "the reactions of ChemistryGroup '$group->{id}' at $group->{location} must share one parent species\n"
        unless keys(%parent)==1;
  }
  for my $table (@{$self->{tables}}) {
    die "ChemistryRateTable '$table->{id}' at $table->{location} is not used by any reaction\n"
        unless $table_used{$table->{id}};
  }
}

sub _parse_options {
  my ($fields,$what,$location,$allowed)=@_;
  my %option;
  for my $field (@$fields) {
    die "invalid option '$field' of $what at $location: expected KEY=VALUE\n"
        unless $field=~/^([A-Za-z_]+)\s*=\s*(.*)$/;
    my ($key,$value)=(uc($1),_trim($2));
    die "unknown option '$key' of $what at $location\n"
        unless $allowed->{$key};
    die "duplicate option '$key' of $what at $location\n"
        if exists $option{$key};
    die "empty option '$key' of $what at $location\n" if $value eq '';
    $option{$key}=$value;
  }
  return %option;
}

sub add_reaction {
  my ($self,$text,$location)=@_;
  $location='' unless defined $location;
  $text=_trim($text);

  my @field=map {_trim($_)} split(/\|/,$text,-1);
  die "invalid ChemistryReaction at $location: expected " .
      "ID | REACTANTS | PRODUCTS | LAW | A | BETA | TA_K | RREF_AU | EXCESS_EV " .
      "[| GROUP=<GROUP>:<VARIANT>] [| RATE=TABLE:<TABLE>:<ROLE>]\n"
      unless @field>=9 && @field<=11;

  my ($id,$reactants,$products,$law,$a,$beta,$ta,$rref,$energy_ev)=@field[0..8];
  my %extra=_parse_options([@field[9..$#field]],"reaction '$id'",$location,
      {GROUP=>1,RATE=>1});
  my ($group,$variant,$rate_table,$rate_role);
  if (defined $extra{GROUP}) {
    die "GROUP of reaction '$id' at $location must be <GROUP>:<VARIANT>\n"
        unless $extra{GROUP}=~/^([A-Za-z][A-Za-z0-9_]*):([A-Za-z][A-Za-z0-9_]*)$/;
    ($group,$variant)=(uc($1),uc($2));
  }
  if (defined $extra{RATE}) {
    die "RATE of reaction '$id' at $location must be TABLE:<TABLE>:<ROLE>\n"
        unless $extra{RATE}=~/^TABLE:([A-Za-z][A-Za-z0-9_]*):([A-Za-z][A-Za-z0-9_]*)$/i;
    ($rate_table,$rate_role)=(uc($1),uc($2));
  }
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
    group=>$group,
    variant=>$variant,
    rate_table=>$rate_table,
    rate_role=>$rate_role,
    location=>$location,
  };

  push @{$self->{reactions}},$reaction;
  $self->{ids}->{$id}=1;
}

# True when any reaction takes its rate from a ChemistryRateTable.
sub has_table_rates {
  my ($self)=@_;
  for my $reaction (@{$self->{reactions}}) {
    return 1 if defined $reaction->{rate_table};
  }
  return 0;
}

sub reaction_count {
  my ($self)=@_;
  return scalar @{$self->{reactions}};
}

# Return a defensive, read-only-style description for input policies that are
# resolved after named databases and explicit reactions have both loaded.
# Stage-6 packet selection needs stable IDs, nominal rates, and concrete
# product stoichiometry but must not depend on this module's private hashes or
# alter the reaction order later emitted by write_header().
sub packet_policy_reactions {
  my ($self)=@_;
  my @result;
  for my $reaction (@{$self->{reactions}}) {
    push @result,{
      id=>$reaction->{id},
      rate=>$reaction->{a},
      products=>[
        map {{species=>$_->{species},coefficient=>$_->{coefficient}}}
            @{$reaction->{products}}
      ],
    };
  }
  return \@result;
}

# Return the stable IDs, nominal rates, and reactant stoichiometry needed by
# the Stage-7 channel-importance resolver.  Keep this separate from the packet
# descriptor above: packet policy is product-oriented, whereas importance
# policy groups complete positive-rate unit-unimolecular sibling families by
# parent.  The returned hashes/arrays are fresh defensive copies, so a policy
# module cannot mutate the database or generated reaction order.
sub channel_importance_policy_reactions {
  my ($self)=@_;
  my @result;
  for my $reaction (@{$self->{reactions}}) {
    push @result,{
      id=>$reaction->{id},
      rate=>$reaction->{a},
      reactants=>[
        map {{species=>$_->{species},coefficient=>$_->{coefficient}}}
            @{$reaction->{reactants}}
      ],
    };
  }
  return \@result;
}

# Return a boolean species table reachable from the supplied source indices.
# A general reaction becomes reachable only when every reactant species is
# reachable; this reduces to ordinary parent-to-daughter traversal for the
# unit-unimolecular network used by the current particle executor while
# remaining meaningful for future continuum chemistry backends.  A zero-rate
# channel is excluded because it cannot create a chemical-only population.
sub reachable_species_from {
  my ($self,$roots)=@_;
  die "reachable_species_from requires a root-species array\n"
      unless ref($roots) eq 'ARRAY';

  my $species_count=scalar keys %{$self->{species_index}};
  my @reachable=(0)x$species_count;
  for my $root (@$roots) {
    die "invalid reachable-species root '$root'\n"
        unless defined($root) && $root=~/^\d+$/ &&
            $root>=0 && $root<$species_count;
    $reachable[$root]=1;
  }

  my $changed=1;
  while ($changed) {
    $changed=0;
    for my $reaction (@{$self->{reactions}}) {
      next unless $reaction->{a}>0.0;
      my $enabled=1;
      for my $reactant (@{$reaction->{reactants}}) {
        if (!$reachable[$reactant->{species}]) {
          $enabled=0;
          last;
        }
      }
      next unless $enabled;
      for my $product (@{$reaction->{products}}) {
        if (!$reachable[$product->{species}]) {
          $reachable[$product->{species}]=1;
          $changed=1;
        }
      }
    }
  }
  return \@reachable;
}

sub write_header {
  my ($self,$path)=@_;
  $self->validate_variants();
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
  print $out "};\n\n";

  # Chemistry variants (srcHartley2/CHEMISTRY_VARIANTS.md).  One entry per
  # reaction in the order above; Group/Variant/RateTable are indices into the
  # arrays below (-1 = none), RateRole is Comet::Chemistry::RateRole.
  my @groups=@{$self->{groups}};
  my @tables=@{$self->{tables}};
  print $out "static const std::size_t GroupCount=".(scalar @groups).";\n";
  for (my $g=0;$g<@groups;$g++) {
    print $out "static const char *const GroupVariants$g\[\] = {".
        join(',',map {"\"$_\""} @{$groups[$g]->{variants}})."};\n";
  }
  print $out "static const Comet::Chemistry::cCompiledGroup Groups[] = {\n";
  if (@groups) {
    for (my $g=0;$g<@groups;$g++) {
      my $group=$groups[$g];
      my ($default)=grep {$group->{variants}->[$_] eq $group->{default}}
          0..$#{$group->{variants}};
      print $out "  {\"$group->{id}\",$default,".
          scalar(@{$group->{variants}}).",GroupVariants$g},\n";
    }
  }
  else {
    print $out "  {\"none\",0,0,0}\n";
  }
  print $out "};\n";
  print $out "static const std::size_t RateTableCount=".(scalar @tables).";\n";
  print $out "static const Comet::Chemistry::cCompiledRateTable RateTables[] = {\n";
  if (@tables) {
    for my $table (@tables) {
      print $out "  {\"$table->{id}\",\"$table->{model}\",".
          _number($table->{scale_s}).","._number($table->{ion_fraction})."},\n";
    }
  }
  else {
    print $out "  {\"none\",\"none\",0.0,0.0}\n";
  }
  print $out "};\n";
  print $out "static const Comet::Chemistry::cCompiledReactionVariant ReactionVariants[] = {\n";
  if (@reactions) {
    for my $reaction (@reactions) {
      my ($group,$variant,$table,$role)=(-1,-1,-1,0);
      if (defined $reaction->{group}) {
        $group=$self->{group_index}->{$reaction->{group}};
        ($variant)=grep {$groups[$group]->{variants}->[$_] eq $reaction->{variant}}
            0..$#{$groups[$group]->{variants}};
      }
      if (defined $reaction->{rate_table}) {
        $table=$self->{table_index}->{$reaction->{rate_table}};
        $role=$TABLE_MODEL_ROLES{$tables[$table]->{model}}->{$reaction->{rate_role}};
      }
      print $out "  {$group,$variant,$table,$role},\n";
    }
  }
  else {
    print $out "  {-1,-1,-1,0}\n";
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
