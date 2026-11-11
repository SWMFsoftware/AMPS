package CometDaughterPacketPolicy;

use strict;
use warnings;
use POSIX qw(LONG_MAX);

# Resolve Stage-6 daughter-packet requests into one immutable count per
# generated reaction.  A complete packet contains every concrete
# stoichiometric product, so the policy is necessarily reaction-level even
# when the input selects a daughter species.  The approved precedence is:
#
#   explicit REACTION:id > maximum matching SPECIES:name > CHEMICAL_ONLY
#       > GLOBAL > implicit one packet.
#
# Taking the maximum of matching product-species requests is intentional: one
# reaction cannot use different packet counts for its coproducts without
# breaking per-packet momentum and excess-energy conservation.

use constant {
  POLICY_IMPLICIT      => 0,
  POLICY_GLOBAL        => 1,
  POLICY_CHEMICAL_ONLY => 2,
  POLICY_SPECIES       => 3,
  POLICY_REACTION      => 4,
};

sub new {
  my ($class,%args)=@_;
  die "CometDaughterPacketPolicy requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';

  my (%species_index,@species);
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=uc(_trim($args{species}->[$i]));
    die "empty species in daughter-packet policy\n" if $name eq '';
    die "duplicate species '$name' in daughter-packet policy\n"
        if exists $species_index{$name};
    push @species,$name;
    $species_index{$name}=$i;
  }

  return bless {
    species=>\@species,
    species_index=>\%species_index,
    enabled=>0,
    mode_is_set=>0,
    maximum=>undef,
    global=>undef,
    chemical_only=>undef,
    species_rules=>{},
    reaction_rules=>{},
    record_count=>0,
  },$class;
}

sub set_mode {
  my ($self,$value,$location)=@_;
  $location='' unless defined $location;
  $value=uc(_trim($value));
  die "ChemistryDaughterPacketMode must be ON or OFF at $location\n"
      unless $value eq 'ON' || $value eq 'OFF';
  die "duplicate ChemistryDaughterPacketMode at $location\n"
      if $self->{mode_is_set};
  $self->{enabled}=($value eq 'ON') ? 1 : 0;
  $self->{mode_is_set}=1;
}

sub is_enabled {
  my ($self)=@_;
  return $self->{enabled};
}

sub set_maximum {
  my ($self,$value,$location)=@_;
  $location='' unless defined $location;
  $value=_trim($value);
  die "ChemistryDaughterPacketMaximum must be an integer of at least two at $location\n"
      unless $value=~/^\d+$/ && $value>=2 && $value<=LONG_MAX;
  die "duplicate ChemistryDaughterPacketMaximum at $location\n"
      if defined $self->{maximum};
  $self->{maximum}=0+$value;
}

sub maximum {
  my ($self)=@_;
  die "daughter-packet maximum is not configured\n"
      unless defined $self->{maximum};
  return $self->{maximum};
}

sub add_count_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid ChemistryDaughterPacketCount at $location\n"
      unless $record=~/^\s*ChemistryDaughterPacketCount\s*(?:\(\s*([^)]+?)\s*\))?\s*[=:]\s*(.*?)\s*$/i;
  my ($selector,$count)=($1,_trim($2));
  die "ChemistryDaughterPacketCount must be a positive integer at $location\n"
      unless $count=~/^\d+$/ && $count>=1 && $count<=LONG_MAX;
  $count=0+$count;

  my ($kind,$key);
  if (!defined($selector) || uc(_trim($selector)) eq 'GLOBAL') {
    ($kind,$key)=('global','GLOBAL');
  }
  else {
    $selector=_trim($selector);
    if (uc($selector) eq 'CHEMICAL_ONLY') {
      ($kind,$key)=('chemical_only','CHEMICAL_ONLY');
    }
    elsif ($selector=~/^SPECIES\s*:\s*(.+)$/i) {
      $key=uc(_trim($1));
      die "unknown daughter-packet species '$key' at $location\n"
          unless exists $self->{species_index}->{$key};
      $kind='species';
    }
    elsif ($selector=~/^REACTION\s*:\s*(.+)$/i) {
      $key=_trim($1);
      die "empty daughter-packet reaction ID at $location\n" if $key eq '';
      $kind='reaction';
    }
    else {
      die "daughter-packet selector must be GLOBAL, CHEMICAL_ONLY, SPECIES:name, or REACTION:id at $location\n";
    }
  }

  my $entry={count=>$count,location=>$location,selector=>$key};
  if ($kind eq 'global') {
    _store_unique($self,'global',$key,$entry);
  }
  elsif ($kind eq 'chemical_only') {
    _store_unique($self,'chemical_only',$key,$entry);
  }
  elsif ($kind eq 'species') {
    die "duplicate daughter-packet policy for '$key' at $location\n"
        if exists $self->{species_rules}->{$key};
    $self->{species_rules}->{$key}=$entry;
  }
  else {
    die "duplicate daughter-packet policy for reaction '$key' at $location\n"
        if exists $self->{reaction_rules}->{$key};
    $self->{reaction_rules}->{$key}=$entry;
  }
  $self->{record_count}++;
}

sub validate_configuration_and_resolve {
  my ($self,%args)=@_;
  my $reactions=$args{reactions};
  my $chemical_only=$args{chemical_only};
  my $transaction_limit=$args{transaction_limit};
  die "daughter-packet resolution requires reaction descriptors\n"
      unless ref($reactions) eq 'ARRAY';
  die "daughter-packet resolution requires chemical-only flags\n"
      unless ref($chemical_only) eq 'ARRAY' &&
          @$chemical_only==@{$self->{species}};

  # Treat the resolver's descriptor API as a checked boundary even though the
  # production caller currently obtains descriptors from the trusted reaction
  # database.  This keeps malformed tests and future database implementations
  # from turning an invalid species index into an unrelated Perl array lookup.
  my %reaction_index;
  for (my $i=0;$i<@$reactions;$i++) {
    my $reaction=$reactions->[$i];
    my $id=$reaction->{id};
    die "invalid reaction descriptor at index $i\n"
        unless defined($id) && $id ne '' &&
            ref($reaction->{products}) eq 'ARRAY';
    die "duplicate reaction ID '$id' in packet-policy descriptors\n"
        if exists $reaction_index{$id};
    $reaction_index{$id}=$i;
    for my $product (@{$reaction->{products}}) {
      my $species=$product->{species};
      die "reaction '$id' has an invalid product species index\n"
          unless defined($species) && $species=~/^\d+$/ &&
              $species<@{$self->{species}};
    }
    _concrete_product_count($reaction);
  }

  if (!$self->is_enabled()) {
    die "daughter-packet records and maximum require ChemistryDaughterPacketMode=ON\n"
        if $self->{record_count}!=0 || defined $self->{maximum};
    return {
      packet_count=>[(1)x scalar(@$reactions)],
      policy_source=>[(POLICY_IMPLICIT)x scalar(@$reactions)],
      selector=>[('IMPLICIT')x scalar(@$reactions)],
      transaction_count=>[map {_concrete_product_count($_)} @$reactions],
    };
  }

  die "ChemistryDaughterPacketMode=ON requires ChemistryDaughterPacketMaximum\n"
      unless defined $self->{maximum};
  die "ChemistryDaughterPacketMode=ON requires a positive transaction limit\n"
      unless defined($transaction_limit) && $transaction_limit=~/^\d+$/ &&
          $transaction_limit>0;

  for my $id (keys %{$self->{reaction_rules}}) {
    die "unknown daughter-packet reaction ID '$id'\n"
        unless exists $reaction_index{$id};
  }

  my (@count,@source,@selector,@transaction);
  my $has_split=0;
  for my $reaction (@$reactions) {
    my $concrete=_concrete_product_count($reaction);
    my ($selected_count,$selected_source,$selected_selector)=
        (1,POLICY_IMPLICIT,'IMPLICIT');

    if ($concrete>0) {
      if (defined $self->{global}) {
        ($selected_count,$selected_source,$selected_selector)=
            ($self->{global}->{count},POLICY_GLOBAL,'GLOBAL');
      }
      my $matches_chemical=0;
      for my $product (@{$reaction->{products}}) {
        $matches_chemical=1 if $chemical_only->[$product->{species}];
      }
      if ($matches_chemical && defined $self->{chemical_only}) {
        ($selected_count,$selected_source,$selected_selector)=
            ($self->{chemical_only}->{count},POLICY_CHEMICAL_ONLY,
             'CHEMICAL_ONLY');
      }

      my @species_matches;
      for my $product (@{$reaction->{products}}) {
        my $name=$self->{species}->[$product->{species}];
        push @species_matches,$self->{species_rules}->{$name}
            if exists $self->{species_rules}->{$name};
      }
      if (@species_matches) {
        my $maximum=1;
        my @names;
        for my $entry (@species_matches) {
          $maximum=$entry->{count} if $entry->{count}>$maximum;
          push @names,$entry->{selector};
        }
        my %seen;
        @names=sort grep {!$seen{$_}++} @names;
        ($selected_count,$selected_source,$selected_selector)=
            ($maximum,POLICY_SPECIES,'SPECIES:'.join(',',@names));
      }
    }

    if (exists $self->{reaction_rules}->{$reaction->{id}}) {
      my $entry=$self->{reaction_rules}->{$reaction->{id}};
      if ($concrete==0) {
        die "daughter-packet policy targets product-free sink reaction '$reaction->{id}'\n";
      }
      if (!(defined($reaction->{rate}) && $reaction->{rate}>0.0) &&
          $entry->{count}>1) {
        die "daughter-packet count greater than one targets zero-rate reaction '$reaction->{id}'\n";
      }
      ($selected_count,$selected_source,$selected_selector)=
          ($entry->{count},POLICY_REACTION,'REACTION:'.$reaction->{id});
    }

    die "reaction '$reaction->{id}' packet count $selected_count exceeds configured maximum $self->{maximum}\n"
        if $selected_count>$self->{maximum};
    die "reaction '$reaction->{id}' daughter transaction overflows integer range\n"
        if $concrete>0 && $selected_count>LONG_MAX/$concrete;
    my $total=$concrete*$selected_count;
    die "reaction '$reaction->{id}' requests $total daughters above transaction limit $transaction_limit\n"
        if $total>$transaction_limit;
    # A broad group/species rule may also match a database channel whose rate
    # is deliberately zero for this case.  Preserve its deterministic table
    # entry, but do not let a reaction that can never execute satisfy the
    # required effective-split gate.
    $has_split=1 if $selected_count>1 &&
        defined($reaction->{rate}) && $reaction->{rate}>0.0;
    push @count,$selected_count;
    push @source,$selected_source;
    push @selector,$selected_selector;
    push @transaction,$total;
  }
  die "ChemistryDaughterPacketMode=ON has no effective packet count above one\n"
      unless $has_split;

  return {
    packet_count=>\@count,policy_source=>\@source,
    selector=>\@selector,transaction_count=>\@transaction,
  };
}

sub write_header {
  my ($self,$path,$resolved)=@_;
  die "daughter-packet header requires a resolved policy\n"
      unless ref($resolved) eq 'HASH';
  my $count=$resolved->{packet_count};
  my $source=$resolved->{policy_source};
  die "daughter-packet generated arrays differ in length\n"
      unless ref($count) eq 'ARRAY' && ref($source) eq 'ARRAY' &&
          @$count==@$source;
  for (my $i=0;$i<@$count;$i++) {
    die "invalid generated daughter-packet count at index $i\n"
        unless defined($count->[$i]) && $count->[$i]=~/^\d+$/ &&
            $count->[$i]>=1 && $count->[$i]<=LONG_MAX;
    die "invalid generated daughter-packet policy source at index $i\n"
        unless defined($source->[$i]) && $source->[$i]=~/^\d+$/ &&
            $source->[$i]>=POLICY_IMPLICIT &&
            $source->[$i]<=POLICY_REACTION;
  }
  open(my $out,'>',$path) or die "cannot write '$path': $!\n";
  print $out "// Generated from the comet application input. Do not edit.\n";
  print $out "#ifndef _COMET_DAUGHTER_PACKET_POLICY_GENERATED_H_\n";
  print $out "#define _COMET_DAUGHTER_PACKET_POLICY_GENERATED_H_\n\n";
  print $out "#include <cstddef>\n\n";
  print $out "namespace Comet { namespace Chemistry { namespace DaughterPacketGenerated {\n";
  print $out "static const std::size_t ReactionCount=".scalar(@$count).";\n";
  my @generated_count=@$count ? @$count : (1);
  my @generated_source=@$source ? @$source : (POLICY_IMPLICIT);
  print $out "static const long int PacketCount[]={".
      join(',',@generated_count)."};\n";
  print $out "static const int PolicySource[]={".
      join(',',@generated_source)."};\n";
  print $out "}}} // namespace Comet::Chemistry::DaughterPacketGenerated\n\n";
  print $out "#endif\n";
  close($out) or die "cannot close '$path': $!\n";
}

sub _store_unique {
  my ($self,$slot,$key,$entry)=@_;
  die "duplicate daughter-packet policy for '$key' at $entry->{location}\n"
      if defined $self->{$slot};
  $self->{$slot}=$entry;
}

sub _concrete_product_count {
  my ($reaction)=@_;
  die "reaction descriptor has no product array\n"
      unless ref($reaction->{products}) eq 'ARRAY';
  my $count=0;
  for my $product (@{$reaction->{products}}) {
    my $coefficient=$product->{coefficient};
    die "packet policy requires positive integer product coefficients\n"
        unless defined($coefficient) && $coefficient>0 &&
            abs($coefficient-int($coefficient+0.5))<=1.0e-12;
    my $integer=int($coefficient+0.5);
    die "concrete product count overflows integer range\n"
        if $count>LONG_MAX-$integer;
    $count+=$integer;
  }
  return $count;
}

sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
