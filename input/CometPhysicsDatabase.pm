package CometPhysicsDatabase;

use strict;
use warnings;

# This module is the input-side catalog for comet physics that is independent
# of a particular nucleus shape or source map.  It deliberately stores names,
# not AMPS species numbers: the selected SpeciesList is the authority for the
# compile-time indices, and this module resolves the names only while writing
# generated C++ tables.
#
# Rubin et al. (2011), Icarus 213, 655--677, is currently the first catalog.
# Its Table 3 contains only parent-parent neutral collision cross sections;
# missing daughter pairs remain missing and must never be replaced by guessed
# values.  Table 4 contains photodissociation and photoionization frequencies
# at 1 AU.  Ion products are outside the neutral model, so those channels are
# represented by the explicit product token NONE.

sub new {
  my ($class,%args)=@_;
  die "CometPhysicsDatabase requires a species list\n"
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
    collisions=>[],
    database_name=>'NONE',
  },$class;
}

# Add one input-defined collision without attaching it to a comet name.  The
# compact constant form is
#   SPECIES0 | SPECIES1 | CONSTANT | SIGMA_M2 | SOURCE
# and a bounded power law is
#   SPECIES0 | SPECIES1 | TEMPERATURE_POWER | SIGMA_REF_M2 |
#   TREF_K | EXPONENT | TMIN_K | TMAX_K | SOURCE
# The validity bounds are mandatory for a temperature law so input cannot
# request an undocumented extrapolation.
sub add_collision {
  my ($self,$text,$location)=@_;
  $location='' unless defined $location;
  my @field=map {_trim($_)} split(/\|/,$text,-1);
  die "NeutralCollision at $location requires at least five fields\n"
      unless @field>=5;

  my ($name0,$name1,$law)=@field[0,1,2];
  $law=uc($law);
  if ($law eq 'CONSTANT') {
    die "constant NeutralCollision at $location expects SPECIES0 | SPECIES1 | CONSTANT | SIGMA_M2 | SOURCE\n"
        unless @field==5;
    die "invalid collision cross section '$field[3]' at $location\n"
        unless _is_number($field[3]) && $field[3]>0.0;
    die "collision source is empty at $location\n" if $field[4] eq '';
    $self->_add_collision($name0,$name1,'Constant',$field[3]+0.0,
        1.0,0.0,0.0,0.0,$field[4]);
  }
  elsif ($law eq 'TEMPERATURE_POWER') {
    die "temperature-power NeutralCollision at $location expects SPECIES0 | SPECIES1 | TEMPERATURE_POWER | SIGMA_REF_M2 | TREF_K | EXPONENT | TMIN_K | TMAX_K | SOURCE\n"
        unless @field==9;
    for my $index (3..7) {
      die "invalid collision number '$field[$index]' at $location\n"
          unless _is_number($field[$index]);
    }
    die "temperature-power collision at $location requires positive sigma, reference temperature, and minimum temperature, with TMAX>=TMIN\n"
        unless $field[3]>0.0 && $field[4]>0.0 && $field[6]>0.0 &&
            $field[7]>=$field[6];
    die "collision source is empty at $location\n" if $field[8] eq '';
    $self->_add_collision($name0,$name1,'TemperaturePowerLaw',
        map({$_+0.0} @field[3..7]),$field[8]);
  }
  else {
    die "unknown NeutralCollision law '$law' at $location\n";
  }

  $self->{database_name}='CUSTOM'
      if $self->{database_name} eq 'NONE';
}

sub load_rubin2011 {
  my ($self,%args)=@_;
  my $network=$args{network};
  die "Rubin2011 requires a CometReactionNetwork instance\n"
      unless defined $network && $network->can('add_reaction');

  # Table 4 prints this OH branch as an interval rather than a single value.
  # Requiring the caller to choose a value prevents the database layer from
  # silently selecting an endpoint or midpoint.  The admissible interval is
  # copied verbatim from the paper (all rates are at 1 AU).
  my $oh_ground_rate=$args{oh_ground_rate};
  die "Rubin2011OHGroundRate is required; Table 4 gives 3e-6 to 6e-6 s^-1\n"
      unless _is_number($oh_ground_rate);
  $oh_ground_rate+=0.0;
  die "Rubin2011OHGroundRate must be in the published [3e-6,6e-6] s^-1 interval\n"
      unless $oh_ground_rate>=3.0e-6 && $oh_ground_rate<=6.0e-6;

  # Footnote (a) says the authors also ran a case with no H2CO-forming
  # methanol branch.  Both variants are paper cases, so the input must select
  # one explicitly rather than inheriting a library default.
  die "Rubin2011MethanolH2COBranch must be ON or OFF\n"
      unless defined $args{methanol_h2co_branch} &&
          ($args{methanol_h2co_branch} eq 'ON' ||
           $args{methanol_h2co_branch} eq 'OFF');

  # User records may supplement the paper's parent-only matrix with a cited
  # daughter dataset.  Keep them and let C++ validation reject duplicates.
  $self->{database_name}=($self->{database_name} eq 'CUSTOM') ?
      'RUBIN2011+CUSTOM' : 'RUBIN2011';
  $self->_load_rubin2011_parent_collisions();

  # Each row below is one physical channel in Table 4.  Excited electronic
  # states are collapsed onto the corresponding neutral AMPS species because
  # the paper follows the chemical species but not separate electronic-state
  # particles.  Their state-dependent energy is retained as channel excess
  # energy.  HELIOCENTRIC means inverse-square scaling from the stated 1 AU
  # reference distance.
  my @reaction=(
    ['CO2_CO_O3P','CO2','CO + O',9.20e-7,4.30],
    ['CO2_CO_O1D','CO2','CO + O',2.80e-7,2.00],
    ['CO2_ION','CO2','NONE',8.03e-7,0.0],

    ['CH3OH_CH3_OH','CH3OH','CH3 + OH',5.60e-7,5.00],
    ['CH3OH_ION','CH3OH','NONE',7.30e-7,0.0],

    ['H2CO_HCO_H','H2CO','HCO + H',6.60e-5,0.39],
    ['H2CO_CO_H2','H2CO','CO + H2',1.40e-4,2.10],
    ['H2CO_CO_2H','H2CO','CO + 2*H',3.20e-5,3.00],
    ['H2CO_ION','H2CO','NONE',7.20e-7,0.0],

    ['C2H6_C2H4_H2','C2H6','C2H4 + H2',3.70e-6,9.00],
    ['C2H6_2CH3','C2H6','2*CH3',8.80e-7,7.40],
    ['C2H6_C2H5_H','C2H6','C2H5 + H',3.30e-6,6.80],
    ['C2H6_CH4_CH2','C2H6','CH4 + CH2',2.20e-6,6.10],
    ['C2H6_ION','C2H6','NONE',4.90e-7,0.0],

    ['HCO_CH_O','HCO','CH + O',3.20e-5,3.00],
    ['HCO_OH_C','HCO','OH + C',3.20e-5,3.00],
    ['HCO_ION','HCO','NONE',1.00e-6,0.0],
    ['C2H5_C2H4_H','C2H5','C2H4 + H',3.30e-6,6.80],
    ['C2H5_ION','C2H5','NONE',4.90e-7,0.0],

    ['CO_C3P_O3P','CO','C + O',7.49e-7,2.60],
    ['CO_C1D_O1D','CO','C + O',9.36e-8,2.30],
    ['CO_ION','CO','NONE',1.16e-6,0.0],
    ['C2H4_C2H2_H2','C2H4','C2H2 + H2',2.40e-5,6.20],
    ['C2H4_C2H2_2H','C2H4','C2H2 + 2*H',2.30e-5,1.70],
    ['C2H4_ION','C2H4','NONE',1.00e-6,0.0],

    ['HCN_CN_H','HCN','CN + H',1.30e-5,3.80],
    ['HCN_ION','HCN','NONE',4.50e-7,0.0],
    ['CN_C_N','CN','C + N',3.20e-6,7.40],
    ['C2H2_C2H_H','C2H2','C2H + H',1.00e-5,3.20],
    ['C2H2_C2_H2','C2H2','C2 + H2',2.70e-6,3.10],
    ['C2H2_ION','C2H2','NONE',8.50e-7,0.0],
    ['C2H_C2_H','C2H','C2 + H',7.00e-6,3.00],
    ['C2H_ION','C2H','NONE',1.00e-6,0.0],
    ['C2_2C','C2','2*C',1.00e-7,3.60],
    ['C2_ION','C2','NONE',9.10e-7,0.0],

    ['H2O_OH_H','H2O','OH + H',1.00e-5,1.70],
    ['H2O_O_H2','H2O','O + H2',6.00e-7,3.80],
    ['H2O_O_2H','H2O','O + 2*H',7.50e-7,0.70],
    ['H2O_ION','H2O','NONE',4.04e-7,0.0],
    ['OH_O_H_GROUND','OH','O + H',$oh_ground_rate,0.36],
    ['OH_O_H_067','OH','O + H',5.00e-7,0.67],
    ['OH_O_H_269','OH','O + H',1.40e-6,2.69],
    ['OH_O_H_384','OH','O + H',3.00e-7,3.84],
    ['OH_ION','OH','NONE',1.20e-5,0.0],

    ['NH3_NH2_H','NH3','NH2 + H',1.70e-4,1.80],
    ['NH3_NH_H2','NH3','NH + H2',4.00e-6,1.70],
    ['NH3_NH_2H','NH3','NH + 2*H',2.00e-6,2.10],
    ['NH3_ION','NH3','NONE',8.00e-7,0.0],
    ['O_ION','O','NONE',2.00e-7,0.0],
    ['NH2_NH_H','NH2','NH + H',2.10e-6,6.40],

    ['CH4_CH3_H','CH4','CH3 + H',2.60e-7,6.50],
    ['CH4_CH2_H2','CH4','CH2 + H2',4.00e-6,5.30],
    ['CH4_CH_H2_H','CH4','CH + H2 + H',6.40e-7,1.70],
    ['CH4_CH2_2H','CH4','CH2 + 2*H',2.10e-6,0.83],
    ['CH4_ION','CH4','NONE',5.94e-7,0.0],
    ['NH_N_H','NH','N + H',5.75e-5,1.25],
    ['CH3_CH2_H','CH3','CH2 + H',7.00e-4,4.00],
    ['CH3_CH_H2','CH3','CH + H2',5.00e-5,3.00],
    ['CH3_ION','CH3','NONE',9.38e-7,0.0],
    ['N_ION','N','NONE',1.90e-7,0.0],
    ['CH2_CH_H','CH2','CH + H',1.00e-4,4.00],
    ['CH2_C_H2','CH2','C + H2',3.00e-5,3.00],
    ['CH2_ION','CH2','NONE',1.05e-6,0.0],
    ['CH_C_H','CH','C + H',9.20e-3,0.45],
    ['C_ION','C','NONE',4.10e-7,0.0],
    ['H2_2H_820','H2','2*H',4.80e-8,8.20],
    ['H2_2H_044','H2','2*H',3.40e-8,0.44],
    ['H2_ION','H2','NONE',6.35e-8,0.0],
    ['H_ION','H','NONE',7.30e-8,0.0],
  );

  if ($args{methanol_h2co_branch} eq 'ON') {
    # Keep this branch next to the other methanol channels in the generated
    # network.  Its placement in the flat table has no physical significance.
    unshift @reaction,
        ['CH3OH_H2CO_H2','CH3OH','H2CO + H2',1.00e-5,7.50];
  }

  for my $entry (@reaction) {
    my ($id,$reactants,$products,$rate,$energy_ev)=@$entry;
    $network->add_reaction(join(' | ',
        $id,$reactants,$products,'HELIOCENTRIC',$rate,0,0,1,$energy_ev),
        'Rubin et al. (2011), Table 4');
  }
}

sub write_collision_header {
  my ($self,$path)=@_;
  open(my $out,'>',$path) or die "cannot write '$path': $!\n";

  print $out "// Generated from the selected comet physics database. Do not edit.\n";
  print $out "// Missing species pairs are intentionally absent; the runtime must report them.\n";
  print $out "#ifndef _COMET_PHYSICS_GENERATED_H_\n";
  print $out "#define _COMET_PHYSICS_GENERATED_H_\n\n";
  print $out "namespace Comet { namespace Physics { namespace Generated {\n";
  print $out "static const char DatabaseName[]=\"$self->{database_name}\";\n";
  print $out "static const std::size_t CollisionCount=".
      scalar(@{$self->{collisions}}).";\n";
  print $out "static const CometaryPhysics::cCollision Collisions[] = {\n";

  if (@{$self->{collisions}}) {
    for my $entry (@{$self->{collisions}}) {
      my $source=$entry->{source};
      $source=~s/\\/\\\\/g;
      $source=~s/"/\\"/g;
      print $out '  {'.join(',',
          $entry->{species0},$entry->{species1},
          'CometaryPhysics::CollisionLaw::'.$entry->{law},
          map({_number($_)} ($entry->{sigma},$entry->{reference_temperature},
              $entry->{temperature_exponent},$entry->{minimum_temperature},
              $entry->{maximum_temperature}))).",\"$source\"},\n";
    }
  }
  else {
    # C++ forbids a standard zero-length array.  CollisionCount remains zero,
    # so the sentinel is never searched or exposed as physical data.
    print $out "  {0,0,CometaryPhysics::CollisionLaw::Constant,1.0,1.0,0.0,0.0,0.0,\"sentinel\"}\n";
  }
  print $out "};\n";
  print $out "}}} // namespace Comet::Physics::Generated\n\n";
  print $out "#endif\n";
  close($out) or die "cannot close '$path': $!\n";
}

sub collision_count {
  my ($self)=@_;
  return scalar @{$self->{collisions}};
}

sub _load_rubin2011_parent_collisions {
  my ($self)=@_;

  # Table 3 is upper triangular in this exact species order.  The values are
  # printed in 10^-15 cm^2, which equals 10^-19 m^2.  The diagonal H2O entry
  # is handled separately because it follows Crifo's temperature law.
  my @species=qw(H2O CO H2CO CO2 CH3OH C2H6 C2H4 C2H2 HCN NH3 CH4);
  my @upper=(
    [undef,3.2,3.6,3.4,3.8,3.8,3.8,3.6,3.4,3.6,3.8],
    [3.2,3.4,3.2,3.6,3.6,3.6,3.4,3.2,3.4,3.6],
    [3.8,3.6,4.0,4.0,4.0,3.8,3.6,3.8,4.0],
    [3.4,3.8,3.8,3.8,3.6,3.4,3.6,3.8],
    [4.2,4.2,4.2,4.0,3.8,4.0,4.2],
    [4.2,4.2,4.0,3.8,4.0,4.2],
    [4.2,4.0,3.8,4.0,4.2],
    [4.0,3.6,3.8,4.0],
    [3.4,3.6,3.8],
    [3.8,4.0],
    [4.2],
  );

  $self->_add_collision('H2O','H2O','TemperaturePowerLaw',
      1.66e-19,300.0,-0.6,50.0,300.0,
      'Crifo (1989), as listed in Rubin et al. (2011), Table 3');

  for (my $row=0;$row<@species;$row++) {
    for (my $offset=0;$offset<@{$upper[$row]};$offset++) {
      next if $row==0 && $offset==0;
      my $column=$row+$offset;
      $self->_add_collision($species[$row],$species[$column],'Constant',
          $upper[$row]->[$offset]*1.0e-19,1.0,0.0,0.0,0.0,
          'Rubin et al. (2011), Table 3');
    }
  }
}

sub _add_collision {
  my ($self,$name0,$name1,$law,$sigma,$reference_temperature,
      $temperature_exponent,$minimum_temperature,$maximum_temperature,
      $source)=@_;
  my $species0=$self->_require_species($name0,'collision database');
  my $species1=$self->_require_species($name1,'collision database');
  push @{$self->{collisions}}, {
    species0=>$species0,
    species1=>$species1,
    law=>$law,
    sigma=>$sigma,
    reference_temperature=>$reference_temperature,
    temperature_exponent=>$temperature_exponent,
    minimum_temperature=>$minimum_temperature,
    maximum_temperature=>$maximum_temperature,
    source=>$source,
  };
}

sub _require_species {
  my ($self,$name,$context)=@_;
  $name=uc($name);
  die "species '$name' required by $context is absent from SpeciesList\n"
      unless exists $self->{species_index}->{$name};
  return $self->{species_index}->{$name};
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
  return defined $value &&
      $value=~/^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?$/;
}

sub _number {
  my ($value)=@_;
  return sprintf('%.17g',$value+0.0);
}

1;
