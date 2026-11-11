package CometBrightnessConfig;

use strict;
use warnings;

# Strict input-side representation of DIM1 resonance-fluorescence output.
# The module deliberately stores no comet or transition constants: every
# fixed photon g factor or named core model and every optional line-center
# optical-depth diagnostic is a case record with provenance.  Runtime output
# is optically thin only; tau is reported and flagged, never converted into an
# invented attenuation law.

use constant MAX_EMISSION_LINES => 16;

sub new {
  my ($class,%args)=@_;
  die "CometBrightnessConfig requires a species list\n"
      unless ref($args{species}) eq 'ARRAY';
  my (@species,%index);
  for (my $i=0;$i<@{$args{species}};$i++) {
    my $name=uc(_trim($args{species}->[$i]));
    die "empty species in brightness species list\n" if $name eq '';
    die "duplicate species '$name' in brightness species list\n"
        if exists $index{$name};
    push @species,$name;
    $index{$name}=$i;
  }
  return bless {
    species=>\@species,species_index=>\%index,mode=>undef,grid=>undef,
    emissions=>[],keys=>{},radial_velocity_km_s=>undef,
    comet_velocity_km_s=>undef,observation_geometry=>undef,
    velocity_domain_policy=>undef,solar_lya_flux_1au_ph_cm2_s=>undef,
  },$class;
}

sub set_comet_velocity_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessCometHeliocentricVelocityKMS at $location\n"
      unless $record=~/^\s*BrightnessCometHeliocentricVelocityKMS\s*[=:]\s*(.*?)\s*$/i;
  die "duplicate comet heliocentric velocity at $location\n"
      if defined($self->{comet_velocity_km_s}) ||
          defined($self->{radial_velocity_km_s});
  $self->{comet_velocity_km_s}=
      _vector3($1,'BrightnessCometHeliocentricVelocityKMS',$location);
}

sub set_observation_geometry_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessObservationGeometry at $location\n"
      unless $record=~/^\s*BrightnessObservationGeometry\s*[=:]\s*(.*?)\s*$/i;
  die "duplicate BrightnessObservationGeometry at $location\n"
      if defined $self->{observation_geometry};
  my @field=map {_trim($_)} split(/\|/,$1,-1);
  die "BrightnessObservationGeometry expects ER_X | ER_Y | ER_Z | LOS_X | LOS_Y | LOS_Z | IMPACT_X | IMPACT_Y | IMPACT_Z at $location\n"
      unless @field==9 && !grep {!_number($_)} @field;
  my @radial=_normalize3([map {$_+0.0} @field[0..2]],
      'heliocentric radial direction',$location);
  my @los=_normalize3([map {$_+0.0} @field[3..5]],
      'line-of-sight direction',$location);
  my @impact=_normalize3([map {$_+0.0} @field[6..8]],
      'impact-parameter direction',$location);
  my $dot=$los[0]*$impact[0]+$los[1]*$impact[1]+$los[2]*$impact[2];
  die "BrightnessObservationGeometry LOS and impact directions must be orthogonal at $location\n"
      if abs($dot)>1.0e-12;
  $self->{observation_geometry}={radial=>\@radial,los=>\@los,
      impact=>\@impact};
}

sub set_velocity_domain_policy_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessVelocityDomainPolicy at $location\n"
      unless $record=~/^\s*BrightnessVelocityDomainPolicy\s*[=:]\s*(.*?)\s*$/i;
  my $policy=uc(_trim($1));
  die "BrightnessVelocityDomainPolicy must be ABORT or EXCLUDE_AND_REPORT at $location\n"
      unless $policy eq 'ABORT' || $policy eq 'EXCLUDE_AND_REPORT';
  die "duplicate BrightnessVelocityDomainPolicy at $location\n"
      if defined $self->{velocity_domain_policy};
  $self->{velocity_domain_policy}=$policy;
}

sub set_radial_velocity_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessHeliocentricRadialVelocityKMS at $location\n"
      unless $record=~/^\s*BrightnessHeliocentricRadialVelocityKMS\s*[=:]\s*(.*?)\s*$/i;
  my $value=_trim($1);
  die "BrightnessHeliocentricRadialVelocityKMS must be a finite number at $location\n"
      unless _number($value);
  die "duplicate BrightnessHeliocentricRadialVelocityKMS at $location\n"
      if defined($self->{radial_velocity_km_s}) ||
          defined($self->{comet_velocity_km_s});
  $self->{radial_velocity_km_s}=$value+0.0;
}

sub set_solar_lya_flux_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessSolarLymanAlphaFlux1AUPhotonCM2S at $location\n"
      unless $record=~/^\s*BrightnessSolarLymanAlphaFlux1AUPhotonCM2S\s*[=:]\s*(.*?)\s*$/i;
  my $value=_trim($1);
  die "BrightnessSolarLymanAlphaFlux1AUPhotonCM2S must be positive at $location\n"
      unless _positive_number($value);
  die "duplicate BrightnessSolarLymanAlphaFlux1AUPhotonCM2S at $location\n"
      if defined $self->{solar_lya_flux_1au_ph_cm2_s};
  $self->{solar_lya_flux_1au_ph_cm2_s}=$value+0.0;
}

sub set_mode_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessMode at $location\n"
      unless $record=~/^\s*BrightnessMode\s*[=:]\s*(.*?)\s*$/i;
  my $mode=uc(_trim($1));
  die "BrightnessMode must be ON or OFF at $location\n"
      unless $mode eq 'ON' || $mode eq 'OFF';
  die "duplicate BrightnessMode at $location\n"
      if defined $self->{mode};
  $self->{mode}=$mode;
}

sub set_grid_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessImpactParameterGrid at $location\n"
      unless $record=~/^\s*BrightnessImpactParameterGrid\s*[=:]\s*(.*?)\s*$/i;
  my @field=map {_trim($_)} split(/\|/,$1,-1);
  die "BrightnessImpactParameterGrid expects LINEAR|LOG | MIN_M | MAX_M | COUNT at $location\n"
      unless @field==4;
  my $mode=uc($field[0]);
  die "BrightnessImpactParameterGrid mode must be LINEAR or LOG at $location\n"
      unless $mode eq 'LINEAR' || $mode eq 'LOG';
  die "BrightnessImpactParameterGrid limits must be non-negative numbers at $location\n"
      unless _nonnegative_number($field[1]) &&
          _nonnegative_number($field[2]);
  die "BrightnessImpactParameterGrid count must be a positive integer at $location\n"
      unless $field[3]=~/^\d+$/ && $field[3]>0;
  my ($minimum,$maximum)=($field[1]+0.0,$field[2]+0.0);
  die "BrightnessImpactParameterGrid maximum must exceed minimum for multiple points at $location\n"
      if $field[3]>1 && !($maximum>$minimum);
  die "BrightnessImpactParameterGrid single-point limits must be equal at $location\n"
      if $field[3]==1 && $maximum!=$minimum;
  die "BrightnessImpactParameterGrid LOG requires a positive minimum at $location\n"
      if $mode eq 'LOG' && !($minimum>0.0);
  die "duplicate BrightnessImpactParameterGrid at $location\n"
      if defined $self->{grid};
  $self->{grid}={mode=>$mode,minimum=>$minimum,maximum=>$maximum,
      count=>$field[3]+0};
}

sub add_emission_record {
  my ($self,$line,$location)=@_;
  $location='' unless defined $location;
  my $record=$line;
  $record=~s/!.*$//;
  die "invalid BrightnessEmission at $location\n"
      unless $record=~/^\s*BrightnessEmission\s*[=:]\s*(.*?)\s*$/i;
  my @field=map {_trim($_)} split(/\|/,$1,-1);
  die "BrightnessEmission expects SPECIES | TRANSITION | PHOTON_G_FACTOR_S-1_OR_MODEL | G_PROVENANCE | THIN | ISOTROPIC_4PI | NONE|LINE_CENTER | CROSS_SECTION_M2 | TAU_THRESHOLD | TAU_PROVENANCE [| EFFECTIVE_WAVELENGTH_NM] at $location\n"
      unless @field==10 || @field==11;
  my ($species,$transition,$g_factor_or_model,$provenance,$emission_mode,
      $phase_convention,$tau_mode,$cross_section,$threshold,
      $tau_provenance,$effective_wavelength_nm)=@field;
  $effective_wavelength_nm=0.0 unless defined $effective_wavelength_nm;
  $species=uc($species);
  die "unknown brightness species '$species' at $location\n"
      unless exists $self->{species_index}->{$species};
  _validate_metadata($transition,'transition',$location);
  _validate_metadata($provenance,'g-factor provenance',$location);
  _validate_metadata($tau_provenance,'tau provenance',$location);
  my %model=(
    "H_LYA_IKL2020"=>1,
    "OH_AX00_SA1988_LOWELL_UNQUENCHED"=>2,
  );
  my $model_name="FIXED";
  my $model_code=0;
  my $g_factor=0.0;
  if (_positive_number($g_factor_or_model)) {
    $g_factor=$g_factor_or_model+0.0;
    die "a fixed photon g factor cannot specify an effective wavelength at $location\n"
        unless _nonnegative_number($effective_wavelength_nm) &&
            $effective_wavelength_nm==0.0;
  }
  else {
    $model_name=uc($g_factor_or_model);
    die "unknown BrightnessEmission g-factor model '$g_factor_or_model' at $location\n"
        unless exists $model{$model_name};
    $model_code=$model{$model_name};
    if ($model_code==1) {
      die "H_LYA_IKL2020 can only be assigned to species H at $location\n"
          unless $species eq 'H';
      die "H_LYA_IKL2020 does not accept an effective wavelength at $location\n"
          unless _nonnegative_number($effective_wavelength_nm) &&
              $effective_wavelength_nm==0.0;
    }
    else {
      die "OH_AX00_SA1988_LOWELL_UNQUENCHED can only be assigned to species OH at $location\n"
          unless $species eq 'OH';
      die "OH energy g-factor conversion requires a positive explicit effective wavelength in nm at $location\n"
          unless _positive_number($effective_wavelength_nm);
    }
  }
  $emission_mode=uc($emission_mode);
  die "BrightnessEmission mode must be THIN at $location; attenuation is not implemented at $location\n"
      unless $emission_mode eq 'THIN';
  $phase_convention=uc($phase_convention);
  die "BrightnessEmission phase convention must be ISOTROPIC_4PI at $location\n"
      unless $phase_convention eq 'ISOTROPIC_4PI';
  $tau_mode=uc($tau_mode);
  die "BrightnessEmission tau diagnostic must be NONE or LINE_CENTER at $location\n"
      unless $tau_mode eq 'NONE' || $tau_mode eq 'LINE_CENTER';
  die "BrightnessEmission cross section and threshold must be non-negative numbers at $location\n"
      unless _nonnegative_number($cross_section) &&
          _nonnegative_number($threshold);
  if ($tau_mode eq 'NONE') {
    die "BrightnessEmission NONE requires zero cross section, zero threshold, and tau provenance NONE at $location\n"
        unless $cross_section==0.0 && $threshold==0.0 &&
            uc($tau_provenance) eq 'NONE';
  }
  else {
    die "BrightnessEmission LINE_CENTER requires positive cross section, positive threshold, and non-NONE tau provenance at $location\n"
        unless $cross_section>0.0 && $threshold>0.0 &&
            uc($tau_provenance) ne 'NONE';
  }
  # Any atomic-H fluorescence line is potentially resonance-thick in a large
  # comet domain.  Requiring the conservative diagnostic is safer than using
  # a transition-name spelling heuristic and still permits arbitrary labels.
  die "BrightnessEmission for H requires LINE_CENTER tau diagnostics at $location\n"
      if $species eq 'H' && $tau_mode ne 'LINE_CENTER';
  my $key=$species."\0".uc($transition);
  die "duplicate BrightnessEmission for $species | $transition at $location\n"
      if exists $self->{keys}->{$key};
  die "BrightnessEmission exceeds the maximum of ".MAX_EMISSION_LINES.
      " configured lines at $location\n"
      if @{$self->{emissions}}>=MAX_EMISSION_LINES;
  $self->{keys}->{$key}=1;
  push @{$self->{emissions}}, {
    species=>$self->{species_index}->{$species},species_name=>$species,
    transition=>$transition,g_factor=>$g_factor,model=>$model_code,
    model_name=>$model_name,
    effective_wavelength_m=>($effective_wavelength_nm+0.0)*1.0e-9,
    provenance=>$provenance,tau_mode=>$tau_mode,
    cross_section=>$cross_section+0.0,threshold=>$threshold+0.0,
    tau_provenance=>$tau_provenance,
  };
}

sub validate_configuration {
  my ($self,%args)=@_;
  my $enabled=defined($self->{mode}) && $self->{mode} eq 'ON';
  if (!$enabled) {
    die "brightness records require BrightnessMode=ON\n"
        if defined($self->{grid}) || @{$self->{emissions}} ||
            defined($self->{radial_velocity_km_s}) ||
            defined($self->{comet_velocity_km_s}) ||
            defined($self->{observation_geometry}) ||
            defined($self->{velocity_domain_policy}) ||
            defined($self->{solar_lya_flux_1au_ph_cm2_s});
    return 1;
  }
  die "BrightnessMode=ON requires Dimension=1\n"
      unless defined($args{dimension}) && $args{dimension}==1;
  die "BrightnessMode=ON requires RadialSimulationMode=ON\n"
      unless $args{radial_enabled};
  die "BrightnessMode=ON requires BrightnessImpactParameterGrid\n"
      unless defined $self->{grid};
  die "BrightnessMode=ON requires at least one BrightnessEmission\n"
      unless @{$self->{emissions}};
  die "BrightnessMode=ON requires BrightnessObservationGeometry\n"
      unless defined $self->{observation_geometry};
  $self->{velocity_domain_policy}='ABORT'
      unless defined $self->{velocity_domain_policy};
  my $uses_core=grep {$_->{model}!=0} @{$self->{emissions}};
  my $uses_h=grep {$_->{model}==1} @{$self->{emissions}};
  die "a named brightness g-factor model requires BrightnessHeliocentricRadialVelocityKMS\n"
      if $uses_core && !defined($self->{radial_velocity_km_s}) &&
          !defined($self->{comet_velocity_km_s});
  die "H_LYA_IKL2020 requires BrightnessSolarLymanAlphaFlux1AUPhotonCM2S\n"
      if $uses_h && !defined($self->{solar_lya_flux_1au_ph_cm2_s});
  return 1;
}

sub enabled {
  my ($self)=@_;
  return defined($self->{mode}) && $self->{mode} eq 'ON';
}
sub grid { return $_[0]->{grid}; }
sub emissions { return $_[0]->{emissions}; }
sub radial_velocity_km_s { return $_[0]->{radial_velocity_km_s}; }
sub observation_geometry { return $_[0]->{observation_geometry}; }
sub velocity_domain_policy { return $_[0]->{velocity_domain_policy}; }
sub resolved_comet_velocity_km_s {
  my ($self)=@_;
  return $self->{comet_velocity_km_s}
      if defined $self->{comet_velocity_km_s};
  return [0.0,0.0,0.0] unless defined $self->{radial_velocity_km_s};
  die "cannot resolve radial comet velocity without observation geometry\n"
      unless defined $self->{observation_geometry};
  return [map {$self->{radial_velocity_km_s}*$_}
      @{$self->{observation_geometry}->{radial}}];
}
sub solar_lya_flux_1au_ph_cm2_s {
  return $_[0]->{solar_lya_flux_1au_ph_cm2_s};
}

sub generated_metadata_calls {
  my ($self)=@_;
  my $result='';
  for (my $line=0;$line<@{$self->{emissions}};$line++) {
    my $entry=$self->{emissions}->[$line];
    $result.='  SetGeneratedMetadata('.$line.',"'.
        _cpp_escape($entry->{transition}).'","'.
        _cpp_escape($entry->{provenance}).'","'.
        _cpp_escape($entry->{tau_provenance}).'");' . "\n";
  }
  return $result;
}

sub _cpp_escape {
  my ($value)=@_;
  $value=~s/\\/\\\\/g;
  $value=~s/"/\\"/g;
  return $value;
}
sub _validate_metadata {
  my ($value,$description,$location)=@_;
  die "BrightnessEmission $description cannot be empty at $location\n"
      if !defined($value) || $value eq '';
  die "BrightnessEmission $description exceeds 159 characters at $location\n"
      if length($value)>=160;
  die "BrightnessEmission $description contains a control character at $location\n"
      if $value=~/[\x00-\x1f\x7f]/;
}
sub _number {
  my ($value)=@_;
  return defined($value) &&
      $value=~/^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?$/;
}
sub _vector3 {
  my ($value,$description,$location)=@_;
  my @field=map {_trim($_)} split(/\|/,$value,-1);
  die "$description expects X | Y | Z at $location\n"
      unless @field==3 && !grep {!_number($_)} @field;
  return [map {$_+0.0} @field];
}
sub _normalize3 {
  my ($value,$description,$location)=@_;
  my $norm=sqrt($value->[0]*$value->[0]+$value->[1]*$value->[1]+
      $value->[2]*$value->[2]);
  die "$description must be nonzero at $location\n" unless $norm>0.0;
  return map {$_/$norm} @$value;
}
sub _nonnegative_number { return _number($_[0]) && $_[0]>=0.0; }
sub _positive_number { return _number($_[0]) && $_[0]>0.0; }
sub _trim {
  my ($value)=@_;
  $value='' unless defined $value;
  $value=~s/^\s+//;
  $value=~s/\s+$//;
  return $value;
}

1;
