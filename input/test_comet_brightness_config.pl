#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use CometBrightnessConfig;

sub expect_failure {
  my ($name,$pattern,$callback)=@_;
  my $ok=eval {$callback->(); 1};
  die "FAIL $name was accepted\n" if $ok;
  die "FAIL $name emitted unexpected error: $@\n" unless $@=~/$pattern/;
}

my @species=qw(H OH H2O CO O);
my $valid=CometBrightnessConfig->new(species=>\@species);
$valid->set_mode_record('BrightnessMode=ON','valid');
$valid->set_grid_record(
    'BrightnessImpactParameterGrid=LOG | 2500 | 1.0e9 | 41','valid');
$valid->set_observation_geometry_record(
    'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','valid');
$valid->add_emission_record(
    'BrightnessEmission=H | H_LYMAN_ALPHA | 2.5e-3 | TEST_H_SOURCE | THIN | ISOTROPIC_4PI | LINE_CENTER | 1.0e-17 | 0.1 | TEST_H_SIGMA',
    'valid');
$valid->add_emission_record(
    'BrightnessEmission=OH | OH_A2SIGMA_X2PI_00 | 2.0e-4 | TEST_OH_SOURCE | THIN | ISOTROPIC_4PI | NONE | 0 | 0 | NONE',
    'valid');
$valid->validate_configuration(dimension=>1,radial_enabled=>1);
die "FAIL valid brightness record count\n"
    unless @{$valid->emissions()}==2 && $valid->enabled();
my $code=$valid->generated_metadata_calls();
die "FAIL generated metadata escaping/content\n"
    unless $code=~/H_LYMAN_ALPHA/ && $code=~/TEST_OH_SOURCE/;

my $modeled=CometBrightnessConfig->new(species=>\@species);
$modeled->set_mode_record('BrightnessMode=ON','modeled');
$modeled->set_grid_record(
    'BrightnessImpactParameterGrid=LOG | 2500 | 1.0e9 | 41','modeled');
$modeled->set_observation_geometry_record(
    'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','modeled');
$modeled->set_radial_velocity_record(
    'BrightnessHeliocentricRadialVelocityKMS=-39.7','modeled');
$modeled->set_solar_lya_flux_record(
    'BrightnessSolarLymanAlphaFlux1AUPhotonCM2S=5.1040518742e11','modeled');
$modeled->add_emission_record(
    'BrightnessEmission=H | H_LYMAN_ALPHA | H_LYA_IKL2020 | CORE_IKL2020 | THIN | ISOTROPIC_4PI | LINE_CENTER | 1e-17 | 0.1 | TEST_H_SIGMA | 0',
    'modeled');
$modeled->add_emission_record(
    'BrightnessEmission=OH | OH_A2SIGMA_X2PI_00 | OH_AX00_SA1988_LOWELL_UNQUENCHED | CORE_LOWELL | THIN | ISOTROPIC_4PI | NONE | 0 | 0 | NONE | 308.5',
    'modeled');
$modeled->validate_configuration(dimension=>1,radial_enabled=>1);
die "FAIL named model representation\n"
    unless $modeled->emissions()->[0]->{model}==1 &&
        $modeled->emissions()->[1]->{model}==2 &&
        abs($modeled->emissions()->[1]->{effective_wavelength_m}-308.5e-9)<1e-20;

my $oxygen=CometBrightnessConfig->new(species=>\@species);
$oxygen->set_mode_record('BrightnessMode=ON','oxygen');
$oxygen->set_grid_record(
    'BrightnessImpactParameterGrid=LOG | 2500 | 1.0e9 | 41','oxygen');
$oxygen->set_observation_geometry_record(
    'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','oxygen');
$oxygen->set_radial_velocity_record(
    'BrightnessHeliocentricRadialVelocityKMS=0.4','oxygen');
$oxygen->set_solar_oi1304_flux_record(
    'BrightnessSolarOxygen1304TripletFlux1AUPhotonCM2S=1.0657e10','oxygen');
$oxygen->set_solar_lyb_flux_record(
    'BrightnessSolarLymanBetaFlux1AUPhotonCM2S=4.961e9','oxygen');
my @oxygen_models=qw(O_I_1304_KVB2022_DIRECT_J2 O_I_1304_SUMER2001_DIRECT_J2
    O_I_1304_SUMER2001_LYB_BOWEN_J2);
for (my $i=0;$i<@oxygen_models;$i++) {
  $oxygen->add_emission_record(
      "BrightnessEmission=O | OI_1304_COMPONENT_$i | $oxygen_models[$i] | CORE_O | THIN | ISOTROPIC_4PI | NONE | 0 | 0 | NONE | 0",
      'oxygen');
}
$oxygen->validate_configuration(dimension=>1,radial_enabled=>1);
die "FAIL O I model representation\n"
    unless join(',',map {$_->{model}} @{$oxygen->emissions()}) eq '3,4,5' &&
        abs($oxygen->solar_oi1304_flux_1au_ph_cm2_s()-1.0657e10)<1 &&
        abs($oxygen->solar_lyb_flux_1au_ph_cm2_s()-4.961e9)<1;

my $omitted=CometBrightnessConfig->new(species=>\@species);
$omitted->validate_configuration(dimension=>1,radial_enabled=>1);
die "FAIL omitted brightness is not disabled\n" if $omitted->enabled();

expect_failure('H without tau',qr/H requires LINE_CENTER/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=H|LYA|1e-3|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
});
expect_failure('fake attenuation mode',qr/mode must be THIN/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|1e-4|SRC|ESCAPE_FACTOR|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
});
expect_failure('unknown species',qr/unknown brightness species/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=C2|SWAN|1e-4|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
});
expect_failure('records while off',qr/require BrightnessMode=ON/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=OFF','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->validate_configuration(dimension=>1,radial_enabled=>1);
});
expect_failure('DIM3 enable',qr/requires Dimension=1/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=ON','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|1e-4|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
  $c->validate_configuration(dimension=>3,radial_enabled=>1);
});
expect_failure('missing radial backend',qr/requires RadialSimulationMode=ON/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=ON','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|1e-4|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
  $c->validate_configuration(dimension=>1,radial_enabled=>0);
});

expect_failure('unsupported phase convention',qr/phase convention must be ISOTROPIC_4PI/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|1e-4|SRC|THIN|DIPOLE|NONE|0|0|NONE','bad');
});
expect_failure('missing tau provenance',qr/non-NONE tau provenance/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=H|LYA|1e-3|SRC|THIN|ISOTROPIC_4PI|LINE_CENTER|1e-17|0.1|NONE','bad');
});
expect_failure('named model missing velocity',qr/requires BrightnessHeliocentricRadialVelocityKMS/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=ON','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->set_observation_geometry_record(
      'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','bad');
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|OH_AX00_SA1988_LOWELL_UNQUENCHED|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE|308.5','bad');
  $c->validate_configuration(dimension=>1,radial_enabled=>1);
});
expect_failure('H model missing irradiance',qr/requires BrightnessSolarLymanAlphaFlux1AUPhotonCM2S/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=ON','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->set_observation_geometry_record(
      'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','bad');
  $c->set_radial_velocity_record('BrightnessHeliocentricRadialVelocityKMS=1','bad');
  $c->add_emission_record(
      'BrightnessEmission=H|LYA|H_LYA_IKL2020|SRC|THIN|ISOTROPIC_4PI|LINE_CENTER|1e-17|0.1|SIGMA|0','bad');
  $c->validate_configuration(dimension=>1,radial_enabled=>1);
});
expect_failure('OH model missing wavelength',qr/requires a positive explicit effective wavelength/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=OH|308NM|OH_AX00_SA1988_LOWELL_UNQUENCHED|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
});

expect_failure('O model on another species',qr/can only be assigned to species O/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=OH|OI|O_I_1304_KVB2022_DIRECT_J2|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
});
expect_failure('O model with wavelength',qr/does not accept an effective wavelength/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->add_emission_record(
      'BrightnessEmission=O|OI|O_I_1304_SUMER2001_DIRECT_J2|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE|130.4','bad');
});
for my $case (['direct',qr/require BrightnessSolarOxygen1304TripletFlux1AUPhotonCM2S/,
    'O_I_1304_SUMER2001_DIRECT_J2'],
    ['Bowen',qr/requires BrightnessSolarLymanBetaFlux1AUPhotonCM2S/,
    'O_I_1304_SUMER2001_LYB_BOWEN_J2']) {
  expect_failure("O $case->[0] model missing irradiance",$case->[1],sub {
    my $c=CometBrightnessConfig->new(species=>\@species);
    $c->set_mode_record('BrightnessMode=ON','bad');
    $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
    $c->set_observation_geometry_record(
        'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','bad');
    $c->set_radial_velocity_record('BrightnessHeliocentricRadialVelocityKMS=1','bad');
    $c->add_emission_record(
        "BrightnessEmission=O|OI|$case->[2]|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE",'bad');
    $c->validate_configuration(dimension=>1,radial_enabled=>1);
  });
}
expect_failure('unused Lyman-beta irradiance',qr/requires O_I_1304_SUMER2001_LYB_BOWEN_J2/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_mode_record('BrightnessMode=ON','bad');
  $c->set_grid_record('BrightnessImpactParameterGrid=LINEAR|0|10|2','bad');
  $c->set_observation_geometry_record(
      'BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0','bad');
  $c->set_solar_lyb_flux_record('BrightnessSolarLymanBetaFlux1AUPhotonCM2S=4e9','bad');
  $c->add_emission_record(
      'BrightnessEmission=OH|OH|2e-4|SRC|THIN|ISOTROPIC_4PI|NONE|0|0|NONE','bad');
  $c->validate_configuration(dimension=>1,radial_enabled=>1);
});

expect_failure('nonorthogonal observation geometry',qr/must be orthogonal/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_observation_geometry_record(
      'BrightnessObservationGeometry=1|0|0|1|0|0|1|0|0','bad');
});
expect_failure('invalid velocity policy',qr/must be ABORT or EXCLUDE_AND_REPORT/,sub {
  my $c=CometBrightnessConfig->new(species=>\@species);
  $c->set_velocity_domain_policy_record(
      'BrightnessVelocityDomainPolicy=CLAMP','bad');
});

print "PASS comet-brightness input policy fixed=1 modeled=1 oxygen=1 default-off=1 invalid=18\n";
