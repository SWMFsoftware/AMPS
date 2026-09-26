#include "../../util/SWMFCoupledProductsContract.h"
#include "../../util/BoundaryProducts.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace CP=Earth::SWMFCoupledProducts;
namespace BP=Earth::BoundaryProducts;
namespace FN=Earth::FluxNumerics;

namespace {

void Require(bool condition,const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void RequireClose(double actual,double expected,double relative,
                  const std::string& message) {
  const double tolerance=relative*std::max(1.0,std::max(std::fabs(actual),
                                                        std::fabs(expected)));
  if (!std::isfinite(actual) || std::fabs(actual-expected)>tolerance)
    throw std::runtime_error(message+": actual="+std::to_string(actual)+
                             ", expected="+std::to_string(expected));
}

template<class Function>
void RequireThrows(Function function,const std::string& message) {
  bool rejected=false;
  try { function(); }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,message);
}

CP::ProductControl ReferenceControl() {
  CP::ProductControl control;
  control.outputMode="POINTS";
  control.speciesName="proton";
  control.charge_e=1.0;
  control.mass_amu=1.0;
  control.boundaryMode="ISOTROPIC";
  control.transmissionMode="DIRECT";
  control.minimumEnergy_MeV=10.0;
  control.maximumEnergy_MeV=30.0;
  control.energyIntervals=2;
  control.transmissionScanPoints=0;
  control.maximumParticlesPerPoint=3456;
  control.energySpacing="LINEAR";
  control.spectrumType="POWER_LAW";
  control.energyBasis="PER_PARTICLE";
  control.spectrumMassNumber=1.0;
  control.intensityUnit="m^-2 s^-1 sr^-1 MeV^-1";
  control.spectrumRelativeUncertainty=0.0;
  // Deliberately insert spectrum keys out of lexical order.  The production
  // fingerprint sorts this map-shaped data, making parser insertion order irrelevant.
  control.spectrumKeyValues.push_back(std::make_pair("SPEC_GAMMA","2"));
  control.spectrumKeyValues.push_back(std::make_pair("SPECTRUM_TYPE","POWER_LAW"));
  CP::EnergyChannelDefinition channel;
  channel.name="P15_25";
  channel.lower_MeV=15.0;
  channel.upper_MeV=25.0;
  control.channels.push_back(channel);
  CP::DetectorResponseDefinition response;
  response.name="TOPHAT";
  response.lower_MeV=15.0;
  response.upper_MeV=25.0;
  response.geometricFactor_m2_sr=0.01;
  control.detectorResponses.push_back(response);
  control.coordinateFrame="GSM";
  CP::ObservationDefinition observation;
  observation.x_km=7000.0;
  observation.y_km=0.0;
  observation.z_km=100.0;
  control.observations.push_back(observation);
  control.shellResolution_deg=15.0; // unused for POINTS but canonically retained
  control.shellGeometry="SPHERICAL";
  CP::ValidateProductControl(control);
  return control;
}

CP::ProductRunSummary ReferenceSummary(const std::string& suffix) {
  CP::ProductRunSummary summary;
  summary.control=ReferenceControl();
  summary.spectrumEvaluationEpochUTC="2024-05-10T12:00:00.000000000Z";
  summary.activeSpectrumTableEpochUTC=summary.spectrumEvaluationEpochUTC;
  summary.spectrumTemporalStatus="EXACT";
  summary.spectrumTemporalFraction=0.0;
  summary.locationCount=1;
  summary.energyCount=3;
  summary.directionCount=1152;
  summary.sampled=100;
  summary.retried=1;
  summary.resolved=99;
  summary.allowed=80;
  summary.terminationCounts={80,19,1};
  summary.maximumUnresolvedFraction=0.01;
  summary.unresolvedTolerance=0.01;
  summary.artifacts={
      "mode3d_points_density"+suffix+".dat",
      "mode3d_points_spectrum"+suffix+".dat",
      "mode3d_points_flux"+suffix+".dat",
      "mode3d_termination_summary"+suffix+".dat"};
  summary.valid=true;
  return summary;
}

void TestControlIdentity() {
  CP::ProductControl control=ReferenceControl();
  const std::string fingerprint=CP::ProductControlFingerprint(control);
  Require(fingerprint=="products-v1-da0e0ed31f85dabe",
          "fixed Step-11 product-control fingerprint changed: "+fingerprint);
  Require(fingerprint==CP::ProductControlFingerprint(control),
          "identical product control must reproduce its fingerprint");

  CP::ProductControl reordered=control;
  std::reverse(reordered.spectrumKeyValues.begin(),
               reordered.spectrumKeyValues.end());
  Require(CP::ProductControlFingerprint(reordered)==fingerprint,
          "raw spectrum map order must not change identity");

  CP::ProductControl changed=control;
  changed.detectorResponses[0].geometricFactor_m2_sr=0.02;
  Require(CP::ProductControlFingerprint(changed)!=fingerprint,
          "detector response must participate in product identity");
  changed=control;
  changed.observations[0].x_km+=1.0;
  Require(CP::ProductControlFingerprint(changed)!=fingerprint,
          "spacecraft state must participate in product identity");
  changed=control;
  changed.channels[0].upper_MeV=26.0;
  Require(CP::ProductControlFingerprint(changed)!=fingerprint,
          "flux channel must participate in product identity");
  CP::ProductControl table=control;
  table.spectrumType="TABLE";
  table.spectrumTableEnergy_MeV={10.0,20.0};
  table.spectrumTableIntensityPerMeV={5.0,2.5};
  const std::string tableFingerprint=CP::BoundarySpectrumFingerprint(table);
  table.spectrumTableIntensityPerMeV[1]=2.500001;
  Require(CP::BoundarySpectrumFingerprint(table)!=tableFingerprint,
          "selected spectrum-table content must participate in identity");
  std::cout << "PASS S11-U01 synchronized physics/observation identity\n";
}

void TestManifestAndNegativeGates() {
  const std::string snapshot="field-v1-reference";
  const double simulationTime_s=3600.0;
  const std::string suffix=Earth::SWMFCoupledAccess::BuildProductSuffix(
      simulationTime_s,snapshot);
  const CP::ProductRunSummary summary=ReferenceSummary(suffix);
  const std::string manifest=CP::BuildProductsManifestJson(
      "PASS",snapshot,"swmf-state-v1-reference","swmf-mesh-v1-reference",
      "2024-05-10T12:00:00.000000000Z",simulationTime_s,"BOX",suffix,summary,
      "complete");
  Require(manifest.find("sep-in-geospace/swmf-coupled-products/v1")!=
              std::string::npos,
          "Step-11 schema missing");
  Require(manifest.find("INSTANTANEOUS_QUASI_STATIC")!=std::string::npos,
          "Phase-1 limitation missing");
  Require(manifest.find("mode3d_points_spectrum")!=std::string::npos &&
          manifest.find("TOPHAT")!=std::string::npos,
          "manifest lost spectrum or detector schema");

  CP::ProductRunSummary invalid=summary;
  invalid.maximumUnresolvedFraction=0.0100001;
  RequireThrows([&]() {
      CP::BuildProductsManifestJson(
          "PASS",snapshot,"state","mesh",summary.spectrumEvaluationEpochUTC,
          simulationTime_s,"BOX",suffix,invalid,"invalid");
    },"unresolved gate must not be relaxed");
  invalid=summary;
  invalid.terminationCounts[0]-=1;
  RequireThrows([&]() { CP::ValidateRunSummary(invalid,true); },
                "termination counts must close");
  invalid=summary;
  invalid.artifacts.pop_back();
  RequireThrows([&]() { CP::ValidateRunSummary(invalid,true); },
                "termination artifact is mandatory");
  RequireThrows([&]() {
      CP::BuildProductsManifestJson(
          "PASS",snapshot,"state","mesh","2024-05-10T12:00:01Z",
          simulationTime_s,"BOX",suffix,summary,"invalid");
    },"field and boundary-spectrum epoch mismatch must fail");
  std::cout << "PASS S11-U02 fail-closed manifest and retained unresolved gate\n";
}

void TestOpenAndBlockedReferenceProducts() {
  const std::vector<double> energy={10.0,20.0,30.0};
  const std::vector<double> open(energy.size(),1.0);
  const std::vector<double> blocked(energy.size(),0.0);
  const double mass=FN::kAtomicMassUnit_kg;
  const auto constantSpectrum=[](double) {
    // Exactly 2 particles m^-2 s^-1 sr^-1 MeV^-1, converted to the per-joule
    // callable required by BoundaryProducts.
    return 2.0/FN::kMeVToJ;
  };
  std::vector<BP::EnergyChannel> channels={
      BP::EnergyChannel("P15_25",15.0,25.0)};
  BP::DetectorResponse response;
  response.name="TOPHAT";
  response.energy_MeV={15.0,25.0};
  response.relativeResponse={1.0,1.0};
  response.geometricFactor_m2_sr=0.01;

  const BP::ProductSet product=BP::EvaluateIsotropicProducts(
      energy,open,open,open,mass,constantSpectrum,channels,{response});
  RequireClose(product.spectrum[1].localPerMeV.nominal,2.0,2.0e-14,
               "open local differential spectrum");
  RequireClose(product.omnidirectionalFlux_m2_s.nominal,160.0*FN::kPi,
               3.0e-14,"open omnidirectional flux analytic reference");
  RequireClose(product.oneWayPlanarFlux_m2_s.nominal,40.0*FN::kPi,
               3.0e-14,"open planar flux analytic reference");
  RequireClose(product.channelFlux_m2_s[0].nominal,80.0*FN::kPi,
               3.0e-14,"clipped energy-channel analytic reference");
  RequireClose(product.detectorRate_s[0].nominal,0.2,3.0e-14,
               "top-hat detector analytic reference");

  const BP::ProductSet zero=BP::EvaluateIsotropicProducts(
      energy,blocked,blocked,blocked,mass,constantSpectrum,channels,{response});
  Require(zero.numberDensity_m3.nominal==0.0 &&
          zero.omnidirectionalFlux_m2_s.nominal==0.0 &&
          zero.channelFlux_m2_s[0].nominal==0.0 &&
          zero.detectorRate_s[0].nominal==0.0,
          "blocked reference must produce physical zero products");
  std::cout << "PASS S11-U03 open/blocked analytic flux-spectrum references\n";
}

void TestUnresolvedBoundsReference() {
  const std::vector<double> energy={10.0,20.0,30.0};
  const std::vector<double> nominal(3,0.5);
  const std::vector<double> lower(3,0.49);
  const std::vector<double> upper(3,0.51);
  const auto constantSpectrum=[](double) { return 1.0/FN::kMeVToJ; };
  const BP::ProductSet product=BP::EvaluateIsotropicProducts(
      energy,nominal,lower,upper,FN::kAtomicMassUnit_kg,constantSpectrum);
  RequireClose(product.omnidirectionalFlux_m2_s.nominal,40.0*FN::kPi,
               3.0e-14,"bounded nominal flux");
  RequireClose(product.omnidirectionalFlux_m2_s.lower,39.2*FN::kPi,
               3.0e-14,"bounded lower flux");
  RequireClose(product.omnidirectionalFlux_m2_s.upper,40.8*FN::kPi,
               3.0e-14,"bounded upper flux");
  Require(product.omnidirectionalFlux_m2_s.lower<
              product.omnidirectionalFlux_m2_s.nominal &&
          product.omnidirectionalFlux_m2_s.nominal<
              product.omnidirectionalFlux_m2_s.upper,
          "unresolved access must remain an interval, not be replaced by zero");
  std::cout << "PASS S11-U04 conservative unresolved-product bounds\n";
}

} // namespace

int main() {
  try {
    TestControlIdentity();
    TestManifestAndNegativeGates();
    TestOpenAndBlockedReferenceProducts();
    TestUnresolvedBoundsReference();
    std::cout << "RESULT: PASS\n";
    return 0;
  }
  catch (const std::exception& error) {
    std::cerr << "RESULT: FAIL: " << error.what() << "\n";
    return 1;
  }
}
