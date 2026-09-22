#include "../../util/BoundaryProducts.h"
#include "../../boundary/spectrum.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace BP=Earth::BoundaryProducts;
namespace FN=Earth::FluxNumerics;

static int failures=0;

static bool Near(double actual,double expected,double relative=1.0e-10,
                 double absolute=1.0e-12) {
  return std::fabs(actual-expected)<=absolute+
      relative*std::max(std::fabs(actual),std::fabs(expected));
}

static void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

template<class Callable>
static void CheckThrows(Callable action,const std::string& message) {
  try { action(); }
  catch (...) { return; }
  Check(false,message);
}

int main() {
  // U-F01 -- unit metadata and all five production spectrum families.  These are
  // comparisons with closed-form definitions, not self-consistency-only checks.
  BP::SpectrumUnits ionUnits;
  ionUnits.energyBasis=BP::EnergyBasis::PerNucleon;
  ionUnits.massNumber=4.0;
  const double inputPerMeVn=37.25;
  Check(Near(ionUnits.PerParticleJToPerCoordinateMeV(
                 ionUnits.PerCoordinateMeVToPerParticleJ(inputPerMeVn)),inputPerMeVn,
             1.0e-13),
        "U-F01 per-MeV/nucleon and per-particle-J round trip");
  Check(Near(ionUnits.ParticleEnergyJ(10.0),40.0*FN::kMeVToJ,1.0e-14),
        "U-F01 per-nucleon coordinate maps to total particle energy");

  const cSpectrum power=cSpectrum::MakePowerLaw(100.0,2.0,10.0,1.0,1000.0);
  Check(Near(power.GetSpectrumPerMeV(20.0),25.0,2.0e-14),
        "U-F01 POWER_LAW closed form");
  const cSpectrum cutoff=cSpectrum::MakePowerLawCutoff(100.0,2.0,10.0,40.0,1.0,1000.0);
  Check(Near(cutoff.GetSpectrumPerMeV(20.0),25.0*std::exp(-0.5),2.0e-14),
        "U-F01 POWER_LAW_CUTOFF closed form");
  const cSpectrum lis=cSpectrum::MakeLisForceField(100.0,2.0,10.0,50.0,1.0,1000.0);
  const double lisEnergy=20.0, shifted=70.0, mass=938.2720813;
  const double lisReference=100.0*std::pow(shifted/10.0,-2.0)*
      (lisEnergy*(lisEnergy+2.0*mass))/(shifted*(shifted+2.0*mass));
  Check(Near(lis.GetSpectrumPerMeV(lisEnergy),lisReference,2.0e-14),
        "U-F01 LIS_FORCE_FIELD closed form");
  const cSpectrum band=cSpectrum::MakeBand(5.0,1.0,3.0,10.0,1.0,1000.0);
  Check(Near(band.GetSpectrumPerMeV(10.0),5.0*std::exp(-1.0),2.0e-14),
        "U-F01 BAND low branch closed form");

  const std::string tableFile="/tmp/amps_u_boundary_products_table.dat";
  {
    std::ofstream out(tableFile.c_str());
    out << "1 100\n10 10\n100 1\n";
  }
  const cSpectrum table=cSpectrum::MakeTable(tableFile,1.0,100.0);
  Check(Near(table.GetSpectrumPerMeV(std::sqrt(10.0)),std::sqrt(1000.0),2.0e-13),
        "U-F01 TABLE log-log interpolation reference");

  // U-F07 -- analytic full-sphere means for raw and normalized PAD/spatial models.
  // Midpoint quadrature in mu converges rapidly; compare the discrete mean with the
  // exact beta-function normalizer and explicitly check normalized integral = 1.
  const std::vector<FN::DirectionSample> sky=FN::BuildEqualSolidAngleDirections(1000,2);
  double sinRawMean=0.0,cosRawMean=0.0,sinNormalizedMean=0.0;
  for (const FN::DirectionSample& d:sky) {
    sinRawMean+=d.solidAngleWeight_sr*BP::PadWeight(
        BP::PadModel::SinAlphaN,d.z,2.0,BP::NormalizationMode::Raw)/(4.0*FN::kPi);
    cosRawMean+=d.solidAngleWeight_sr*BP::PadWeight(
        BP::PadModel::CosAlphaN,d.z,2.0,BP::NormalizationMode::Raw)/(4.0*FN::kPi);
    sinNormalizedMean+=d.solidAngleWeight_sr*BP::PadWeight(
        BP::PadModel::SinAlphaN,d.z,2.0,BP::NormalizationMode::UnitMean)/(4.0*FN::kPi);
  }
  Check(Near(sinRawMean,2.0/3.0,2.0e-6) &&
        Near(cosRawMean,1.0/3.0,2.0e-6),
        "U-F07 sin^2 and cos^2 raw means match analytic references");
  Check(Near(sinNormalizedMean,1.0,2.0e-6),
        "U-F07 unit-mean PAD integrates to one");
  Check(Near(BP::SpatialWeight(BP::SpatialModel::DaysideNightside,1.0,
                 3.0,1.0,BP::NormalizationMode::UnitMean),1.5) &&
        Near(BP::SpatialWeight(BP::SpatialModel::DaysideNightside,-1.0,
                 3.0,1.0,BP::NormalizationMode::UnitMean),0.5),
        "U-F07 normalized day/night factors preserve unit hemispheric mean");

  // U-F08 -- exact row, geometric (log-intensity) midpoint, deterministic gap and
  // out-of-range policies.  The 10 -> 100 midpoint reference is sqrt(1000).
  std::vector<BP::TemporalSpectrumRow> rows(2);
  rows[0].time_s=0.0; rows[0].epochUTC="2000-01-01T00:00:00";
  rows[0].intensity=std::vector<double>({10.0,100.0});
  rows[1].time_s=10.0; rows[1].epochUTC="2000-01-01T00:00:10";
  rows[1].intensity=std::vector<double>({100.0,10.0});
  const BP::TemporalSpectrumSelection exact=BP::SelectTemporalSpectrum(
      rows,0.0,0.0,BP::OutOfRangePolicy::Clamp,BP::GapPolicy::InterpolateAndFlag);
  const BP::TemporalSpectrumSelection middle=BP::SelectTemporalSpectrum(
      rows,5.0,0.0,BP::OutOfRangePolicy::Clamp,BP::GapPolicy::InterpolateAndFlag);
  Check(exact.status==BP::TemporalStatus::Exact && exact.intensity==rows[0].intensity,
        "U-F08 identical timestamp selects exact row");
  Check(Near(middle.intensity[0],std::sqrt(1000.0),1.0e-13) &&
        Near(middle.intensity[1],std::sqrt(1000.0),1.0e-13),
        "U-F08 log-intensity midpoint reference");
  const BP::TemporalSpectrumSelection held=BP::SelectTemporalSpectrum(
      rows,5.0,2.0,BP::OutOfRangePolicy::Clamp,BP::GapPolicy::HoldNearest);
  Check(held.status==BP::TemporalStatus::GapHeld && held.gapFlag &&
        held.intensity==rows[0].intensity,
        "U-F08 gap hold uses deterministic earlier-row tie break");
  const BP::TemporalSpectrumSelection before=BP::SelectTemporalSpectrum(
      rows,-1.0,0.0,BP::OutOfRangePolicy::Zero,BP::GapPolicy::InterpolateAndFlag);
  Check(before.status==BP::TemporalStatus::ZeroBefore &&
        Near(before.intensity[0],0.0),
        "U-F08 out-of-range ZERO policy and status");
  CheckThrows([&]() { BP::SelectTemporalSpectrum(
      rows,5.0,2.0,BP::OutOfRangePolicy::Clamp,BP::GapPolicy::Fail); },
      "U-F08 gap FAIL policy");

  // U-F05/U-F09 plus F1/F15/F16-style numerical references.  A constant boundary
  // intensity and constant access make trapezoidal values exact.  The top-hat response
  // rate is J*G*DeltaE, the triangular response is J*G*base/2.
  const std::vector<double> energy{1.0,2.0,3.0};
  const std::vector<double> half{0.5,0.5,0.5};
  const std::vector<double> zero{0.0,0.0,0.0};
  const auto unitPerJ=[](double) { return 1.0; };
  std::vector<BP::EnergyChannel> channels;
  channels.push_back(BP::EnergyChannel{"edge",0.5,2.5});
  BP::DetectorResponse top;
  top.name="top"; top.energy_MeV={1.0,3.0}; top.relativeResponse={1.0,1.0};
  top.geometricFactor_m2_sr=2.0;
  BP::DetectorResponse triangle;
  triangle.name="triangle"; triangle.energy_MeV={1.0,2.0,3.0};
  triangle.relativeResponse={0.0,1.0,0.0}; triangle.geometricFactor_m2_sr=2.0;
  const BP::ProductSet products=BP::EvaluateIsotropicProducts(
      energy,half,half,half,FN::kAtomicMassUnit_kg,unitPerJ,channels,
      std::vector<BP::DetectorResponse>({top,triangle}));
  const double expectedOmni=4.0*FN::kPi*0.5*2.0*FN::kMeVToJ;
  Check(Near(products.omnidirectionalFlux_m2_s.nominal,expectedOmni) &&
        Near(products.oneWayPlanarFlux_m2_s.nominal,0.25*expectedOmni),
        "U-F09 analytic isotropic omni/one-way planar relation");
  Check(Near(products.channelFlux_m2_s[0].nominal,
             4.0*FN::kPi*0.5*1.5*FN::kMeVToJ),
        "U-F05 channel edge clipping reference");
  Check(Near(products.detectorRate_s[0].nominal,2.0*0.5*2.0*FN::kMeVToJ) &&
        Near(products.detectorRate_s[1].nominal,2.0*0.5*1.0*FN::kMeVToJ),
        "U-F09 top-hat and triangular detector response references");
  const BP::ProductSet blocked=BP::EvaluateIsotropicProducts(
      energy,zero,zero,zero,FN::kAtomicMassUnit_kg,unitPerJ);
  Check(Near(blocked.omnidirectionalFlux_m2_s.nominal,0.0) &&
        Near(blocked.numberDensity_m3.nominal,0.0),
        "F16 all-blocked access produces zero products");

  // Written-spectrum closure: trapezoid the stored omnidirectional differential
  // spectrum in per-MeV units and recover the reported integral exactly.
  double reconstructed=0.0;
  for (std::size_t i=0;i+1<products.spectrum.size();++i)
    reconstructed+=0.5*(products.spectrum[i].omnidirectionalPerMeV.nominal+
                        products.spectrum[i+1].omnidirectionalPerMeV.nominal)*
                   (products.spectrum[i+1].energy_MeV-products.spectrum[i].energy_MeV);
  Check(Near(reconstructed,products.omnidirectionalFlux_m2_s.nominal,2.0e-14),
        "F4 written differential spectrum reconstructs integral");

  // General phase-space branch reference: p_local/p_boundary=2 multiplies j by 4.
  const BP::Bounds phase=BP::MapBoundaryIntensity(
      BP::Bounds(0.5,0.4,0.6),BP::Bounds(10.0,9.0,11.0),
      BP::CharacteristicMapping::GeneralPhaseSpace,2.0,1.0);
  Check(Near(phase.nominal,20.0) && Near(phase.lower,14.4) && Near(phase.upper,26.4),
        "Step-6 j/p^2 phase-space mapping reference");

  std::remove(tableFile.c_str());
  if (failures!=0) {
    std::cerr << failures << " BoundaryProducts test(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: U-F01 U-F05 U-F07 U-F08 U-F09 and F1/F4/F15/F16 kernels\n";
  return EXIT_SUCCESS;
}
