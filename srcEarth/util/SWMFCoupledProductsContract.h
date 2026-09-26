#ifndef _SRC_EARTH_UTIL_SWMF_COUPLED_PRODUCTS_CONTRACT_H_
#define _SRC_EARTH_UTIL_SWMF_COUPLED_PRODUCTS_CONTRACT_H_

//======================================================================================
// SWMFCoupledProductsContract.h
//======================================================================================
//
// Dependency-free Roadmap Step-11 contract for SWMF-coupled energetic-particle
// density, differential-spectrum, integral-flux, channel, and detector-rate products.
//
// The expensive production path lives in 3d/DensityMode3D.cpp and depends on AMPS,
// MPI, and the live SWMF coupler.  The scientific identity and release checks below
// deliberately do not.  Keeping them here lets a small numerical test exercise the
// exact manifest/fingerprint rules without substituting a forward Monte-Carlo sampler
// for the backward-characteristic calculation.
//
// Step 11 has four important invariants:
//
//   * one immutable B/u snapshot, one absolute UTC, and one boundary-spectrum
//     evaluation epoch feed every requested product;
//   * the complete physics control surface (species, energy/access grid, boundary
//     spectrum, channels, and responses) has a deterministic, layout-independent
//     fingerprint;
//   * a PASS manifest contains a closed spectrum artifact, density and flux artifacts
//     (or a shell artifact containing both), and a trajectory-termination summary;
//   * unresolved trajectories retain conservative bounds and may not be hidden as
//     physical zero.  A PASS manifest requires the existing configured unresolved
//     tolerance; this file neither changes nor invents that tolerance.
//
// The fingerprint is an identity guard, not a cryptographic signature.  It uses the
// same canonical 64-bit FNV-1a implementation as the Step-9 snapshot contract and
// hashes floating-point values by their IEEE bit pattern after normalizing signed zero.
//======================================================================================

#include "SWMFCoupledAccessContract.h"
#include "SWMFSnapshotContract.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Earth {
namespace SWMFCoupledProducts {

struct EnergyChannelDefinition {
  std::string name;
  double lower_MeV{0.0};
  double upper_MeV{0.0};
};

struct DetectorResponseDefinition {
  std::string name;
  double lower_MeV{0.0};
  double upper_MeV{0.0};
  double geometricFactor_m2_sr{0.0};
};

struct ObservationDefinition {
  std::string epochUTC; // empty only for untimed POINTS inputs
  double x_km{0.0};
  double y_km{0.0};
  double z_km{0.0};
};

// Value-owned description of every input that changes the Step-6 fold.  Raw spectrum
// keys are included because spectrum families have different parameter sets; sorting
// those key/value pairs before hashing makes map iteration and parser insertion order
// irrelevant.  Channel and response order is intentionally retained because it is
// part of the emitted VARIABLES schema.
struct ProductControl {
  std::string outputMode;
  std::string speciesName;
  double charge_e{0.0};
  double mass_amu{0.0};
  std::string boundaryMode;
  std::string transmissionMode;
  double minimumEnergy_MeV{0.0};
  double maximumEnergy_MeV{0.0};
  int energyIntervals{0};
  int transmissionScanPoints{0};
  int maximumParticlesPerPoint{0};
  std::string energySpacing;
  std::string spectrumType;
  std::string energyBasis;
  double spectrumMassNumber{1.0};
  std::string intensityUnit;
  double spectrumRelativeUncertainty{0.0};
  std::vector<std::pair<std::string,std::string> > spectrumKeyValues;
  // For TABLE spectra these are the actual post-time-selection nodes used by the
  // fold, not merely the input filename.  Hashing them prevents a modified table file
  // at the same path from masquerading as the archived boundary distribution.
  std::vector<double> spectrumTableEnergy_MeV;
  std::vector<double> spectrumTableIntensityPerMeV;
  std::vector<EnergyChannelDefinition> channels;
  std::vector<DetectorResponseDefinition> detectorResponses;
  std::string coordinateFrame;
  std::vector<ObservationDefinition> observations;
  std::vector<double> shellAltitude_km;
  double shellResolution_deg{0.0};
  std::string shellGeometry;
};

inline void RequireFinite_(double value,const char* label) {
  if (!std::isfinite(value))
    throw std::invalid_argument(std::string("non-finite coupled product ")+label);
}

inline void ValidateProductControl(const ProductControl& control) {
  const std::string outputMode=Earth::SWMFCoupledAccess::UpperTrimmed(
      control.outputMode);
  if (outputMode!="POINTS" && outputMode!="TRAJECTORY" && outputMode!="SHELLS")
    throw std::invalid_argument(
        "coupled product output mode must be POINTS, TRAJECTORY, or SHELLS");
  if (control.speciesName.empty() || control.boundaryMode.empty() ||
      control.transmissionMode.empty() || control.energySpacing.empty() ||
      control.spectrumType.empty() || control.energyBasis.empty() ||
      control.intensityUnit.empty())
    throw std::invalid_argument("coupled product control is missing required metadata");
  RequireFinite_(control.charge_e,"charge");
  RequireFinite_(control.mass_amu,"mass");
  RequireFinite_(control.minimumEnergy_MeV,"minimum energy");
  RequireFinite_(control.maximumEnergy_MeV,"maximum energy");
  RequireFinite_(control.spectrumMassNumber,"spectrum mass number");
  RequireFinite_(control.spectrumRelativeUncertainty,"spectrum uncertainty");
  if (control.charge_e==0.0 || !(control.mass_amu>0.0) ||
      !(control.minimumEnergy_MeV>0.0) ||
      !(control.maximumEnergy_MeV>control.minimumEnergy_MeV) ||
      control.energyIntervals<1 || !(control.spectrumMassNumber>0.0) ||
      control.spectrumRelativeUncertainty<0.0)
    throw std::invalid_argument("coupled product control contains invalid physics values");
  if (control.spectrumTableEnergy_MeV.size()!=
      control.spectrumTableIntensityPerMeV.size())
    throw std::invalid_argument("coupled spectrum table arrays have different sizes");
  for (std::size_t i=0;i<control.spectrumTableEnergy_MeV.size();++i) {
    RequireFinite_(control.spectrumTableEnergy_MeV[i],"spectrum table energy");
    RequireFinite_(control.spectrumTableIntensityPerMeV[i],
                   "spectrum table intensity");
    if (!(control.spectrumTableEnergy_MeV[i]>0.0) ||
        control.spectrumTableIntensityPerMeV[i]<0.0 ||
        (i>0 && !(control.spectrumTableEnergy_MeV[i]>
                   control.spectrumTableEnergy_MeV[i-1])))
      throw std::invalid_argument("coupled spectrum table nodes are invalid");
  }

  std::vector<std::string> names;
  for (std::size_t i=0;i<control.channels.size();++i) {
    const EnergyChannelDefinition& channel=control.channels[i];
    RequireFinite_(channel.lower_MeV,"channel lower energy");
    RequireFinite_(channel.upper_MeV,"channel upper energy");
    if (channel.name.empty() || !(channel.lower_MeV>0.0) ||
        !(channel.upper_MeV>channel.lower_MeV) ||
        std::find(names.begin(),names.end(),channel.name)!=names.end())
      throw std::invalid_argument("invalid or duplicate coupled flux channel");
    names.push_back(channel.name);
  }
  names.clear();
  for (std::size_t i=0;i<control.detectorResponses.size();++i) {
    const DetectorResponseDefinition& response=control.detectorResponses[i];
    RequireFinite_(response.lower_MeV,"response lower energy");
    RequireFinite_(response.upper_MeV,"response upper energy");
    RequireFinite_(response.geometricFactor_m2_sr,"response geometric factor");
    if (response.name.empty() || !(response.lower_MeV>0.0) ||
        !(response.upper_MeV>response.lower_MeV) ||
        !(response.geometricFactor_m2_sr>=0.0) ||
        std::find(names.begin(),names.end(),response.name)!=names.end())
      throw std::invalid_argument("invalid or duplicate coupled detector response");
    names.push_back(response.name);
  }
  if (control.coordinateFrame.empty())
    throw std::invalid_argument("coupled product observation frame is empty");
  if (outputMode=="POINTS" || outputMode=="TRAJECTORY") {
    if (control.observations.empty())
      throw std::invalid_argument("coupled point/trajectory product has no observations");
    for (std::size_t i=0;i<control.observations.size();++i) {
      RequireFinite_(control.observations[i].x_km,"observation X");
      RequireFinite_(control.observations[i].y_km,"observation Y");
      RequireFinite_(control.observations[i].z_km,"observation Z");
      if (outputMode=="TRAJECTORY" && control.observations[i].epochUTC.empty())
        throw std::invalid_argument(
            "coupled trajectory observation is missing its ephemeris epoch");
    }
  }
  else {
    RequireFinite_(control.shellResolution_deg,"shell resolution");
    if (control.shellAltitude_km.empty() || !(control.shellResolution_deg>0.0) ||
        control.shellGeometry.empty())
      throw std::invalid_argument("coupled shell observation geometry is incomplete");
    for (std::size_t i=0;i<control.shellAltitude_km.size();++i) {
      RequireFinite_(control.shellAltitude_km[i],"shell altitude");
      if (control.shellAltitude_km[i]<0.0)
        throw std::invalid_argument("coupled shell altitude must be non-negative");
    }
  }
}

inline void AddControlCommon_(Earth::SWMFSnapshot::Detail::StableHash64& hash,
                              const ProductControl& control) {
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.outputMode));
  hash.AddString(control.speciesName);
  hash.AddDouble(control.charge_e);
  hash.AddDouble(control.mass_amu);
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.boundaryMode));
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.transmissionMode));
  hash.AddDouble(control.minimumEnergy_MeV);
  hash.AddDouble(control.maximumEnergy_MeV);
  hash.AddSigned(control.energyIntervals);
  hash.AddSigned(control.transmissionScanPoints);
  hash.AddSigned(control.maximumParticlesPerPoint);
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.energySpacing));
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.spectrumType));
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.energyBasis));
  hash.AddDouble(control.spectrumMassNumber);
  hash.AddString(control.intensityUnit);
  hash.AddDouble(control.spectrumRelativeUncertainty);

  std::vector<std::pair<std::string,std::string> > spectrum=
      control.spectrumKeyValues;
  std::sort(spectrum.begin(),spectrum.end());
  hash.AddUnsigned(static_cast<std::uint64_t>(spectrum.size()));
  for (std::size_t i=0;i<spectrum.size();++i) {
    hash.AddString(spectrum[i].first);
    hash.AddString(spectrum[i].second);
  }
  hash.AddUnsigned(static_cast<std::uint64_t>(
      control.spectrumTableEnergy_MeV.size()));
  for (std::size_t i=0;i<control.spectrumTableEnergy_MeV.size();++i) {
    hash.AddDouble(control.spectrumTableEnergy_MeV[i]);
    hash.AddDouble(control.spectrumTableIntensityPerMeV[i]);
  }
}

inline std::string BoundarySpectrumFingerprint(const ProductControl& control) {
  ValidateProductControl(control);
  Earth::SWMFSnapshot::Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/boundary-spectrum-state/v1");
  // Reuse the complete common block rather than attempting to guess which raw keys
  // matter for each spectrum family.  It includes the selected TABLE node values.
  AddControlCommon_(hash,control);
  return hash.Hex("spectrum-v1-");
}

inline std::string ChannelSchemaFingerprint(const ProductControl& control) {
  ValidateProductControl(control);
  Earth::SWMFSnapshot::Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/channel-schema/v1");
  hash.AddUnsigned(static_cast<std::uint64_t>(control.channels.size()));
  for (std::size_t i=0;i<control.channels.size();++i) {
    hash.AddString(control.channels[i].name);
    hash.AddDouble(control.channels[i].lower_MeV);
    hash.AddDouble(control.channels[i].upper_MeV);
  }
  return hash.Hex("channel-v1-");
}

inline std::string DetectorResponseFingerprint(const ProductControl& control) {
  ValidateProductControl(control);
  Earth::SWMFSnapshot::Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/detector-response-schema/v1");
  hash.AddUnsigned(static_cast<std::uint64_t>(control.detectorResponses.size()));
  for (std::size_t i=0;i<control.detectorResponses.size();++i) {
    hash.AddString(control.detectorResponses[i].name);
    hash.AddDouble(control.detectorResponses[i].lower_MeV);
    hash.AddDouble(control.detectorResponses[i].upper_MeV);
    hash.AddDouble(control.detectorResponses[i].geometricFactor_m2_sr);
  }
  return hash.Hex("response-v1-");
}

inline std::string ObservationStateFingerprint(const ProductControl& control) {
  ValidateProductControl(control);
  Earth::SWMFSnapshot::Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/observation-state/v1");
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.outputMode));
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.coordinateFrame));
  hash.AddUnsigned(static_cast<std::uint64_t>(control.observations.size()));
  for (std::size_t i=0;i<control.observations.size();++i) {
    hash.AddString(control.observations[i].epochUTC);
    hash.AddDouble(control.observations[i].x_km);
    hash.AddDouble(control.observations[i].y_km);
    hash.AddDouble(control.observations[i].z_km);
  }
  hash.AddUnsigned(static_cast<std::uint64_t>(control.shellAltitude_km.size()));
  for (std::size_t i=0;i<control.shellAltitude_km.size();++i)
    hash.AddDouble(control.shellAltitude_km[i]);
  hash.AddDouble(control.shellResolution_deg);
  hash.AddString(Earth::SWMFCoupledAccess::UpperTrimmed(control.shellGeometry));
  return hash.Hex("observation-v1-");
}

inline std::string ProductControlFingerprint(const ProductControl& control) {
  ValidateProductControl(control);
  Earth::SWMFSnapshot::Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/swmf-coupled-products-control/v1");
  AddControlCommon_(hash,control);
  hash.AddString(BoundarySpectrumFingerprint(control));
  hash.AddString(ChannelSchemaFingerprint(control));
  hash.AddString(DetectorResponseFingerprint(control));
  hash.AddString(ObservationStateFingerprint(control));
  return hash.Hex("products-v1-");
}

struct ProductRunSummary {
  ProductControl control;
  std::string spectrumEvaluationEpochUTC;
  std::string activeSpectrumTableEpochUTC;
  std::string spectrumTemporalStatus;
  bool spectrumTemporalGap{false};
  double spectrumTemporalFraction{0.0};
  int locationCount{0};
  int energyCount{0};
  int directionCount{0};
  long long sampled{0};
  long long retried{0};
  long long resolved{0};
  long long allowed{0};
  std::vector<long long> terminationCounts;
  double maximumUnresolvedFraction{0.0};
  double unresolvedTolerance{0.0};
  std::vector<std::string> artifacts;
  bool valid{false};
};

inline std::string ArtifactRole(const std::string& path) {
  // Classify the emitted filename, not parent directories.  A campaign directory
  // named "spectrum" must not turn its density artifact into a SPECTRUM role.
  const std::size_t slash=path.find_last_of("/\\");
  const std::string name=(slash==std::string::npos) ? path : path.substr(slash+1);
  if (name.find("termination_summary")!=std::string::npos) return "TERMINATION";
  if (name.find("spectrum")!=std::string::npos) return "SPECTRUM";
  if (name.find("density_flux")!=std::string::npos) return "DENSITY_FLUX";
  if (name.find("density")!=std::string::npos) return "DENSITY";
  if (name.find("flux")!=std::string::npos) return "FLUX";
  return "UNKNOWN";
}

inline void ValidateRunSummary(const ProductRunSummary& summary,
                               bool requireCompleteArtifacts) {
  ValidateProductControl(summary.control);
  if (!summary.valid || summary.spectrumEvaluationEpochUTC.empty() ||
      summary.spectrumTemporalStatus.empty() || summary.locationCount<1 ||
      summary.energyCount<2 || summary.directionCount<1)
    throw std::invalid_argument("coupled product run summary is incomplete");
  RequireFinite_(summary.spectrumTemporalFraction,"temporal fraction");
  RequireFinite_(summary.maximumUnresolvedFraction,"maximum unresolved fraction");
  RequireFinite_(summary.unresolvedTolerance,"unresolved tolerance");
  if (summary.spectrumTemporalFraction<0.0 || summary.spectrumTemporalFraction>1.0 ||
      summary.maximumUnresolvedFraction<0.0 ||
      summary.maximumUnresolvedFraction>1.0 || summary.unresolvedTolerance<0.0 ||
      summary.unresolvedTolerance>1.0)
    throw std::invalid_argument("coupled product run summary contains invalid fractions");
  if (summary.sampled<=0 || summary.retried<0 || summary.resolved<0 ||
      summary.allowed<0 || summary.resolved>summary.sampled ||
      summary.allowed>summary.resolved)
    throw std::invalid_argument("coupled product termination counts are invalid");
  long long terminationTotal=0;
  for (std::size_t i=0;i<summary.terminationCounts.size();++i) {
    if (summary.terminationCounts[i]<0)
      throw std::invalid_argument("negative coupled product termination count");
    terminationTotal+=summary.terminationCounts[i];
  }
  if (terminationTotal!=summary.sampled)
    throw std::invalid_argument(
        "coupled product termination counts do not close to sampled trajectories");
  if (requireCompleteArtifacts &&
      summary.maximumUnresolvedFraction>summary.unresolvedTolerance)
    throw std::invalid_argument(
        "coupled product unresolved fraction exceeds the configured release gate");

  bool haveDensity=false,haveFlux=false,haveCombined=false;
  bool haveSpectrum=false,haveTermination=false;
  for (std::size_t i=0;i<summary.artifacts.size();++i) {
    if (summary.artifacts[i].empty() ||
        std::find(summary.artifacts.begin(),summary.artifacts.begin()+i,
                  summary.artifacts[i])!=summary.artifacts.begin()+i)
      throw std::invalid_argument("coupled product artifact list is empty or duplicated");
    const std::string role=ArtifactRole(summary.artifacts[i]);
    haveDensity|=(role=="DENSITY");
    haveFlux|=(role=="FLUX");
    haveCombined|=(role=="DENSITY_FLUX");
    haveSpectrum|=(role=="SPECTRUM");
    haveTermination|=(role=="TERMINATION");
    if (role=="UNKNOWN")
      throw std::invalid_argument("unclassified coupled product artifact");
  }
  if (requireCompleteArtifacts &&
      (!haveSpectrum || !haveTermination ||
       !(haveCombined || (haveDensity && haveFlux))))
    throw std::invalid_argument(
        "PASS coupled products require spectrum, density, flux, and termination artifacts");
}

inline std::string BuildProductsManifestJson(
    const std::string& result,const std::string& snapshotId,
    const std::string& contentFingerprint,const std::string& meshRevision,
    const std::string& epochUTC,double simulationTime_s,
    const std::string& boundaryPolicy,const std::string& outputSuffix,
    const ProductRunSummary& summary,const std::string& message) {
  if (result!="PASS" && result!="FAILED")
    throw std::invalid_argument("coupled product manifest result must be PASS or FAILED");
  if (snapshotId.empty() || contentFingerprint.empty() || meshRevision.empty() ||
      epochUTC.empty())
    throw std::invalid_argument(
        "coupled product manifest requires complete snapshot provenance");
  RequireFinite_(simulationTime_s,"simulation time");
  if (simulationTime_s<0.0)
    throw std::invalid_argument("coupled product simulation time must be non-negative");
  const std::string canonicalBoundary=
      Earth::SWMFCoupledAccess::UpperTrimmed(boundaryPolicy);
  if (canonicalBoundary!="BOX" && canonicalBoundary!="SHUE")
    throw std::invalid_argument("coupled product boundary policy must be BOX or SHUE");
  if (outputSuffix!=Earth::SWMFCoupledAccess::BuildProductSuffix(
          simulationTime_s,snapshotId))
    throw std::invalid_argument(
        "coupled product manifest suffix does not match time/snapshot identity");
  if (summary.spectrumEvaluationEpochUTC!=epochUTC)
    throw std::invalid_argument(
        "field and boundary-spectrum evaluation epochs are not synchronized");
  ValidateRunSummary(summary,result=="PASS");

  const std::string controlFingerprint=ProductControlFingerprint(summary.control);
  const std::string spectrumFingerprint=BoundarySpectrumFingerprint(summary.control);
  const std::string channelFingerprint=ChannelSchemaFingerprint(summary.control);
  const std::string responseFingerprint=DetectorResponseFingerprint(summary.control);
  const std::string observationFingerprint=ObservationStateFingerprint(summary.control);
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << "{\n"
      << "  \"schema\": \"sep-in-geospace/swmf-coupled-products/v1\",\n"
      << "  \"RESULT\": \"" << result << "\",\n"
      << "  \"phase_1_interpretation\": \"INSTANTANEOUS_QUASI_STATIC\",\n"
      << "  \"characteristic_mapping\": \"STATIC_MAGNETIC\",\n"
      << "  \"snapshot_id\": \"" << Earth::SWMFCoupledAccess::JsonEscape(snapshotId) << "\",\n"
      << "  \"content_fingerprint\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(contentFingerprint) << "\",\n"
      << "  \"mesh_revision\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(meshRevision) << "\",\n"
      << "  \"epoch_utc\": \"" << Earth::SWMFCoupledAccess::JsonEscape(epochUTC) << "\",\n"
      << "  \"boundary_spectrum_evaluation_epoch_utc\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(summary.spectrumEvaluationEpochUTC)
      << "\",\n"
      << "  \"active_spectrum_table_epoch_utc\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(summary.activeSpectrumTableEpochUTC)
      << "\",\n"
      << "  \"simulation_time_s\": " << std::setprecision(17)
      << simulationTime_s << ",\n"
      << "  \"output_mode\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(
             Earth::SWMFCoupledAccess::UpperTrimmed(summary.control.outputMode))
      << "\",\n"
      << "  \"outer_boundary_policy\": \"" << canonicalBoundary << "\",\n"
      << "  \"output_suffix\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(outputSuffix) << "\",\n"
      << "  \"product_control_fingerprint\": \"" << controlFingerprint << "\",\n"
      << "  \"boundary_spectrum_fingerprint\": \"" << spectrumFingerprint << "\",\n"
      << "  \"channel_schema_fingerprint\": \"" << channelFingerprint << "\",\n"
      << "  \"detector_response_fingerprint\": \"" << responseFingerprint << "\",\n"
      << "  \"observation_state_fingerprint\": \"" << observationFingerprint << "\",\n"
      << "  \"spectrum_energy_basis\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(summary.control.energyBasis) << "\",\n"
      << "  \"spectrum_mass_number\": " << summary.control.spectrumMassNumber << ",\n"
      << "  \"spectrum_intensity_unit\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(summary.control.intensityUnit) << "\",\n"
      << "  \"spectrum_relative_uncertainty\": "
      << summary.control.spectrumRelativeUncertainty << ",\n"
      << "  \"spectrum_temporal_status\": \""
      << Earth::SWMFCoupledAccess::JsonEscape(summary.spectrumTemporalStatus) << "\",\n"
      << "  \"spectrum_temporal_gap\": "
      << (summary.spectrumTemporalGap ? "true" : "false") << ",\n"
      << "  \"spectrum_temporal_fraction\": "
      << summary.spectrumTemporalFraction << ",\n"
      << "  \"location_count\": " << summary.locationCount << ",\n"
      << "  \"energy_count\": " << summary.energyCount << ",\n"
      << "  \"direction_count\": " << summary.directionCount << ",\n"
      << "  \"termination\": {\"sampled\": " << summary.sampled
      << ", \"retried\": " << summary.retried
      << ", \"resolved\": " << summary.resolved
      << ", \"allowed\": " << summary.allowed
      << ", \"maximum_unresolved_fraction\": "
      << summary.maximumUnresolvedFraction
      // Channel and detector response weights are non-negative. Their weighted
      // unresolved fraction therefore cannot exceed the maximum contributing
      // energy/location bin. Publishing that conservative upper bound makes the
      // response-support gate auditable without recomputing a second product fold.
      << ", \"response_weighted_unresolved_upper_bound\": "
      << summary.maximumUnresolvedFraction
      << ", \"unresolved_tolerance\": " << summary.unresolvedTolerance
      << ", \"counts\": [";
  for (std::size_t i=0;i<summary.terminationCounts.size();++i) {
    if (i!=0) out << ", ";
    out << summary.terminationCounts[i];
  }
  out << "]},\n"
      << "  \"channels\": [";
  for (std::size_t i=0;i<summary.control.channels.size();++i) {
    if (i!=0) out << ", ";
    out << "\"" << Earth::SWMFCoupledAccess::JsonEscape(
        summary.control.channels[i].name) << "\"";
  }
  out << "],\n  \"detector_responses\": [";
  for (std::size_t i=0;i<summary.control.detectorResponses.size();++i) {
    if (i!=0) out << ", ";
    out << "\"" << Earth::SWMFCoupledAccess::JsonEscape(
        summary.control.detectorResponses[i].name) << "\"";
  }
  out << "],\n  \"artifacts\": [";
  for (std::size_t i=0;i<summary.artifacts.size();++i) {
    if (i!=0) out << ", ";
    out << "{\"role\": \"" << ArtifactRole(summary.artifacts[i])
        << "\", \"path\": \""
        << Earth::SWMFCoupledAccess::JsonEscape(summary.artifacts[i]) << "\"}";
  }
  out << "],\n"
      << "  \"message\": \"" << Earth::SWMFCoupledAccess::JsonEscape(message)
      << "\"\n}\n";
  return out.str();
}

} // namespace SWMFCoupledProducts
} // namespace Earth

#endif // _SRC_EARTH_UTIL_SWMF_COUPLED_PRODUCTS_CONTRACT_H_
