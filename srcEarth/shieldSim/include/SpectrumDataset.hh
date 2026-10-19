#ifndef SHIELDSIM_SPECTRUM_DATASET_HH
#define SHIELDSIM_SPECTRUM_DATASET_HH

/*
 * Canonical spectrum container shared by the Stage-1 legacy source, the
 * Stage-2 published-data adapters, and the Stage-11 calculated-spectrum
 * exchange contract.
 *
 * The only Stage-1 external adapter is the explicitly labeled legacy
 * three-column adapter.  A sidecar is mandatory because the old numeric file
 * alone cannot say whether energy means total energy or energy/nucleon, nor
 * whether values are rates, fluence, intensity, or arbitrary test weights.
 */

#include "StrictJson.hh"

#include <string>
#include <vector>

namespace shield_sim {

/*
 * A component is one species/direction/time aggregate on its own energy grid.
 * The legacy three-column format happened to put protons and alpha particles
 * on one grid, but published products do not: AMS-02 has different proton and
 * helium grids and GOES has two directional sensors.  Keeping components
 * separate prevents an adapter from manufacturing common bins or silently
 * averaging directions merely to fit the historical representation.
 *
 * Energy bounds and representative coordinates are in the canonical unit
 * named by energyUnit.  A representative value is optional because a bin
 * average does not generally have a unique effective energy.  Published
 * adapters therefore retain bounds even when no representative coordinate is
 * defensible.
 */
struct SpectrumComponent {
  std::string componentId;
  std::string speciesId;
  std::string geant4Particle;
  int atomicNumber = 0;
  // Zero means that the published population mixes isotopes and must not be
  // relabeled as a single transport ion (for example AMS-02 inclusive helium).
  int massNumber = 0;
  std::string energyVariable;
  std::string energyUnit;
  std::string valueType;
  std::string valueUnit;
  std::string angularRepresentation;
  std::string directionId;
  std::string uncertaintySemantics;
  std::vector<double> binLower;
  std::vector<double> binUpper;
  std::vector<double> representative;
  std::vector<double> value;
  // Absolute lower/upper curves are distinct from uncertainty magnitudes.
  // AMPS unresolved-trajectory bounds are conservative physical envelopes,
  // not standard deviations, so they must never be placed in the statistical
  // uncertainty arrays below.
  std::vector<double> lowerBound;
  std::vector<double> upperBound;
  std::vector<double> totalUncertainty;
  // Published products often separate uncertainty components because they
  // have different correlations.  In particular, the AMS-02 time-independent
  // term must not be treated as an independent random error at every Bartels
  // rotation.  Keeping these vectors separate prevents a later comparison
  // routine from manufacturing an unsupported diagonal covariance matrix.
  std::vector<double> statisticalUncertainty;
  std::vector<double> timeDependentSystematicUncertainty;
  std::vector<double> timeIndependentSystematicUncertainty;

  // Some count-derived products contain legitimate zero-count bins.  A zero
  // flux in such a bin does not have zero statistical uncertainty, so the
  // count and validity mask are retained explicitly.  Invalid bins must be
  // excluded from quantitative comparisons rather than interpreted as exact
  // zero measurements.
  std::vector<double> counts;
  std::vector<bool> valid;

  // Calculated AMPS spectra carry trajectory-resolution diagnostics rather
  // than statistical counting errors.  Keep those diagnostics next to the
  // values so a consumer never has to recover quality from a plot-oriented
  // Tecplot file or assume that an unresolved trajectory was forbidden.
  std::vector<double> unresolvedFraction;
  std::vector<double> sampledDirections;
  std::vector<double> resolvedDirections;
  std::vector<double> allowedDirections;
};

struct SpectrumDataset {
  std::string id;
  std::string energyVariable;
  std::string energyUnit;
  std::string valueType;
  std::string valueUnit;
  std::string angularModel;
  std::string scienceStatus;
  std::string qualityStatus;
  double maximumUnresolvedFraction = 0.0;
  double unresolvedFractionAcceptanceLimit = 0.0;
  std::string spectrumPath;
  std::string metadataPath;
  std::string contentHash;
  Json canonicalMetadata;
  std::vector<double> energyMeV;
  std::vector<double> binWidthMeV;
  std::vector<double> protonValue;
  std::vector<double> alphaValue;

  // A published differential spectrum is sampled as a piecewise-constant
  // density over its documented bin edges.  These arrays are populated only
  // after a component has passed the explicit transport-admission checks in
  // CaseDefinition; ingestion alone never makes a component transportable.
  // The legacy center-only adapter intentionally leaves them empty and keeps
  // its historical representative-energy sampling behavior.
  std::vector<double> transportBinLowerMeV;
  std::vector<double> transportBinUpperMeV;
  // Empty endpoint arrays mean the bins above are piecewise constant.  When
  // both arrays match the bin count, each interval is a point-spectrum
  // segment: positive endpoints are joined in log(E)-log(J), while a segment
  // touching zero is joined linearly.  PrimaryGeneratorAction integrates and
  // samples that declared curve exactly, including clipped energy limits.
  std::vector<double> transportValueLowerPerMeV;
  std::vector<double> transportValueUpperPerMeV;
  std::string transportInterpolation;
  std::string transportComponentId;

  // Versioned adapters populate components. The older vectors above remain a
  // compatibility projection consumed by the direct generator. Published
  // datasets are not projected during ingestion; CaseDefinition may create a
  // transport projection only after a selected component passes its explicit
  // particle/unit/angular-purpose admission contract.
  std::vector<SpectrumComponent> components;

  static SpectrumDataset LoadLegacyThreeColumn(const std::string& spectrumPath,
                                               const std::string& metadataPath);

  // Read the strict, self-describing Stage-11 JSON exchange document.  This
  // method performs ingestion only.  CaseDefinition separately decides
  // whether a selected component is physically admissible for transport.
  static SpectrumDataset LoadSpectrumExchangeV1(const std::string& spectrumPath);
};

} // namespace shield_sim

#endif
