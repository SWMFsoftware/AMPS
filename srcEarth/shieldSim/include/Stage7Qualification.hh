#ifndef SHIELDSIM_STAGE7_QUALIFICATION_HH
#define SHIELDSIM_STAGE7_QUALIFICATION_HH

/*
 * Stage-7 published-spectrum qualification utilities.
 *
 * These algorithms deliberately sit outside Geant4.  Transport produces a
 * response at declared incident nodes; this module performs the independently
 * testable operations that connect those responses to changing spectra and
 * time intervals.  None of the routines invents a particle isotope, angular
 * boundary, missing cadence, covariance, or response outside its tabulated
 * domain.
 */

#include "StrictJson.hh"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace shield_sim {

struct AdaptiveResponseFoldResult {
  double value = 0.0;
  double previousValue = 0.0;
  double finalRelativeChange = 0.0;
  std::size_t subdivisionsPerSourceBin = 0;
  std::size_t refinements = 0;
  std::size_t responseEvaluations = 0;
  bool converged = false;
};

// Linear interpolation is applied to the response ordinate as a function of
// log(E).  Queries outside the closed node domain are errors: Stage 7 has no
// extrapolation mode.
double InterpolateResponseInLogEnergy(const std::vector<double>& energyMeV,
                                      const std::vector<double>& response,
                                      double queryEnergyMeV);

// Integrate a piecewise-constant differential source over its exact energy
// edges.  Composite midpoint quadrature is repeatedly doubled until the fold
// changes by no more than relativeTolerance.  Empty/zero folds use an absolute
// machine-scale comparison and never divide by zero.
AdaptiveResponseFoldResult FoldPiecewiseConstantSpectrum(
    const std::vector<double>& sourceEdgesMeV,
    const std::vector<double>& differentialSource,
    const std::vector<double>& responseEnergyMeV,
    const std::vector<double>& response,
    double relativeTolerance = 0.02,
    std::size_t maximumRefinements = 12);

struct NeymanStratum {
  std::string id;
  double sourceWeight = 0.0;
  double pilotStandardDeviation = 0.0;
};

struct NeymanAllocation {
  std::string id;
  std::size_t histories = 0;
};

// Allocate independent histories across source/response strata.  The optimum
// continuous allocation is proportional to |weight|*pilot_sigma.  A declared
// minimum protects rare or zero-pilot strata; deterministic largest-remainder
// rounding preserves the exact total and makes reruns reproducible.
std::vector<NeymanAllocation> AllocateNeymanHistories(
    const std::vector<NeymanStratum>& strata,
    std::size_t totalHistories,
    std::size_t minimumPerStratum);

// Parse exactly YYYY-MM-DDTHH:MM:SSZ in the proleptic Gregorian calendar.
// Numeric offsets, local time, fractional seconds, and leap-second spellings
// are rejected because the Stage-7 exchange time system is explicitly UTC.
std::int64_t ParseUtcTimestamp(const std::string& utc);

// Strict JSON command kernels used by shieldSimExposureTool and core tests.
// Each returns a canonical, provenance-bearing result document.
Json IntegrateTimeResolvedExposure(const Json& input);
Json FoldSpectrumDocument(const Json& input);
Json AllocateHistoriesDocument(const Json& input);

} // namespace shield_sim

#endif
