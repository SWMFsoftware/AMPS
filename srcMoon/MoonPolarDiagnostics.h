#ifndef SRCMOON_MOON_POLAR_DIAGNOSTICS_H
#define SRCMOON_MOON_POLAR_DIAGNOSTICS_H

#include <map>
#include <string>
#include <vector>

namespace Moon {
namespace PolarDiagnostics {

// Piecewise-linear evaluation of the MODEL boundary/layer-centre state [K]
// at explicit depth [m]. There is no observational interpolation. The surface
// is z=0; layers must be strictly increasing positive depths. A positive-depth
// query outside the last centre is rejected, never extrapolated. Failure is
// transactional, leaving *temperatureK untouched.
bool InterpolateThermalProfileK(double depthM,double boundaryTemperatureK,
    const std::vector<double>& layerDepthM,
    const std::vector<double>& layerTemperatureK,
    double* temperatureK,std::string* error);

// Thin I20 adapter around the NORMAL production surface/thermal services.
// The configured geometry-only CSV declares site, logical model support,
// explicit UTC/TDB epoch, depth [m], sample weight and temporal-support ID.
// Output records individual samples, not a guessed annual climatology. A
// separately qualified operator must aggregate them under the frozen native
// product's time/quantity definition before observational scoring.
//
// Rank zero evaluates this deterministic thermal history; other ranks return
// without work. Failure preserves partial output and returns an explanation;
// the owning linked runner assigns ERROR/SKIPPED according to its contract.
bool RunPolarDiagnostics(const std::string& outputPrefix,
    std::map<std::string,std::string>* metrics,std::string* error);

}  // namespace PolarDiagnostics
}  // namespace Moon

#endif
