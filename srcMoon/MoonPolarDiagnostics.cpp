// Production-path I20 thermal sampling. This file contains no heat equation,
// source law, observation correction, or native-data interpolation.
#include "pic.h"
#include "Moon.h"
#include "MoonPolarDiagnostics.h"
#include "MoonInput.h"
#include "MoonSurfaceTemperature.h"
#include "MoonDivinerObservation.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <new>
#include <set>
#include <sstream>
#include <vector>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

namespace Moon {
namespace PolarDiagnostics {
namespace {

struct Query {
  std::string id,site,timeSupport;
  std::size_t element=0;
  double epoch=0,depthM=0,weight=0;
};

// This is the same restricted, unquoted numeric/identifier CSV convention as
// the native integration diagnostics. Commas/quotes in an identity are errors,
// not a reason to truncate or silently change a native site identifier.
std::vector<std::string> Fields(const std::string& line) {
  std::vector<std::string> fields;
  std::istringstream stream(line);
  std::string field;
  while (std::getline(stream,field,',')) fields.push_back(field);
  if (!line.empty() && line.back()==',') fields.emplace_back();
  return fields;
}

bool Scalar(const std::string& text,double* value) {
  if (text.empty()) return false;
  char* end=nullptr;
  errno=0;
  *value=std::strtod(text.c_str(),&end);
  return errno!=ERANGE && end!=text.c_str() && *end=='\0' && std::isfinite(*value);
}

std::string Number(double value) {
  std::ostringstream stream;
  stream << std::setprecision(17) << value;
  return stream.str();
}

// Only the MODEL state is interpolated: z=0 is its actual boundary value;
// positive z values interpolate between that boundary and finite-volume layer
// centres. No observational temperature is ever supplied to this operation.
// The deepest layer centre is the last supported point; no extrapolation to
// the numerical bottom or an unobserved depth is performed.
bool AtDepth(std::size_t element,double depth,double surface,double* value,
             std::string* error) {
  std::vector<double> depths,temperatures;
  if (!SurfaceTemperature::ProductionThermalProfile(element,&depths,&temperatures,error))
    return false;
  return InterpolateThermalProfileK(depth,surface,depths,temperatures,value,error);
}

}  // namespace

bool InterpolateThermalProfileK(double depth,double boundary,
    const std::vector<double>& depths,const std::vector<double>& temperatures,
    double* value,std::string* error) {
  if (value==nullptr || error==nullptr) return false;
  error->clear();
  if (!std::isfinite(depth) || depth<0 || !std::isfinite(boundary) || boundary<=0 ||
      depths.size()!=temperatures.size()) {
    *error="invalid thermal model boundary/depth/profile dimensions";return false;
  }
  double previous=0;
  for (std::size_t i=0;i<depths.size();++i) {
    if (!std::isfinite(depths[i]) || depths[i]<=previous ||
        !std::isfinite(temperatures[i]) || temperatures[i]<=0) {
      *error="thermal model layer depths/temperatures must be finite, ordered and positive";return false;
    }
    previous=depths[i];
  }
  if (depth==0) { *value=boundary;return true; }
  if (depths.empty() || depth>depths.back()) {
    *error="I20 query depth exceeds the production layer-centre domain; "
        "supply a supported depth or explicitly refine/deepen the thermal model";return false;
  }
  const auto upper=std::lower_bound(depths.begin(),depths.end(),depth);
  const std::size_t i=static_cast<std::size_t>(upper-depths.begin());
  const double lowDepth=i ? depths[i-1]:0;
  const double lowTemperature=i ? temperatures[i-1]:boundary;
  const double temperature=(*upper==depth) ? temperatures[i]:
      lowTemperature+(temperatures[i]-lowTemperature)*(depth-lowDepth)/(*upper-lowDepth);
  if (!std::isfinite(temperature) || temperature<=0) {
    *error="nonphysical interpolated production thermal temperature";return false;
  }
  *value=temperature;return true;
}

bool RunPolarDiagnostics(const std::string& outputPrefix,
    std::map<std::string,std::string>* metrics,std::string* error) {
  if (metrics==nullptr || error==nullptr) return false;
  error->clear();
  if (PIC::ThisThread!=0) return true;
#ifdef _NO_SPICE_CALLS_
  *error="I20 thermal history requires the actual SPICE-enabled production model";
  return false;
#else
  if (!Runtime::HasConfiguration()) { *error="I20 requires an explicit application receipt"; return false; }
  const auto& configuration=Runtime::GetConfiguration();
  if (!configuration.m2EnvironmentEnabled || configuration.surfaceTemperatureMode!="thermal-model" ||
      configuration.m2PolarQueriesFile.empty() || !SurfaceTemperature::ProductionStateInitialized()) {
    *error="I20 needs thermal-model production state and an explicit m2_polar_queries_file";
    return false;
  }
  const std::filesystem::path target=outputPrefix+".I20_thermal_model.csv";
  const std::filesystem::path partial=target.string()+".partial";
  if (std::filesystem::exists(target) || std::filesystem::exists(partial)) {
    *error="I20 refuses to overwrite an earlier complete or partial thermal export";
    return false;
  }
  try {
    std::ifstream input(configuration.m2PolarQueriesFile);
    std::string line;
    if (!input || !std::getline(input,line)) {
      *error="I20 geometry-only query CSV is unavailable or empty";
      return false;
    }
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line!="query_id,site_id,logical_surface_element,epoch,time_system,depth_m,sample_weight,time_support_id") {
      *error="I20 geometry-only query CSV has a missing or unsupported exact header";
      return false;
    }
    std::vector<Query> queries;
    std::set<std::string> identities;
    while (std::getline(input,line)) {
      if (!line.empty() && line.back()=='\r') line.pop_back();
      if (line.empty()) continue;
      const auto fields=Fields(line);
      Query query;
      double element=0;
      if (fields.size()!=8 || line.find('"')!=std::string::npos ||
          fields[0].empty() || fields[1].empty() || fields[7].empty() ||
          !Scalar(fields[2],&element) || element<0 || std::floor(element)!=element ||
          element>std::numeric_limits<int>::max() ||
          !Scalar(fields[5],&query.depthM) || query.depthM<0 ||
          !Scalar(fields[6],&query.weight) || query.weight<0) {
        *error="invalid I20 query identity, integer element, nonnegative depth/weight or CSV shape";
        return false;
      }
      query.id=fields[0];query.site=fields[1];query.timeSupport=fields[7];
      query.element=static_cast<std::size_t>(element);
      if (!identities.insert(query.id).second) { *error="duplicate I20 query_id"; return false; }
      if (fields[4]=="UTC") {
        str2et_c(fields[3].c_str(),&query.epoch);
        if (failed_c()) { *error="SPICE cannot parse the explicitly UTC I20 query epoch"; return false; }
      } else if (fields[4]=="TDB") {
        if (!Scalar(fields[3],&query.epoch)) { *error="invalid explicit TDB query seconds"; return false; }
      } else { *error="I20 query time_system must explicitly be UTC or TDB"; return false; }
      queries.push_back(query);
    }
    if (!input.eof() || queries.empty()) { *error="I20 query CSV is empty or unreadable"; return false; }
    std::sort(queries.begin(),queries.end(),[](const Query& a,const Query& b) {
      return a.epoch!=b.epoch ? a.epoch<b.epoch:a.id<b.id;
    });
    // Input geometry is validated before the potentially long spin-up. Its
    // end epoch is explicit in .in; there is no inference from PRP averages.
    SpiceDouble start=0,end=0;
    str2et_c(configuration.thermalSpinupCycleStartUtc.c_str(),&start);
    str2et_c(configuration.thermalSpinupCycleEndUtc.c_str(),&end);
    if (failed_c() || !(end>start) || queries.front().epoch<end) {
      *error="I20 history precedes its declared independent spin-up end, or the interval is invalid";
      return false;
    }
    std::vector<SurfaceTemperature::SurfaceEnvironmentSample> samples;
    if (!SurfaceTemperature::SampleProductionSurfaceEnvironment(end,&samples,error)) return false;
    for (const auto& query:queries) if (query.element>=samples.size()) {
      *error="I20 site mapping names an unavailable production logical surface element";return false;
    }
    SurfaceTemperature::SpinUpMetrics spinup;
    if (!SurfaceTemperature::SpinUpProductionState(&spinup,error)) return false;
    // Preserve the ACTUAL complete model initial condition, after the
    // existing periodic solver has been applied and before any query advance.
    // This is model-state provenance, never an independent observation or an
    // automatically qualified thermal reference. The selected I20 diagnostic
    // explicitly requests it, independently of particle restart checkpoints.
    const std::string initialStatePath=outputPrefix+".I20_initial_thermal_state.csv";
    if (!SurfaceTemperature::WriteProductionInitialThermalStateCsv(initialStatePath,error)) return false;
    SurfaceTemperature::ProductionThermalInitializationInfo initialization;
    if (!SurfaceTemperature::GetProductionThermalInitialization(&initialization,error)) return false;
    (*metrics)["initial_thermal_state_file"]=initialStatePath;
    (*metrics)["initial_thermal_surface_elements"]=std::to_string(samples.size());
    (*metrics)["initial_thermal_layer_count"]=std::to_string(configuration.thermal.layerCount);
    (*metrics)["initial_thermal_state_tdb_s"]=Number(initialization.ephemerisTimeS);
    // The dispatcher records the shared initialization mode/metrics once
    // after this operator returns; avoid duplicate keys in the native receipt.
    std::ofstream output(partial);
    if (!output) { *error="cannot create I20 thermal CSV; inspect output permissions/free disk"; return false; }
    output << std::setprecision(17)
        << "query_id,site_id,logical_surface_element,ephemeris_time_tdb_s,depth_m,sample_weight,time_support_id,"
           "model_latitude_deg,model_longitude_east_deg,model_physical_temperature_K,"
           "model_bolometric_brightness_temperature_K,absorbed_solar_flux_W_m-2,terrain_sun_visible,"
           "earth_visible_disc_fraction,frame,surface_retention_timescale_s,"
           "na_surface_residence_time_s,na_surface_thermal_retention,"
           "ar_surface_residence_time_s,ar_surface_thermal_retention\n";
    double epoch=end,maximumChange=0;
    bool haveSamples=false;
    for (const auto& query:queries) {
      if (query.epoch>epoch) {
        if (!SurfaceTemperature::AdvanceProductionState(epoch,query.epoch,error)) return false;
        maximumChange=std::max(maximumChange,SurfaceTemperature::LastAdvanceMaximumSurfaceTemperatureChangeK());
        epoch=query.epoch;haveSamples=false;
      }
      if (!haveSamples) {
        if (!SurfaceTemperature::SampleProductionSurfaceEnvironment(epoch,&samples,error)) return false;
        haveSamples=true;
      }
      const auto& sample=samples[query.element];
      double temperature=0;
      if (!AtDepth(query.element,query.depthM,sample.surfaceTemperatureK,&temperature,error)) return false;
      output << query.id << ',' << query.site << ',' << query.element << ',' << epoch << ','
          << query.depthM << ',' << query.weight << ',' << query.timeSupport << ','
          << sample.latitudeDeg << ',' << sample.longitudeEastDeg << ',' << temperature << ',';
      // Subsurface kinetic temperature is not a brightness temperature. An
      // empty field preserves the distinction instead of applying epsilon
      // sigma T^4 to a layer that the instrument does not directly observe.
      if (query.depthM==0)
        output << DivinerObservation::GreyBolometricBrightnessTemperatureK(
            temperature,configuration.thermal.infraredEmissivity);
      output << ',' << sample.absorbedSolarFluxWm2 << ',' << (sample.terrainSunVisible?1:0)
          << ',' << sample.solarDiscVisibleFraction << ',' << Frames::BodyFixed;
      // Retention is a kinetic surface-state diagnostic at this exact epoch,
      // not an annual statistic or a permanent/geometric PSR label. Read the
      // same species kernel as the actual desorption source, independently
      // of checkpoint writes. Never assign surface kinetics to a 2-cm layer.
      if (query.depthM==0) {
        double sodiumResidence=0,argonResidence=0;
        bool sodiumRetained=false,argonRetained=false;
        if (!SurfaceTemperature::SampleProductionElementKinetics(
                static_cast<int>(query.element),_NA_SPEC_,&sodiumResidence,&sodiumRetained,error) ||
            !SurfaceTemperature::SampleProductionElementKinetics(
                static_cast<int>(query.element),_AR_SPEC_,&argonResidence,&argonRetained,error))
          return false;
        output << ',' << SurfaceTemperature::ProductionRetentionTimescaleS()
            << ',' << sodiumResidence << ',' << (sodiumRetained?1:0)
            << ',' << argonResidence << ',' << (argonRetained?1:0);
      } else output << ",,,,,";
      output << '\n';
      if (!output) { *error="I20 thermal CSV write failed; preserve partial file and free disk"; return false; }
    }
    output.close();
    if (!output) { *error="I20 thermal CSV close failed"; return false; }
    std::filesystem::rename(partial,target);
    (*metrics)["thermal_model_rows"]=std::to_string(queries.size());
    (*metrics)["thermal_model_csv"]=target.string();
    (*metrics)["surface_kinetics_export_schema"]="native-species-residence-retention-v1";
    (*metrics)["thermal_export_scope"]="production-per-epoch-state; annual aggregation requires independent qualification";
    (*metrics)["thermal_spinup_completed_cycles"]=std::to_string(spinup.completedCycles);
    (*metrics)["thermal_spinup_end_tdb_s"]=Number(spinup.cycleEndEphemerisTimeS);
    (*metrics)["maximum_production_update_temperature_change_k"]=Number(maximumChange);
    (*metrics)["deterministic_thermal_worker_ranks"]="1";
    return true;
  } catch (const std::bad_alloc&) {
    *error="I20 thermal export exhausted memory while retaining query/production state; "
        "reduce the explicitly requested query count or thermal/surface resolution, "
        "or request a larger allocation; no partial result is qualified";
  } catch (const std::exception& failure) {
    *error=std::string("I20 thermal export infrastructure failure: ")+failure.what();
  }
  return false;
#endif
}

}  // namespace PolarDiagnostics
}  // namespace Moon
