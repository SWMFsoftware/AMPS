#ifndef _SRC_EARTH_UTIL_STANDALONE_PRODUCT_CONTRACT_H_
#define _SRC_EARTH_UTIL_STANDALONE_PRODUCT_CONTRACT_H_

//======================================================================================
// StandaloneProductContract.h -- Roadmap Step 7
//======================================================================================
// This header is the dependency-free startup contract shared by the GRIDLESS and
// Mode3D standalone paths.  Step 7 deliberately treats those paths as two field
// representations of one calculation: model aliases, requested products, output
// domains, driver requirements, and epoch coherence are decided here, before either
// backend starts trajectories.
//
// The file intentionally contains no PIC, MPI, SPICE, Geopack, or SWMF types.  Besides
// making the rules unit-testable, that separation prevents standalone empirical-field
// selection from leaking into the SWMF component.  Live SWMF ingestion starts at
// Roadmap Step 9 and is not represented by this contract.
//======================================================================================

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Earth {
namespace StandaloneProducts {

inline std::string UpperTrim(std::string value) {
  while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front())))
    value.erase(value.begin());
  while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back())))
    value.pop_back();
  for (char& c:value) c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return value;
}

// Map historical CCMC, ViRBO, and Fortran entry-point spellings to one public name.
// All dispatch and provenance records use the returned spelling.
inline std::string CanonicalFieldModel(const std::string& input) {
  const std::string model=UpperTrim(input);
  if (model=="TS96" || model=="T96S") return "T96";
  if (model=="TS01" || model=="T01S") return "T01";
  if (model=="TS05" || model=="T05S" || model=="T04S" || model=="TS04")
    return "T05";
  if (model=="TA16RBF") return "TA16";
  return model;
}

// Production standalone models named by the Step-7 release plan.
inline bool IsReleasedFieldModel(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  return model=="DIPOLE" || model=="IGRF" || model=="T96" || model=="T01" ||
         model=="T05" || model=="TA15N" || model=="TA15B" || model=="TA16";
}

// NONE is retained only as the exact zero-field reference used by F1/F2/F12/F15/F16
// and I-F01/I-F02.  Keeping it explicit avoids weakening those analytic validation
// tests while preventing manifests from presenting NONE as a physical field model.
inline bool IsValidationOnlyFieldModel(const std::string& input) {
  return CanonicalFieldModel(input)=="NONE";
}

inline bool IsSupportedFieldModel(const std::string& input) {
  return IsReleasedFieldModel(input) || IsValidationOnlyFieldModel(input);
}

inline bool IsExternalFieldModel(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  return model=="T96" || model=="T01" || model=="T05" || model=="TA15N" ||
         model=="TA15B" || model=="TA16";
}

inline bool RequiresGeopackInitialization(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  return model!="DIPOLE" && model!="NONE";
}

// The returned names are the canonical quantities consumed by each wrapper.  Do not
// require unrelated solar-wind columns merely because a particular source table happens
// to contain them: accepting a table with the right columns is different from silently
// inventing a missing physical driver.
inline std::vector<std::string> RequiredDriverColumns(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  if (model=="DIPOLE" || model=="IGRF" || model=="NONE") return {};
  if (model=="T96") return {"BYIMF","BZIMF","PDYN","DST"};
  if (model=="T01") return {"BYIMF","BZIMF","PDYN","DST","G1","G2","G3"};
  if (model=="T05")
    return {"BYIMF","BZIMF","PDYN","DST","W1","W2","W3","W4","W5","W6"};
  if (model=="TA15N" || model=="TA15B")
    return {"BYIMF","BZIMF","PDYN","XIND"};
  if (model=="TA16") return {"BYIMF","PDYN","DST"};
  throw std::invalid_argument("No standalone driver contract for FIELD_MODEL='"+
                              input+"'");
}

// Canonical native units expected by the AMPS-side PARMOD construction.  Empty means
// that the quantity is a model-defined coefficient for which this layer performs no
// dimensional conversion (for example W1..W6).  Explicit units on such coefficients
// remain in provenance, but are not guessed here.
inline std::string ExpectedDriverUnit(const std::string& columnInput) {
  const std::string column=UpperTrim(columnInput);
  if (column=="BYIMF" || column=="BZIMF" || column=="DST" || column=="SYMHC" ||
      (column.size()==3 && column.substr(0,2)=="BZ" &&
       column[2]>='1' && column[2]<='6')) return "NT";
  if (column=="VSW") return "KM/S";
  if (column=="DEN_P") return "CM^-3";
  if (column=="PDYN") return "NPA";
  if (column=="G1" || column=="G2" || column=="G3" || column=="XIND") return "1";
  return "";
}

inline std::string CanonicalDriverUnit(std::string unit) {
  unit=UpperTrim(unit);
  std::string compact;
  for (char c:unit) {
    if (!std::isspace(static_cast<unsigned char>(c)) && c!='{' && c!='}')
      compact.push_back(c);
  }
  if (compact=="NANOTESLA" || compact=="NANOTESLAS") return "NT";
  if (compact=="NANOPASCAL" || compact=="NANOPASCALS") return "NPA";
  if (compact=="KM/SEC" || compact=="KMS^-1" || compact=="KM*S^-1") return "KM/S";
  if (compact=="CM-3" || compact=="CM**-3" || compact=="1/CM3") return "CM^-3";
  if (compact=="DIMENSIONLESS" || compact=="UNITLESS" || compact=="NONE") return "1";
  return compact;
}

inline bool DriverUnitMatches(const std::string& column,const std::string& unit) {
  const std::string expected=ExpectedDriverUnit(column);
  return expected.empty() || CanonicalDriverUnit(unit)==expected;
}

// Parse a driver cell without the permissive std::stod prefix behavior.  A token such
// as "2.0bad", NaN, or infinity must never reach PARMOD as a plausible physical zero
// or finite driver.  Optional, model-unused columns remain the loader's concern; this
// helper is applied to every column returned by RequiredDriverColumns().
inline double ParseFiniteDriverValue(const std::string& column,
                                     const std::string& token) {
  std::size_t used=0;
  double value=0.0;
  try {
    value=std::stod(token,&used);
  }
  catch (const std::exception&) {
    throw std::invalid_argument("Driver column '"+column+
                                "' has nonnumeric value '"+token+"'");
  }
  if (used!=token.size() || !std::isfinite(value))
    throw std::invalid_argument("Driver column '"+column+
                                "' has invalid finite value '"+token+"'");
  return value;
}

// Validate only the columns actually consumed by the selected model.  A legacy simple
// header without bracketed units is allowed only when the caller explicitly declares
// that it is using the fixed AMPS-wizard schema; an explicit unit is never ignored.
inline void ValidateDriverUnits(const std::string& model,
                                const std::map<std::string,std::string>& units,
                                bool fixedLegacySchema) {
  for (const std::string& column:RequiredDriverColumns(model)) {
    const std::string expected=ExpectedDriverUnit(column);
    if (expected.empty()) continue;
    const auto it=units.find(column);
    if (it==units.end() || UpperTrim(it->second).empty()) {
      if (fixedLegacySchema) continue;
      throw std::invalid_argument("Driver column '"+column+
          "' has no declared unit; expected "+expected);
    }
    if (!DriverUnitMatches(column,it->second)) {
      throw std::invalid_argument("Driver column '"+column+"' has unit '"+
          it->second+"'; expected "+expected);
    }
  }
}

struct ProductSelection {
  bool cutoff{false};
  bool fluxSpectrum{false};
};

inline ProductSelection ParseProductSelection(const std::string& target) {
  std::string normalized=UpperTrim(target);
  for (char& c:normalized)
    if (c=='+' || c==',' || c=='|' || c==';') c=' ';

  ProductSelection selection;
  std::istringstream input(normalized);
  std::string token;
  while (input >> token) {
    if (token=="ALL" || token=="BOTH") {
      selection.cutoff=true;
      selection.fluxSpectrum=true;
    }
    else if (token=="CUTOFF" || token=="CUTOFF_RIGIDITY") selection.cutoff=true;
    else if (token=="DENSITY" || token=="DENSITY_SPECTRUM" || token=="FLUX" ||
             token=="FLUX_SPECTRUM" || token=="SPECTRUM")
      selection.fluxSpectrum=true;
    else
      throw std::invalid_argument("Unknown CALC_TARGET component '"+token+"'");
  }
  if (!selection.cutoff && !selection.fluxSpectrum)
    throw std::invalid_argument(
        "CALC_TARGET must request CUTOFF_RIGIDITY, DENSITY_SPECTRUM/FLUX, BOTH, or ALL");
  return selection;
}

inline std::string CanonicalOutputMode(const std::string& mode) {
  const std::string result=UpperTrim(mode);
  if (result!="POINTS" && result!="TRAJECTORY" && result!="SHELLS")
    throw std::invalid_argument("OUTPUT_MODE must be POINTS, TRAJECTORY, or SHELLS");
  return result;
}

enum class FieldRepresentation { Gridless, Mesh };

struct SnapshotEpochs {
  std::string field;
  std::string drivers;
  std::string boundarySpectrum;
  std::string ephemeris;
  std::string output;
};

inline void RequireOneEpoch(const SnapshotEpochs& epochs) {
  if (epochs.field.empty())
    throw std::invalid_argument("Standalone snapshot requires a non-empty field epoch");
  const std::string* values[]={&epochs.drivers,&epochs.boundarySpectrum,
                              &epochs.ephemeris,&epochs.output};
  const char* labels[]={"driver","boundary spectrum","ephemeris","output"};
  for (int i=0;i<4;++i) {
    if (values[i]->empty() || *values[i]!=epochs.field) {
      std::ostringstream message;
      message << "Standalone snapshot epoch mismatch: field='" << epochs.field
              << "' but " << labels[i] << "='" << *values[i] << "'";
      throw std::runtime_error(message.str());
    }
  }
}

struct RunPlan {
  std::string fieldModel;
  std::string outputMode;
  ProductSelection products;
  FieldRepresentation representation{FieldRepresentation::Gridless};
  SnapshotEpochs epochs;
  std::string snapshotId;
  std::string driverSource{"INLINE"};
  bool driverColumnsValidated{false};
  bool driverUnitsValidated{false};
  bool driverEpochValidated{false};
  bool geopackInitialized{false};
  bool fieldValidityValidated{false};

  void Validate(bool externalDriversAreInline=false) {
    fieldModel=CanonicalFieldModel(fieldModel);
    outputMode=CanonicalOutputMode(outputMode);
    if (!IsSupportedFieldModel(fieldModel))
      throw std::invalid_argument(
          "Unsupported standalone FIELD_MODEL; supported production models are "
          "DIPOLE, IGRF, T96, T01, T05/TS05, TA15N, TA15B, and TA16 "
          "(NONE is reserved for analytic validation)");
    if (!products.cutoff && !products.fluxSpectrum)
      throw std::invalid_argument("Standalone run plan requests no product");
    RequireOneEpoch(epochs);
    if (snapshotId.empty())
      throw std::runtime_error("Standalone field snapshot has no deterministic identity");
    if (IsExternalFieldModel(fieldModel) && !externalDriversAreInline &&
        (!driverColumnsValidated || !driverUnitsValidated))
      throw std::runtime_error(
          "External field driver columns/units were not validated before snapshot creation");
    if (IsExternalFieldModel(fieldModel) && !driverEpochValidated)
      throw std::runtime_error(
          "External field driver table does not cover the authoritative snapshot epoch");
    if (RequiresGeopackInitialization(fieldModel) && !geopackInitialized)
      throw std::runtime_error("Geopack RECALC/IGRF state was not initialized");
    if (!fieldValidityValidated)
      throw std::runtime_error("Field validity/domain contract was not validated");
  }
};

inline const char* RepresentationName(FieldRepresentation value) {
  return value==FieldRepresentation::Gridless ? "GRIDLESS" : "MESH";
}

inline std::string JsonEscape(const std::string& value) {
  std::string out;
  for (char c:value) {
    if (c=='\\') out+="\\\\";
    else if (c=='\"') out+="\\\"";
    else if (c=='\n') out+="\\n";
    else if (c=='\r') out+="\\r";
    else if (c=='\t') out+="\\t";
    else out.push_back(c);
  }
  return out;
}

inline std::string BuildManifestJson(const RunPlan& plan) {
  std::ostringstream out;
  out << "{\n"
      << "  \"schema\": \"sep-in-geospace/standalone-products/v1\",\n"
      << "  \"field_model\": \"" << JsonEscape(plan.fieldModel) << "\",\n"
      << "  \"validation_only_field\": "
      << (IsValidationOnlyFieldModel(plan.fieldModel) ? "true" : "false") << ",\n"
      << "  \"field_representation\": \"" << RepresentationName(plan.representation) << "\",\n"
      << "  \"snapshot_id\": \"" << JsonEscape(plan.snapshotId) << "\",\n"
      << "  \"epoch_utc\": \"" << JsonEscape(plan.epochs.field) << "\",\n"
      << "  \"output_mode\": \"" << JsonEscape(plan.outputMode) << "\",\n"
      << "  \"cutoff\": " << (plan.products.cutoff ? "true" : "false") << ",\n"
      << "  \"flux_spectrum\": " << (plan.products.fluxSpectrum ? "true" : "false") << ",\n"
      << "  \"driver_source\": \"" << JsonEscape(plan.driverSource) << "\",\n"
      << "  \"snapshot_epoch_coherent\": true,\n"
      << "  \"driver_columns_validated\": "
      << (plan.driverColumnsValidated ? "true" : "false") << ",\n"
      << "  \"driver_units_validated\": "
      << (plan.driverUnitsValidated ? "true" : "false") << ",\n"
      << "  \"driver_epoch_validated\": "
      << (plan.driverEpochValidated ? "true" : "false") << ",\n"
      << "  \"geopack_initialized\": "
      << (plan.geopackInitialized ? "true" : "false") << ",\n"
      << "  \"field_validity_validated\": "
      << (plan.fieldValidityValidated ? "true" : "false") << "\n"
      << "}\n";
  return out.str();
}

} // namespace StandaloneProducts
} // namespace Earth

#endif
