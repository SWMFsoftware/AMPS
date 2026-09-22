#ifndef _SRC_EARTH_UTIL_STANDALONE_PRODUCT_CONTRACT_H_
#define _SRC_EARTH_UTIL_STANDALONE_PRODUCT_CONTRACT_H_

//======================================================================================
// StandaloneProductContract.h -- Roadmap Step 7
//======================================================================================
// A dependency-free, fail-fast run plan shared by the direct/gridless and materialized
// Mode3D standalone paths.  This contract keeps product selection, released field-model
// aliases, output-domain support, and snapshot-time coherence out of backend dispatch
// if-chains.  It does not contain SWMF coupling behavior (roadmap Steps 9--11).
//======================================================================================

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>

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

inline std::string CanonicalFieldModel(const std::string& input) {
  const std::string model=UpperTrim(input);
  if (model=="TS96" || model=="T96S") return "T96";
  if (model=="TS01" || model=="T01S") return "T01";
  if (model=="TS05" || model=="T05S" || model=="T04S" || model=="TS04") return "T05";
  if (model=="TA16RBF") return "TA16";
  return model;
}

inline bool IsReleasedFieldModel(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  return model=="DIPOLE" || model=="IGRF" || model=="T96" || model=="T01" ||
         model=="T05" || model=="TA15N" || model=="TA15B" || model=="TA16";
}

inline bool IsExternalFieldModel(const std::string& input) {
  const std::string model=CanonicalFieldModel(input);
  return model=="T96" || model=="T01" || model=="T05" || model=="TA15N" ||
         model=="TA15B" || model=="TA16";
}

struct ProductSelection {
  bool cutoff{false};
  bool fluxSpectrum{false};
};

inline ProductSelection ParseProductSelection(const std::string& target) {
  const std::string normalized=UpperTrim(target);
  ProductSelection selection;
  selection.cutoff=normalized.find("CUTOFF")!=std::string::npos ||
                   normalized=="ALL" || normalized=="BOTH";
  selection.fluxSpectrum=normalized.find("DENSITY")!=std::string::npos ||
                         normalized.find("FLUX")!=std::string::npos ||
                         normalized.find("SPECTRUM")!=std::string::npos ||
                         normalized=="ALL" || normalized=="BOTH";
  if (!selection.cutoff && !selection.fluxSpectrum) {
    throw std::invalid_argument(
        "CALC_TARGET must request CUTOFF_RIGIDITY, DENSITY_SPECTRUM/FLUX, BOTH, or ALL");
  }
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
    if (!values[i]->empty() && *values[i]!=epochs.field) {
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
  bool driverColumnsValidated{false};
  bool driverUnitsValidated{false};
  bool fieldValidityValidated{false};

  void Validate(bool externalDriversAreInline=false) {
    fieldModel=CanonicalFieldModel(fieldModel);
    outputMode=CanonicalOutputMode(outputMode);
    if (!IsReleasedFieldModel(fieldModel))
      throw std::invalid_argument(
          "Unsupported standalone FIELD_MODEL; released models are "
          "DIPOLE, IGRF, T96, T01, T05/TS05, TA15N, TA15B, and TA16");
    if (!products.cutoff && !products.fluxSpectrum)
      throw std::invalid_argument("Standalone run plan requests no product");
    RequireOneEpoch(epochs);
    if (IsExternalFieldModel(fieldModel) && !externalDriversAreInline &&
        (!driverColumnsValidated || !driverUnitsValidated))
      throw std::runtime_error(
          "External field driver columns/units were not validated before snapshot creation");
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
    if (c=='\\' || c=='\"') out.push_back('\\');
    if (c=='\n') out+="\\n";
    else out.push_back(c);
  }
  return out;
}

inline std::string BuildManifestJson(const RunPlan& plan) {
  std::ostringstream out;
  out << "{\n"
      << "  \"schema\": \"sep-in-geospace/standalone-products/v1\",\n"
      << "  \"field_model\": \"" << JsonEscape(plan.fieldModel) << "\",\n"
      << "  \"field_representation\": \"" << RepresentationName(plan.representation) << "\",\n"
      << "  \"epoch_utc\": \"" << JsonEscape(plan.epochs.field) << "\",\n"
      << "  \"output_mode\": \"" << JsonEscape(plan.outputMode) << "\",\n"
      << "  \"cutoff\": " << (plan.products.cutoff ? "true" : "false") << ",\n"
      << "  \"flux_spectrum\": " << (plan.products.fluxSpectrum ? "true" : "false") << ",\n"
      << "  \"snapshot_epoch_coherent\": true,\n"
      << "  \"driver_columns_validated\": "
      << (plan.driverColumnsValidated ? "true" : "false") << ",\n"
      << "  \"driver_units_validated\": "
      << (plan.driverUnitsValidated ? "true" : "false") << ",\n"
      << "  \"field_validity_validated\": "
      << (plan.fieldValidityValidated ? "true" : "false") << "\n"
      << "}\n";
  return out.str();
}

} // namespace StandaloneProducts
} // namespace Earth

#endif
