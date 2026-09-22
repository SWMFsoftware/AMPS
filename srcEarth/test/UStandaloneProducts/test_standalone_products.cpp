#include "../../util/StandaloneProductContract.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace SP=Earth::StandaloneProducts;
static int failures=0;

static void Check(bool condition,const std::string& message) {
  if (!condition) { ++failures; std::cerr << "FAIL: " << message << "\n"; }
}

template<class Callable>
static void CheckThrows(Callable callable,const std::string& message) {
  try { callable(); }
  catch (...) { return; }
  Check(false,message);
}

int main() {
  // I-F01/I-F02: every documented alias resolves to one released canonical model.
  Check(SP::CanonicalFieldModel("TS96")=="T96" &&
        SP::CanonicalFieldModel("T01S")=="T01" &&
        SP::CanonicalFieldModel("T04S")=="T05" &&
        SP::CanonicalFieldModel("TA16RBF")=="TA16",
        "field aliases canonicalize deterministically");
  const char* released[]={"DIPOLE","IGRF","T96","T01","T05","TA15N","TA15B","TA16"};
  for (const char* model:released)
    Check(SP::IsReleasedFieldModel(model),std::string("released model ")+model);
  Check(!SP::IsReleasedFieldModel("T89"),"unreleased model rejected");

  // I-F03: compact combined targets and all three standalone output domains.
  const SP::ProductSelection both=SP::ParseProductSelection(
      "CUTOFF_RIGIDITY+DENSITY_SPECTRUM");
  Check(both.cutoff && both.fluxSpectrum,"combined target selects both products");
  Check(SP::CanonicalOutputMode("trajectory")=="TRAJECTORY" &&
        SP::CanonicalOutputMode("shells")=="SHELLS","standalone output modes");
  CheckThrows([](){ SP::ParseProductSelection("DO_NOTHING"); },
              "unknown target fails before execution");

  // I-F04/I-F05: field, driver, spectrum, ephemeris and output all use one epoch.
  SP::RunPlan plan;
  plan.fieldModel="TS05";
  plan.outputMode="POINTS";
  plan.products=both;
  plan.representation=SP::FieldRepresentation::Gridless;
  plan.epochs.field="2017-09-10T00:00:00";
  plan.epochs.drivers=plan.epochs.field;
  plan.epochs.boundarySpectrum=plan.epochs.field;
  plan.epochs.ephemeris=plan.epochs.field;
  plan.epochs.output=plan.epochs.field;
  plan.driverColumnsValidated=true;
  plan.driverUnitsValidated=true;
  plan.fieldValidityValidated=true;
  plan.Validate();
  Check(plan.fieldModel=="T05","validated run plan retains canonical model");
  const std::string manifest=SP::BuildManifestJson(plan);
  Check(manifest.find("\"cutoff\": true")!=std::string::npos &&
        manifest.find("\"flux_spectrum\": true")!=std::string::npos &&
        manifest.find("2017-09-10T00:00:00")!=std::string::npos,
        "manifest records both products and authoritative epoch");

  SP::RunPlan mismatch=plan;
  mismatch.epochs.boundarySpectrum="2017-09-10T00:05:00";
  CheckThrows([&](){ mismatch.Validate(); },"mismatched spectrum epoch fails fast");
  SP::RunPlan badDrivers=plan;
  badDrivers.driverColumnsValidated=false;
  CheckThrows([&](){ badDrivers.Validate(); },"unvalidated external driver mapping fails fast");
  SP::RunPlan dipole=plan;
  dipole.fieldModel="DIPOLE";
  dipole.driverColumnsValidated=false;
  dipole.driverUnitsValidated=false;
  dipole.Validate();
  Check(true,"analytic field does not require external driver table");

  if (failures) {
    std::cerr << failures << " StandaloneProductContract test(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: Step-7 aliases, products, domains, snapshot epochs, and fail-fast gates\n";
  return EXIT_SUCCESS;
}
