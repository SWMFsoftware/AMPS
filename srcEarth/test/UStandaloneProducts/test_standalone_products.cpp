#include "../../util/StandaloneProductContract.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace SP=Earth::StandaloneProducts;

namespace {
int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

template<class Callable>
void CheckThrows(Callable callable,const std::string& message) {
  try { callable(); }
  catch (...) { return; }
  Check(false,message);
}

bool Contains(const std::vector<std::string>& values,const std::string& value) {
  for (const std::string& item:values) if (item==value) return true;
  return false;
}

SP::RunPlan ValidPlan(const std::string& model,
                      const SP::ProductSelection& products) {
  SP::RunPlan plan;
  plan.fieldModel=model;
  plan.outputMode="POINTS";
  plan.products=products;
  plan.representation=SP::FieldRepresentation::Gridless;
  plan.epochs.field="2017-09-10T00:00:00";
  plan.epochs.drivers=plan.epochs.field;
  plan.epochs.boundarySpectrum=plan.epochs.field;
  plan.epochs.ephemeris=plan.epochs.field;
  plan.epochs.output=plan.epochs.field;
  plan.snapshotId="field-v1-reference";
  plan.driverSource="driver.dat";
  plan.driverColumnsValidated=true;
  plan.driverUnitsValidated=true;
  plan.driverEpochValidated=true;
  plan.geopackInitialized=true;
  plan.fieldValidityValidated=true;
  return plan;
}
} // namespace

int main() {
  // S7-C01: aliases and released/validation-only roles.  This is an exact table,
  // not a nonempty/smoke assertion.
  Check(SP::CanonicalFieldModel("TS96")=="T96" &&
        SP::CanonicalFieldModel("T01S")=="T01" &&
        SP::CanonicalFieldModel("T04S")=="T05" &&
        SP::CanonicalFieldModel("TA16RBF")=="TA16",
        "documented aliases must canonicalize deterministically");
  const char* released[]={"DIPOLE","IGRF","T96","T01","T05","TA15N","TA15B","TA16"};
  for (const char* model:released)
    Check(SP::IsReleasedFieldModel(model),std::string("released model ")+model);
  Check(SP::IsSupportedFieldModel("NONE") &&
        SP::IsValidationOnlyFieldModel("NONE") &&
        !SP::IsReleasedFieldModel("NONE"),
        "NONE must remain available only for analytic validation");
  Check(!SP::IsSupportedFieldModel("T89"),"unreleased T89 must fail startup");

  // S7-C02: model-to-driver mapping checks the exact quantities consumed by each
  // wrapper.  The TA15 assertion prevents a regression to the unrelated T05 W/BZ
  // history requirements.
  const std::vector<std::string> t01=SP::RequiredDriverColumns("T01");
  Check(t01.size()==7 && Contains(t01,"G1") && Contains(t01,"G3") &&
        !Contains(t01,"W1"),"T01 driver contract");
  const std::vector<std::string> t05=SP::RequiredDriverColumns("TS05");
  Check(t05.size()==10 && Contains(t05,"W1") && Contains(t05,"W6"),
        "T05 driver contract");
  const std::vector<std::string> ta15=SP::RequiredDriverColumns("TA15N");
  Check(ta15.size()==4 && Contains(ta15,"XIND") && !Contains(ta15,"W1") &&
        !Contains(ta15,"BZ1"),"TA15 four-parameter contract");

  // S7-C03: unit normalization and fail-fast unit validation.  The negative Pa-vs-
  // nPa case is important because accepting it would change pressure by 1e9.
  Check(SP::DriverUnitMatches("BYIMF","nanotesla") &&
        SP::DriverUnitMatches("PDYN","nPa") &&
        SP::DriverUnitMatches("VSW","km/sec") &&
        SP::DriverUnitMatches("DEN_P","cm-3") &&
        !SP::DriverUnitMatches("PDYN","Pa"),
        "native driver units and aliases");
  std::map<std::string,std::string> t96Units{
      {"BYIMF","nT"},{"BZIMF","nT"},{"PDYN","nPa"},{"DST","nT"}};
  SP::ValidateDriverUnits("T96",t96Units,false);
  std::map<std::string,std::string> wrongUnits=t96Units;
  wrongUnits["PDYN"]="Pa";
  CheckThrows([&](){ SP::ValidateDriverUnits("T96",wrongUnits,false); },
              "wrong explicit pressure unit must fail");
  t96Units.erase("DST");
  CheckThrows([&](){ SP::ValidateDriverUnits("T96",t96Units,false); },
              "missing JSON unit must fail");
  SP::ValidateDriverUnits("T96",std::map<std::string,std::string>(),true);
  Check(SP::ParseFiniteDriverValue("PDYN","2.75")==2.75,
        "finite driver value parsing");
  CheckThrows([](){ SP::ParseFiniteDriverValue("PDYN","2.75nPa"); },
              "partially numeric driver token must fail");
  CheckThrows([](){ SP::ParseFiniteDriverValue("PDYN","nan"); },
              "NaN driver value must fail");
  CheckThrows([](){ SP::ParseFiniteDriverValue("PDYN","inf"); },
              "infinite driver value must fail");

  // S7-C04: cutoff-only, flux-only, and combined selection plus all standalone
  // output domains.  Unknown suffixes are rejected instead of substring-matching.
  const SP::ProductSelection cutoff=SP::ParseProductSelection("CUTOFF_RIGIDITY");
  const SP::ProductSelection flux=SP::ParseProductSelection("DENSITY_SPECTRUM");
  const SP::ProductSelection both=SP::ParseProductSelection(
      "CUTOFF_RIGIDITY+DENSITY_SPECTRUM");
  Check(cutoff.cutoff && !cutoff.fluxSpectrum,"cutoff-only selection");
  Check(!flux.cutoff && flux.fluxSpectrum,"flux-only selection");
  Check(both.cutoff && both.fluxSpectrum,"combined selection");
  Check(SP::CanonicalOutputMode("points")=="POINTS" &&
        SP::CanonicalOutputMode("trajectory")=="TRAJECTORY" &&
        SP::CanonicalOutputMode("shells")=="SHELLS","three output domains");
  CheckThrows([](){ SP::ParseProductSelection("CUTOFF_BOGUS"); },
              "unknown product must fail before dispatch");
  CheckThrows([](){ SP::CanonicalOutputMode("VOLUME"); },
              "unsupported output domain must fail before dispatch");

  // S7-C05/I-F07 contract: field, driver, spectrum, ephemeris, and output share
  // exactly one epoch and immutable snapshot identity.
  SP::RunPlan plan=ValidPlan("TS05",both);
  plan.Validate(false);
  Check(plan.fieldModel=="T05","validated plan stores canonical model");
  SP::RunPlan mismatch=plan;
  mismatch.epochs.boundarySpectrum="2017-09-10T00:05:00";
  CheckThrows([&](){ mismatch.Validate(false); },
              "mismatched spectrum epoch must fail");
  SP::RunPlan noIdentity=plan;
  noIdentity.snapshotId.clear();
  CheckThrows([&](){ noIdentity.Validate(false); },
              "missing immutable snapshot identity must fail");
  SP::RunPlan badCoverage=plan;
  badCoverage.driverEpochValidated=false;
  CheckThrows([&](){ badCoverage.Validate(false); },
              "out-of-range driver epoch must fail");
  SP::RunPlan badGeopack=plan;
  badGeopack.geopackInitialized=false;
  CheckThrows([&](){ badGeopack.Validate(false); },
              "missing Geopack initialization must fail");
  SP::RunPlan badSchema=plan;
  badSchema.driverColumnsValidated=false;
  CheckThrows([&](){ badSchema.Validate(false); },
              "unvalidated file-backed driver columns must fail");
  SP::RunPlan badField=plan;
  badField.fieldValidityValidated=false;
  CheckThrows([&](){ badField.Validate(false); },
              "invalid field snapshot must fail");
  SP::RunPlan badModel=plan;
  badModel.fieldModel="T89";
  CheckThrows([&](){ badModel.Validate(false); },
              "unsupported standalone model must fail");
  SP::RunPlan inlineDrivers=plan;
  inlineDrivers.driverColumnsValidated=false;
  inlineDrivers.driverUnitsValidated=false;
  inlineDrivers.Validate(true);

  // S7-C06: exact manifest fields make the run auditable and distinguish the NONE
  // reference backend from released physics.  JSON escaping is checked explicitly.
  plan.driverSource="drivers/quoted_\"name\".dat";
  const std::string manifest=SP::BuildManifestJson(plan);
  Check(manifest.find("\"snapshot_id\": \"field-v1-reference\"")!=std::string::npos &&
        manifest.find("\"cutoff\": true")!=std::string::npos &&
        manifest.find("\"flux_spectrum\": true")!=std::string::npos &&
        manifest.find("quoted_\\\"name\\\"")!=std::string::npos,
        "manifest schema and escaping");
  SP::RunPlan zero=ValidPlan("NONE",flux);
  zero.geopackInitialized=false;
  zero.driverEpochValidated=true;
  zero.Validate(true);
  Check(SP::BuildManifestJson(zero).find("\"validation_only_field\": true")!=
        std::string::npos,"NONE manifest must be marked validation-only");

  if (failures) {
    std::cerr << failures << " Step-7 standalone contract test(s) failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "PASS: Step-7 aliases, driver contracts, units, products, epochs, and manifest\n";
  return EXIT_SUCCESS;
}
