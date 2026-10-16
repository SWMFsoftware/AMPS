/*
 * Dependency-light Stage-1 unit tests.
 *
 * These tests exercise the parser, material resolver, spectrum adapter,
 * canonical hashes, layer conversion, and adaptive-accuracy mathematics
 * without running particle transport.  Failures print the exact contract that
 * was violated and return nonzero for CTest/CI.
 */

#include "AccuracyController.hh"
#include "CaseDefinition.hh"
#include "ContentHash.hh"
#include "MaterialLibrary.hh"
#include "SpectrumDataset.hh"
#include "StrictJson.hh"

#include <cmath>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int failures=0;

void Check(bool condition,const std::string& description){
  if(!condition){
    ++failures;
    std::cerr<<"FAIL: "<<description<<'\n';
  }
}

template<class Function>
void CheckThrows(Function function,const std::string& description){
  try { function(); }
  catch(const std::exception&) { return; }
  ++failures;
  std::cerr<<"FAIL: expected rejection: "<<description<<'\n';
}

bool Near(double a,double b,double tolerance=1.0e-12){
  return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

} // namespace

int main(int argc,char** argv){
  if(argc>2){
    std::cerr<<"usage: stage1_core_tests [validation-evidence-directory]\n";
    return 2;
  }
  const std::string root=SHIELDSIM_SOURCE_DIR;

  Check(shield_sim::Sha256("abc")==
        "sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "SHA-256 matches the published 'abc' test vector");
  CheckThrows([]{ shield_sim::Json::Parse("{\"a\":1,\"a\":2}","duplicate-test"); },
              "duplicate JSON keys");
  CheckThrows([]{ shield_sim::Json::Parse("NaN","nan-test"); },"non-JSON NaN token");

  const auto library=shield_sim::MaterialLibrary::Load(root+"/data/material_library.v1.json");
  const auto& al=library.Resolve("Al","shield");
  Check(Near(al.densityGPerCm3,2.699),"Al density is the explicit library value");
  Check(library.Resolve("SiO2","absorber").densityGPerCm3==2.32,
        "SiO2 ambiguity is repaired to the 2.32 g/cm3 reference");
  Check(library.Resolve("InGaAs","absorber").densityGPerCm3==5.50,
        "Ga0.47In0.53As uses the approved 5.50 g/cm3 density");
  CheckThrows([&]{ (void)library.Resolve("CFRP","shield"); },"disabled CFRP template cannot transport");
  const auto hdpe=library.ElementalMassFractions(library.Resolve("HDPE","shield"));
  Check(Near(hdpe.at("C")+hdpe.at("H"),1.0),"atom-count recipe converts to normalized mass fractions");

  const shield_sim::Json override=shield_sim::Json::Parse(R"json(
    {"material_id":"Al_project_density","revision":"case-1","extends":"Al_natural_G4_11_4@1",
     "status":"project_specific",
     "density":{"value":2.70,"unit":"g/cm3","uncertainty":{"type":"standard","value":0.01,"unit":"g/cm3"},"reference_condition":"project measurement"},
     "provenance":[{"properties":["density"],"source":"project density record","version":"1"}]}
  )json");
  const auto overridden=shield_sim::MaterialLibrary::Load(
      root+"/data/material_library.v1.json",shield_sim::Json::Array{override});
  Check(Near(overridden.Resolve("Al_project_density","shield").densityGPerCm3,2.70),
        "case material density override wins exactly");
  Check(Near(overridden.Resolve("Al","shield").densityGPerCm3,2.699),
        "case override does not mutate its parent");
  Check(overridden.Resolve("Al_project_density","shield").contentHash!=al.contentHash,
        "transport-property override changes content hash");

  const shield_sim::Json isotopeOverride=shield_sim::Json::Parse(R"json(
    {"material_id":"B_enriched_test","revision":"case-1","extends":"B_natural_G4_11_4@1",
     "status":"project_specific",
     "composition":{"basis":"elemental_mass_fraction","isotopes":"explicit","residual_policy":"exact",
       "components":[{"isotope":"B-10","value":0.95},{"isotope":"B-11","value":0.05}]},
     "provenance":[{"properties":["composition","isotopes"],"source":"synthetic isotope test record","version":"1"}]}
  )json");
  const auto isotopeLibrary=shield_sim::MaterialLibrary::Load(
      root+"/data/material_library.v1.json",shield_sim::Json::Array{isotopeOverride});
  const auto isotopeFractions=isotopeLibrary.ElementalMassFractions(
      isotopeLibrary.Resolve("B_enriched_test","shield"));
  Check(Near(isotopeFractions.at("isotope:B-10"),0.95) &&
        Near(isotopeFractions.at("isotope:B-11"),0.05),
        "case override can explicitly replace natural isotope abundance");

  const auto spectrum=shield_sim::SpectrumDataset::LoadLegacyThreeColumn(
      root+"/tests/data/mono_100MeV_proton.dat",
      root+"/tests/data/mono_spectrum_metadata.v1.json");
  Check(spectrum.energyMeV.size()==1 && spectrum.binWidthMeV[0]==1.0,
        "single-bin width comes from mandatory metadata");
  CheckThrows([&]{ shield_sim::SpectrumDataset::LoadLegacyThreeColumn(
      root+"/tests/data/mono_100MeV_proton.dat",""); },
      "legacy spectrum without metadata");

  const auto example=shield_sim::CaseDefinition::Load(root+"/examples/stage1_case.json");
  Check(example->shieldLayers.size()==2,"ordered multilayer case retains both layers");
  Check(Near(example->shieldLayers[0].thicknessMm,1.0),
        "Al areal density converts analytically to 1 mm");
  Check(Near(example->ShieldTotalThicknessMm(),3.0),"resolved shield thickness sums layers");
  Check(Near(example->accuracy.global.target,0.10) && example->accuracy.targetWasDefaulted,
        "omitted accuracy resolves and archives exact 10% default");
  Check(example->absorbers[0].geometryHash!=example->absorbers[1].geometryHash,
        "absorber identity participates in geometry content hash");
  const std::string serialized=example->resolvedJson.DumpCanonical();
  Check(shield_sim::Json::Parse(serialized,"resolved-roundtrip").DumpCanonical()==serialized,
        "resolved case canonical serialization round-trips exactly");

  shield_sim::AccuracyPlan statistical;
  statistical.global.scope="statistical";
  statistical.global.target=0.10;
  statistical.minHistories=4;
  statistical.maxHistories=1000;
  statistical.batchSize=4;
  shield_sim::AccuracyController controller(statistical);
  for(int i=0;i<100;++i) controller.AddSample("TID",i%2?1.1:0.9);
  const auto assessment100=controller.Assess("TID");
  Check(assessment100.status==shield_sim::AccuracyStatus::Met,
        "synthetic estimator meets a 10% statistical target");
  const double u100=assessment100.statisticalStandardUncertainty;
  for(int i=0;i<300;++i) controller.AddSample("TID",i%2?1.1:0.9);
  const double u400=controller.Assess("TID").statisticalStandardUncertainty;
  Check(u400<u100 && Near(u400/u100,0.5,0.03),
        "standard uncertainty follows expected 1/sqrt(N) behavior");

  shield_sim::AccuracyPlan combined=statistical;
  combined.global.scope="combined";
  shield_sim::AccuracyController incomplete(combined);
  for(int i=0;i<100;++i) incomplete.AddSample("TID",i%2?1.1:0.9);
  incomplete.MarkSystematicBudgetIncomplete("TID");
  Check(incomplete.Assess("TID").status==shield_sim::AccuracyStatus::NotAssessable,
        "combined request with incomplete systematics is not_assessable");
  incomplete.SetKnownSystematicFloor("TID",0.20);
  Check(incomplete.Assess("TID").status==shield_sim::AccuracyStatus::NotMet,
        "known systematic floor above target is reproducibly not_met");

  if(failures){
    std::cerr<<failures<<" Stage-1 core test(s) failed\n";
    return 1;
  }

  // Write machine-readable values only after all assertions pass.  The
  // dependent figure test visualizes these exact production-test values; it
  // never substitutes image comparison for the numerical acceptance checks.
  if(argc==2){
    namespace fs=std::filesystem;
    shield_sim::Json::Array comparisons;
    auto add=[&](const std::string& id,double expected,double actual,
                 const std::string& unit){
      comparisons.emplace_back(shield_sim::Json::Object{{"actual",actual},{"expected",expected},
        {"id",id},{"unit",unit}});
    };
    add("al_areal_to_geometric_thickness",1.0,example->shieldLayers[0].thicknessMm,"mm");
    add("multilayer_total_thickness",3.0,example->ShieldTotalThicknessMm(),"mm");
    add("legacy_single_bin_width",1.0,spectrum.binWidthMeV[0],"MeV");
    add("hdpe_mass_fraction_closure",1.0,hdpe.at("C")+hdpe.at("H"),"fraction");
    add("statistical_uncertainty_ratio_N400_to_N100",0.5,u400/u100,"ratio");
    shield_sim::Json evidence(shield_sim::Json::Object{
      {"comparisons",shield_sim::Json(std::move(comparisons))},
      {"schema_version","shieldSim.validation_evidence.v1"},
      {"test_case","stage1_core"}});
    fs::create_directories(argv[1]);
    std::ofstream output(fs::path(argv[1])/"stage1_core.json",std::ios::binary|std::ios::trunc);
    if(!output){ std::cerr<<"cannot write Stage-1 core validation evidence\n"; return 2; }
    output<<evidence.DumpCanonical()<<'\n';
    if(!output){ std::cerr<<"failed writing Stage-1 core validation evidence\n"; return 2; }
  }
  std::cout<<"All Stage-1 core tests passed\n";
  return 0;
}
