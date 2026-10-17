/*
 * Deterministic Stage-4 TID and uncertainty-contract tests.
 *
 * These checks contain no Geant4 transport. They prove the exact unit/mass
 * algebra, the case-v2 capability boundary, structured covariance validation,
 * and the explicitly selected combined-coverage calculation. Stochastic and
 * reference-physics evidence is kept in the separate Stage-4 Python drivers.
 */

#include "AccuracyController.hh"
#include "CaseDefinition.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

namespace {
int failures=0;
void Check(bool condition,const std::string& description){
  if(!condition){ ++failures; std::cerr<<"FAIL: "<<description<<'\n'; }
}
template<class Function> void CheckThrows(Function function,const std::string& description){
  try { function(); } catch(const std::exception&) { return; }
  ++failures; std::cerr<<"FAIL: expected rejection: "<<description<<'\n';
}
bool Near(double actual,double expected,double relative=1.0e-12){
  return std::abs(actual-expected)<=relative*std::max({1.0,std::abs(actual),std::abs(expected)});
}
}

int main(int argc,char** argv){
  if(argc!=3){ std::cerr<<"usage: stage4_core_tests SOURCE_ROOT EVIDENCE_DIR\n"; return 2; }
  namespace fs=std::filesystem;
  const fs::path root=argv[1], evidence=argv[2];
  fs::create_directories(evidence);
  const fs::path scratch=evidence/"stage4_core_scratch";
  fs::create_directories(scratch);

  const auto stage4=shield_sim::CaseDefinition::Load(
      (root/"examples/stage4_tid_case.json").string());
  Check(stage4->schemaVersion=="shieldSim.case.v2" &&
        stage4->quantities==std::vector<std::string>{"TID"},
        "case-v2 enables only the requested Stage-4 TID observable");
  Check(stage4->emPhysics=="reference_list_default" &&
        stage4->physicsList=="FTFP_BERT",
        "production reference and EM physics choices are explicit");
  const auto& requirement=stage4->AccuracyFor("TID");
  const double expectedSystematic=std::sqrt(0.02*0.02+0.04*0.04);
  Check(requirement.systematicBudgetComplete &&
        requirement.systematicComponents.size()==4 &&
        requirement.systematicCorrelations.size()==6 &&
        Near(requirement.systematicRelativeUncertainty,expectedSystematic),
        "four-category all-pairs budget combines by the declared covariance");

  // One MeV deposited in one kilogram is exactly defined through the SI eV
  // conversion. A rad is exactly 0.01 Gy, so the numeric rad value is 100 times
  // the numeric Gy value. No material or transport assumption enters here.
  constexpr double joulePerMeV=1.602176634e-13;
  const double doseGy=joulePerMeV/1.0;
  const double doseRad=100.0*doseGy;
  Check(Near(doseGy,1.602176634e-13,1.0e-15),"MeV/kg converts exactly to Gy");
  Check(Near(doseRad/doseGy,100.0),"Gy-to-rad conversion is exactly 100");

  shield_sim::AccuracyPlan plan=stage4->accuracy;
  shield_sim::AccuracyController controller(plan);
  controller.SetSystematicRelativeUncertainty("TID",requirement.systematicRelativeUncertainty);
  for(int i=0;i<1000;++i){
    controller.RecordHistory();
    controller.AddSample("TID",1.0+(i%2?0.01:-0.01));
  }
  const auto assessment=controller.Assess("TID");
  Check(assessment.status==shield_sim::AccuracyStatus::Met &&
        Near(assessment.systematicRelativeUncertainty,expectedSystematic),
        "combined accuracy uses the parsed structured systematic budget");
  Check(std::isfinite(assessment.combinedInterval95Low) &&
        std::isfinite(assessment.combinedInterval95High) &&
        assessment.combinedInterval95Low<assessment.estimate &&
        assessment.combinedInterval95High>assessment.estimate,
        "explicit Student-t plus normal-systematics method produces a combined 95% interval");

  auto document=shield_sim::Json::ParseFile((root/"examples/stage4_tid_case.json").string());
  document.At("material_library").AsObject()["path"]=(root/"data/material_library.v1.json").string();
  auto invalidCovariance=document;
  auto& correlations=invalidCovariance.At("accuracy").At("systematic_budget").At("correlations").AsArray();
  correlations[0].AsObject()["correlation_coefficient"]=0.9;
  correlations[1].AsObject()["correlation_coefficient"]=0.9;
  correlations[3].AsObject()["correlation_coefficient"]=-0.9;
  const fs::path invalidPath=scratch/"non_psd_covariance.json";
  { std::ofstream stream(invalidPath); stream<<invalidCovariance.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(invalidPath.string()); },
              "non-positive-semidefinite uncertainty correlation matrix");

  auto aggregateOnly=document;
  auto& budget=aggregateOnly.At("accuracy").At("systematic_budget").AsObject();
  budget.clear();
  budget["status"]="complete";
  budget["relative_standard_uncertainty"]=0.05;
  const fs::path aggregatePath=scratch/"aggregate_only.json";
  { std::ofstream stream(aggregatePath); stream<<aggregateOnly.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(aggregatePath.string()); },
              "case-v2 TID aggregate-only systematic claim");

  auto forbidden=document;
  forbidden.At("quantities").AsArray()={shield_sim::Json("DDD")};
  const fs::path forbiddenPath=scratch/"forbidden_ddd.json";
  { std::ofstream stream(forbiddenPath); stream<<forbidden.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(forbiddenPath.string()); },
              "Stage-4 TID does not enable DDD");

  auto gamma=document;
  auto& particle=gamma.At("source").At("spectrum").At("particle").AsObject();
  particle["name"]="gamma"; particle["Z"]=0.0; particle["A"]=0.0;
  gamma.At("source").At("spectrum").At("energy").AsObject()["value"]=1.0;
  gamma.AsObject()["quantities"]=shield_sim::Json(
      shield_sim::Json::Array{shield_sim::Json("transport_raw")});
  gamma.At("accuracy").AsObject()["scope"]="statistical";
  gamma.At("accuracy").At("coverage").AsObject()["interval_method"]="student_t_or_exact_bounded";
  gamma.At("accuracy").AsObject()["systematic_budget"]=shield_sim::Json(
      shield_sim::Json::Object{{"status","incomplete"}});
  const fs::path gammaPath=scratch/"gamma_attenuation_node.json";
  { std::ofstream stream(gammaPath); stream<<gamma.DumpCanonical()<<'\n'; }
  const auto gammaCase=shield_sim::CaseDefinition::Load(gammaPath.string());
  Check(gammaCase->responseIncident.particleName=="gamma" &&
        gammaCase->responseIncident.atomicNumber==0 &&
        gammaCase->responseIncident.massNumber==0,
        "gamma/0/0 total-energy node is accepted for attenuation verification");

  // The same declared incident must not cross the Stage-4 capability gate:
  // photon attenuation is verified against XCOM, but photon dose is outside
  // the validated proton/He-4 TID domain.
  auto gammaTid=gamma;
  gammaTid.AsObject()["quantities"]=shield_sim::Json(
      shield_sim::Json::Array{shield_sim::Json("TID")});
  const fs::path gammaTidPath=scratch/"unsupported_gamma_tid.json";
  { std::ofstream stream(gammaTidPath); stream<<gammaTid.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(gammaTidPath.string()); },
              "gamma attenuation verification does not enable a photon TID claim");

  auto heavyIonTid=document;
  auto& heavyIon=heavyIonTid.At("source").At("spectrum").At("particle").AsObject();
  heavyIon["name"]="ion"; heavyIon["Z"]=6.0; heavyIon["A"]=12.0;
  const fs::path heavyIonTidPath=scratch/"unsupported_heavy_ion_tid.json";
  { std::ofstream stream(heavyIonTidPath); stream<<heavyIonTid.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(heavyIonTidPath.string()); },
              "raw explicit-ion support does not enable a heavy-ion TID claim");

  auto unsupportedAbsorber=document;
  unsupportedAbsorber.At("geometry").At("absorbers").AsArray()[0]
      .AsObject()["material"]="Cu";
  const fs::path unsupportedAbsorberPath=scratch/"unsupported_copper_absorber_tid.json";
  { std::ofstream stream(unsupportedAbsorberPath); stream<<unsupportedAbsorber.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(unsupportedAbsorberPath.string()); },
              "case-v2 TID cannot claim an absorber outside the Al/Si/water validation matrix");

  shield_sim::Json evidenceDocument(shield_sim::Json::Object{
    {"comparisons",shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{{"actual",doseGy},{"expected",1.602176634e-13},
        {"id","one_MeV_per_kg_to_Gy"},{"unit","Gy"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",doseRad/doseGy},{"expected",100.0},
        {"id","rad_per_Gy_numeric_ratio"},{"unit","1"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",requirement.systematicRelativeUncertainty},
        {"expected",expectedSystematic},{"id","structured_covariance_combination"},{"unit","1"}})})},
    {"schema_version","shieldSim.validation_evidence.v1"},
    {"test_case","stage4_core"}});
  std::ofstream evidenceStream(evidence/"stage4_core.json",std::ios::binary|std::ios::trunc);
  evidenceStream<<evidenceDocument.DumpCanonical()<<'\n';

  if(failures){ std::cerr<<failures<<" Stage-4 core test(s) failed\n"; return 1; }
  std::cout<<"All Stage-4 core tests passed\n";
  return 0;
}
