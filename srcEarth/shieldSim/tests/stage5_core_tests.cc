/*
 * Deterministic Stage-5 displacement-damage tests.
 *
 * These tests exercise the authenticated table loader, exact reference nodes,
 * log-log interpolation, analytic bin integration, ASTM group intersection,
 * the 95 MeV-mb silicon convention, and strict case capability boundaries.
 * Geant4 boundary scoring is exercised separately by test_stage5_transport.py.
 */

#include "CaseDefinition.hh"
#include "DamageResponseLibrary.hh"
#include "AccuracyController.hh"
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
bool Near(double actual,double expected,double relative=1.0e-11){
  return std::abs(actual-expected)<=relative*std::max({1.0e-300,std::abs(actual),
                                                       std::abs(expected)});
}
}

int main(int argc,char** argv){
  if(argc!=3){ std::cerr<<"usage: stage5_core_tests SOURCE_ROOT EVIDENCE_DIR\n"; return 2; }
  namespace fs=std::filesystem;
  const fs::path root=argv[1], evidence=argv[2];
  fs::create_directories(evidence);
  const fs::path scratch=evidence/"stage5_core_scratch";
  fs::create_directories(scratch);

  const auto library=shield_sim::DamageResponseLibrary::Load(
      (root/"data/niel/damage_response_library.v1.json").string());
  Check(library.LibraryId()=="shieldSim.silicon_displacement_damage" &&
        library.LibraryVersion()=="1.0.0" &&
        library.TargetMaterialId()=="Si_natural_G4_11_4",
        "manifest resolves the approved natural-silicon response identity");

  const auto proton100=library.Evaluate("proton",100.0);
  const auto alpha100=library.Evaluate("alpha",100.0);
  Check(proton100.status==shield_sim::DamageResponseStatus::Available &&
        Near(proton100.nielMeVCm2PerG,2.9757e-3),
        "100 MeV proton exact SR-NIEL node");
  Check(alpha100.status==shield_sim::DamageResponseStatus::Available &&
        Near(alpha100.nielMeVCm2PerG,1.5904e-2),
        "100 MeV He-4 exact SR-NIEL node");

  // Halfway in log(E) must be halfway in log(S), which is an exact geometric
  // mean identity independent of the local power-law slope.
  const double midpoint=std::sqrt(100.0*150.0);
  const auto protonMid=library.Evaluate("proton",midpoint);
  const double expectedMid=std::sqrt(2.9757e-3*2.2803e-3);
  Check(protonMid.status==shield_sim::DamageResponseStatus::Available &&
        Near(protonMid.nielMeVCm2PerG,expectedMid),
        "charged-particle log-log interpolation geometric-mean identity");
  Check(library.Evaluate("proton",1.0e-4).status==
        shield_sim::DamageResponseStatus::OutOfDomain,
        "proton response is not extrapolated below its first positive node");
  Check(library.Evaluate("electron",1.0).status==
        shield_sim::DamageResponseStatus::UnsupportedSpecies,
        "electron DDD remains explicitly unsupported");

  const auto neutron105=library.Evaluate("neutron",1.05);
  const double reference=95.0e-27*6.02214076e23/28.0855;
  Check(neutron105.status==shield_sim::DamageResponseStatus::Available &&
        Near(neutron105.nielMeVCm2PerG,1.6573e-3) &&
        Near(neutron105.neutronHardnessFactor,1.6573e-3/reference,5.0e-12),
        "ASTM [1.0,1.1) MeV neutron group and 95 MeV-mb normalization");
  Check(library.Evaluate("neutron",20.0).status==
        shield_sim::DamageResponseStatus::Available &&
        library.Evaluate("neutron",20.0+1.0e-12).status==
        shield_sim::DamageResponseStatus::OutOfDomain,
        "20 MeV final ASTM edge is inclusive and never extrapolated");

  // A constant differential spectrum over one E722 group has the closed-form
  // fold phi*(E1-E0)*S_group. This directly tests explicit source-bin edges.
  const double phi=10.0;
  const auto neutronFold=library.FoldPiecewiseConstantDifferential("neutron",
      {{1.0,1.1,phi}});
  const double expectedNeutronDdd=phi*0.1*1.6573e-3;
  Check(neutronFold.status==shield_sim::DamageResponseStatus::Available &&
        Near(neutronFold.fluencePerCm2,1.0) &&
        Near(neutronFold.dddMeVPerG,expectedNeutronDdd) &&
        Near(neutronFold.neutronEquivalentFluencePerCm2,
             expectedNeutronDdd/reference,5.0e-12),
        "analytic single-group ASTM fold");

  // Integrate the local proton power law analytically over one complete table
  // interval. This checks integration rather than merely point interpolation.
  const double e0=100.0,e1=150.0,s0=2.9757e-3,s1=2.2803e-3;
  const double slope=std::log(s1/s0)/std::log(e1/e0);
  const double expectedCharged=s0*e0*(std::pow(e1/e0,slope+1.0)-1.0)/(slope+1.0);
  const auto protonFold=library.FoldPiecewiseConstantDifferential("proton",
      {{e0,e1,1.0}});
  Check(protonFold.status==shield_sim::DamageResponseStatus::Available &&
        Near(protonFold.dddMeVPerG,expectedCharged),
        "analytic log-log charged-particle interval fold");
  Check(library.FoldPiecewiseConstantDifferential("neutron",{{19.9,20.1,1.0}}).status==
        shield_sim::DamageResponseStatus::OutOfDomain,
        "nonzero partially out-of-domain bin is not clipped");
  Check(library.FoldPiecewiseConstantDifferential("neutron",{{20.0,21.0,0.0}}).status==
        shield_sim::DamageResponseStatus::Available,
        "zero out-of-domain fluence does not invalidate a fold");

  const auto stage5=shield_sim::CaseDefinition::Load(
      (root/"examples/stage5_ddd_neq_case.json").string());
  Check(stage5->damageResponse &&
        stage5->quantities==std::vector<std::string>({"DDD","n_eq"}) &&
        stage5->responseIncident.particleName=="neutron",
        "case-v2 Stage-5 example resolves DDD, neutron-only n_eq, and its library");

  auto document=shield_sim::Json::ParseFile((root/"examples/stage5_ddd_neq_case.json").string());
  document.At("material_library").AsObject()["path"]=
      (root/"data/material_library.v1.json").string();
  document.At("damage_response").AsObject()["path"]=
      (root/"data/niel/damage_response_library.v1.json").string();
  auto missing=document;
  missing.AsObject().erase("damage_response");
  const fs::path missingPath=scratch/"missing_damage_response.json";
  { std::ofstream out(missingPath); out<<missing.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(missingPath.string()); },
              "DDD/n_eq case without an authenticated response library");

  auto copper=document;
  copper.At("geometry").At("absorbers").AsArray()[0].AsObject()["material"]="Cu";
  const fs::path copperPath=scratch/"copper_damage.json";
  { std::ofstream out(copperPath); out<<copper.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(copperPath.string()); },
              "non-silicon Stage-5 absorber");

  auto electron=document;
  auto& particle=electron.At("source").At("spectrum").At("particle").AsObject();
  particle["name"]="gamma"; particle["Z"]=0.0; particle["A"]=0.0;
  const fs::path electronPath=scratch/"unsupported_incident.json";
  { std::ofstream out(electronPath); out<<electron.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(electronPath.string()); },
              "unsupported Stage-5 incident species");

  // A combined request with no defensible response-function uncertainty stays
  // runnable but must be not assessable. This is distinct from claiming a
  // complete budget while omitting the required physics_reference component,
  // which is a schema error.
  auto combinedIncomplete=document;
  combinedIncomplete.At("accuracy").AsObject()["scope"]="combined";
  combinedIncomplete.At("accuracy").At("coverage").AsObject()["interval_method"]=
      "student_t_mc_plus_normal_systematics_quadrature";
  const fs::path combinedPath=scratch/"combined_incomplete_damage.json";
  { std::ofstream out(combinedPath); out<<combinedIncomplete.DumpCanonical()<<'\n'; }
  const auto combinedCase=shield_sim::CaseDefinition::Load(combinedPath.string());
  shield_sim::AccuracyController damageAccuracy(combinedCase->accuracy);
  damageAccuracy.MarkSystematicBudgetIncomplete("DDD");
  for(int i=0;i<200;++i){
    damageAccuracy.RecordHistory();
    damageAccuracy.AddSample("DDD",1.0+(i%2?0.01:-0.01));
  }
  Check(damageAccuracy.Assess("DDD").status==shield_sim::AccuracyStatus::NotAssessable,
        "combined Stage-5 accuracy is not assessable without physics-reference uncertainty");

  auto missingPhysics=document;
  auto stage4Document=shield_sim::Json::ParseFile(
      (root/"examples/stage4_tid_case.json").string());
  auto budgetWithoutPhysics=stage4Document.At("accuracy").At("systematic_budget");
  auto& components=budgetWithoutPhysics.At("components").AsArray();
  components.erase(std::remove_if(components.begin(),components.end(),[](const auto& item){
    return item.At("category").AsString()=="physics_reference";
  }),components.end());
  auto& correlations=budgetWithoutPhysics.At("correlations").AsArray();
  correlations.erase(std::remove_if(correlations.begin(),correlations.end(),[](const auto& item){
    return item.At("component_i").AsString()=="registered_physics_reference_envelope" ||
           item.At("component_j").AsString()=="registered_physics_reference_envelope";
  }),correlations.end());
  missingPhysics.At("accuracy").AsObject()["systematic_budget"]=budgetWithoutPhysics;
  const fs::path missingPhysicsPath=scratch/"complete_without_physics_reference.json";
  { std::ofstream out(missingPhysicsPath); out<<missingPhysics.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(missingPhysicsPath.string()); },
              "complete Stage-5 uncertainty budget without physics_reference provenance");

  shield_sim::Json evidenceDocument(shield_sim::Json::Object{
    {"comparisons",shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{{"actual",proton100.nielMeVCm2PerG},
        {"expected",2.9757e-3},{"id","SR_NIEL_proton_100MeV_node"},{"unit","MeV cm2/g"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",alpha100.nielMeVCm2PerG},
        {"expected",1.5904e-2},{"id","SR_NIEL_alpha_100MeV_node"},{"unit","MeV cm2/g"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",protonMid.nielMeVCm2PerG},
        {"expected",expectedMid},{"id","charged_log_log_geometric_midpoint"},{"unit","MeV cm2/g"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",neutronFold.dddMeVPerG},
        {"expected",expectedNeutronDdd},{"id","ASTM_single_group_DDD_fold"},{"unit","MeV/g"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",neutronFold.neutronEquivalentFluencePerCm2},
        {"expected",expectedNeutronDdd/reference},{"id","ASTM_95MeVmb_neutron_equivalent"},{"unit","cm-2"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",protonFold.dddMeVPerG},
        {"expected",expectedCharged},{"id","analytic_log_log_interval_fold"},{"unit","MeV/g"}})
    })},
    {"schema_version","shieldSim.validation_evidence.v1"},
    {"test_case","stage5_core"}});
  std::ofstream evidenceStream(evidence/"stage5_core.json",std::ios::binary|std::ios::trunc);
  evidenceStream<<evidenceDocument.DumpCanonical()<<'\n';

  if(failures){ std::cerr<<failures<<" Stage-5 core test(s) failed\n"; return 1; }
  std::cout<<"All Stage-5 core tests passed\n";
  return 0;
}
