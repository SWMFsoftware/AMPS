/*
 * Deterministic Stage-3 geometry/response tests.
 *
 * These tests prove the algebra around Geant4: rigid transforms, v2 case
 * semantics, response identity/zero/linearity, uncertainty propagation, and
 * fatal cache-boundary checks.  Stochastic transport equivalence belongs to
 * test_stage3_transport.py; keeping the two layers separate prevents a noisy
 * Monte Carlo pass from hiding an exact software error.
 */

#include "CaseDefinition.hh"
#include "ContentHash.hh"
#include "GeometryTransform.hh"
#include "ResponseMatrix.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
int failures=0;
void Check(bool value,const std::string& description){
  if(!value){ ++failures; std::cerr<<"FAIL: "<<description<<'\n'; }
}
template<class F> void CheckThrows(F function,const std::string& description){
  try { function(); } catch(const std::exception&) { return; }
  ++failures; std::cerr<<"FAIL: expected rejection: "<<description<<'\n';
}
bool Near(double a,double b,double relative=1.0e-12){
  return std::abs(a-b)<=relative*std::max({1.0,std::abs(a),std::abs(b)});
}

shield_sim::Json Cell(const std::string& id,double response,double error,
                      const std::string& geometry="fixed"){
  const shield_sim::Json keyDocument(
      shield_sim::Json::Object{{"geometry",geometry}});
  const std::string responseKey=shield_sim::Sha256(keyDocument.DumpCanonical());
  const double energy=id=="b"?200.0:id=="c"?300.0:id=="zero"?400.0:100.0;
  const double sampleVariance=100.0*error*error;
  const double sum=100.0*response;
  const double sumSquares=99.0*sampleVariance+sum*sum/100.0;
  return shield_sim::Json(shield_sim::Json::Object{
    {"absorber_id","test"},{"case_hash","case"},{"histories",100.0},
    {"incident",shield_sim::Json(shield_sim::Json::Object{
      {"A",1},{"Z",1},{"energy_coordinate","total_kinetic_energy"},
      {"energy_coordinate_unit","MeV"},{"energy_coordinate_value",energy},
      {"incident_id",id},{"mu",1.0},{"particle","proton"},{"phi_rad",0.0},
      {"total_kinetic_energy_MeV",energy}})},
    {"normalization","per primary"},{"random_seed",1.0},
    {"response_key",responseKey},
    {"response_key_document",keyDocument},
    {"schema_version","shieldSim.response_cell.v1"},
    {"score_covariance",shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{
        {"covariance_of_mean",error*error},{"histories",100.0},
        {"sample_covariance",sampleVariance},{"score_id_i","score"},
        {"score_id_j","score"},{"sum_products",sumSquares},
        {"unit","MeV/primary * MeV/primary"}})})},
    {"scores",shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{{"histories",100.0},
        {"mean_per_primary",response},{"sample_variance",sampleVariance},
        {"score_id","score"},{"standard_error",error},{"sum",sum},
        {"sum_squares",sumSquares},{"unit","MeV/primary"}})})},
    {"secondary_production_by_layer",shield_sim::Json(shield_sim::Json::Array{})},
    {"validation_status","test"}});
}
}

int main(int argc,char** argv){
  if(argc!=3){ std::cerr<<"usage: stage3_core_tests SOURCE_ROOT EVIDENCE_DIR\n"; return 2; }
  namespace fs=std::filesystem;
  const fs::path root=argv[1], evidence=argv[2];
  fs::create_directories(evidence);
  const fs::path scratch=evidence/"stage3_core_scratch";
  fs::create_directories(scratch);

  const double half=std::sqrt(0.5);
  shield_sim::RigidTransform transform{{10.0,-2.0,3.0},{half,0.0,half,0.0}};
  const auto world=transform.TransformPoint({1.0,2.0,3.0});
  const auto local=transform.InverseTransformPoint(world);
  Check(Near(local.x,1.0)&&Near(local.y,2.0)&&Near(local.z,3.0),
        "rigid transform and inverse close exactly");
  const auto q=transform.rotation.Canonical();
  const auto negative=shield_sim::Quaterniond{-q.w,-q.x,-q.y,-q.z}.Canonical();
  Check(Near(q.w,negative.w)&&Near(q.x,negative.x)&&Near(q.y,negative.y)&&Near(q.z,negative.z),
        "q and -q have one canonical representation");
  CheckThrows([]{ (void)shield_sim::Quaterniond{2.0,0.0,0.0,0.0}.Canonical(); },
              "non-unit transform quaternion");

  const auto caseV2=shield_sim::CaseDefinition::Load(
      (root/"examples/stage3_response_cell_case.json").string());
  Check(caseV2->schemaVersion=="shieldSim.case.v2" && caseV2->responseIncident.enabled,
        "case-v2 response cell resolves explicitly");
  Check(caseV2->responseIncident.atomicNumber==1 && caseV2->responseIncident.massNumber==1 &&
        Near(caseV2->responseIncident.totalKineticEnergyMeV,100.0),
        "response species and total energy resolve exactly");
  Check(caseV2->shieldLayers.size()==2 && caseV2->shieldLayers[0].id=="outer_al" &&
        caseV2->shieldLayers[1].id=="inner_hdpe",
        "case-v2 preserves the declared physical layer order");
  Check(Near(caseV2->ShieldTotalThicknessMm(),3.0) &&
        Near(caseV2->ShieldTotalArealDensityGPerCm2(),2.699*0.1+0.950*0.2),
        "multilayer thickness and normal areal-density ledgers are exact");
  // The shield transform origin is its upstream reference surface.  Layer
  // centers are therefore -19.5 and -18.0 mm and the final interface is
  // -17.0 mm for this hand-calculated identity-transform fixture.
  Check(Near(caseV2->shieldTransform.TransformPoint({0.0,0.0,0.5}).z,-19.5) &&
        Near(caseV2->shieldTransform.TransformPoint({0.0,0.0,2.0}).z,-18.0) &&
        Near(caseV2->shieldTransform.TransformPoint({0.0,0.0,3.0}).z,-17.0),
        "multilayer interface and center positions follow the declared local frame");

  // Reversing unlike layers must change the configuration identity and the
  // parser must not sort them back into a material-based order.
  auto reversedDocument=shield_sim::Json::ParseFile(
      (root/"examples/stage3_response_cell_case.json").string());
  reversedDocument.At("material_library").AsObject()["path"]=
      (root/"data/material_library.v1.json").string();
  auto& reversedLayers=reversedDocument.At("geometry").At("shield_stack").At("layers").AsArray();
  std::reverse(reversedLayers.begin(),reversedLayers.end());
  const fs::path reversedPath=scratch/"reversed_layers.json";
  { std::ofstream out(reversedPath); out<<reversedDocument.DumpCanonical()<<'\n'; }
  const auto reversedCase=shield_sim::CaseDefinition::Load(reversedPath.string());
  Check(reversedCase->shieldLayers.front().id=="inner_hdpe" &&
        reversedCase->geometryHash!=caseV2->geometryHash,
        "A|B and B|A retain distinct order and geometry hashes");
  auto forbiddenQuantity=reversedDocument;
  forbiddenQuantity.At("quantities").AsArray()={shield_sim::Json("DDD")};
  const fs::path forbiddenPath=scratch/"forbidden_stage3_quantity.json";
  { std::ofstream out(forbiddenPath); out<<forbiddenQuantity.DumpCanonical()<<'\n'; }
  CheckThrows([&]{ (void)shield_sim::CaseDefinition::Load(forbiddenPath.string()); },
              "case-v2 DDD remains forbidden after the Stage-4 TID gate");
  const auto& absorber=caseV2->absorbers.front();
  const double volumeCm3=(absorber.sizeXmm/10.0)*(absorber.sizeYmm/10.0)*
                         (absorber.thicknessMm/10.0);
  const double massKg=volumeCm3*absorber.densityGPerCm3/1000.0;
  Check(Near(volumeCm3,2.5)&&Near(massKg,0.005825),
        "slab analytic volume and mass match the resolved dimensions/density");
  Check(Near(absorber.thicknessMm/0.5,2.0),
        "non-interacting slab chord follows t/mu");

  const fs::path cellA=scratch/"a.json",cellB=scratch/"b.json",bad=scratch/"bad.json";
  { std::ofstream out(cellA); out<<Cell("a",2.0,0.1).DumpCanonical()<<'\n'; }
  { std::ofstream out(cellB); out<<Cell("b",6.0,0.2).DumpCanonical()<<'\n'; }
  { std::ofstream out(bad); out<<Cell("c",1.0,0.1,"different").DumpCanonical()<<'\n'; }
  auto matrix=shield_sim::ResponseMatrix::Assemble({cellA.string(),cellB.string()});
  const fs::path manifest=scratch/"matrix.json";
  matrix.Write(manifest.string());
  matrix=shield_sim::ResponseMatrix::Load(manifest.string());
  const auto selectedMatrix=shield_sim::ResponseMatrix::AssembleSelected(
      {cellA.string(),cellB.string()},{"score"});
  Check(selectedMatrix.Cells().size()==2 &&
        selectedMatrix.Cells().front().scores.size()==1 &&
        selectedMatrix.Cells().front().meanCovariance.size()==1,
        "selected-score assembly retains a complete requested covariance submatrix");
  CheckThrows([&]{ (void)shield_sim::ResponseMatrix::AssembleSelected(
      {cellA.string(),cellB.string()},{"absent"}); },
      "selected-score assembly rejects an absent observable");
  CheckThrows([&]{ (void)shield_sim::ResponseMatrix::Assemble({cellA.string(),bad.string()}); },
              "incompatible response keys");

  const auto weights=shield_sim::Json::Parse(R"json(
    {"schema_version":"shieldSim.response_weights.v1",
     "matrix_hash":"MATRIX_HASH",
     "normalization":"unit test normalized weights",
     "quadrature_policy":"exact hand weights",
     "interpolation_policy":"exact_nodes",
     "weights":[{"incident_id":"a","weight":0.25},{"incident_id":"b","weight":0.75}],
     "source_weight_covariance":[[0.0001,0.0],[0.0,0.0004]]}
  )json");
  auto resolvedWeights=weights;
  resolvedWeights.AsObject()["matrix_hash"]=matrix.MatrixHash();
  const auto fold=matrix.Fold(resolvedWeights);
  const auto& result=fold.At("results").AsArray().front().AsObject();
  const double expectedMean=0.25*2.0+0.75*6.0;
  const double expectedMatrix=std::sqrt(std::pow(0.25*0.1,2)+std::pow(0.75*0.2,2));
  const double expectedSource=std::sqrt(2.0*2.0*0.0001+6.0*6.0*0.0004);
  Check(Near(result.at("value").AsNumber(),expectedMean),"response superposition is exact");
  Check(Near(result.at("matrix_statistical_standard_uncertainty").AsNumber(),expectedMatrix),
        "independent matrix Monte Carlo uncertainty propagates analytically");
  Check(Near(result.at("source_standard_uncertainty").AsNumber(),expectedSource),
        "correlated-source quadratic form propagates analytically");
  const auto& foldedCovariance=fold.At("result_covariance").AsArray().front().AsObject();
  Check(Near(foldedCovariance.at("combined_covariance").AsNumber(),
             expectedMatrix*expectedMatrix+expectedSource*expectedSource),
        "folded primitive-score covariance equals analytic matrix plus source covariance");
  auto unknown=resolvedWeights;
  unknown.At("weights").AsArray()[0].AsObject()["incident_id"]="outside_grid";
  CheckThrows([&]{ (void)matrix.Fold(unknown); },"out-of-domain response node/extrapolation");
  auto indefinite=resolvedWeights;
  indefinite.AsObject()["source_weight_covariance"]=shield_sim::Json(
      shield_sim::Json::Array{
        shield_sim::Json(shield_sim::Json::Array{1.0,2.0}),
        shield_sim::Json(shield_sim::Json::Array{2.0,1.0})});
  CheckThrows([&]{ (void)matrix.Fold(indefinite); },
              "symmetric but non-positive-semidefinite source covariance");

  // Zero-response cells use an absolute comparison: zero mean, variance, and
  // propagated uncertainty must remain exact rather than becoming 0/0 in a
  // relative-tolerance check.
  const fs::path zeroPath=scratch/"zero.json";
  { std::ofstream out(zeroPath); out<<Cell("zero",0.0,0.0).DumpCanonical()<<'\n'; }
  auto zeroMatrix=shield_sim::ResponseMatrix::Assemble({zeroPath.string()});
  auto zeroWeights=shield_sim::Json::Parse(R"json(
    {"schema_version":"shieldSim.response_weights.v1",
     "matrix_hash":"HASH","normalization":"one","quadrature_policy":"exact",
     "interpolation_policy":"exact_nodes",
     "weights":[{"incident_id":"zero","weight":1.0}]}
  )json");
  zeroWeights.AsObject()["matrix_hash"]=zeroMatrix.MatrixHash();
  const auto zeroFold=zeroMatrix.Fold(zeroWeights);
  const auto& zeroResult=zeroFold.At("results").AsArray().front().AsObject();
  Check(zeroResult.at("value").AsNumber()==0.0 &&
        zeroResult.at("combined_standard_uncertainty").AsNumber()==0.0,
        "exact zero response and uncertainty survive folding");

  shield_sim::Json evidenceDocument(shield_sim::Json::Object{
    {"comparisons",shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{{"actual",volumeCm3},{"expected",2.5},
        {"id","slab_volume"},{"unit","cm3"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",massKg},{"expected",0.005825},
        {"id","slab_mass"},{"unit","kg"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",result.at("value")},{"expected",expectedMean},
        {"id","response_linear_fold"},{"unit","MeV/primary"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",result.at("matrix_statistical_standard_uncertainty")},
        {"expected",expectedMatrix},{"id","response_statistical_uncertainty"},{"unit","MeV/primary"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",result.at("source_standard_uncertainty")},
        {"expected",expectedSource},{"id","response_source_covariance"},{"unit","MeV/primary"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",caseV2->ShieldTotalThicknessMm()},
        {"expected",3.0},{"id","multilayer_total_thickness"},{"unit","mm"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",caseV2->ShieldTotalArealDensityGPerCm2()},
        {"expected",2.699*0.1+0.950*0.2},{"id","multilayer_total_areal_density"},
        {"unit","g/cm2"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",zeroResult.at("value")},
        {"expected",0.0},{"id","response_exact_zero"},{"unit","MeV/primary"}})})},
    {"schema_version","shieldSim.validation_evidence.v1"},{"test_case","stage3_core"}});
  std::ofstream evidenceOut(evidence/"stage3_core.json",std::ios::binary|std::ios::trunc);
  evidenceOut<<evidenceDocument.DumpCanonical()<<'\n';

  if(failures){ std::cerr<<failures<<" Stage-3 core test(s) failed\n"; return 1; }
  std::cout<<"All Stage-3 core tests passed\n";
  return 0;
}
