#include "../../util/SWMFCoupledAccessContract.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace CA=Earth::SWMFCoupledAccess;

namespace {

void Require(bool condition,const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void RequireClose(double actual,double expected,double tolerance,
                  const std::string& message) {
  if (!std::isfinite(actual) || std::fabs(actual-expected)>tolerance) {
    throw std::runtime_error(message+": actual="+std::to_string(actual)+
                             ", expected="+std::to_string(expected));
  }
}

CA::BoundaryPolicy BoxPolicy(double earthRadius_m) {
  CA::BoundaryPolicy policy;
  policy.kind=CA::OuterBoundaryKind::Box;
  policy.computationalBox.minimum_m={{-40.0*earthRadius_m,-20.0*earthRadius_m,
                                      -20.0*earthRadius_m}};
  policy.computationalBox.maximum_m={{20.0*earthRadius_m,20.0*earthRadius_m,
                                      20.0*earthRadius_m}};
  CA::ValidateBoundaryPolicy(policy);
  return policy;
}

CA::BoundaryPolicy ShuePolicy(double earthRadius_m,double r0_Re=10.0,
                              double alpha=0.5) {
  CA::BoundaryPolicy policy=BoxPolicy(earthRadius_m);
  policy.kind=CA::OuterBoundaryKind::Shue;
  policy.shue.r0_Re=r0_Re;
  policy.shue.alpha=alpha;
  policy.shue.earthRadius_m=earthRadius_m;
  policy.shue.tailCapX_m=policy.computationalBox.minimum_m[0];
  CA::ValidateBoundaryPolicy(policy);
  return policy;
}

void TestCadence() {
  CA::CadenceGate gate;
  Require(gate.Evaluate(0.0,60.0).action==CA::CadenceAction::Run,
          "first complete snapshot must run");
  gate.CommitCompleted(0.0);
  const CA::CadenceDecision early=gate.Evaluate(30.0,60.0);
  Require(early.action==CA::CadenceAction::Skip,"early callback must skip");
  RequireClose(early.nextDueTime_s,60.0,0.0,"next due time");
  Require(gate.Evaluate(60.0,60.0).action==CA::CadenceAction::Run,
          "exact cadence boundary must run");
  gate.CommitCompleted(60.0);
  Require(gate.Evaluate(60.0,60.0).action==CA::CadenceAction::Skip,
          "duplicate epoch in one process must skip");
  Require(gate.Evaluate(59.0,60.0).action==CA::CadenceAction::RejectStaleTime,
          "in-process clock rollback must fail closed");

  // Restart is represented by a new scheduler instance.  It must reproduce the exact
  // restart epoch rather than inherit hidden state from the old process.
  CA::CadenceGate restarted;
  Require(restarted.Evaluate(60.0,60.0).action==CA::CadenceAction::Run,
          "fresh restart gate must evaluate restart epoch");
  restarted.CommitCompleted(60.0);
  Require(restarted.Evaluate(61.0,0.0).action==CA::CadenceAction::Run,
          "zero cadence must run each distinct callback");

  bool rejected=false;
  try { (void)gate.Evaluate(std::numeric_limits<double>::quiet_NaN(),60.0); }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"non-finite cadence time must be rejected");
  std::cout << "PASS S10-U01 collective-cadence reference\n";
}

void TestProductIdentity() {
  const std::string suffix=CA::BuildProductSuffix(600.125,"field-v1-abc:def");
  Require(suffix==
              ".swmf_t0000000600.125000000s_sidfield-v1-abc_def",
          "product suffix fixed reference changed: "+suffix);
  Require(CA::BuildProductSuffix(600.125,"field-v1-abc:def")==suffix,
          "restart must reproduce exact suffix");
  Require(CA::BuildProductSuffix(600.125,"field-v1-other")!=suffix,
          "different content ID must change suffix");
  Require(CA::BuildProductSuffix(601.125,"field-v1-abc:def")!=suffix,
          "different PT time must change suffix");
  std::cout << "PASS S10-U02 deterministic time/snapshot product identity\n";
}

void TestShueReference() {
  constexpr double re=6371200.0;
  const CA::ShueParameters autoParameters=CA::ResolveShueParameters(
      "AUTO","AUTO",2.0,-5.0,re,-40.0*re);

  // Fixed published-formula references for Pdyn=2 nPa and IMF Bz=-5 nT.  These
  // constants are intentionally stored rather than calculated by another copy of the
  // production formula inside the test.
  RequireClose(autoParameters.r0_Re,9.8062382137514579,2.0e-14,
               "Shue AUTO r0 reference");
  RequireClose(autoParameters.alpha,0.62523085238506471,2.0e-15,
               "Shue AUTO alpha reference");
  RequireClose(CA::ShueRadiusMeters(autoParameters,1.0),
               62477504.907453291,2.0e-8,
               "Shue subsolar radius reference");
  RequireClose(CA::ShueRadiusMeters(autoParameters,0.0),
               96368903.61832805,2.0e-8,
               "Shue terminator radius reference");

  bool rejected=false;
  try { (void)CA::ResolveShueParameters("AUTO","AUTO",0.0,-5.0,re,-40.0*re); }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"non-positive PDYN must reject SHUE AUTO");
  std::cout << "PASS S10-U03 Shue analytic reference solution\n";
}

void TestBoundaryClassification() {
  constexpr double re=6371200.0;
  const CA::BoundaryPolicy box=BoxPolicy(re);
  Require(CA::ClassifyPosition(box,{{0.0,0.0,0.0}},true)==
              CA::PositionDisposition::Interior,
          "BOX origin must be interior");
  Require(CA::ClassifyPosition(box,{{21.0*re,0.0,0.0}},true)==
              CA::PositionDisposition::PhysicalEscape,
          "BOX face exit must be physical");
  Require(CA::ClassifyPosition(box,{{0.0,0.0,0.0}},false)==
              CA::PositionDisposition::MeshUnavailable,
          "in-box data hole must not be physical escape");

  const CA::BoundaryPolicy shue=ShuePolicy(re);
  Require(CA::ClassifyPosition(shue,{{9.0*re,0.0,0.0}},true)==
              CA::PositionDisposition::Interior,
          "point inside Shue nose");
  Require(CA::ClassifyPosition(shue,{{11.0*re,0.0,0.0}},true)==
              CA::PositionDisposition::PhysicalEscape,
          "point beyond Shue nose");

  const CA::SegmentResult nose=CA::ClassifySegment(
      shue,{{9.0*re,0.0,0.0}},{{11.0*re,0.0,0.0}});
  Require(nose.disposition==CA::SegmentDisposition::PhysicalEscape &&
          nose.surface=="SHUE","Shue nose segment must physically escape");
  RequireClose(nose.fraction,0.5,2.0e-14,"Shue nose crossing fraction");

  // At X=0 and alpha=0.5, the magnetopause is 10*sqrt(2) Re.  A Y box face at
  // 10 Re is therefore a numerical-data boundary, not physical access.
  CA::BoundaryPolicy sideLimited=shue;
  sideLimited.computationalBox.maximum_m[1]=10.0*re;
  CA::ValidateBoundaryPolicy(sideLimited);
  const CA::SegmentResult side=CA::ClassifySegment(
      sideLimited,{{0.0,9.0*re,0.0}},{{0.0,11.0*re,0.0}});
  Require(side.disposition==CA::SegmentDisposition::MeshUnavailable &&
          side.surface=="YMAX",
          "computational side face inside Shue surface must be unavailable");

  const CA::SegmentResult tail=CA::ClassifySegment(
      shue,{{-39.0*re,0.0,0.0}},{{-41.0*re,0.0,0.0}});
  Require(tail.disposition==CA::SegmentDisposition::PhysicalEscape &&
          tail.surface=="XMIN_TAIL_CAP",
          "configured nightside cap must be a physical escape");
  RequireClose(tail.fraction,0.5,2.0e-14,"tail-cap crossing fraction");
  std::cout << "PASS S10-U04 physical escape versus unavailable mesh\n";
}

void TestManifest() {
  const std::vector<std::string> artifacts={"cutoff_3d_points.dat",
                                             "cutoff_3d_dir_access_loc_000000.dat"};
  const std::string suffix=CA::BuildProductSuffix(3600.0,"field-v1-reference");
  const std::string manifest=CA::BuildAccessManifestJson(
      "PASS","field-v1-reference","swmf-state-v1-reference",
      "swmf-mesh-v1-reference","2024-05-10T12:00:00.000000000Z",3600.0,
      "POINTS","SHUE",suffix,artifacts,"complete");
  Require(manifest.find("\"RESULT\": \"PASS\"")!=std::string::npos,
          "manifest RESULT missing");
  Require(manifest.find("INSTANTANEOUS_QUASI_STATIC")!=std::string::npos,
          "Phase-1 interpretation missing");
  Require(manifest.find("field-v1-reference")!=std::string::npos,
          "snapshot identity missing");
  Require(manifest.find("cutoff_3d_dir_access_loc_000000.dat")!=std::string::npos,
          "artifact list missing");
  bool rejected=false;
  try {
    (void)CA::BuildAccessManifestJson(
        "UNKNOWN","id","content","mesh","epoch",0.0,"POINTS","BOX","",{},"");
  }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"invalid manifest result must be rejected");

  rejected=false;
  try {
    (void)CA::BuildAccessManifestJson(
        "PASS","field-v1-reference","","mesh","epoch",3600.0,
        "POINTS","BOX",suffix,artifacts,"complete");
  }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"missing manifest provenance must be rejected");

  rejected=false;
  try {
    (void)CA::BuildAccessManifestJson(
        "PASS","field-v1-reference","content","mesh","epoch",3600.0,
        "POINTS","BOX",suffix,{},"complete");
  }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"PASS manifest without an artifact must be rejected");

  rejected=false;
  try {
    (void)CA::BuildAccessManifestJson(
        "PASS","field-v1-reference","content","mesh","epoch",3600.0,
        "POINTS","BOX",suffix+"_wrong",artifacts,"complete");
  }
  catch (const std::invalid_argument&) { rejected=true; }
  Require(rejected,"manifest suffix/provenance mismatch must be rejected");
  std::cout << "PASS S10-U05 fail-closed artifact-manifest contract\n";
}

} // namespace

int main() {
  try {
    TestCadence();
    TestProductIdentity();
    TestShueReference();
    TestBoundaryClassification();
    TestManifest();
    std::cout << "RESULT: PASS\n";
    return EXIT_SUCCESS;
  }
  catch (const std::exception& error) {
    std::cerr << "RESULT: FAIL: " << error.what() << "\n";
    return EXIT_FAILURE;
  }
}
