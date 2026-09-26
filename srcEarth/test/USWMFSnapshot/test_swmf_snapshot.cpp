#include "../../util/SWMFSnapshotContract.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

namespace SS=Earth::SWMFSnapshot;

namespace {

int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

template<typename F>
void CheckThrows(F function,const std::string& message) {
  bool threw=false;
  try { function(); }
  catch (const std::exception&) { threw=true; }
  Check(threw,message);
}

SS::Snapshot ReferenceSnapshot() {
  SS::Snapshot snapshot;
  snapshot.epochUTC="2024-05-10T12:00:00Z";
  snapshot.simulationTime_s=43200.0;
  snapshot.usedLeafBlocks=2;
  snapshot.blockCellsX=2;
  snapshot.blockCellsY=2;
  snapshot.blockCellsZ=1;
  snapshot.domain.enabled=true;
  snapshot.domain.minimum_m[0]=0.0;
  snapshot.domain.minimum_m[1]=0.0;
  snapshot.domain.minimum_m[2]=-1.0;
  snapshot.domain.maximum_m[0]=4.0;
  snapshot.domain.maximum_m[1]=2.0;
  snapshot.domain.maximum_m[2]=1.0;

  for (long int block=0;block<2;++block) {
    for (int i=0;i<2;++i) {
      for (int j=0;j<2;++j) {
        SS::Cell cell;
        cell.blockId=block;
        cell.i=i; cell.j=j; cell.k=0;
        cell.position_m[0]=2.0*block+0.5+i;
        cell.position_m[1]=0.5+j;
        cell.position_m[2]=0.0;
        const double n=static_cast<double>(1+block*4+i*2+j);
        cell.magneticField_T[0]=n*1.0e-9;
        cell.magneticField_T[1]=-n*2.0e-9;
        cell.magneticField_T[2]=n*3.0e-9;
        cell.plasmaVelocity_m_s[0]=350000.0+100.0*n;
        cell.plasmaVelocity_m_s[1]=-2500.0+10.0*n;
        cell.plasmaVelocity_m_s[2]=50.0*n;
        snapshot.cells.push_back(cell);
      }
    }
  }
  SS::FinalizeIdentity(snapshot);
  return snapshot;
}

std::string ReadBytes(const std::string& path) {
  std::ifstream input(path.c_str(),std::ios::in|std::ios::binary);
  if (!input) throw std::runtime_error("cannot read test artifact: "+path);
  std::ostringstream bytes;
  bytes << input.rdbuf();
  return bytes.str();
}

} // namespace

int main(int argc,char** argv) {
  if (argc!=3) {
    std::cerr << "usage: test_swmf_snapshot OUTPUT_A.csv OUTPUT_B.csv\n";
    return 2;
  }

  SS::Snapshot reference=ReferenceSnapshot();
  SS::Validate(reference,true);

  // S9-U01: deterministic physical identity is a numerical reference, not merely a
  // nonempty-string check.  Any change to canonical ordering, dimensions, units, epoch,
  // or B/u values must be reviewed explicitly by updating this known answer.
  Check(reference.meshRevision=="swmf-mesh-v1-8f6e2aeb23c9ab4f",
        "mesh revision differs from the fixed Step-9 reference: "+
        reference.meshRevision);
  Check(reference.contentFingerprint=="swmf-state-v1-f9fca28cd353b89b",
        "content fingerprint differs from the fixed Step-9 reference: "+
        reference.contentFingerprint);
  Check(reference.snapshotId=="field-v1-18530916da122c7a",
        "snapshot ID differs from the fixed Step-9 reference: "+reference.snapshotId);

  // S9-U02: exact ideal-MHD sign and units reference.
  const double velocity[3]={100000.0,0.0,0.0};
  const double magnetic[3]={0.0,0.0,2.0e-5};
  double electric[3]={0.0,0.0,0.0};
  SS::DeriveElectricField(velocity,magnetic,electric);
  Check(electric[0]==0.0 && electric[1]==2.0 && electric[2]==0.0,
        "E=-u x B analytic reference failed");
  Check(reference.electricFieldMode==SS::kMagneticOnlyMode &&
        !SS::UsesExperimentalIdealMhdElectricField(reference),
        "Phase-1 snapshot default is not magnetic-only");
  SS::Snapshot experimental=reference;
  experimental.electricFieldMode=SS::kExperimentalIdealMhdMode;
  SS::FinalizeIdentity(experimental);
  Check(SS::UsesExperimentalIdealMhdElectricField(experimental) &&
        experimental.snapshotId!=reference.snapshotId,
        "experimental derived-E mode is not explicit in physical identity");

  // S9-U03: portable write/read must preserve every value and identity. A second write
  // of the parsed result must be byte-for-byte identical, which catches unstable
  // ordering, precision, locale, or metadata serialization—not just numeric equality.
  SS::Write(reference,argv[1]);
  const SS::Snapshot replay=SS::Read(argv[1]);
  SS::Write(replay,argv[2]);
  const SS::Comparison exact=SS::Compare(reference,replay,0.0,0.0);
  Check(exact.passed && exact.failedComponents==0 &&
        replay.snapshotId==reference.snapshotId,
        "exact exported-snapshot replay failed");
  Check(ReadBytes(argv[1])==ReadBytes(argv[2]),
        "SWMF serialization is not byte-stable after read/write");

  // S9-U04: canonical identity must not depend on record/rank partition order.  Reverse
  // and interleave the same cells to model 1x1, 2x8, and 8x16 collection layouts.
  SS::Snapshot reversed=reference;
  std::reverse(reversed.cells.begin(),reversed.cells.end());
  SS::FinalizeIdentity(reversed);
  Check(reversed.snapshotId==reference.snapshotId,
        "snapshot identity depends on cell traversal order");
  Check(reversed.meshRevision==reference.meshRevision,
        "mesh revision depends on cell traversal order");
  SS::Snapshot interleaved=reference;
  std::vector<SS::Cell> reordered;
  for (std::size_t offset=0;offset<4;++offset)
    for (std::size_t n=offset;n<interleaved.cells.size();n+=4)
      reordered.push_back(interleaved.cells[n]);
  interleaved.cells=reordered;
  SS::FinalizeIdentity(interleaved);
  Check(interleaved.snapshotId==reference.snapshotId,
        "snapshot identity depends on simulated MPI/thread layout");

  // S9-U05: restart reproducibility and sensitivity to a physical change.
  const SS::Snapshot restarted=ReferenceSnapshot();
  Check(restarted.snapshotId==reference.snapshotId,
        "same-epoch restart did not reproduce the snapshot identity");
  SS::Snapshot changed=reference;
  changed.cells[0].magneticField_T[0]+=1.0e-15;
  SS::FinalizeIdentity(changed);
  Check(changed.snapshotId!=reference.snapshotId,
        "physical B change did not change snapshot identity");
  Check(changed.meshRevision==reference.meshRevision,
        "B-only change incorrectly changed mesh revision");
  SS::Snapshot remeshed=reference;
  remeshed.cells[0].position_m[0]+=1.0e-6;
  SS::FinalizeIdentity(remeshed);
  Check(remeshed.meshRevision!=reference.meshRevision,
        "cell-centre change did not change mesh revision");

  // S9-U06: all malformed/corrupted inputs fail closed.
  SS::Snapshot corrupted=reference;
  corrupted.cells[0].magneticField_T[0]+=1.0e-12; // do not update declared identity
  CheckThrows([&]() { SS::Validate(corrupted,true); },
              "content corruption was accepted under the old identity");
  SS::Snapshot staleMesh=reference;
  staleMesh.meshRevision="swmf-mesh-v1-0000000000000000";
  CheckThrows([&]() { SS::Validate(staleMesh,true); },
              "stale mesh revision was accepted");
  SS::Snapshot badUnit=reference;
  badUnit.magneticFieldUnit="nT";
  CheckThrows([&]() { SS::Validate(badUnit,false); },
              "nT snapshot was accepted as SI tesla");
  SS::Snapshot badFrame=reference;
  badFrame.frame="GSE";
  CheckThrows([&]() { SS::Validate(badFrame,false); },
              "non-GSM snapshot was accepted");
  SS::Snapshot badElectricConvention=reference;
  badElectricConvention.electricFieldConvention="E=+u_cross_B";
  CheckThrows([&]() { SS::Validate(badElectricConvention,false); },
              "wrong ideal-MHD electric-field sign was accepted");
  SS::Snapshot duplicate=reference;
  duplicate.cells[1]=duplicate.cells[0];
  CheckThrows([&]() { SS::Validate(duplicate,false); },
              "duplicate compact cell was accepted");
  SS::Snapshot missing=reference;
  missing.cells.pop_back();
  CheckThrows([&]() { SS::Validate(missing,false); },
              "incomplete snapshot was accepted");
  SS::Snapshot nonfinite=reference;
  nonfinite.cells[0].plasmaVelocity_m_s[0]=
      std::numeric_limits<double>::infinity();
  CheckThrows([&]() { SS::Validate(nonfinite,false); },
              "non-finite plasma velocity was accepted");
  SS::Snapshot badTime=reference;
  badTime.simulationTime_s=-1.0;
  CheckThrows([&]() { SS::Validate(badTime,false); },
              "negative coupled simulation time was accepted");
  SS::Snapshot badMode=reference;
  badMode.electricFieldMode="ON";
  CheckThrows([&]() { SS::Validate(badMode,false); },
              "unlabelled electric-field mode was accepted");

  const std::string truncatedPath=std::string(argv[2])+".truncated";
  {
    std::ofstream truncated(truncatedPath.c_str());
    truncated << "# schema=" << SS::kSchema << "\n";
  }
  CheckThrows([&]() { (void)SS::Read(truncatedPath); },
              "truncated snapshot file was accepted");

  // S9-U07: cross-path tolerance is quantitative and cannot pass a changed state by
  // checking filenames alone.  Re-finalize each synthetic state so both inputs remain
  // internally valid before applying the comparison gate.  Absolute thresholds are
  // quantity-specific so a velocity tolerance can never accidentally relax B or E.
  SS::ComparisonTolerances tolerances;
  tolerances.positionAbsolute_m=1.0e-9;
  tolerances.magneticFieldAbsolute_T=1.0e-30;
  tolerances.plasmaVelocityAbsolute_m_s=1.0e-9;
  tolerances.electricFieldAbsolute_V_m=1.0e-18;
  tolerances.relative=1.0e-8;
  SS::Snapshot near=reference;
  near.cells[0].magneticField_T[0]*=(1.0+5.0e-9);
  SS::FinalizeIdentity(near);
  Check(SS::Compare(reference,near,tolerances).passed,
        "within-tolerance cross-path field was rejected");
  SS::Snapshot far=reference;
  far.cells[0].magneticField_T[0]*=(1.0+5.0e-5);
  SS::FinalizeIdentity(far);
  const SS::Comparison mismatch=SS::Compare(reference,far,tolerances);
  Check(!mismatch.passed && mismatch.failedComponents>0,
        "out-of-tolerance cross-path field was accepted");
  SS::Snapshot nearbyPosition=reference;
  nearbyPosition.cells[0].position_m[0]+=5.0e-10;
  SS::FinalizeIdentity(nearbyPosition);
  Check(SS::Compare(reference,nearbyPosition,tolerances).passed,
        "within-tolerance cell-centre displacement was rejected");
  SS::Snapshot displaced=reference;
  displaced.cells[0].position_m[0]+=2.0e-9;
  SS::FinalizeIdentity(displaced);
  Check(!SS::Compare(reference,displaced,tolerances).passed,
        "out-of-tolerance cell-centre displacement was accepted");
  SS::ComparisonTolerances invalid=tolerances;
  invalid.plasmaVelocityAbsolute_m_s=-1.0;
  CheckThrows([&]() { (void)SS::Compare(reference,replay,invalid); },
              "negative quantity-specific comparison tolerance was accepted");

  // S9-U08: unavailable, stale, corrupt, and queued states fail closed. The queue
  // models a new coupling receive arriving while a batch owns the active generation.
  SS::PublicationQueue queue;
  Check(queue.State()==SS::PublicationState::Unavailable,
        "new publication queue is not unavailable");
  CheckThrows([&]() { (void)queue.Freeze(); },
              "unavailable snapshot was frozen");
  queue.Publish(reference);
  Check(queue.State()==SS::PublicationState::Ready,
        "validated snapshot did not become ready");
  const std::string frozenId=queue.Freeze().snapshotId;
  CheckThrows([&]() { (void)queue.Freeze(); },
              "the same snapshot generation was frozen twice");
  queue.Publish(changed);
  Check(queue.State()==SS::PublicationState::Queued && queue.HasPending() &&
        queue.Active().snapshotId==frozenId,
        "new receive replaced the frozen generation instead of being queued");
  CheckThrows([&]() { queue.Release(changed.snapshotId); },
              "stale batch identity released another generation");
  queue.Release(frozenId);
  Check(queue.State()==SS::PublicationState::Ready &&
        queue.Active().snapshotId==changed.snapshotId,
        "queued generation was not promoted after release");
  CheckThrows([&]() { queue.Publish(corrupted); },
              "corrupt pending snapshot was published");
  Check(queue.State()==SS::PublicationState::Failed,
        "corrupt receive did not block the previous ready generation");
  CheckThrows([&]() { (void)queue.Active(); },
              "corrupt receive silently fell back to the previous product state");

  SS::PublicationQueue corruptBehindFrozen;
  corruptBehindFrozen.Publish(reference);
  const std::string originalFrozenId=corruptBehindFrozen.Freeze().snapshotId;
  CheckThrows([&]() { corruptBehindFrozen.Publish(corrupted); },
              "corrupt receive behind a frozen batch was accepted");
  Check(corruptBehindFrozen.State()==SS::PublicationState::Failed &&
        corruptBehindFrozen.Active().snapshotId==originalFrozenId,
        "corrupt queued receive modified the immutable active batch");
  corruptBehindFrozen.Release(originalFrozenId);
  Check(corruptBehindFrozen.State()==SS::PublicationState::Failed,
        "release fell back to the old snapshot after a corrupt queued receive");
  CheckThrows([&]() { (void)corruptBehindFrozen.Active(); },
              "old snapshot remained available after corrupt queued receive release");

  // Product status is independently fail-closed. STARTED/RUNNING are not passing
  // states; production writes FAILED first and PASS only after all requested products.
  const std::string failed=SS::BuildProductStatusJson(
      "FAILED",reference.snapshotId,reference.epochUTC,".test","snapshot.csv",
      true,true,"not completed");
  const std::string passed=SS::BuildProductStatusJson(
      "PASS",reference.snapshotId,reference.epochUTC,".test","snapshot.csv",
      true,true,"complete");
  Check(failed.find("\"status\": \"FAILED\"")!=std::string::npos &&
        failed.find("\"status\": \"PASS\"")==std::string::npos,
        "precalculation status is not fail-closed");
  Check(passed.find("\"status\": \"PASS\"")!=std::string::npos,
        "completed product status is not PASS");
  CheckThrows([&]() {
      (void)SS::BuildProductStatusJson(
          "RUNNING",reference.snapshotId,reference.epochUTC,"","",true,false,"");
    },"nonterminal RUNNING status was accepted");

  if (failures!=0) {
    std::cerr << failures << " Step-9 SWMF snapshot test(s) failed\n";
    return 1;
  }
  std::cout << "PASS: Step-9 SWMF snapshot contract and numerical references\n";
  return 0;
}
