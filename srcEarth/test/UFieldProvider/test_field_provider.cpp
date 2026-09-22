#include "../../util/FieldProvider.h"
#include "../../gridless/DipoleInterface.h"

#include <cmath>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

class UniformSnapshot : public Earth::Field::IFieldSnapshot {
public:
  explicit UniformSnapshot(const Earth::Field::SnapshotMetadata& metadata)
      : metadata_(metadata) {}

  const Earth::Field::SnapshotMetadata& Metadata() const override {
    return metadata_;
  }

  Earth::Field::FieldSample Sample(
      const Earth::Field::FieldQuery& query) const override {
    Earth::Field::FieldSample sample;
    sample.snapshotId=metadata_.snapshotId;
    sample.interpolation=metadata_.interpolation;
    sample.status=Earth::Field::ValidateQuery(metadata_,query,&sample.message);
    if (sample.status!=Earth::Field::FieldSampleStatus::Valid) return sample;

    sample.magneticField_T[0]=3.0e-5;
    sample.magneticField_T[1]=-2.0e-6;
    sample.magneticField_T[2]=1.0e-5;
    sample.electricField_V_m[0]=0.0;
    sample.electricField_V_m[1]=1.0e-3;
    sample.electricField_V_m[2]=0.0;
    sample.status=Earth::Field::FieldSampleStatus::Valid;
    return sample;
  }

private:
  // A value copy is central to the test: changes to the provider after CreateSnapshot
  // cannot alter the metadata or field state seen by this immutable batch snapshot.
  const Earth::Field::SnapshotMetadata metadata_;
};

class UniformProvider : public Earth::Field::IFieldProvider {
public:
  UniformProvider() : revision_(1) {}

  std::string SourceId() const override { return "TEST:UNIFORM"; }

  std::shared_ptr<const Earth::Field::IFieldSnapshot> CreateSnapshot(
      const Earth::Field::SnapshotRequest& request) override {
    Earth::Field::SnapshotMetadata metadata;
    metadata.sourceId=SourceId();
    metadata.modelName="UNIFORM";
    metadata.epochUTC=request.epochUTC;
    metadata.frame=Earth::Field::CoordinateFrame::GSM;
    metadata.interpolation=Earth::Field::InterpolationMode::DirectAnalytic;
    metadata.magneticFieldAvailable=true;
    metadata.electricFieldAvailable=true;
    metadata.immutableDuringBatch=true;
    metadata.valid=true;
    metadata.domain.enabled=true;
    for (int d=0;d<3;++d) {
      metadata.domain.minimum_m[d]=-10.0;
      metadata.domain.maximum_m[d]=10.0;
    }
    metadata.snapshotId=Earth::Field::MakeSnapshotId(
        metadata.sourceId,metadata.epochUTC,
        request.requestId+"|revision="+std::to_string(revision_));
    return std::shared_ptr<const Earth::Field::IFieldSnapshot>(
        new UniformSnapshot(metadata));
  }

  void SetRevision(int revision) { revision_=revision; }

private:
  int revision_;
};

} // anonymous namespace

int main() {
  using namespace Earth::Field;

  // U-P01: stable snapshot identity and sensitivity to physical state.
  const std::string id1=MakeSnapshotId("DIPOLE","2012-05-17T01:00:00","M=1|tilt=0");
  const std::string id2=MakeSnapshotId("DIPOLE","2012-05-17T01:00:00","M=1|tilt=0");
  const std::string id3=MakeSnapshotId("DIPOLE","2012-05-17T01:00:00","M=1|tilt=1");
  Check(id1==id2,"snapshot IDs must be deterministic");
  Check(id1!=id3,"snapshot IDs must change when a field-defining driver changes");

  UniformProvider provider;
  SnapshotRequest request;
  request.epochUTC="2012-05-17T01:00:00";
  request.requestId="batch-001";
  std::shared_ptr<const IFieldSnapshot> snapshot=provider.CreateSnapshot(request);

  // U-P02: complete epoch/frame/SI/interpolation/validity metadata.
  std::string detail;
  Check(ValidateMetadata(snapshot->Metadata(),&detail)==FieldSampleStatus::Valid,
        "a complete SI/GSM immutable metadata record must validate: "+detail);
  Check(snapshot->Metadata().positionUnit=="m", "position unit must be m");
  Check(snapshot->Metadata().magneticFieldUnit=="T", "B unit must be T");
  Check(snapshot->Metadata().electricFieldUnit=="V/m", "E unit must be V/m");
  Check(snapshot->Metadata().frame==CoordinateFrame::GSM, "frame must be explicit GSM");

  // U-P03: nominal sample and exact snapshot identity propagation.
  FieldQuery query;
  query.position_m[0]=1.0;
  query.position_m[1]=2.0;
  query.position_m[2]=3.0;
  query.epochUTC=request.epochUTC;
  query.requireElectricField=true;
  FieldSample sample=snapshot->Sample(query);
  Check(sample.ok(),"in-domain sample at the frozen epoch must be valid");
  Check(sample.snapshotId==snapshot->Metadata().snapshotId,
        "sample must carry its snapshot identifier");
  Check(std::fabs(sample.magneticField_T[0]-3.0e-5)<1.0e-18,
        "provider must preserve SI magnetic-field values");

  // U-P04: stale-epoch requests fail explicitly.
  query.epochUTC="2012-05-17T01:05:00";
  sample=snapshot->Sample(query);
  Check(sample.status==FieldSampleStatus::StaleEpoch,
        "mismatched query epoch must return STALE_EPOCH");

  // U-P05: domain and invalid-coordinate failures remain distinct.
  query.epochUTC=request.epochUTC;
  query.position_m[0]=11.0;
  sample=snapshot->Sample(query);
  Check(sample.status==FieldSampleStatus::OutsideDomain,
        "out-of-domain point must return OUTSIDE_DOMAIN");
  query.position_m[0]=std::nan("");
  sample=snapshot->Sample(query);
  Check(sample.status==FieldSampleStatus::InvalidRequest,
        "non-finite point must return INVALID_REQUEST");

  // U-P06: the provider may advance, but an already returned snapshot is immutable.
  const std::string frozenId=snapshot->Metadata().snapshotId;
  provider.SetRevision(2);
  std::shared_ptr<const IFieldSnapshot> next=provider.CreateSnapshot(request);
  Check(next->Metadata().snapshotId!=frozenId,
        "new provider revision must produce a new snapshot ID");
  Check(snapshot->Metadata().snapshotId==frozenId,
        "old immutable snapshot identity must not change with provider state");

  // U-P07: cutoff and flux product synchronization gate.
  try {
    RequireSameSnapshot(snapshot->Metadata(),snapshot->Metadata(),"unit test");
  }
  catch (...) {
    Check(false,"matching snapshots must pass the product synchronization gate");
  }
  bool mismatchCaught=false;
  try {
    RequireSameSnapshot(snapshot->Metadata(),next->Metadata(),"unit test");
  }
  catch (const std::runtime_error&) {
    mismatchCaught=true;
  }
  Check(mismatchCaught,"different cutoff/flux snapshots must fail the synchronization gate");

  // U-P08: unit and E-capability claims are validated rather than inferred.
  SnapshotMetadata bad=snapshot->Metadata();
  bad.magneticFieldUnit="nT";
  Check(ValidateMetadata(bad,&detail)==FieldSampleStatus::InvalidRequest,
        "non-SI magnetic-field metadata must fail validation");
  bad=snapshot->Metadata();
  bad.electricFieldAvailable=false;
  query.position_m[0]=0.0;
  query.requireElectricField=true;
  Check(ValidateQuery(bad,query,&detail)==FieldSampleStatus::SourceUnavailable,
        "required but unavailable E must fail explicitly");

  // U-P09: production analytic DIPOLE point values in the declared SI/GSM contract.
  const Earth::GridlessMode::Dipole::Params alignedDipole=
      Earth::GridlessMode::Dipole::MakeParams(1.0,0.0);
  double xEquator[3]={Earth::GridlessMode::Dipole::Re_m,0.0,0.0};
  double xPole[3]={0.0,0.0,Earth::GridlessMode::Dipole::Re_m};
  double B[3]={0.0,0.0,0.0};
  Earth::GridlessMode::Dipole::GetB_Tesla(xEquator,B,alignedDipole);
  Check(std::fabs(B[0])<1.0e-18 && std::fabs(B[1])<1.0e-18 &&
        std::fabs(B[2]+Earth::GridlessMode::Dipole::B_eq_Re)<1.0e-16,
        "aligned dipole equator must be -B_eq along GSM Z");
  Earth::GridlessMode::Dipole::GetB_Tesla(xPole,B,alignedDipole);
  Check(std::fabs(B[0])<1.0e-18 && std::fabs(B[1])<1.0e-18 &&
        std::fabs(B[2]-2.0*Earth::GridlessMode::Dipole::B_eq_Re)<1.0e-16,
        "aligned dipole pole must be +2 B_eq along GSM Z");

  // U-P10: a provider-owned analytic snapshot must not alias mutable global DIPOLE
  // configuration.  F4's threaded execution exposed this separate Step-3 contract
  // defect while the all-unresolved NaN was being diagnosed: trajectory workers must
  // be able to sample one frozen dipole while another provider/global setup is created,
  // and every worker must see bitwise-stable field values.
  const Earth::GridlessMode::Dipole::Params frozenDipole=
      Earth::GridlessMode::Dipole::MakeParams(0.85,-17.0);
  double xGeneral[3]={1.7*Earth::GridlessMode::Dipole::Re_m,
                     0.2*Earth::GridlessMode::Dipole::Re_m,
                    -0.4*Earth::GridlessMode::Dipole::Re_m};
  double frozenReference[3]={0.0,0.0,0.0};
  Earth::GridlessMode::Dipole::GetB_Tesla(xGeneral,frozenReference,frozenDipole);

  // Deliberately move the legacy global field to a very different physical state.
  // An old implementation sampled gParams from every supposedly immutable evaluator,
  // so this operation changed already-created snapshots and raced under F4 threads.
  Earth::GridlessMode::Dipole::SetMomentScale(1.7);
  Earth::GridlessMode::Dipole::SetTiltDeg(41.0);
  double frozenAfterGlobalChange[3]={0.0,0.0,0.0};
  Earth::GridlessMode::Dipole::GetB_Tesla(
      xGeneral,frozenAfterGlobalChange,frozenDipole);
  for (int d=0;d<3;++d) {
    Check(frozenAfterGlobalChange[d]==frozenReference[d],
          "owned DIPOLE parameters must remain unchanged after legacy global reconfiguration");
  }

  // Exercise the same read-only object from the 16-worker count used by F4.  This is
  // intentionally a concurrency contract test, not a performance benchmark.
  std::atomic<bool> concurrentStable(true);
  std::vector<std::thread> workers;
  for (int worker=0;worker<16;++worker) {
    workers.emplace_back([&]() {
      for (int sampleIndex=0;sampleIndex<2000;++sampleIndex) {
        double value[3]={0.0,0.0,0.0};
        Earth::GridlessMode::Dipole::GetB_Tesla(xGeneral,value,frozenDipole);
        for (int d=0;d<3;++d) {
          if (value[d]!=frozenReference[d]) concurrentStable.store(false);
        }
      }
    });
  }
  for (std::thread& worker:workers) worker.join();
  Check(concurrentStable.load(),
        "one immutable DIPOLE snapshot must be stable across 16 concurrent readers");

  // Restore process-global defaults so this test remains safe if embedded in a larger
  // in-process test runner in the future.
  Earth::GridlessMode::Dipole::SetMomentScale(1.0);
  Earth::GridlessMode::Dipole::SetTiltDeg(0.0);

  if (failures!=0) {
    std::cerr << "UFieldProvider: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UFieldProvider: PASS (10 contracts)\n";
  return EXIT_SUCCESS;
}
