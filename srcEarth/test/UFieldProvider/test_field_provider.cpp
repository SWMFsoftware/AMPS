#include "../../util/FieldProvider.h"
#include "../../gridless/DipoleInterface.h"

#include <algorithm>
#include <atomic>
#include <cmath>
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

bool Near(double actual,double expected,double relative=1.0e-13,
          double absolute=1.0e-18) {
  return std::fabs(actual-expected)<=absolute+
      relative*std::max(std::fabs(actual),std::fabs(expected));
}

class UniformSnapshot : public Earth::Field::IFieldSnapshot {
public:
  UniformSnapshot(const Earth::Field::SnapshotMetadata& metadata,
                  const double magnetic_T[3],const double electric_V_m[3])
      : metadata_(metadata) {
    for (int d=0;d<3;++d) {
      magnetic_T_[d]=magnetic_T[d];
      electric_V_m_[d]=electric_V_m[d];
    }
  }

  const Earth::Field::SnapshotMetadata& Metadata() const override {
    return metadata_;
  }

  Earth::Field::FieldSample Sample(
      const Earth::Field::FieldQuery& query) const override {
    Earth::Field::FieldSample sample;
    sample.snapshotId=metadata_.snapshotId;
    sample.interpolation=metadata_.interpolation;
    sample.status=Earth::Field::ValidateQuery(metadata_,query,&sample.message);
    if (!sample.ok()) return sample;
    for (int d=0;d<3;++d) {
      sample.magneticField_T[d]=magnetic_T_[d];
      sample.electricField_V_m[d]=electric_V_m_[d];
    }
    return sample;
  }

private:
  const Earth::Field::SnapshotMetadata metadata_;
  double magnetic_T_[3];
  double electric_V_m_[3];
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

    // requestId is intentionally absent. It labels a consumer, not field state.
    metadata.snapshotId=Earth::Field::MakeSnapshotId(
        metadata.sourceId,metadata.epochUTC,
        "revision="+std::to_string(revision_));
    const double magnetic[3]={3.0e-5*revision_,-2.0e-6,1.0e-5};
    const double electric[3]={0.0,1.0e-3*revision_,0.0};
    return std::shared_ptr<const Earth::Field::IFieldSnapshot>(
        new UniformSnapshot(metadata,magnetic,electric));
  }

  void SetRevision(int revision) { revision_=revision; }

private:
  int revision_;
};

} // anonymous namespace

int main() {
  using namespace Earth::Field;
  namespace Dipole=Earth::GridlessMode::Dipole;

  // U-P01: compare the deterministic ID with an independently precomputed FNV-1a
  // reference, then prove sensitivity to a physical driver.
  const std::string id1=MakeSnapshotId(
      "DIPOLE","2012-05-17T01:00:00","M=1|tilt=0");
  const std::string id2=MakeSnapshotId(
      "DIPOLE","2012-05-17T01:00:00","M=1|tilt=0");
  const std::string id3=MakeSnapshotId(
      "DIPOLE","2012-05-17T01:00:00","M=1|tilt=1");
  Check(id1=="field-v1-9a27e878e9de572e",
        "snapshot ID must match the fixed FNV-1a reference");
  Check(id1==id2,"identical physical state must have deterministic identity");
  Check(id1!=id3,"field-defining state change must change identity");

  UniformProvider provider;
  SnapshotRequest requestA;
  requestA.epochUTC="2012-05-17T01:00:00";
  requestA.requestId="cutoff";
  std::shared_ptr<const IFieldSnapshot> snapshot=provider.CreateSnapshot(requestA);
  SnapshotRequest requestB=requestA;
  requestB.requestId="flux";
  const std::shared_ptr<const IFieldSnapshot> samePhysical=
      provider.CreateSnapshot(requestB);
  Check(snapshot->Metadata().snapshotId==samePhysical->Metadata().snapshotId,
        "consumer request labels must not change physical snapshot identity");

  // U-P02: complete frame/unit/interpolation/domain metadata is mandatory.
  std::string detail;
  Check(ValidateMetadata(snapshot->Metadata(),&detail)==FieldSampleStatus::Valid,
        "complete immutable SI/GSM metadata must validate: "+detail);
  Check(snapshot->Metadata().positionUnit=="m" &&
        snapshot->Metadata().magneticFieldUnit=="T" &&
        snapshot->Metadata().electricFieldUnit=="V/m",
        "snapshot units must be exactly m, T, and V/m");

  // U-P03: nominal sample values and identity propagation are checked against exact
  // constants, not merely finite/nonzero conditions.
  FieldQuery query;
  query.position_m[0]=1.0;
  query.position_m[1]=2.0;
  query.position_m[2]=3.0;
  query.epochUTC=requestA.epochUTC;
  query.requireElectricField=true;
  FieldSample sample=snapshot->Sample(query);
  Check(sample.ok(),"in-domain sample at frozen epoch must be valid");
  Check(sample.snapshotId==snapshot->Metadata().snapshotId,
        "sample must propagate exact snapshot identity");
  Check(sample.magneticField_T[0]==3.0e-5 &&
        sample.magneticField_T[1]==-2.0e-6 &&
        sample.magneticField_T[2]==1.0e-5 &&
        sample.electricField_V_m[1]==1.0e-3,
        "uniform provider must preserve exact SI reference values");

  // U-P04: stale epoch is explicit and cannot become a valid/zero sample.
  query.epochUTC="2012-05-17T01:05:00";
  sample=snapshot->Sample(query);
  Check(sample.status==FieldSampleStatus::StaleEpoch,
        "mismatched query epoch must return STALE_EPOCH");

  // U-P05: domain and malformed-coordinate failures remain distinct. Domain bounds
  // are inclusive by contract.
  query.epochUTC=requestA.epochUTC;
  query.position_m[0]=10.0;
  Check(snapshot->Sample(query).ok(),"domain maximum must be inclusive");
  query.position_m[0]=10.0+1.0e-12;
  Check(snapshot->Sample(query).status==FieldSampleStatus::OutsideDomain,
        "out-of-domain point must return OUTSIDE_DOMAIN");
  query.position_m[0]=std::nan("");
  Check(snapshot->Sample(query).status==FieldSampleStatus::InvalidRequest,
        "non-finite position must return INVALID_REQUEST");

  // U-P06: advancing a provider creates a new object and cannot mutate the old field
  // values or metadata.
  const std::string frozenId=snapshot->Metadata().snapshotId;
  provider.SetRevision(2);
  std::shared_ptr<const IFieldSnapshot> next=provider.CreateSnapshot(requestA);
  query.position_m[0]=0.0;
  const FieldSample oldValue=snapshot->Sample(query);
  const FieldSample newValue=next->Sample(query);
  Check(next->Metadata().snapshotId!=frozenId,
        "new physical revision must produce a new snapshot ID");
  Check(snapshot->Metadata().snapshotId==frozenId &&
        oldValue.magneticField_T[0]==3.0e-5 &&
        newValue.magneticField_T[0]==6.0e-5,
        "old snapshot must preserve owned field state after provider advances");

  // U-P07: product synchronization accepts exact identity and rejects replacement.
  try {
    RequireSameSnapshot(snapshot->Metadata(),samePhysical->Metadata(),"unit test");
  }
  catch (...) {
    Check(false,"matching physical snapshots must pass synchronization gate");
  }
  bool mismatchCaught=false;
  try {
    RequireSameSnapshot(snapshot->Metadata(),next->Metadata(),"unit test");
  }
  catch (const std::runtime_error&) { mismatchCaught=true; }
  Check(mismatchCaught,"different snapshots must fail synchronization gate");
  SnapshotMetadata reusedIdentity=snapshot->Metadata();
  reusedIdentity.frame=CoordinateFrame::GEO;
  mismatchCaught=false;
  try {
    RequireSameSnapshot(snapshot->Metadata(),reusedIdentity,"unit test");
  }
  catch (const std::runtime_error&) { mismatchCaught=true; }
  Check(mismatchCaught,
        "reused identity with a changed metadata contract must fail synchronization");

  // U-P08: false schema/unit/interpolation/E claims fail instead of being tolerated.
  SnapshotMetadata bad=snapshot->Metadata();
  bad.schemaVersion=2;
  Check(ValidateMetadata(bad,&detail)==FieldSampleStatus::InvalidRequest,
        "unsupported metadata schema must fail");
  bad=snapshot->Metadata();
  bad.magneticFieldUnit="nT";
  Check(ValidateMetadata(bad,&detail)==FieldSampleStatus::InvalidRequest,
        "non-SI B unit must fail");
  bad=snapshot->Metadata();
  bad.interpolation=InterpolationMode::Unknown;
  Check(ValidateMetadata(bad,&detail)==FieldSampleStatus::InvalidRequest,
        "unknown interpolation mode must fail");
  bad=snapshot->Metadata();
  bad.electricFieldAvailable=false;
  query.requireElectricField=true;
  Check(ValidateQuery(bad,query,&detail)==FieldSampleStatus::SourceUnavailable,
        "required but unavailable E must fail explicitly");
  query.requireElectricField=false;

  // U-P09: production dipole values against closed-form independent references.
  const Dipole::Params aligned=Dipole::MakeParams(1.0,0.0);
  double xEquator[3]={Dipole::Re_m,0.0,0.0};
  double xPole[3]={0.0,0.0,Dipole::Re_m};
  double field[3]={0.0,0.0,0.0};
  Dipole::GetB_Tesla(xEquator,field,aligned);
  Check(Near(field[0],0.0) && Near(field[1],0.0) &&
        Near(field[2],-Dipole::B_eq_Re),
        "aligned dipole equator must be exactly -B_eq along GSM Z");
  Dipole::GetB_Tesla(xPole,field,aligned);
  Check(Near(field[0],0.0) && Near(field[1],0.0) &&
        Near(field[2],2.0*Dipole::B_eq_Re),
        "aligned dipole pole must be exactly +2 B_eq along GSM Z");

  const Dipole::Params tilted=Dipole::MakeParams(0.8,30.0);
  const double xAxis[3]={2.0*Dipole::Re_m*tilted.m_hat[0],0.0,
                         2.0*Dipole::Re_m*tilted.m_hat[2]};
  Dipole::GetB_Tesla(xAxis,field,tilted);
  const double axialMagnitude=0.25*0.8*Dipole::B_eq_Re;
  Check(Near(field[0],axialMagnitude*tilted.m_hat[0]) &&
        Near(field[1],0.0) &&
        Near(field[2],axialMagnitude*tilted.m_hat[2]),
        "tilted dipole axis value must follow the closed-form 2M/r^3 reference");

  // U-P10: provider-owned dipole parameters remain stable after deliberate legacy
  // global reconfiguration and under the 16 concurrent readers used by F4.
  const Dipole::Params frozenDipole=Dipole::MakeParams(0.85,-17.0);
  const double xGeneral[3]={1.7*Dipole::Re_m,0.2*Dipole::Re_m,-0.4*Dipole::Re_m};
  double frozenReference[3]={0.0,0.0,0.0};
  Dipole::GetB_Tesla(xGeneral,frozenReference,frozenDipole);
  Dipole::SetMomentScale(1.7);
  Dipole::SetTiltDeg(41.0);
  double afterGlobalChange[3]={0.0,0.0,0.0};
  Dipole::GetB_Tesla(xGeneral,afterGlobalChange,frozenDipole);
  for (int d=0;d<3;++d)
    Check(afterGlobalChange[d]==frozenReference[d],
          "owned dipole state changed after legacy global reconfiguration");

  std::atomic<bool> concurrentStable(true);
  std::vector<std::thread> workers;
  for (int worker=0;worker<16;++worker) {
    workers.emplace_back([&]() {
      for (int i=0;i<2000;++i) {
        double value[3]={0.0,0.0,0.0};
        Dipole::GetB_Tesla(xGeneral,value,frozenDipole);
        for (int d=0;d<3;++d)
          if (value[d]!=frozenReference[d]) concurrentStable.store(false);
      }
    });
  }
  for (std::thread& worker:workers) worker.join();
  Check(concurrentStable.load(),
        "immutable dipole snapshot must be stable across concurrent readers");

  Dipole::SetMomentScale(1.0);
  Dipole::SetTiltDeg(0.0);

  if (failures!=0) {
    std::cerr << "UFieldProvider: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UFieldProvider: PASS (10 strict contracts)\n";
  return EXIT_SUCCESS;
}
