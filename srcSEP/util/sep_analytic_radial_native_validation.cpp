#include "sep_analytic_radial_native_validation.h"
#include "../sep.h"
#include "sep_background_runtime.h"
#include "sep_initialization.h"
#include "sep_analytic_radial_background.h"
#include "sep_run_configuration.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <list>
#include <string>

namespace SEP { namespace Testing { namespace {
namespace FL=PIC::FieldLine;

struct Checks {
  unsigned failures=0;
  std::string message;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    ++failures;
    if(message.size()<1800) {
      if(!message.empty())message+="; ";
      message+=reason;
    }
  }
};

Result Finish(const Checks& check,const std::string& message) {
  Result result;
  result.status=check.failures?Status::Fail:Status::Pass;
  result.message=check.failures?check.message:message;
  result.metrics.push_back({"assertion_failures",static_cast<double>(check.failures),0,"<=","count"});
  result.configuration.push_back("scope=fixed-analytic-radial-startup-and-list-registration;production_closure_qualified=false");
  result.configuration.push_back("thermal_profile=not-created;imported_velocity=declared-list-value;no-particle-or-wave-mutation");
  return result;
}

Result ActualStartup() {
  Result missing;
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||
      (DomainType!=DomainType_ParkerSpiral&&DomainType!=DomainType_StraitLine)) {
    missing.status=Status::Skip;
    missing.message="requires actual fixed analytic generated-domain startup; imported/time-dependent providers are excluded";
    return missing;
  }
  // This is an independent expected value from the declared configuration,
  // not a call to the radial helper whose native wiring is under test. The
  // documented legacy declaration applies only when no runtime input exists.
  std::array<double,3> origin{{0,0,0}};
  double speed=400000.0;
  if(Initialization::HasActive()) {
    const auto& configuration=Initialization::Active();
    origin={{configuration.parkerOriginM.x,configuration.parkerOriginM.y,configuration.parkerOriginM.z}};
    speed=configuration.solarWindSpeedMPerS;
  }
  Checks check;
  std::uint64_t vertices=0;
  double maximumVelocityRelativeError=0;
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* vertex=FL::FieldLinesAll[line].GetFirstVertex();vertex;vertex=vertex->GetNext()) {
      ++vertices;
      const double* position=vertex->GetX();
      double magnetic[3]={},velocity[3]={},previousMagnetic[3]={},previousVelocity[3]={};
      vertex->GetMagneticField(magnetic);vertex->GetPlasmaVelocity(velocity);
      vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexMagneticField,previousMagnetic);
      vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity,previousVelocity);
      const double dx=position[0]-origin[0],dy=position[1]-origin[1],dz=position[2]-origin[2];
      const double radius=std::hypot(std::hypot(dx,dy),dz);
      check.Require(std::isfinite(radius)&&radius>0,"startup vertex has no finite configured-origin radius");
      const double relative[3]={dx,dy,dz};
      for(int dimension=0;dimension<3;++dimension) {
        check.Require(std::isfinite(magnetic[dimension])&&std::isfinite(velocity[dimension]),
                      "startup current B/U is not finite");
        const double expected=speed*(relative[dimension]/radius);
        maximumVelocityRelativeError=std::max(maximumVelocityRelativeError,
            std::fabs(velocity[dimension]-expected)/speed);
        check.Require(previousMagnetic[dimension]==magnetic[dimension],
                      "fixed startup previous B differs from the actual current magnetic field");
        check.Require(previousVelocity[dimension]==velocity[dimension],
                      "fixed startup previous U differs from the actual declared current velocity");
      }
      double density=0,previousDensity=0;
      vertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&density);
      vertex->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&previousDensity);
      const double expectedDensity=5.0e6*std::pow(1.495978707e11/radius,2);
      check.Require(std::isfinite(density)&&density>0&&
                    std::fabs(density/expectedDensity-1)<=2.0e-14,
                    "startup density differs from the existing declared inverse-square SI profile");
      check.Require(previousDensity==density,"fixed startup previous n differs from current n");
    }
  check.Require(vertices>0,"actual startup contains no native field-line vertices");
  check.Require(maximumVelocityRelativeError<=2.0e-14,
                "native startup velocity does not use the configured radial origin and speed");
  auto result=Finish(check,"actual startup carries declared radial U and identical finite fixed-background current/previous B,U,n");
  result.metrics.push_back({"native_vertices",static_cast<double>(vertices),0,">","count"});
  result.metrics.push_back({"velocity_relative_error",maximumVelocityRelativeError,2.0e-14,"<=",""});
  return result;
}

class AppendedLines {
  const long originalCount_=FL::nFieldLine;
 public:
  ~AppendedLines() {
    // Test lines contain no particles and are appended only to unused slots.
    // Delete releases their native vertex/segment allocations, then restores
    // the original registry count. Existing geometry/decomposition is untouched.
    for(long line=originalCount_;line<FL::nFieldLine;++line)FL::FieldLinesAll[line].Delete();
    FL::nFieldLine=originalCount_;
  }
};

Result GeneratedAndImportedVelocities() {
  Checks check;
  if(FL::nFieldLine+2>FL::nFieldLineMax)
    return Finish(Checks{1,"native velocity fixture requires two unused field-line slots"},"");
  const long initialCount=FL::nFieldLine;
  {
    AppendedLines owned;
    const double origin[3]={2.0e10,-1.0e10,3.0e9};
    const double initial[3]={5.0e10,3.0e10,3.0e9};
    constexpr double speed=725000.0;
    // The public native adapter constructs this geometry-only record before
    // calling BuildParkerLine. Diagnose its prerequisite nonfatally first:
    // the legacy builder incorrectly requires unrelated time-step/source/
    // output fields and otherwise terminates inside PIC::exit. Do not invent
    // those fields merely to get past that real production contract defect.
    // Once the builder validates geometry independently, the actual native
    // adapter below must still pass every velocity/ownership assertion.
    Initialization::Configuration geometry;
    geometry.parkerOriginM={origin[0],origin[1],origin[2]};
    geometry.parkerInitialPointM={initial[0],initial[1],initial[2]};
    geometry.parkerLengthM=1.0e9;geometry.parkerPointCount=7;
    geometry.solarWindSpeedMPerS=speed;geometry.solarRotationRateRadPerS=2.86533e-6;
    geometry.innerRadiusM=5.0e10;geometry.outerRadiusM=5.1e10;
    geometry.minimumCellSizeM=1;geometry.globalCellSizeM=1;geometry.maximumMeshLevel=0;
    geometry.solarRefinementEnabled=false;geometry.solarSurfaceCellSizeM=1;
    geometry.solarTransitionOuterRadiusM=geometry.outerRadiusM;
    geometry.tubeRefinementEnabled=false;geometry.tubeReferenceRadiusM=geometry.outerRadiusM;
    geometry.tubeRadiusAtReferenceM=1;geometry.tubeCenterCellSizeM=1;
    std::vector<Initialization::Vec3> preflightPoints;
    const auto preflight=Initialization::BuildParkerLine(geometry,&preflightPoints);
    if(!preflight.ok()) {
      check.Require(false,"configured native generator geometry is blocked by unrelated full-configuration validation: "+preflight.message);
      return Finish(check,"");
    }
    std::list<cFieldLine> generated;
    ParkerSpiral::CreateFileLine(&generated,origin,initial,1.0e9,7,speed,2.86533e-6,5.0e-9);
    check.Require(generated.size()==7,"configured generator changed the requested vertex count");
    for(const auto& point:generated) {
      const double dx=point.x[0]-origin[0],dy=point.x[1]-origin[1],dz=point.x[2]-origin[2];
      const double radius=std::hypot(std::hypot(dx,dy),dz);
      const double expected[3]={speed*dx/radius,speed*dy/radius,speed*dz/radius};
      for(int dimension=0;dimension<3;++dimension)
        check.Require(std::isfinite(point.U[dimension])&&
                      std::fabs(point.U[dimension]-expected[dimension])<=2.0e-14*speed,
                      "configured list generator does not declare translated-origin nondefault-speed U");
    }
    Mesh::InitFieldLineAMPS(&generated);
    auto* vertex=FL::FieldLinesAll[initialCount].GetFirstVertex();
    for(const auto& point:generated) {
      check.Require(vertex!=nullptr,"native generated-line registration lost a vertex");
      if(!vertex)break;
      double velocity[3];vertex->GetPlasmaVelocity(velocity);
      for(int dimension=0;dimension<3;++dimension)
        check.Require(velocity[dimension]==point.U[dimension],"native registration overwrote declared generator U");
      vertex=vertex->GetNext();
    }
    // An arbitrary imported velocity need not be parallel to B or the line.
    // This is a storage/ownership assertion only; no analytic derivative or
    // transport qualification is inferred for this deliberately nonradial list.
    std::list<cFieldLine> imported=generated;
    unsigned index=0;
    for(auto& point:imported) {
      point.U[0]=-123000.0-17.0*index;point.U[1]=93000.0+29.0*index;
      point.U[2]=31000.0-7.0*index;++index;
    }
    Mesh::InitFieldLineAMPS(&imported);
    vertex=FL::FieldLinesAll[initialCount+1].GetFirstVertex();
    for(const auto& point:imported) {
      check.Require(vertex!=nullptr,"native imported-line registration lost a vertex");
      if(!vertex)break;
      double velocity[3];vertex->GetPlasmaVelocity(velocity);
      for(int dimension=0;dimension<3;++dimension)
        check.Require(velocity[dimension]==point.U[dimension],"native registration overwrote an imported declared U");
      vertex=vertex->GetNext();
    }
  }
  check.Require(FL::nFieldLine==initialCount,"native fixture did not restore the original geometry registry");
  return Finish(check,"actual configured generator declares complete SI U and PIC list registration preserves generated/imported velocities exactly");
}

Result AtomicStartupRejection() {
  Result missing;
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||
      (DomainType!=DomainType_ParkerSpiral&&DomainType!=DomainType_StraitLine)) {
    missing.status=Status::Skip;missing.message="requires fixed analytic generated-domain startup for native history transaction";return missing;
  }
  Checks check;
  struct Row {double* data;std::vector<double> original,requested;};
  struct Restore {
    std::vector<Row> rows;
    ~Restore(){for(const auto& row:rows)std::copy(row.original.begin(),row.original.end(),row.data);}
  } owned;
  FL::cFieldLineVertex* first=nullptr;FL::cFieldLineVertex* second=nullptr;
  for(int line=0;line<FL::nFieldLine;++line)
    for(auto* vertex=FL::FieldLinesAll[line].GetFirstVertex();vertex;vertex=vertex->GetNext()) {
      if(!first)first=vertex;else if(!second)second=vertex;
      PIC::Datum::cDatumStored* datums[]={&FL::DatumAtVertexMagneticField,&FL::DatumAtVertexPlasmaVelocity,
          &FL::DatumAtVertexPlasmaDensity,&FL::DatumAtVertexPrevious::DatumAtVertexMagneticField,
          &FL::DatumAtVertexPrevious::DatumAtVertexPlasmaVelocity,&FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity};
      for(auto* datum:datums) {
        auto* data=vertex->GetDatum_ptr(*datum);
        if(!data) {check.Require(false,"native atomic history fixture lacks required datum storage");return Finish(check,"");}
        Row row;row.data=data;row.original.assign(data,data+datum->length);owned.rows.push_back(row);
      }
    }
  if(!first||!second) {check.Require(false,"atomic history fixture needs two actual generated vertices");return Finish(check,"");}
  // The first history is a finite sentinel, while the later current endpoint
  // violates the declared radial field. A streaming partial initialization
  // would overwrite that sentinel before noticing the bad later vertex. Save
  // all requested bytes after injecting the defect, then require exact atomic
  // rejection; RAII restores every original datum even on an early failure.
  first->GetDatum_ptr(FL::DatumAtVertexPrevious::DatumAtVertexMagneticField)[0]=19;
  const int invalidRank=PIC::nTotalThreads>1?1:0;
  if(PIC::ThisThread==invalidRank)
    second->GetDatum_ptr(FL::DatumAtVertexPlasmaVelocity)[0]+=0.25*Background::AnalyticRadial::DeclaredNativeWind().speedMPerS;
  for(auto& row:owned.rows)row.requested.assign(row.data,row.data+row.original.size());
  const auto status=Background::AnalyticRadial::InitializeFixedNativeHistory();
  check.Require(!status.ok()&&status.code==Transport::StatusCode::UnsupportedConfiguration,
                "bad declared endpoint did not produce explicit unsupported analytic startup status");
  check.Require(status.message.find("rank "+std::to_string(invalidRank))!=std::string::npos,
                "collective startup rejection lost the originating failing-rank context");
  for(const auto& row:owned.rows)
    check.Require(std::memcmp(row.data,row.requested.data(),row.requested.size()*sizeof(double))==0,
                  "rejected analytic history transaction partially modified native primitives/history");
  return Finish(check,"a later invalid endpoint rejects the entire native fixed-history transaction without any earlier B/U/n write");
}
Result MoverIndependentIdentity() {
  Checks check;
  {
  struct Restore {
    Checks& check;
    Run::Configuration original=Run::Active().get();
    explicit Restore(Checks& checks):check(checks){}
    ~Restore(){const auto status=Run::InstallActive(original);check.Require(status.ok(),"cannot restore original frozen run configuration: "+status.message);}
  } owned(check);
  // Change actual frozen run values/provenance, not just the dispatch pointer:
  // the historical background hash leaked both through Run::Fingerprint.
  // Hold all physical background values fixed, restore the original frozen
  // run on every return, and do not select a runtime mover (that would reset
  // its completed-dispatch diagnostic). No particle operation is performed.
  const auto original=Background::CurrentConfigurationFingerprint();
  const auto runtimeMover=Mover::CurrentProductionMover();
  const auto completedDispatches=Mover::CompletedDispatchCount();
  for(const auto mover:{Mover::ProductionMover::Parker,Mover::ProductionMover::FocusedTransportDiffusion,
      Mover::ProductionMover::FocusedTransportMeanFreePath}) {
    auto variant=owned.original;variant.mover.value=mover;
    variant.mover.source=Run::ValueSource::CommandLine;
    const auto installed=Run::InstallActive(variant);
    check.Require(installed.ok(),"identity fixture cannot install explicit frozen mover variant: "+installed.message);
    if(!installed.ok())continue;
    check.Require(Background::CurrentConfigurationFingerprint()==original,
                  "frozen mover value/provenance changed the same physical background identity");
  }
  check.Require(Mover::CurrentProductionMover()==runtimeMover&&Mover::CompletedDispatchCount()==completedDispatches,
                "background identity fixture altered runtime mover/dispatch diagnostics");
  }
  return Finish(check,"actual background fingerprint excludes frozen mover value/provenance while the full run/checkpoint identity and runtime dispatch remain separate");
}
} // namespace

std::vector<Descriptor> AnalyticRadialNativeDescriptors() {
  std::vector<Descriptor> result;
  const char* names[]={"Fixed analytic startup primitive histories","Configured/generated and imported native velocity ownership","Atomic native analytic history rejection","Mover-independent native background identity"};
  const TestCallback callbacks[]={ActualStartup,GeneratedAndImportedVelocities,AtomicStartupRejection,MoverIndependentIdentity};
  for(unsigned index=0;index<4;++index) {
    Descriptor descriptor;descriptor.id="ARBN0"+std::to_string(index+1);
    descriptor.name=names[index];descriptor.description=names[index];
    descriptor.group="analytic-radial-native";descriptor.initialization=InitializationLevel::FieldLineModel;
    descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Deterministic;
    descriptor.repeatSupported=true;descriptor.seedPolicy="independent closed-form radial values; no random draws";
    descriptor.stateIsolation=index==0?"read-only actual native startup":index==1?
        "private appended lines deleted and registry count restored; existing native geometry untouched":index==2?
        "all current/previous B/U/n saved and restored around a controlled invalid-endpoint transaction":
        "frozen run mover value/provenance saved and restored; runtime dispatch and all background values held fixed";
    descriptor.supportedBuildModes="native PIC field-line storage; fixed analytic generated-domain startup for ARBN01";
    descriptor.callback=callbacks[index];result.push_back(descriptor);
  }
  return result;
}
} }
