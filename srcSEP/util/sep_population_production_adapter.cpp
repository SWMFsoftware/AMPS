#include "sep_population_production_adapter.h"
#include "sep_source_lifecycle.h"
#include "sep_run_configuration.h"
#include "../transport_common.h"
#include "../reduced_shock_particles.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace SEP { namespace PopulationControl { namespace {
namespace PB=PIC::ParticleBuffer;
namespace FL=PIC::FieldLine;
using EventKey=std::pair<int,int>;
using BinKey=std::tuple<int,int,int>;
std::map<EventKey,SourceLifecycle::Event> events;
std::vector<ConservativeMergeReport> reports;
bool passReady=false;

std::uint64_t Identity(long int pointer) {
  std::uint64_t schema=0,id=0;
  const auto* data=PB::GetParticleDataPointer(pointer);
  std::memcpy(&schema,data+SEP::Offset::TransportSchema,sizeof(schema));
  std::memcpy(&id,data+SEP::Offset::StableParticleId,sizeof(id));
  if(schema!=Transport::PICAdapter::ParticleTransportSchema||id==0)
    throw std::runtime_error("conservative merger requires initialized persisted particle history");
  return id;
}
int VelocityBin(double velocity,double maximum,int count) {
  // A zero-width velocity dimension contains one physical value and therefore
  // one populated bin. Other dimensions retain the configured bin count; no
  // spatial/momentum coarsening is performed to force a requested population.
  if(maximum==0.0)return 0;
  const double index=count*((velocity/maximum+1.0)*0.5);
  if(index<=0.0)return 0;
  if(index>=count)return count-1;
  return static_cast<int>(std::floor(index));
}
void ScaleDiagnostics(Conservative::UnpreservedMoments& moments,double weight) {
  moments.parallelCurrentCMPerS*=weight;
  moments.parallelStressJ*=weight;
  moments.perpendicularStressJ*=weight;
}
bool Finite(const Conservative::UnpreservedMoments& moments) {
  return std::isfinite(moments.parallelCurrentCMPerS)&&
      std::isfinite(moments.parallelStressJ)&&std::isfinite(moments.perpendicularStressJ);
}
void Publish(ConservativeMergeReport report) {
  const char* status=report.status==ConservativeMergeStatus::Reduced?"reduced":
      report.status==ConservativeMergeStatus::SparseNoReduction?"sparse-no-reduction":"numerical-rejection";
  std::ostringstream text;
  text<<std::setprecision(17)<<"SEP_CONSERVATIVE_MERGE line="<<report.fieldLine
      <<" species="<<report.species<<" segment="<<report.segment<<" status="<<status
      <<" before="<<report.before<<" after="<<report.after<<" target="<<report.target
      <<" sparse_bins="<<report.sparseBins<<" max_relative_residual="<<report.maximumRelativeResidual
      <<" numerical_stops="<<report.numericalStops
      <<" delta_current_C_m_per_s="<<report.afterUnpreserved.parallelCurrentCMPerS-report.beforeUnpreserved.parallelCurrentCMPerS
      <<" delta_parallel_stress_J="<<report.afterUnpreserved.parallelStressJ-report.beforeUnpreserved.parallelStressJ
      <<" delta_perpendicular_stress_J="<<report.afterUnpreserved.perpendicularStressJ-report.beforeUnpreserved.perpendicularStressJ
      <<" detail="<<report.detail<<'\n';
  std::cout<<text.str();reports.push_back(report);
}
void Merge(int species,long int& head,int spatialBins,int parallelBins,int normalBins,
           double spatialRange,int minimum,int maximum,int threshold) {
  if(species<0||species>=PIC::nTotalSpecies||minimum<0||maximum<minimum||maximum<1||
     spatialBins<1||parallelBins<1||normalBins<1||threshold<2||spatialRange!=1.0)
    throw std::runtime_error("conservative merger has an invalid existing population/bin contract");
  std::vector<long int> pointers;
  for(long int particle=head;particle!=-1;particle=PB::GetNext(particle))
    if(PB::GetI(particle)==species)pointers.push_back(particle);
  if(pointers.size()<static_cast<std::size_t>(maximum))return;
  if(pointers.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("conservative merger population exceeds controller count range");
  if(!passReady)throw std::runtime_error("conservative merger lacks a collective population pass");
  const int line=PB::GetFieldLineId(pointers.front());
  const double firstCoordinate=PB::GetFieldLineCoord(pointers.front());
  if(line<0||line>=FL::nFieldLine||!std::isfinite(firstCoordinate)||firstCoordinate<0||
     firstCoordinate>=FL::FieldLinesAll[line].GetTotalSegmentNumber())
    throw std::runtime_error("conservative merger particle has an invalid physical segment");
  const int segment=static_cast<int>(std::floor(firstCoordinate));
  auto* nativeSegment=FL::FieldLinesAll[line].GetSegment(segment);
  if(!nativeSegment||&head!=&nativeSegment->FirstParticleIndex||nativeSegment->Thread!=PIC::ThisThread)
    throw std::runtime_error("conservative merger list does not belong to this physical segment owner");
  const auto event=events.find(EventKey(line,species));
  if(event==events.end())throw std::runtime_error("conservative merger has no persisted physical source event");
  // Front insertion changes array indices without changing the surviving
  // material segment. Key the numerical population event by that segment's
  // stable label, not by its transient index or MPI owner. Conversion of a
  // negative signed label to uint64_t is defined modulo 2^64 and injective
  // over int64_t, so inward extensions retain distinct reproducible streams.
  const auto material=FL::FieldLinesAll[line].ExternallyManagedGeometry?
      FL::FieldLinesAll[line].FirstMaterialLabel+segment:segment;
  auto random=SourceLifecycle::Random(event->second,static_cast<std::uint64_t>(material),UINT64_C(0x4d45524745));
  const double tolerance=Run::Active().get().populationControl.relativeInvariantTolerance;
  std::sort(pointers.begin(),pointers.end(),[](long int a,long int b){return Identity(a)<Identity(b);});
  std::vector<Conservative::Particle> particles;
  double parallelMaximum=0,normalMaximum=0;
  for(long int pointer:pointers) {
    const double coordinate=PB::GetFieldLineCoord(pointer);
    if(PB::GetFieldLineId(pointer)!=line||!std::isfinite(coordinate)||
       std::floor(coordinate)!=segment)
      throw std::runtime_error("conservative merger cannot mix physical field-line segments");
    Conservative::Particle particle;
    particle.identity=Identity(pointer);particle.species=species;
    particle.weight=PB::GetIndividualStatWeightCorrection(pointer);
    particle.localFraction=coordinate-segment;
    particle.parallelMPerS=PB::GetVParallel(pointer);particle.normalMPerS=PB::GetVNormal(pointer);
    particle.massKg=PIC::MolecularData::GetMass(species);
    particle.chargeC=PIC::MolecularData::GetElectricCharge(species);
    if(ReducedParticles::Native::Enabled()) {
      // Anonymous native particle extensions migrate and split with the whole
      // record. Read the immutable source record; never manufacture a missing
      // birth energy from an evolved velocity or from a neighbouring particle.
      ReducedParticles::Native::SourceRecord source;
      std::memcpy(&source,PB::GetParticleDataPointer(pointer)+SEP::Offset::ReducedSourceRecord,sizeof(source));
      if(source.schema!=1)throw std::runtime_error("conservative reduced merger requires source birth-energy history");
      particle.hasBirthEnergy=true;particle.birthKineticEnergyJ=source.birthKineticEnergyJ;
    }
    parallelMaximum=std::max(parallelMaximum,std::fabs(particle.parallelMPerS));
    normalMaximum=std::max(normalMaximum,std::fabs(particle.normalMPerS));
    particles.push_back(particle);
  }
  std::vector<std::vector<double>> constraints;
  std::vector<double> before,after;
  const auto evaluation=Conservative::EvaluateConstraints(particles,SpeedOfLight,&constraints,&before);
  if(!evaluation.ok())throw std::runtime_error(evaluation.message);
  after=before;
  std::map<BinKey,std::vector<std::size_t>> bins;
  for(std::size_t i=0;i<particles.size();++i) {
    const auto& p=particles[i];
    const int spatial=std::min(spatialBins-1,static_cast<int>(std::floor(p.localFraction*spatialBins)));
    bins[BinKey(spatial,VelocityBin(p.parallelMPerS,parallelMaximum,parallelBins),
        VelocityBin(p.normalMPerS,normalMaximum,normalBins))].push_back(i);
  }
  ConservativeMergeReport report;
  report.fieldLine=line;report.species=species;report.segment=segment;
  report.before=report.after=static_cast<int>(pointers.size());
  report.target=minimum+(maximum-minimum)/2;
  const double baseWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
  if(!std::isfinite(baseWeight)||baseWeight<=0)
    throw std::runtime_error("conservative merger requires a finite positive physical base weight");
  report.beforeUnpreserved=Conservative::EvaluateUnpreserved(particles,before,SpeedOfLight);
  ScaleDiagnostics(report.beforeUnpreserved,baseWeight);
  if(!Finite(report.beforeUnpreserved))
    throw std::runtime_error("conservative merger input physical diagnostics overflow finite SI storage");
  report.afterUnpreserved=report.beforeUnpreserved;
  bool stopped=false;
  for(const auto& bin:bins) {
    std::vector<std::size_t> active=bin.second;
    while(report.after>report.target) {
      if(active.size()<static_cast<std::size_t>(threshold)) {++report.sparseBins;break;}
      std::vector<Conservative::Particle> candidates;
      for(std::size_t index:active) {
        auto p=particles[index];p.weight=after[index];candidates.push_back(p);
      }
      const auto endpoints=Conservative::Prepare(candidates,tolerance,SpeedOfLight);
      if(endpoints.weights.status==PIC::ConservativeResampling::Status::NoNullDirection) {
        ++report.sparseBins;break;
      }
      if(endpoints.weights.status!=PIC::ConservativeResampling::Status::Success) {
        // Stop BEFORE another random choice. Earlier validated choices form a
        // weight martingale and remain accepted. Selectively rolling them all
        // back on an outcome-dependent later failure would bias the ensemble.
        ++report.numericalStops;stopped=true;
        report.detail="stopped before invalid candidate: "+endpoints.weights.message;
        break;
      }
      // Preflight BOTH possible global endpoints against the original physical
      // moments, including cumulative rounding from earlier accepted deletions.
      // The stop decision cannot depend on which next outcome would be sampled.
      // A bounded number of such conditional-unbiased updates followed by this
      // predictable stop preserves E[w_after|complete initial records]=w_before.
      bool endpointsValid=true;
      for(const auto* endpoint:{&endpoints.weights.positive,&endpoints.weights.negative}) {
        auto candidate=after;
        for(std::size_t i=0;i<active.size();++i)candidate[active[i]]=(*endpoint)[i];
        const double residual=Conservative::MaximumRelativeResidual(constraints,before,candidate);
        auto diagnostic=Conservative::EvaluateUnpreserved(particles,candidate,SpeedOfLight);
        ScaleDiagnostics(diagnostic,baseWeight);
        endpointsValid=endpointsValid&&std::isfinite(residual)&&residual<=tolerance&&Finite(diagnostic);
      }
      if(!endpointsValid) {
        ++report.numericalStops;stopped=true;
        report.detail="stopped before candidate: both global endpoints must pass SI and finite-diagnostic gates";
        break;
      }
      std::vector<double> selected;
      const auto selection=Conservative::Select(endpoints,random.UniformOpen01(),&selected);
      if(!selection.ok())throw std::runtime_error(selection.message);
      std::vector<std::size_t> survivors;
      for(std::size_t i=0;i<active.size();++i) {
        after[active[i]]=selected[i];
        if(selected[i]>0)survivors.push_back(active[i]);
        else --report.after;
      }
      active.swap(survivors);
    }
    if(stopped||report.after<=report.target)break;
  }
  report.maximumRelativeResidual=Conservative::MaximumRelativeResidual(constraints,before,after);
  if(!std::isfinite(report.maximumRelativeResidual)||report.maximumRelativeResidual>tolerance) {
    throw std::runtime_error("internal conservative merger inconsistency: an accepted global endpoint failed its repeated SI gate");
  }
  report.afterUnpreserved=Conservative::EvaluateUnpreserved(particles,after,SpeedOfLight);
  ScaleDiagnostics(report.afterUnpreserved,baseWeight);
  if(!Finite(report.afterUnpreserved)) {
    throw std::runtime_error("internal conservative merger inconsistency: accepted diagnostics are nonfinite");
  }
  // Commit only weights/deletion. Every surviving phase-space and anonymous
  // extension byte (including ID, RNG/cache and residual hazard) stays intact.
  for(std::size_t i=0;i<pointers.size();++i) {
    if(after[i]==0.0)PB::DeleteParticle(pointers[i],head);
    else if(after[i]!=before[i])PB::SetIndividualStatWeightCorrection(after[i],pointers[i]);
  }
  report.status=report.after<report.before?ConservativeMergeStatus::Reduced:
      report.numericalStops?ConservativeMergeStatus::NumericalRejection:
                            ConservativeMergeStatus::SparseNoReduction;
  if(!stopped)report.detail=report.after>report.target?
      "existing bins retained; sparse/rank-limited remainder":
      "existing bins retained; requested target reached";
  Publish(report);
}
} // namespace
void RegisterConservativeMerger() {PIC::ParticleSplitting::FledLine::MergeController=Merge;}
void BeginConservativePopulationPass() {
  events.clear();reports.clear();passReady=false;
  int minimum=0,maximum=0;
  MPI_Allreduce(&FL::nFieldLine,&minimum,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&FL::nFieldLine,&maximum,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  const int localAvailable=FL::nFieldLine==0||FL::FieldLinesAll?1:0;
  int available=0;
  MPI_Allreduce(&localAvailable,&available,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(minimum!=maximum||minimum<0||!available)
    throw std::runtime_error("population scheduler requires identical physical line catalogs on all MPI owners");
  for(int line=0;line<FL::nFieldLine;++line) {
    const int local=FL::FieldLinesAll[line].IsInitialized()?1:0;
    int all=0,any=0;
    MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local,&any,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    if(all!=any)throw std::runtime_error("population scheduler physical line initialization differs across MPI ranks");
    if(all)for(int species=0;species<PIC::nTotalSpecies;++species)
      events.emplace(EventKey(line,species),SourceLifecycle::BeginNativeEvent(
          SourceLifecycle::Kind::PopulationControl,species,line));
  }
  passReady=true;
}
const std::vector<ConservativeMergeReport>& ConservativeMergeReports() {return reports;}
}} // namespace SEP::PopulationControl
