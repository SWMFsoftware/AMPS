#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "models/chemical_solver/ChemistryDiagnostics.h"

namespace {

using ChemicalSolver::cChemistryDiagnostics;

bool Equal(double left,double right) {
  return std::fabs(left-right)<=1.0e-13*
      std::max(1.0,std::max(std::fabs(left),std::fabs(right)));
}

void RecordCanonicalEvent(cChemistryDiagnostics &diagnostics,int event) {
  const std::size_t reaction=static_cast<std::size_t>(event%4);
  diagnostics.RecordEligible(reaction);

  // Every fourth visit is deliberately a no-reaction visit.  This verifies
  // that eligible-parent exposure is distinct from a selected physical
  // reaction/channel.
  if (event%4==0 && event%8==0) return;

  const bool sink=(reaction==3);
  diagnostics.RecordSelection(reaction,sink);
  if (sink) return;

  const std::size_t productCount=(reaction==0) ? 1 : 2;
  for (std::size_t product=0;product<productCount;product++) {
    // Include zero, one, and more-than-one model multiplicities.  Reaction 2
    // has two separate terms for the same conceptual product species; the
    // diagnostic must not collapse them.
    const std::uint64_t multiplicity=
        static_cast<std::uint64_t>((event+static_cast<int>(product))%3);
    const std::uint64_t reused=multiplicity>0 && product==0 ? 1 : 0;
    const std::uint64_t allocated=multiplicity-reused;
    const double number=10.0*event+product+0.25;
    const double mass=2.0*number;
    const double momentum[3]={number,-2.0*number,3.0*number};
    diagnostics.RecordProduct(reaction,product,multiplicity,allocated,reused,
        number,mass,momentum,4.0*number,
        number,mass,momentum,4.0*number,0.5*number);
  }
}

int TestProductShapes() {
  cChemistryDiagnostics diagnostics;
  diagnostics.Configure({1,2,2,0},true);

  const double p0[3]={1.0,2.0,3.0};
  diagnostics.RecordEligible(0,2);
  diagnostics.RecordSelection(0,false);
  diagnostics.RecordProduct(0,0,1,0,1,4.0,8.0,p0,16.0,
      5.0,10.0,p0,20.0,0.5);

  diagnostics.RecordEligible(1);
  diagnostics.RecordSelection(1,false);
  diagnostics.RecordProduct(1,0,2,1,1,6.0,12.0,p0,24.0,
      7.0,14.0,p0,28.0,0.7);
  diagnostics.RecordProduct(1,1,0,0,0,1.0,2.0,p0,4.0,
      0.0,0.0,p0,0.0,0.1);

  diagnostics.RecordEligible(2);
  diagnostics.RecordSelection(2,false);
  diagnostics.RecordProduct(2,0,1,0,1,2.0,4.0,p0,8.0,
      3.0,6.0,p0,12.0,0.3);
  diagnostics.RecordProduct(2,1,1,1,0,2.0,4.0,p0,8.0,
      4.0,8.0,p0,16.0,0.4);

  diagnostics.RecordEligible(3);
  diagnostics.RecordSelection(3,true);
  diagnostics.RecordDeferral(1,2);
  diagnostics.RecordEmergency(1);

  const cChemistryDiagnostics::cReactionRecord &one=diagnostics.GetReaction(0);
  if (one.Eligible!=2 || one.Reactions!=1 || one.Products.size()!=1 ||
      one.Products[0].MultiplicityOne!=1 || one.Products[0].Reused!=1 ||
      !Equal(one.Products[0].Number,5.0)) return 10;

  const cChemistryDiagnostics::cReactionRecord &two=diagnostics.GetReaction(1);
  if (two.Products.size()!=2 || two.Products[0].MultiplicityMany!=1 ||
      two.Products[1].MultiplicityZero!=1 || two.Deferrals!=2 ||
      two.Emergencies!=1) return 11;

  const cChemistryDiagnostics::cReactionRecord &repeated=
      diagnostics.GetReaction(2);
  if (repeated.Products.size()!=2 ||
      !Equal(repeated.Products[0].Number,3.0) ||
      !Equal(repeated.Products[1].Number,4.0)) return 12;

  const cChemistryDiagnostics::cReactionRecord &sink=diagnostics.GetReaction(3);
  if (sink.Products.size()!=0 || sink.Sinks!=1 || sink.Destructions!=1)
    return 13;
  if (diagnostics.CountValues().size()!=4*7+5*6 ||
      diagnostics.MomentValues().size()!=5*13) return 14;

  bool rejectedInconsistentSlots=false;
  try {
    diagnostics.RecordProduct(0,0,2,0,0,0.0,0.0,p0,0.0,
        0.0,0.0,p0,0.0,0.0);
  }
  catch (const std::invalid_argument&) {
    rejectedInconsistentSlots=true;
  }
  if (!rejectedInconsistentSlots) return 15;
  return 0;
}

int TestDisabledContract() {
  cChemistryDiagnostics diagnostics;
  diagnostics.Configure({1,2,0},false);

  // Disabled methods return before index validation.  Instrumentation can be
  // compiled into a chemistry hot path without requiring dummy reaction
  // tables and without adding serialized fields to model output.
  const double momentum[3]={0.0,0.0,0.0};
  diagnostics.RecordEligible(99);
  diagnostics.RecordSelection(99,true);
  diagnostics.RecordProduct(99,99,2,1,1,1.0,2.0,momentum,3.0,
      1.0,2.0,momentum,3.0,0.0);
  diagnostics.RecordDeferral(99);
  diagnostics.RecordEmergency(99);
  if (diagnostics.IsEnabled() || diagnostics.GetReactionCount()!=0 ||
      !diagnostics.CountValues().empty() ||
      !diagnostics.MomentValues().empty()) return 20;
  return 0;
}

int TestCollectiveReduction(int rank,int size) {
  const int eventCount=240;
  cChemistryDiagnostics local,expected;
  local.Configure({1,2,2,0},true);
  expected.Configure({1,2,2,0},true);

  for (int event=0;event<eventCount;event++) {
    if (event%size==rank) RecordCanonicalEvent(local,event);
    RecordCanonicalEvent(expected,event);
  }

  std::vector<std::uint64_t> localCounts=local.CountValues();
  std::vector<std::uint64_t> reducedCounts(localCounts.size(),0);
  std::vector<double> localMoments=local.MomentValues();
  std::vector<double> reducedMoments(localMoments.size(),0.0);
  MPI_Allreduce(localCounts.data(),reducedCounts.data(),
      static_cast<int>(localCounts.size()),MPI_UINT64_T,MPI_SUM,
      MPI_COMM_WORLD);
  MPI_Allreduce(localMoments.data(),reducedMoments.data(),
      static_cast<int>(localMoments.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  const std::vector<std::uint64_t> expectedCounts=expected.CountValues();
  const std::vector<double> expectedMoments=expected.MomentValues();
  if (reducedCounts!=expectedCounts) return 30;
  if (reducedMoments.size()!=expectedMoments.size()) return 31;
  for (std::size_t value=0;value<reducedMoments.size();value++)
    if (!Equal(reducedMoments[value],expectedMoments[value])) return 32;

  // This checksum is intentionally independent of the number of ranks.  The
  // shell runner compares it across the required 1-, 4-, and 12-rank cases.
  if (rank==0) {
    std::uint64_t countChecksum=0;
    double momentChecksum=0.0;
    for (std::size_t value=0;value<reducedCounts.size();value++)
      countChecksum+=static_cast<std::uint64_t>(value+1)*reducedCounts[value];
    for (std::size_t value=0;value<reducedMoments.size();value++)
      momentChecksum+=(value+1)*reducedMoments[value];
    std::printf("stage0-mpi-signature counts=%llu moments=%.17e events=%d\n",
        static_cast<unsigned long long>(countChecksum),momentChecksum,
        eventCount);
  }
  return 0;
}

} // namespace

int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int result=TestProductShapes();
  if (result==0) result=TestDisabledContract();
  if (result==0) result=TestCollectiveReduction(rank,size);
  int globalResult=0;
  MPI_Allreduce(&result,&globalResult,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if (rank==0 && globalResult!=0)
    std::fprintf(stderr,"stage0 chemistry diagnostic test failed: %d\n",
        globalResult);
  MPI_Finalize();
  return globalResult;
}
