#ifndef SRCMOON_TEST_COLLECTIVE_SOURCE_ADMISSION_H
#define SRCMOON_TEST_COLLECTIVE_SOURCE_ADMISSION_H
#include <mpi.h>
#include <iostream>
#include <string>
namespace Moon { namespace SodiumCohortInput {
struct Admission {int severity=0;std::string diagnostic;};
// Severity0=available,1=authoritative prerequisite unavailable,2=ERROR.
// All ranks agree BEFORE anybody enters amps_init or finalizes MPI. Rank is
// diagnostic metadata only and never enters a physical RNG/parent key.
inline Admission CollectiveAdmission(MPI_Comm communicator,int local,
    const std::string& localDiagnostic,const std::string& phase){
 int rank=0;MPI_Comm_rank(communicator,&rank);
 if(local<0||local>2)local=2;
 if(local)std::cerr<<"cohort "<<phase<<" rank="<<rank<<" severity="<<local<<": "<<localDiagnostic<<'\n';
 Admission result;MPI_Allreduce(&local,&result.severity,1,MPI_INT,MPI_MAX,communicator);
 if(result.severity){
  if(local==result.severity&&!localDiagnostic.empty())result.diagnostic=localDiagnostic;
  else result.diagnostic=phase+" admission refused on another rank ("+(result.severity==2?"ERROR":"missing prerequisite")+"); inspect retained per-rank native diagnostics. No rank entered physical startup.";
 }
 return result;
}
}}
#endif
