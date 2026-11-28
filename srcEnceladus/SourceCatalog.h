#ifndef ENCELADUS_SOURCE_CATALOG_H
#define ENCELADUS_SOURCE_CATALOG_H
#include "CaseSelection.h"
#include <map>

namespace Enceladus { namespace ForwardModel {
// A prepared source run binds precisely the bytes that will be consumed.
// Absolute local paths are retained as provenance; the physical fingerprint
// uses logical bundle paths/content so identical data on two machines agree.
struct ResolvedRun {
  Configuration configuration;
  SourceDistribution distribution;
  std::string caseId,measurement,basis,publication,sourceSet,dataRoot;
  std::string canonicalInput,sourceBytes,fingerprint;
  std::string executableSha256;
  int mpiRanks=1;
  long completedIterations=0;
  std::map<std::string,std::string> assetBytes,assetHashes,assetLogicalPaths;
};
struct CaseInventory {
  std::string id,status,detail;
};
Status ResolveRun(const RuntimeSelection&,const std::string& configuredRoot,ResolvedRun*);
Status ListRuntimeCases(const RuntimeSelection&,const std::string& configuredRoot,std::vector<CaseInventory>*);
std::string DescribeRun(const ResolvedRun&);
Status WriteRunSnapshot(const ResolvedRun&,const std::string& directory,
                        const NumericalPlan* plan=nullptr);
} }
#endif
