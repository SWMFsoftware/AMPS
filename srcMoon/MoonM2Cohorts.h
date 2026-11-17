#ifndef MOON_M2_COHORTS_H
#define MOON_M2_COHORTS_H
#include <filesystem>
#include <vector>
namespace Moon { namespace M2Cohorts {
// Allocate passive fraction fields before the real particle buffer exists.
// Fields are absent in normal builds/runs; enabled fields migrate with the
// complete production particle record rather than a rank-local pointer map.
void RegisterStorage(int cohorts);
void Begin(const std::vector<int>& elements,const std::vector<double>& initialAtoms,
           const std::filesystem::path& directory);
void EndStep(double elapsedS);
void Snapshot(double elapsedS);
void End();
} }
#endif
