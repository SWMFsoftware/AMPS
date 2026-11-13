#ifndef MOON_PARTICLE_BUDGET_H_
#define MOON_PARTICLE_BUDGET_H_

#include "pic.h"
#include <filesystem>

namespace Moon { namespace ParticleBudget {
// Supporting I11 diagnostics for the current neutral-Na sphere/global-weight
// contract. Begin/Snapshot are collective; End disables only this observer.
// Snapshots compare actual list inventories with cumulative accepted events,
// independently of cell sampling, output cadence and source-rate estimates.
void Begin(const std::filesystem::path& directory);
void Snapshot(long long completedParticleSteps,const char* checkpoint);
void End();
bool Valid();
long long MaximumCountResidual();
long long Count(PIC::ParticleEvents::Event event);
long long GasCount();
} }
#endif
