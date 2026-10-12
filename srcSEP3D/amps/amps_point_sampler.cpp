// AMPS binding of srcSEP3D point sampling.  The sampling cycle, interval
// bookkeeping, ownership and MPI contract are documented in
// amps_point_sampler.h; the physics and the product in
// ../output/point_sampling.h.

#include "amps_point_sampler.h"

#include "amps_particle_adapter.h"
#include "../output/point_sampling.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace SEP3D {
namespace AMPS {
namespace PointSampler {
namespace {

namespace PS = Output::PointSampling;

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

// Stop every rank with the sampler's diagnosis.  AMPS' exit() reports the
// source location and terminates the MPI job.
[[noreturn]] void Stop(const std::string& message) {
  const std::string text = "srcSEP3D point sampling: " + message;
  ::exit(__LINE__, __FILE__, text.c_str());
  std::abort();  // not reached; AMPS exit() does not return
}

// Complete sampler state of this process.
struct State {
  bool installed = false;
  std::vector<PS::PointDefinition> points;
  std::vector<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> nodes;  // per point
  std::vector<PS::SpeciesChannels> species;                   // per slot
  std::vector<int> slotOfAmpsSpecies;  // AMPS species index -> slot, -1 = not sampled
  PS::ProductSelection products;
  PS::Layout layout;
  std::vector<double> sums;        // this rank's interval sums
  std::vector<double> cumulative;  // cumulative fluence (meaningful on rank 0)
  std::string filePrefix;
  CellBackgroundResolver resolver = nullptr;
  // Interval bookkeeping (identical on every rank: all ranks call Sample()).
  long sampleCalls = 0;
  double firstSampleTimeS = 0.0;
  double lastSampleTimeS = 0.0;
  bool intervalPrinted = false;
};

State gState;

// Reduce the interval sums to rank 0 (collective over all ranks).
std::vector<double> ReduceToRoot() {
  std::vector<double> global(gState.sums.size(), 0.0);
  MPI_Reduce(gState.sums.data(), global.data(),
             static_cast<int>(gState.sums.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  return global;
}

// Start a new AMPS sampling interval.  If the previous interval ended
// without a Print() (AMPS suppressed that output), its fluence is still
// added to the cumulative fluence so the time integral has no gap.
void BeginInterval() {
  if (gState.sampleCalls > 0 && !gState.intervalPrinted) {
    const std::vector<double> global = ReduceToRoot();
    if (PIC::ThisThread == 0)
      PS::FoldFluence(gState.layout, gState.species, global, &gState.cumulative);
  }
  std::fill(gState.sums.begin(), gState.sums.end(), 0.0);
  gState.sampleCalls = 0;
  gState.intervalPrinted = false;
}

// Sample the macroparticles of the cell containing point 'p' (owner only).
void SamplePoint(std::size_t p) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = gState.nodes[p];
  if (node->Thread != PIC::ThisThread || node->block == nullptr) return;
  const PS::PointDefinition& point = gState.points[p];

  double x[3];
  point.positionM.CopyTo(x);
  int i = 0, j = 0, k = 0;
  const long cell = PIC::Mesh::mesh->FindCellIndex(x, i, j, k, node, false);
  if (cell < 0) Stop("point '" + point.name + "' has no AMR cell in its block");
  PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(cell);
  const double volume = center == nullptr ? -1.0 : center->Measure;
  if (!std::isfinite(volume) || volume <= 0.0)
    Stop("point '" + point.name + "' lies in a cell without a positive volume");

  Core::Vec3 bulkVelocity, magneticDirection;
  const Core::Status background =
      gState.resolver(point.positionM, node, &bulkVelocity, &magneticDirection);
  if (!background.ok())
    Stop("point '" + point.name + "': cell background unavailable: " +
         background.message);

  long ptr = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
  while (ptr != -1) {
    Adapters::ParticleRecord record;
    const Core::Status read = Movers::ReadParticle(ptr, &record);
    if (!read.ok())
      Stop("point '" + point.name + "': particle read failed: " + read.message);
    const int slot = record.species >= 0 &&
        record.species < static_cast<int>(gState.slotOfAmpsSpecies.size())
        ? gState.slotOfAmpsSpecies[record.species] : -1;
    if (slot >= 0) {
      const PS::SpeciesChannels& species = gState.species[slot];
      // The AMPS velocity is the plasma-frame velocity vector the mover
      // rebuilt from p, mu and the gyrophase at the particle's local field;
      // only its direction is used (|v| follows from p).
      double v[3];
      PIC::ParticleBuffer::GetV(v, PIC::ParticleBuffer::GetParticleDataPointer(ptr));
      PS::ParticleInput input;
      input.weight = record.statisticalWeight;
      input.massKg = species.massKg;
      input.momentumKgMPerS = record.momentumKgMPerS;
      input.mu = record.mu;
      input.plasmaVelocityDirection = Core::Vec3(v);
      PS::ParticleKinematics kinematics;
      const Core::Status evaluated = PS::EvaluateKinematics(
          input, bulkVelocity, magneticDirection, &kinematics);
      if (!evaluated.ok())
        Stop("point '" + point.name + "': " + evaluated.message);
      const double dt =
          PIC::ParticleWeightTimeStep::GlobalTimeStep[record.species];
      PS::Accumulate(gState.layout, species, p, static_cast<std::size_t>(slot),
                     input.weight, volume, dt, kinematics, &gState.sums);
    }
    ptr = PIC::ParticleBuffer::GetNext(ptr);
  }
}

}  // namespace

Core::Status Install(
    const RuntimeModel::PointSamplingOptions& options,
    const std::vector<RuntimeModel::CompiledSpeciesRecord>& compiled,
    const std::string& filePrefix, CellBackgroundResolver resolver) {
  if (!options.enabled) return Core::Status::OK();
  if (gState.installed) return Invalid("point sampling is already installed");
  if (resolver == nullptr) return Invalid("point sampling needs a background resolver");

  State state;
  state.filePrefix = filePrefix;
  state.resolver = resolver;
  const Core::Status products = PS::ParseProducts(options.products, &state.products);
  if (!products.ok()) return products;

  // Species: bind each configured symbol to the compiled AMPS table.
  state.slotOfAmpsSpecies.assign(static_cast<std::size_t>(PIC::nTotalSpecies), -1);
  for (const RuntimeModel::PointSamplingSpeciesOptions& configured : options.species) {
    const RuntimeModel::CompiledSpeciesRecord* match = nullptr;
    std::string available;
    for (const RuntimeModel::CompiledSpeciesRecord& record : compiled) {
      available += (available.empty() ? "" : ", ") + record.symbol;
      if (record.symbol == configured.symbol) match = &record;
    }
    if (match == nullptr)
      return Invalid("sampling species '" + configured.symbol +
                     "' is not a compiled AMPS species (compiled: " +
                     available + ")");
    PS::SpeciesChannels channels;
    channels.ampsIndex = match->ampsIndex;
    channels.symbol = match->symbol;
    channels.massKg = match->massKg;
    channels.scale = configured.energyScale == RuntimeModel::EnergyChannelSpacing::Linear
        ? PS::EnergyScale::Linear : PS::EnergyScale::Logarithmic;
    const Core::Status grid = PS::BuildEnergyEdges(
        configured.minimumEnergyMeV, configured.maximumEnergyMeV,
        configured.energyIntervals, channels.scale, &channels.edgesMeV);
    if (!grid.ok())
      return Invalid("sampling species '" + configured.symbol + "': " + grid.message);
    state.slotOfAmpsSpecies[static_cast<std::size_t>(match->ampsIndex)] =
        static_cast<int>(state.species.size());
    state.species.push_back(channels);
  }

  // Points: bind each to its AMR tree node.  The tree is replicated on every
  // rank, so every rank reaches the same verdict without communication.
  for (const RuntimeModel::PointSamplingPointOptions& configured : options.points) {
    PS::PointDefinition point;
    point.name = configured.name;
    point.positionM = configured.positionM;
    point.origin = configured.origin;
    double x[3];
    point.positionM.CopyTo(x);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = PIC::Mesh::mesh->findTreeNode(x);
    const std::string where = point.origin.empty() ? "" : point.origin + ": ";
    if (node == nullptr)
      return Invalid(where + "sampling point '" + point.name +
                     "' lies outside the AMPS domain");
    if (!node->IsUsedInCalculationFlag)
      return Invalid(where + "sampling point '" + point.name +
                     "' lies in an inactive (unallocated) mesh region");
    // The cell index follows from the node geometry alone, so every rank
    // records the sampled cell (written into the product) without owning it.
    int i = 0, j = 0, k = 0;
    if (PIC::Mesh::mesh->FindCellIndex(x, i, j, k, node, false) < 0)
      return Invalid(where + "sampling point '" + point.name +
                     "' has no AMR cell in its block");
    point.cellSizeM = Core::Vec3((node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_,
                                 (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_,
                                 (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_);
    point.cellCenterM = Core::Vec3(node->xmin[0] + (i + 0.5) * point.cellSizeM.x,
                                   node->xmin[1] + (j + 0.5) * point.cellSizeM.y,
                                   node->xmin[2] + (k + 0.5) * point.cellSizeM.z);
    point.cellKnown = true;
    state.points.push_back(point);
    state.nodes.push_back(node);
  }

  state.layout = PS::Layout(state.points.size(), state.species);
  state.sums.assign(state.layout.size(), 0.0);
  state.cumulative.assign(state.layout.size(), 0.0);
  state.installed = true;
  gState = std::move(state);
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sample, Print);
  return Core::Status::OK();
}

bool Installed() { return gState.installed; }

void Sample() {
  if (!gState.installed) return;
  if (PIC::CollectingSampleCounter == 0) BeginInterval();
  const double time = PIC::SimulationTime::Get();
  if (gState.sampleCalls == 0) gState.firstSampleTimeS = time;
  gState.lastSampleTimeS = time;
  ++gState.sampleCalls;
  for (std::size_t p = 0; p < gState.points.size(); ++p) SamplePoint(p);
}

void Print(int dataOutputFileNumber) {
  if (!gState.installed) return;
  const std::vector<double> global = ReduceToRoot();
  gState.intervalPrinted = true;
  if (PIC::ThisThread != 0) return;

  // The AMPS sampling length of this output must be exactly the number of
  // sampling calls this module received; anything else means the interval
  // bookkeeping and the AMPS cycle disagree.
  if (gState.sampleCalls != PIC::LastSampleLength)
    Stop("interval has " + std::to_string(gState.sampleCalls) +
         " sampling calls but PIC::LastSampleLength = " +
         std::to_string(PIC::LastSampleLength));
  if (gState.sampleCalls <= 0) return;

  PS::FoldFluence(gState.layout, gState.species, global, &gState.cumulative);
  PS::IntervalInfo info;
  info.outputFileNumber = dataOutputFileNumber;
  info.sampleCalls = gState.sampleCalls;
  info.firstSampleTimeS = gState.firstSampleTimeS;
  info.lastSampleTimeS = gState.lastSampleTimeS;
  const std::string text = PS::FormatTecplot(
      gState.layout, gState.points, gState.species, gState.products, info,
      global, gState.cumulative);

  const std::string path = std::string(PIC::OutputDataFileDirectory) + "/" +
      gState.filePrefix + ".point-sampling.out=" +
      std::to_string(dataOutputFileNumber) + ".dat";
  const std::string temporary = path + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    output << text;
    output.flush();
    if (!output) Stop("cannot write '" + temporary + "'");
  }
  if (std::rename(temporary.c_str(), path.c_str()) != 0)
    Stop("cannot rename '" + temporary + "' to '" + path + "'");
  std::printf("$PREFIX: srcSEP3D point sampling: wrote %s (%ld sampling calls)\n",
              path.c_str(), gState.sampleCalls);
}

}  // namespace PointSampler
}  // namespace AMPS
}  // namespace SEP3D
