/*
 * Native one-rank qualification of ElectricallyChargedDust::DustInjection__Sphere.
 *
 * The fixture intentionally builds a real 3-D AMR mesh, registers a spherical
 * internal boundary, allocates the native particle buffer, and lets the
 * production boundary-injection mover attach the new grain to its owning
 * cell.  A deterministic application callback supplies one physically valid
 * grain state so the acceptance oracle is mass/weight/ownership accounting,
 * not a statistical tolerance on a pseudorandom size distribution.
 */
#include "pic.h"
#include "constants.h"
#include "Exosphere.h"
#include "Dust.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

const double kSphereCentreM[3] = {0.20, -0.10, 0.05};
const double kSphereRadiusM = 0.40;
const double kGrainRadiusM = 2.0e-7;
const double kGrainSpeedMPerS = 0.10;
const double kIndividualWeightCorrection = 0.25;
const double kTimeStepS = 0.10;
const double kBasePhysicalWeight = 10.0;

double Resolution(double*) { return 0.25; }
double UnitLoad(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) { return 1.0; }

void ConfigureLeafBlocks(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (node == NULL) return;
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->block != NULL) {
      for (int species = 0; species < PIC::nTotalSpecies; ++species) {
        node->block->SetLocalTimeStep(kTimeStepS, species);
        node->block->SetLocalParticleWeight(kBasePhysicalWeight, species);
      }
    }
    return;
  }
  for (int child = 0; child < (1 << DIM); ++child)
    ConfigureLeafBlocks(node->downNode[child]);
}

bool DeterministicGrain(
    double* position_m, double* velocity_m_per_s,
    double& grain_radius_m, double& individual_weight_correction,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*& owner) {
  // Launch along a fixed outward unit normal.  The shared production helper
  // applies x=x_centre+R*n and v=v0*n in SI units.  A nonzero centre catches
  // the former origin-only DustInjection__Sphere implementation.
  const double outward[3] = {1.0, 0.0, 0.0};
  ElectricallyChargedDust::Core::SphericalSourceState(
      kSphereCentreM, kSphereRadiusM, kGrainSpeedMPerS, outward,
      position_m, velocity_m_per_s);
  grain_radius_m = kGrainRadiusM;
  individual_weight_correction = kIndividualWeightCorrection;
  owner = PIC::Mesh::mesh->findTreeNode(position_m, NULL);
  return owner != NULL && owner->block != NULL &&
      owner->Thread == PIC::Mesh::mesh->ThisThread;
}

long int FindOnlyParticle() {
  long int found = -1;
  for (long int particle = 0;
       particle < PIC::ParticleBuffer::MaxNPart; ++particle) {
    PIC::ParticleBuffer::byte* data =
        PIC::ParticleBuffer::GetParticleDataPointer(particle);
    if (PIC::ParticleBuffer::IsParticleAllocated(data)) {
      if (found != -1) return -2;
      found = particle;
    }
  }
  return found;
}

}  // namespace

int main(int, char**) {
  PIC::Init_BeforeParser();

  // Dust::Init_BeforeParser() constructs the power-law normalizer and
  // therefore requires the physical radius interval to be defined first.
  // These SI limits also bracket the deterministic 0.2 micrometre fixture.
  ElectricallyChargedDust::minDustRadius = 1.0e-7;
  ElectricallyChargedDust::maxDustRadius = 1.0e-6;
  ElectricallyChargedDust::SizeDistribution::PowerIndex = 4.0;
  ElectricallyChargedDust::Init_BeforeParser();

  // Register sampling storage before allocating blocks or a sphere.  The
  // sphere owns two sample windows for every physical element/species; dust
  // injection writes its represented mass flux into the collecting window.
  PIC::Mesh::initCellSamplingDataBuffer();
  cInternalSphericalData::SetGeneralSurfaceMeshParameters(4, 8);
  PIC::BC::InternalBoundary::Sphere::Init();
  cInternalBoundaryConditionsDescriptor descriptor =
      PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  cInternalSphericalData* sphere =
      static_cast<cInternalSphericalData*>(descriptor.BoundaryElement);
  double centre[3] = {kSphereCentreM[0], kSphereCentreM[1],
                      kSphereCentreM[2]};
  sphere->SetSphereGeometricalParameters(centre, kSphereRadiusM);
  sphere->localResolution = Resolution;

  double xmin[3] = {-1.0, -1.0, -1.0};
  double xmax[3] = { 1.0,  1.0,  1.0};
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::Mesh::mesh->init(xmin, xmax, Resolution);
  PIC::Mesh::mesh->SetNodeNeibResolutionLevelLimit();
  PIC::Mesh::mesh->CreateBottomBranchNodeList(PIC::Mesh::mesh->rootTree);
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  ConfigureLeafBlocks(PIC::Mesh::mesh->rootTree);
  PIC::DomainBlockDecomposition::UpdateBlockTable();

  PIC::ParticleBuffer::Init(128);
  for (int species = 0; species < PIC::nTotalSpecies; ++species) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] =
        kBasePhysicalWeight;
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species] = kTimeStepS;
  }

  ElectricallyChargedDust::MeanDustDensity = 1000.0;
  ElectricallyChargedDust::InitialGrainSpeed = kGrainSpeedMPerS;
  // Init_BeforeParser() has already initialized the size-law normalizer.  The
  // deterministic callback below deliberately bypasses stochastic radius
  // draws, so reinitializing it here would add no coverage and would emit a
  // diagnostic size-distribution file into the checkout.
  ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity = 0.01;
  ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity = 1.0;
  ElectricallyChargedDust::GrainVelocityGroup::nGroups = 1;
  ElectricallyChargedDust::GrainVelocityGroup::Init();
  ElectricallyChargedDust::GenerateNewDustGrainInternalProperties =
      DeterministicGrain;
  ElectricallyChargedDust::Injection::ResetMassReservoir();
  ElectricallyChargedDust::SampledDustMassInjectionRate = 0.0;

  const double grain_mass_kg =
      ElectricallyChargedDust::Core::SphericalGrainMass(
          kGrainRadiusM, ElectricallyChargedDust::MeanDustDensity);
  const double represented_mass_kg = grain_mass_kg *
      kBasePhysicalWeight * kIndividualWeightCorrection;
  // Two calls represent two equal physical time intervals.  Half a numerical
  // grain is requested in each interval.  The signed residual reservoir must
  // emit once on the first call and exactly repay its half-grain debt on the
  // second, closing requested and represented mass across 2*dt.
  ElectricallyChargedDust::TotalMassDustProductionRate =
      represented_mass_kg / (2.0 * kTimeStepS);

  std::cout << "[CASE] RUN  CORE-DUST-INJECTION-INTEGRATION-001 "
               "native spherical source, AMR owner and mass ledger"
            << std::endl;
  const long int first_injected =
      ElectricallyChargedDust::DustInjection__Sphere(
          _INTERNAL_BOUNDARY_TYPE_SPHERE_, sphere);
  const double first_residual =
      ElectricallyChargedDust::Injection::GetMassReservoir();
  const long int second_injected =
      ElectricallyChargedDust::DustInjection__Sphere(
          _INTERNAL_BOUNDARY_TYPE_SPHERE_, sphere);
  const double final_residual =
      ElectricallyChargedDust::Injection::GetMassReservoir();

  bool pass = first_injected == 1 && second_injected == 0 &&
      PIC::ParticleBuffer::GetAllPartNum() == 1;
  const long int particle = FindOnlyParticle();
  if (particle < 0) pass = false;
  if (particle >= 0) {
    PIC::ParticleBuffer::byte* data =
        PIC::ParticleBuffer::GetParticleDataPointer(particle);
    double x[3], velocity[3];
    PIC::ParticleBuffer::GetX(x, data);
    PIC::ParticleBuffer::GetV(velocity, data);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* owner =
        PIC::Mesh::mesh->findTreeNode(x, NULL);
    const double correction =
        PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
    pass = pass && owner != NULL && owner->block != NULL &&
        owner->Thread == PIC::Mesh::mesh->ThisThread &&
        PIC::ParticleBuffer::GetI(data) == _DUST_SPEC_ &&
        std::fabs(ElectricallyChargedDust::GetGrainRadius(data) -
                  kGrainRadiusM) <= 1.0e-15 * kGrainRadiusM &&
        std::fabs(ElectricallyChargedDust::GetGrainMass(data) -
                  grain_mass_kg) <= 1.0e-14 * grain_mass_kg &&
        std::fabs(correction - kIndividualWeightCorrection) <= 1.0e-15 &&
        velocity[0] > 0.0;
  }

  const double mass_tolerance = 64.0 *
      std::numeric_limits<double>::epsilon() * represented_mass_kg;
  pass = pass && first_residual < 0.0 &&
      std::fabs(final_residual) <= mass_tolerance &&
      std::fabs(ElectricallyChargedDust::SampledDustMassInjectionRate -
                represented_mass_kg / kTimeStepS) <=
          2.0e-14 * represented_mass_kg / kTimeStepS;

  std::cout << std::setprecision(17)
            << "dust_native first_injected=" << first_injected
            << " second_injected=" << second_injected
            << " represented_mass_kg=" << represented_mass_kg
            << " first_residual_kg=" << first_residual
            << " final_residual_kg=" << final_residual
            << " sampled_mass_rate_kg_s="
            << ElectricallyChargedDust::SampledDustMassInjectionRate
            << std::endl;
  std::cout << "[CORE-DUST-INJECTION-INTEGRATION-001] "
            << (pass ? "PASS" : "FAIL") << std::endl;

  MPI_Finalize();
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
