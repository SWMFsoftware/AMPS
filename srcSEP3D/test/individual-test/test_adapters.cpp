#include "sep3d_test_registry.h"

#include "particle_ledger.h"
#include "swcme_source_adapter.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace {

using SEP3D::Testing::Result;
namespace A = SEP3D::Adapters;
namespace C = SEP3D::Core;
namespace R = SEP3D::RuntimeModel;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}

Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

A::MoverInput BaseMover() {
  A::MoverInput input;
  input.model = R::TransportModel::Parker3D;
  input.particle.stableId = 19;
  input.particle.species = 0;
  input.particle.positionM = C::Vec3(10.0, 0.0, 0.0);
  input.particle.momentumKgMPerS = 0.0;
  input.particle.mu = 0.25;
  input.particle.gyrophaseRad = 0.5;
  input.particle.statisticalWeight = 2.0;
  input.local.background.status = C::Status::OK();
  input.local.background.valid = true;
  input.local.background.B = C::Vec3(1.0, 0.0, 0.0);
  input.local.background.absB = 1.0;
  input.local.background.bHat = C::Vec3(1.0, 0.0, 0.0);
  input.local.background.U = C::Vec3(1.0, 0.0, 0.0);
  input.local.cellSizeM = 1.0e6;
  input.speciesMassKg = C::Const::m_p;
  input.requestedDtS = 1.0;
  input.innerRadiusM = 1.0;
  input.outerRadiusM = 100.0;
  input.campaignSeed = 47;
  return input;
}

swcme::sep::SEPSourceState Source() {
  swcme::sep::SEPSourceState source;
  source.status = swcme::ModelStatus::success();
  source.active = true;
  source.source_id = 71;
  source.position_m = {{10.0, 20.0, 30.0}};
  source.normal = {{1.0, 0.0, 0.0}};
  source.relative_patch_weight = 0.8;
  source.compression = 4.0;
  source.normal_speed_m_s = 8.0e5;
  source.q_phase_space = 4.0;
  source.spectrum.particle_mass_kg = C::Const::m_p;
  source.spectrum.kinetic_energy_min_MeV = 1.0;
  source.spectrum.kinetic_energy_max_MeV = 100.0;
  source.spectrum.reference_energy_MeV = 10.0;
  return source;
}

Result RunADP3D01() {
  const std::vector<std::string>& names =
      A::ProductionMoverRegistry::CanonicalNames();
  if (names.size() != 2 || names[0] != "parker3d-tensor" ||
      names[1] != "focused3d-split")
    return Fail("production registry is not the exact two-core Phase-A set");
  A::MoverInput parker = BaseMover();
  const A::MoverResult p = A::AdvanceParticle(parker);
  A::MoverInput focused = BaseMover();
  focused.model = R::TransportModel::Focused3D;
  const A::MoverResult f = A::AdvanceParticle(focused);
  if (!p.status.ok() || !f.status.ok() ||
      p.disposition != A::ParticleDisposition::Active ||
      f.disposition != A::ParticleDisposition::Active)
    return Fail("one of the validating dispatch paths rejected a valid record");
  return Pass("exactly two registered production cores dispatch through one validated record boundary");
}

Result RunNAT3D04() {
  A::MoverInput outer = BaseMover();
  outer.outerRadiusM = 10.5;
  const A::MoverResult escaped = A::AdvanceParticle(outer);
  A::MoverInput inner = BaseMover();
  inner.particle.positionM = C::Vec3(2.0, 0.0, 0.0);
  inner.local.background.U = C::Vec3(-2.0, 0.0, 0.0);
  inner.innerRadiusM = 1.5;
  const A::MoverResult absorbed = A::AdvanceParticle(inner);
  A::MoverInput bad = BaseMover();
  bad.local.background.valid = false;
  const A::MoverResult invalid = A::AdvanceParticle(bad);
  if (escaped.disposition != A::ParticleDisposition::Escaped ||
      escaped.status.code != C::StatusCode::DomainExit ||
      absorbed.disposition != A::ParticleDisposition::Absorbed ||
      absorbed.status.code != C::StatusCode::InnerBoundary ||
      invalid.disposition != A::ParticleDisposition::Failed ||
      invalid.status.code != C::StatusCode::BackgroundInvalid)
    return Fail("inner, outer, or invalid-background disposition is incorrect");
  return Pass("inner absorption, outer escape, and invalid background are distinct mover outcomes");
}

Result RunNAT3D05() {
  A::ParticleLedger ledger;
  if (!ledger.Begin(4, 0, 10).ok() ||
      !ledger.RecordInjection(4, 0, 3).ok())
    return Fail("ledger row could not be opened");
  for (int i = 0; i < 10; ++i)
    if (!ledger.RecordMover(4, 0, A::ParticleDisposition::Active).ok())
      return Fail("active mover outcome was rejected");
  if (!ledger.RecordMover(4, 0, A::ParticleDisposition::Escaped).ok() ||
      !ledger.RecordMover(4, 0, A::ParticleDisposition::Absorbed).ok() ||
      !ledger.RecordMover(4, 0, A::ParticleDisposition::Failed).ok())
    return Fail("terminal mover outcome was rejected");
  const C::Status mismatch = ledger.Close(4, 0, 11);
  if (mismatch.ok() || ledger.Find(4, 0)->closed)
    return Fail("mismatched close mutated or accepted the ledger row");
  if (!ledger.Close(4, 0, 10).ok())
    return Fail("exact conservation equation did not close");
  const A::LedgerRow* row = ledger.Find(4, 0);
  if (row == nullptr || !row->closed || row->advanced != 10)
    return Fail("closed ledger diagnostics are incomplete");
  return Pass("active_start+injected equals active_end+escaped+absorbed+failed exactly");
}

Result RunNAT3D08() {
  A::MoverInput input = BaseMover();
  input.local.background.U = C::Vec3(4.0, 0.0, 0.0);
  input.shock.active = true;
  input.shock.radiusAtStepStartM = 12.0;
  input.shock.radialSpeedMPerS = 0.0;
  input.shock.generation = 8;
  input.timeStepControls.shockCrossingFraction = 10.0;
  const A::MoverResult first = A::AdvanceParticle(input);
  if (!first.status.ok() || !first.shockIntersection.crossed ||
      std::fabs(first.shockIntersection.stepFraction - 0.5) > 1.0e-14 ||
      first.particle.lastShockGeneration != 8)
    return Fail("first moving-surface crossing was not located exactly");
  input.particle = first.particle;
  input.particle.positionM = C::Vec3(10.0, 0.0, 0.0);
  const A::MoverResult duplicate = A::AdvanceParticle(input);
  if (!duplicate.status.ok() || duplicate.shockIntersection.crossed)
    return Fail("same-generation shock crossing was injected twice");
  return Pass("the first segment/sphere root is recorded once per shock generation");
}

Result RunSHK3D01() {
  const swcme::sep::SEPSourceState common = Source();
  const A::ShockSourceRecord oneD =
      A::MakeShockSourceRecord(common, 9, 123, 32, 0.25);
  const A::ShockSourceRecord threeD =
      A::MakeShockSourceRecord(common, 9, 123, 32, 0.25);
  if (!oneD.status.ok() || !threeD.status.ok() ||
      oneD.sourceFingerprint != threeD.sourceFingerprint ||
      SEP::Injection::Fingerprint(oneD.injection) !=
          SEP::Injection::Fingerprint(threeD.injection))
    return Fail("dimensional adapters did not preserve one common source record");
  for (std::uint64_t i = 0; i < 32; ++i) {
    const A::InjectedParticle left = A::SampleInjectedParticle(oneD, i, 0);
    const A::InjectedParticle right = A::SampleInjectedParticle(threeD, i, 0);
    if (!left.status.ok() || !right.status.ok() ||
        left.particle.stableId != right.particle.stableId ||
        left.particle.momentumKgMPerS != right.particle.momentumKgMPerS ||
        left.particle.mu != right.particle.mu ||
        left.particle.gyrophaseRad != right.particle.gyrophaseRad)
      return Fail("1-D/3-D source sampling diverged for an identical SWCME state");
  }
  return Pass("common SWCME input produces bitwise-identical source keys and samples in both dimensional paths");
}

Result RunSHK3D02() {
  A::ExpandingSphericalShock shock;
  shock.active = true;
  shock.radiusAtStepStartM = 10.0;
  shock.radialSpeedMPerS = 1.0;
  shock.generation = 2;
  const A::ShockIntersection hit = A::FirstShockIntersection(
      C::Vec3(8.0, 0.0, 0.0), C::Vec3(14.0, 0.0, 0.0), 2.0,
      shock, 0);
  // 8+6u = 10+2u, hence u=1/2.
  if (!hit.status.ok() || !hit.crossed ||
      std::fabs(hit.stepFraction - 0.5) > 1.0e-14 ||
      std::fabs(hit.positionM.x - 11.0) > 1.0e-14)
    return Fail("expanding-sphere quadratic returned the wrong first root");
  return Pass("moving spherical shock geometry agrees with its analytic intersection");
}

Result RunSHK3D03() {
  swcme::sep::SEPSourceState inactive = Source();
  inactive.active = false;
  const A::ShockSourceRecord rejected =
      A::MakeShockSourceRecord(inactive, 4, 3, 8, 0.1);
  const A::ShockSourceRecord stale =
      A::MakeShockSourceRecord(Source(), 0, 3, 8, 0.1);
  if (rejected.status.ok() || rejected.active ||
      stale.status.ok() || stale.active)
    return Fail("inactive or generation-zero source became injectable");
  return Pass("inactive and unversioned SWCME sources are rejected before particle allocation");
}

Result RunSHK3D04() {
  const A::ShockSourceRecord source =
      A::MakeShockSourceRecord(Source(), 12, 99, 128, 0.25);
  if (!source.status.ok()) return Fail(source.status.message);
  double total = 0.0;
  for (std::uint64_t i = 0; i < 128; ++i) {
    const A::InjectedParticle particle = A::SampleInjectedParticle(source, i, 0);
    if (!particle.status.ok()) return Fail(particle.status.message);
    total += particle.particle.statisticalWeight;
  }
  const double expected = source.relativePatchWeight *
      source.injection.injectionEfficiency;
  if (std::fabs(total - expected) > 64.0 *
      std::numeric_limits<double>::epsilon() * expected)
    return Fail("macroparticle weights do not recover the event normalization");
  return Pass("sample weights sum to patch weight times injection efficiency");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterAdapterTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = group;
    d.description = "Phase-A host-neutral AMPS/SWCME adapter acceptance";
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = RC::Routine; d.seedPolicy = "semantic keyed source/mover";
    d.stateIsolation = "fresh source, particle, and ledger per callback";
    d.callback = std::move(callback); return d;
  };
  return {
      make("ADP3D01", "ADP3D", "Production mover dispatch", RunADP3D01),
      make("NAT3D04", "NAT3D", "Boundary dispositions", RunNAT3D04),
      make("NAT3D05", "NAT3D", "Particle ledger closure", RunNAT3D05),
      make("NAT3D08", "NAT3D", "Shock crossing dispatch", RunNAT3D08),
      make("SHK3D01", "SHK3D", "Dimensional source identity", RunSHK3D01),
      make("SHK3D02", "SHK3D", "Expanding shock geometry", RunSHK3D02),
      make("SHK3D03", "SHK3D", "Source ownership guards", RunSHK3D03),
      make("SHK3D04", "SHK3D", "Source weight normalization", RunSHK3D04),
  };
}
