#include "../adapters/swcme3d_spherical_material_adapter.h"
#include "spherical_shock_only_fixture.h"
#include "sep_test_registry.h"

#include <cstddef>
#include <stdexcept>

namespace SEP::Testing {

Result VerifySphericalMaterialAdapter() {
  namespace M = Background::Material;
  namespace S = CoronaSwcme::SphericalShockOnly;
  std::size_t assertions = 0;
  const auto expect = [&](bool check, const char* message) {
    ++assertions;
    if (!check) throw std::runtime_error(message);
  };
  Result result;
  try {
    auto definition = S::Verification::Configuration(-1);
    definition.front.V0_m_s = 1.5e6;
    auto event = S::Event::Create(definition);
    expect(event.ok(), event.status.message.c_str());
    const auto controls = S::Verification::IntegrationControls();
    const auto view = SW3DSphericalAdapter::Capture(event.value, controls, 0, 20, 19, 0);
    expect(view.ok(), view.status.message.c_str());
    const auto& interval = *view.value;
    const auto& descriptor = interval.Identity();
    expect(descriptor.geometryGenerationPresent && descriptor.geometryGeneration == 0 &&
           descriptor.eventIdentity == event.value->Identity() && descriptor.backgroundGeneration == 19,
           "detached interval lost supplied lineage/physical identity");
    expect(M::Require(descriptor.profile, M::Use::SnapshotExchange).ok() &&
           interval.RequireStep(0, 20, 20, M::Region::ParkerExterior).code ==
               Core::StatusCode::UnsupportedCapability,
           "kinematic candidate granted native particle tracking");
    const std::string capturedIdentity = M::Payload(descriptor);

    // An immutable metadata object is not an installed PIC material owner.
    // Exercise the real binding boundary with a legitimate model-owned SWCME
    // snapshot, without publishing it into the process-wide native store.
    // Snapshot exchange and local motion queries must not grant this still
    // unqualified profile source, particle-step or downstream-volume authority.
    const Background::BackgroundSnapshot detachedSnapshot(
        Background::Provider::Swcme, Background::Ownership::ModelOwned,
        0, 0, 20, 19, descriptor.eventIdentity,
        "S13AD01 detached immutable snapshot metadata; no native publication");
    expect(interval.BindSnapshot(detachedSnapshot).code ==
               Core::StatusCode::UnsupportedCapability &&
           M::Require(descriptor.profile, M::Use::PhysicalDownstream,
                      M::Region::Sheath).code == Core::StatusCode::UnsupportedCapability,
           "detached metadata or an RH limit granted native/downstream admission");

    // Translation across the bridge must preserve the canonical primitive,
    // magnetic polarity and both derivative tensors to machine equality. The
    // core suite independently verifies the physical formulas.
    M::Query query;
    query.positionM = {0, 0, S::Verification::Au};
    query.epochS = 0;
    query.region = M::Region::ParkerExterior;
    query.required = M::Bit(M::Capability::Primitives) | M::Bit(M::Capability::Species) |
        M::Bit(M::Capability::Thermodynamics) | M::Bit(M::Capability::VelocityGradient) |
        M::Bit(M::Capability::MagneticGradient);
    const auto departureQuery = query;
    const auto state = interval.Sample(query);
    const auto direct = event.value->Evaluate({0, 0, S::Verification::Au}, 0);
    expect(state.ok(), state.status.message.c_str());
    expect(direct.ok(), direct.status.message.c_str());
    expect(state.value.magneticFieldT[2] < 0 &&
           *state.value.massDensityKgM3 == direct.value.ambient.primitive.plasma.massDensityKgM3 &&
           state.value.velocityGradientPerS->values == direct.value.ambient.gradientU &&
           state.value.magneticGradientTPerM->values == direct.value.ambient.gradientB,
           "detached adapter changed signs, density or a Cartesian derivative");
    expect(state.value.species->size() == 3 && state.value.species->back().identity == "alpha" &&
           state.value.species->back().temperatureK == definition.ambient.composition.alphaTemperatureK,
           "detached adapter lost configured alpha density/temperature");
    const auto& primitive = direct.value.ambient.primitive;
    expect(state.value.velocityMPerS == M::Vector{primitive.velocityMPerS.x,
               primitive.velocityMPerS.y, primitive.velocityMPerS.z} &&
           state.value.magneticFieldT == M::Vector{primitive.magneticFieldT.x,
               primitive.magneticFieldT.y, primitive.magneticFieldT.z} &&
           *state.value.pressurePa == primitive.plasma.pressurePa &&
           *state.value.gammaAdiabatic == definition.ambient.composition.gammaAdiabatic &&
           state.value.species->at(0).numberDensityM3 == primitive.plasma.electronNumberDensityM3 &&
           state.value.species->at(0).temperatureK == primitive.electronTemperatureK &&
           state.value.species->at(1).numberDensityM3 == primitive.plasma.protonNumberDensityM3 &&
           state.value.species->at(1).temperatureK == primitive.protonTemperatureK &&
           state.value.species->at(2).numberDensityM3 ==
               definition.ambient.composition.alphaToProtonNumberRatio * primitive.plasma.protonNumberDensityM3,
           "detached adapter changed the canonical velocity, pressure/EOS or species composition");

    // Capturing a second generation of the SAME physical event must preserve
    // each interval's supplied native lineage. A process-global owner/cache
    // would silently substitute the later geometry generation into old reads.
    // Exact equality is appropriate: this is value translation, not a second
    // numerical integration or an observational-error comparison.
    const auto sibling = SW3DSphericalAdapter::Capture(event.value, controls, 0, 20, 20, 7);
    expect(sibling.ok(), sibling.status.message.c_str());
    const auto siblingState = sibling.value->Sample(query);
    expect(siblingState.ok(), siblingState.status.message.c_str());
    expect(siblingState.value.backgroundGeneration == 20 &&
           siblingState.value.geometryGeneration == 7 &&
           siblingState.value.eventIdentity == state.value.eventIdentity &&
           siblingState.value.velocityMPerS == state.value.velocityMPerS &&
           M::Payload(interval.Identity()) == capturedIdentity,
           "another detached capture replaced physical identity or native generation lineage");

    // Concurrent owners with opposite DECLARED magnetic sectors must not
    // alias a shared last-evaluated field. Reuse the existing fixture for the
    // opposite polarity and compare the bridge with that actual core sample;
    // this does not invent a new background or infer a sign from |B|.
    auto oppositeDefinition = S::Verification::Configuration(1);
    oppositeDefinition.front.V0_m_s = definition.front.V0_m_s;
    const auto oppositeEvent = S::Event::Create(oppositeDefinition);
    expect(oppositeEvent.ok(), oppositeEvent.status.message.c_str());
    const auto oppositeView = SW3DSphericalAdapter::Capture(oppositeEvent.value, controls, 0, 20, 19, 0);
    expect(oppositeView.ok(), oppositeView.status.message.c_str());
    const auto oppositeState = oppositeView.value->Sample(departureQuery);
    const auto oppositeDirect = oppositeEvent.value->Evaluate({0, 0, S::Verification::Au}, 0);
    expect(oppositeState.ok() && oppositeDirect.ok(), "opposite declared sector is unavailable");
    expect(oppositeState.value.magneticFieldT == M::Vector{oppositeDirect.value.ambient.primitive.magneticFieldT.x,
               oppositeDirect.value.ambient.primitive.magneticFieldT.y,
               oppositeDirect.value.ambient.primitive.magneticFieldT.z} &&
           oppositeState.value.magneticFieldT[2] > 0 &&
           oppositeState.value.velocityGradientPerS->values == oppositeDirect.value.ambient.gradientU &&
           oppositeState.value.magneticGradientTPerM->values == oppositeDirect.value.ambient.gradientB &&
           oppositeState.value.magneticGradientTPerM->branchIdentity !=
               state.value.magneticGradientTPerM->branchIdentity &&
           oppositeState.value.eventIdentity != state.value.eventIdentity,
           "opposite physical owners lost signed fields, derivative branches or event identity");

    auto wrongRegion = query;
    wrongRegion.region = M::Region::CoronalOpen;
    expect(interval.Sample(wrongRegion).status.code == Core::StatusCode::OutOfDomain,
           "a query relabeled Parker exterior as coronal open plasma");
    auto wrongSide = query;
    wrongSide.side = M::Side::UpstreamLimit;
    wrongSide.interfaceIdentity = descriptor.eventIdentity + ":leading-sphere";
    expect(interval.Sample(wrongSide).status.code == Core::StatusCode::UnsupportedCapability,
           "a smooth query impersonated a one-sided shock-limit query");
    auto futureQuery = query;
    futureQuery.epochS = 21;
    expect(!interval.Sample(futureQuery).ok(), "a detached sample exceeded its captured closed interval");
    const auto moved = interval.Move(query.positionM, 0, 20);
    expect(moved.ok() && moved.value.positionM[2] > query.positionM[2], "core material motion was not exported");
    expect(!interval.Move(query.positionM, 0, 21).ok(), "closed interval allowed a later endpoint");
    const auto directMove = event.value->Move({query.positionM[0], query.positionM[1],
                                              query.positionM[2]}, 0, 20, controls);
    expect(directMove.ok(), directMove.status.message.c_str());
    expect(moved.value.positionM == M::Vector{directMove.value.positionM.x,
               directMove.value.positionM.y, directMove.value.positionM.z} &&
           moved.value.velocityMPerS == M::Vector{directMove.value.velocityMPerS.x,
               directMove.value.velocityMPerS.y, directMove.value.velocityMPerS.z} &&
           !moved.value.deformation,
           "detached motion changed the core endpoint/velocity or invented a general deformation export");
    expect(interval.Inverse(query.positionM, 0).status.code == Core::StatusCode::UnsupportedCapability &&
           interval.Contact(query.positionM, 0).status.code == Core::StatusCode::UnsupportedCapability,
           "absent inverse/contact capabilities were silently manufactured");

    const auto front = event.value->Surface({.4, .2, 1}, 0);
    expect(front.ok() && front.value.jump.has_value(), "controlled fixture has no canonical oblique RH limit");
    const auto& x = front.value.positionM;
    const auto surface = interval.Shock({x.x, x.y, x.z}, 0);
    expect(surface.ok(), surface.status.message.c_str());
    expect(surface.value.upstream && surface.value.downstream &&
           !surface.value.downstream->species && !surface.value.downstream->velocityGradientPerS,
           "adapter did not preserve the limited one-sided RH primitive");
    expect(!interval.Shock(query.positionM, 0).ok(), "off-front spherical query accepted");
    query.positionM = {0, 0, 2 * S::Verification::Rs};
    query.region = M::Region::CoronalOpen;
    expect(!interval.Sample(query).ok(), "undisturbed ambient was admitted behind the shock as downstream");

    // The existing manufactured fixture also supplies a stationary front.
    // Its genuine canonical sub-fast classification has no downstream jump.
    // A geometrical sphere must never become an acceleration/source record
    // merely because a caller requests a front at a known radius.
    const auto stationary = S::Event::Create(S::Verification::Configuration(-1));
    expect(stationary.ok(), stationary.status.message.c_str());
    const auto stationaryView = SW3DSphericalAdapter::Capture(stationary.value, controls, 0, 20, 19, 0);
    expect(stationaryView.ok(), stationaryView.status.message.c_str());
    const auto stationaryFront = stationary.value->Surface({.4, .2, 1}, 0);
    expect(stationaryFront.ok() && !stationaryFront.value.jump,
           "stationary control did not retain its genuine no-fast-shock state");
    const auto& stationaryPoint = stationaryFront.value.positionM;
    const auto stationarySurface = stationaryView.value->Shock(
        {stationaryPoint.x, stationaryPoint.y, stationaryPoint.z}, 0);
    expect(stationarySurface.ok() && stationarySurface.value.classification == "no-fast-shock" &&
           stationarySurface.value.upstream && !stationarySurface.value.downstream,
           "sub-fast geometric surface acquired an invented shock/downstream source primitive");

    // Drop the caller's event handle after capture. The detached interval must
    // retain the same immutable physical owner throughout its own lifetime;
    // no configured global event, PIC clock or later capture owns these reads.
    event.value.reset();
    const auto retained = interval.Sample(departureQuery);
    expect(retained.ok() && retained.value.backgroundGeneration == 19 &&
           retained.value.geometryGeneration == 0 &&
           retained.value.magneticFieldT == state.value.magneticFieldT &&
           M::Payload(interval.Identity()) == capturedIdentity,
           "detached owner lifetime or physical identity depended on the caller's event handle");
    result.status = Status::Pass;
    result.message = "Linked spherical core/SEP primitive, motion and RH export agree and keep native admission closed";
  } catch (const std::exception& error) {
    result.status = Status::Fail;
    result.message = error.what();
  }
  result.configuration.push_back("scope=detached-canonical-adapter;native-material-admission=closed");
  result.configuration.push_back("physical-inputs=existing-manufactured-open-field-fixture;new-physics=none");
  result.metrics.emplace_back("adapter_assertions", static_cast<double>(assertions), 0, ">", "count");
  return result;
}

}  // namespace SEP::Testing
