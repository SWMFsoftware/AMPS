#include <cmath>
#include <iostream>
#include <stdexcept>

#include "sep_corona_swcme/coronal_event_ambient.h"
#include "sep_coronal_cme/background_event.h"
#include "sep_coronal_cme/constants.h"

namespace C = SEP::CoronalCME;
namespace A = SEP::CoronaSwcme;
namespace {
constexpr double Rs = C::Constants::kSolarRadiusM, Au = 149597870700.;
void Require(bool value, const std::string& text) {
  if (!value)
    throw std::runtime_error(text);
}
A::AmbientDefinition Definition() {
  A::AmbientDefinition d;
  // A complete declared analytic verification plasma, using the maintained
  // canonical ambient factory WITHOUT inactive sheath/ejecta/source inputs.
  // None of these values is inferred from a missing species/field datum.
  d.composition = {5. / 3., 1.5e6, 1.5e6, 6e6, .04, false};
  d.ambient.sourceSurfaceRadiusM = 2.5 * Rs;
  d.ambient.referenceRadiusM = Au;
  d.ambient.electronDensityAtReferenceM3 = 5e6;
  d.ambient.closedBaseElectronDensityM3 = 1e14;
  d.ambient.rotationRateRadPerS = 2.86533e-6;
  d.ambient.minimumMagneticFieldT = 1e-14;
  d.ambient.traceStepM = 1e6;
  d.ambient.gradientRelativeStep = 1e-5;
  d.ambient.traceMaximumSteps = 200000;
  d.ambient.windTablePoints = 513;
  d.ambient.harmonics = "1,0,0.0001,0";
  d.support = {0, 1200, Rs, Rs, 2 * Au, 1e-6, 2000};
  d.coordinateFrame = "HCI-like-inertial";
  d.physicsFingerprint = "S14BG09-actual-canonical-dipole-ambient-v1";
  return d;
}
C::BackgroundEventDefinition Event(const A::AmbientDefinition& ambient, double latitude,
                                   bool exterior) {
  C::BackgroundEventDefinition d;
  d.front.basis = C::BuildRadialPrincipalBasis(latitude, 0, .2).value;
  const double center = exterior ? 4 * Rs : .97 * Rs,
               axis = exterior ? .25 * Rs : .08 * Rs;
  for (unsigned i = 0; i < 4; ++i) {
    auto& h = d.front.components[i];
    h.referenceValueM = i == 0 ? center : axis;
    h.transitionDurationS = 1200;
  }
  d.support = {0,
               1200,
               Rs,
               1.05 * Rs,
               Rs,
               2 * Au,
               ambient.coordinateFrame,
               "synthetic-epoch-v1",
               exterior ? C::AttachmentPolicy::AllowDetachedFront
                        : C::AttachmentPolicy::RequireSolarIntersection,
               C::PlasmaCoveragePolicy::RequireEveryExposedPatch};
  d.numerics = {16, 1e-8, 8, 16};
  d.gammaAdiabatic = ambient.composition.gammaAdiabatic;
  d.includeElectronMass = ambient.composition.includeElectronMass;
  d.inputIdentity = "S14BG09-canonical-ambient-binding";
  d.ambientIdentity = ambient.physicsFingerprint;
  return d;
}
}  // namespace
int main() {
  try {
    const auto definition = Definition();
    const auto actual = A::AmbientModel::Create(definition);
    Require(actual.ok(), actual.status.message);
    const auto binding = A::BindAmbientEpoch(actual.value, 0, 7, "synthetic-epoch-v1",
                                             definition.physicsFingerprint);
    Require(binding.ok(), binding.status.message);
    const auto epoch = binding.value;
    const C::Vec3 direction{std::sqrt(.5), 0, std::sqrt(.5)};
    const auto north = actual.value->Evaluate(Au * direction, 0),
               south =
                   actual.value->Evaluate({Au * direction.x, 0, -Au * direction.z}, 0);
    Require(north.ok() && south.ok() && north.value.magneticSector == 1 &&
                south.value.magneticSector == -1,
            "actual signed Parker sectors missing");
    const double radialSpeed = C::Dot(north.value.velocityMPerS, direction);
    const double radialFlux =
        north.value.plasma.massDensityKgM3 * radialSpeed * Au * Au;
    Require(std::abs(radialFlux - actual.value->MassFluxPerSteradianKgPerS()) <=
                2e-12 * actual.value->MassFluxPerSteradianKgPerS(),
            "independent canonical radial mass-flux invariant failed");
    const auto open = epoch->At({0, 0, 1.2 * Rs}, C::QuerySide::Ordinary),
               closed = epoch->At({1.2 * Rs, 0, 0}, C::QuerySide::Ordinary);
    Require(open.ok() && closed.ok() &&
                open.value.region == C::AmbientRegion::PfssOpen &&
                closed.value.region == C::AmbientRegion::PfssClosed &&
                closed.value.sector == 0,
            "actual open/closed states or one-sided tensors missing");
    Require(C::Norm(closed.value.primitive.velocityMPerS -
                    C::Vec3{0, definition.ambient.rotationRateRadPerS * 1.2 * Rs, 0}) <
                1e-8,
            "independent closed-corotation state mismatch");
    const auto below = epoch->At({0, 0, definition.ambient.sourceSurfaceRadiusM - 1},
                                 C::QuerySide::Minus);
    const auto above = epoch->At({0, 0, definition.ambient.sourceSurfaceRadiusM + 1},
                                 C::QuerySide::Plus);
    Require(below.ok() && above.ok() &&
                below.value.region == C::AmbientRegion::PfssOpen &&
                above.value.region == C::AmbientRegion::ParkerExterior,
            "source-surface one-sided branch/derivative query failed");
    Require(!epoch->At({0, 0, 1.2 * Rs}, C::QuerySide::Plus).ok(),
            "wrong PFSS radial trace was relabeled as Parker");
    Require(!epoch->At({0, 0, Au}, C::QuerySide::Minus).ok(),
            "wrong Parker radial trace was relabeled as PFSS");
    Require(
        !epoch->At({0, 0, definition.ambient.sourceSurfaceRadiusM}, C::QuerySide::Minus)
             .ok(),
        "unavailable exact opposite source-surface trace was fabricated");
    Require(!A::BindAmbientEpoch(actual.value, -1, 7, "synthetic-epoch-v1",
                                 definition.physicsFingerprint)
                 .ok(),
            "uncaptured history epoch accepted");
    Require(
        !A::BindAmbientEpoch(actual.value, 0, 7, "synthetic-epoch-v1", "wrong-event")
             .ok(),
        "foreign actual ambient identity accepted");
    const auto null = epoch->At({Au, 0, 0}, C::QuerySide::Ordinary);
    Require(!null.ok(), "actual HCS/null query was silently filled");
    for (const auto& target :
         std::vector<std::pair<double, bool>>{{C::Constants::kPi / 2, false},
                                              {0, false},
                                              {-C::Constants::kPi / 2, false},
                                              {C::Constants::kPi / 2, true}}) {
      const auto d = Event(definition, target.first, target.second);
      const auto event = C::BackgroundEvent::Create(d);
      Require(event.ok(), event.status.message);
      const auto prepared = event.value->Prepare(0, epoch, C::QuerySide::Ordinary);
      Require(prepared.ok(), prepared.status.message);
      Require(prepared.value->ambient == epoch &&
                  prepared.value->shocks->backgroundGeneration == 7,
              "actual owning ambient epoch/generation was replaced");
      for (const auto& patch : prepared.value->frontPatches) {
        Require(patch.plasmaSupport == C::PatchSupport::PlasmaAvailable &&
                    patch.ambient.derivativesValid,
                "assembled canonical state has missing coverage/tensors");
        Require(
            patch.parameterCellCoverage ==
                C::ParameterCellCoverage::InsidePlasmaSupport,
            "strict canonical coverage did not certify whole exposed parameter cell");
      }
      Require(prepared.value->shocks->measures.sourceActiveAreaM2 == 0,
              "background epoch created a particle source");
    }
    // Resolve an expanding, detached front in the ACTUAL signed canonical
    // Parker ambient, rather than only classifying static geometry. The
    // declared one-million-metre-per-second semiaxis rates are verification
    // inputs, not inferred acceleration or universal CME defaults.
    auto moving = Event(definition, C::Constants::kPi / 4, true);
    for (unsigned i = 1; i < 4; ++i) {
      moving.front.components[i].initialRateMPerS = 1e6;
      moving.front.components[i].finalRateMPerS = 1e6;
    }
    const auto expanding = C::BackgroundEvent::Create(moving);
    Require(expanding.ok(), expanding.status.message);
    const auto laterBinding = A::BindAmbientEpoch(
        actual.value, 100, 8, "synthetic-epoch-v1", definition.physicsFingerprint);
    Require(laterBinding.ok(), laterBinding.status.message);
    const auto advanced = expanding.value->Prepare(100, laterBinding.value);
    Require(advanced.ok(), advanced.status.message);
    unsigned fast = 0, oblique = 0;
    for (const auto& patch : advanced.value->shocks->patches) {
      if (!patch.fast)
        continue;
      ++fast;
      Require(patch.jump.compressionRatio > 1 && patch.jump.entropyLogIncrement > 0 &&
                  patch.jump.residuals.maximum <= 1e-9,
              "actual canonical ambient jump failed the unchanged RH gate");
      if (patch.jump.upstreamCharacteristics.obliquityRad > .1 &&
          patch.jump.upstreamCharacteristics.obliquityRad < C::Constants::kPi / 2 - .1)
        ++oblique;
    }
    Require(fast > 0 && oblique > 0 &&
                advanced.value->shocks->measures.sourceActiveAreaM2 == 0,
            "expanding canonical ambient front lacks resolved fast oblique/no-source "
            "patches");
    std::cout << "[S14BG09] PASS actual immutable PFSS/open/closed/Parker ambient, "
                 "both signed sectors, source-surface one-sided tensors, EOS and full "
                 "background-only epochs\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "[S14BG09] FAIL " << error.what() << '\n';
    return 1;
  }
}
