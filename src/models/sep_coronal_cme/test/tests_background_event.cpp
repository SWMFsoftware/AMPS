#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>

#include "sep_coronal_cme/background_event.h"
#include "sep_coronal_cme/clipped_ellipsoid_quadrature.h"
#include "sep_coronal_cme/constants.h"
#include "test_framework.h"

namespace SCCMTest {
namespace {
using namespace SEP::CoronalCME;
constexpr double Rs = 6.957e8;
bool Near(double a, double b, double relative = 1e-10) {
  return std::abs(a - b) <= relative * std::max({1., std::abs(a), std::abs(b)});
}
ScalarGeometryHistory Law(double value, double first, double last,
                          double duration = 1200) {
  ScalarGeometryHistory h;
  h.referenceValueM = value;
  h.initialRateMPerS = first;
  h.finalRateMPerS = last;
  h.transitionDurationS = duration;
  return h;
}
BackgroundEventDefinition Attached() {
  BackgroundEventDefinition d;
  d.front.basis =
      BuildRadialPrincipalBasis(.7853981633974483, 0, .2617993877991494).value;
  d.front.components = {{Law(.97 * Rs, 0, 0), Law(.08 * Rs, 25000, 600000),
                         Law(.12 * Rs, 5000, 200000), Law(.10 * Rs, 4000, 160000)}};
  d.support = {0,
               1200,
               Rs,
               1.05 * Rs,
               Rs,
               3e10,
               "fixture",
               "fixture",
               AttachmentPolicy::RequireSolarIntersection,
               PlasmaCoveragePolicy::RequireEveryExposedPatch};
  d.numerics = {16, 1e-8, 16, 32};
  d.gammaAdiabatic = 5. / 3.;
  d.inputIdentity = "manufactured-attached-event-v1";
  d.ambientIdentity = "explicit-uniform-verification-v1";
  return d;
}
class UniformEpoch final : public AmbientEpoch {
 public:
  mutable AmbientEpochIdentity id;
  bool badPressure = false, oneSided = false, throwLate = false, mutateIdentity = false;
  double magneticFieldT = 1e-6;
  mutable unsigned calls = 0;
  explicit UniformEpoch(const BackgroundEventDefinition& d, double time,
                        std::uint64_t generation = 1) {
    id = {d.ambientIdentity, d.support.frame, d.support.epochUtc,   time, time, time,
          d.gammaAdiabatic,  generation,      d.includeElectronMass};
  }
  const AmbientEpochIdentity& Identity() const override {
    return id;
  }
  SEP::Core::Result<EventAmbientSample> At(Vec3, QuerySide side) const override {
    if (throwLate && ++calls == 5)
      throw std::runtime_error("injected late ambient acquisition failure");
    if (mutateIdentity)
      ++id.generation;
    // explicit adversarial owner, never production behavior
    EventAmbientSample s;
    s.electronDensityM3 = 1e10;
    s.electronTemperatureK = 1e6;
    s.ions = {{"proton", 1, 1, Constants::kProtonMassKg, 1e6}};
    // Explicit uniform, static verification plasma; this fixture makes no
    // claim to be the maintained PFSS/wind provider tested separately S14BG09.
    s.primitive.magneticFieldT = {magneticFieldT, 0, 0};
    const auto eos = EvaluatePlasmaFromElectronDensity(
        s.electronDensityM3, s.electronTemperatureK, s.ions, id.includeElectronMass,
        id.gammaAdiabatic, magneticFieldT);
    s.primitive.massDensityKgM3 = eos.value.massDensityKgM3;
    s.primitive.pressurePa = eos.value.pressurePa * (badPressure ? 2 : 1);
    s.region = AmbientRegion::PfssOpen;
    s.topology = FieldLineTopology::OpenToOuterBoundary;
    s.sector = 1;
    s.derivativesValid = true;
    s.oneSided = oneSided;
    s.side = side;
    return SEP::Core::Result<EventAmbientSample>::Success(s);
  }
};
std::string CompleteGeometry() {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "  #begin shock_geometry\nmodel = triaxial-ellipsoid\naxis_convention = "
         "radial-principal-axis\n"
      << "initial_parameterization = center-and-principal-axes\nreference_time_s = "
         "0\ndirection_frame = fixture\n"
      << "direction_longitude_rad = 0\ndirection_latitude_rad = "
         ".7853981633974483\nlateral_axis_tilt_rad = .2617993877991494\n"
      << "orientation_evolution = fixed\nevolution = "
         "independent-component-laws\ninitial_apex_radius_m = 0\n"
      << "snapshots_file = none\nsnapshot_interpolation = none\nsnapshot_extrapolation "
         "= reject\n";
  const auto d = Attached();
  const char* names[] = {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                         "lateral_semiaxis_2"};
  for (std::size_t i = 0; i < 4; ++i) {
    const auto& h = d.front.components[i];
    out << "#begin " << names[i] << "\ninitial_value_m = " << h.referenceValueM
        << "\nkinematics = smooth-rate-transition\ninitial_rate_m_per_s = "
        << h.initialRateMPerS << "\nfinal_rate_m_per_s = " << h.finalRateMPerS
        << "\ntransition_start_time_s = 0\ntransition_duration_s = "
           "1200\nafter_transition = constant-final-rate\n#end "
        << names[i] << '\n';
  }
  out << "#end shock_geometry\n#begin cme_piston\nenabled = false\ngeometry_source = "
         "none\ngeometry_record = none\n#end cme_piston\n"
      << "#begin piston_geometry\nmodel = none\nrecord_file = "
         "none\nminimum_front_separation_m = 0\n#end piston_geometry\n";
  return out.str();
}
ModelConfiguration Configuration() {
  auto input = ReplaceOnce(Fixture(), "end_time_s = 86400", "end_time_s = 1200");
  input = ReplaceOnce(input, "qualified_source_inner_radius_m = 7.02657e8",
                      "qualified_source_inner_radius_m = 730485000");
  const auto at = input.rfind("#end physics");
  Require(at != std::string::npos, "fixture wrapper missing");
  input.insert(at, CompleteGeometry());
  return ParseGood(input).schema5;
}
GeometryAsset Asset(const GeometryHistory& history, bool table,
                    const std::string& frame = "fixture") {
  std::ostringstream text;
  text << std::setprecision(17) << "SCCM_FIXED_ELLIPSOID 1\nframe "
       << std::quoted(frame)
       << "\nreference_time_s 0\ndirection_rad 0 .7853981633974483 "
          ".2617993877991494\nevolution "
       << std::quoted(table ? "tabulated-snapshots" : "independent-component-laws")
       << '\n';
  if (table) {
    text << "samples 2\n";
    // A constant, explicitly supplied verification table provides an exact
    // interpolation oracle; this does not infer missing component slopes.
    for (double time : {0., 1200.}) {
      text << time;
      for (const auto& component : history.components)
        text << ' ' << component.referenceValueM << " 0";
      text << '\n';
    }
  } else {
    const char* names[] = {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                           "lateral_semiaxis_2"};
    for (unsigned i = 0; i < 4; ++i) {
      const auto& component = history.components[i];
      text << "component " << std::quoted(names[i]) << ' ' << component.referenceValueM
           << ' ' << component.initialRateMPerS << ' ' << component.finalRateMPerS
           << ' ' << component.transitionStartS << ' ' << component.transitionDurationS
           << '\n';
    }
  }
  GeometryAsset asset;
  asset.bytes = text.str();
  asset.expectedSha256 = ComputeContentChecksum(asset.bytes);
  return asset;
}
void S14BG01() {
  const auto config = Configuration();
  const auto d = Attached();
  const auto materialized =
      MaterializeBackgroundEvent(config, d.support, d.numerics, d.ambientIdentity);
  Require(materialized.ok(), materialized.status.message);
  auto byApex = config;
  byApex.assignments["shock_geometry.initial_parameterization"] =
      "apex-and-principal-axes";
  byApex.assignments["shock_geometry.initial_apex_radius_m"] = "730485000";
  byApex.assignments["shock_geometry.center_distance.initial_value_m"] = "0";
  byApex.physicsFingerprint = ComputePhysicsFingerprint(byApex.assignments);
  const auto apex =
      MaterializeBackgroundEvent(byApex, d.support, d.numerics, d.ambientIdentity);
  Require(apex.ok(), apex.status.message);
  for (double t : {0., 600., 1200.}) {
    const auto a = materialized.value.front.At(t), b = apex.value.front.At(t);
    Require(a.ok() && b.ok(), "complete histories unavailable");
    Require(a.value.centerDistanceM.value == b.value.centerDistanceM.value &&
                a.value.centerDistanceM.firstDerivative ==
                    b.value.centerDistanceM.firstDerivative,
            "apex parameterization changed center history/rates");
  }
  byApex.assignments["shock_geometry.center_distance.initial_value_m"] = "1";
  byApex.physicsFingerprint = ComputePhysicsFingerprint(byApex.assignments);
  Require(!MaterializeBackgroundEvent(byApex, d.support, d.numerics, d.ambientIdentity)
               .ok(),
          "conflicting radial authorities accepted");
  auto missing = config;
  missing.assignments.erase("shock_geometry.lateral_semiaxis_2.final_rate_m_per_s");
  missing.physicsFingerprint = ComputePhysicsFingerprint(missing.assignments);
  Require(!MaterializeBackgroundEvent(missing, d.support, d.numerics, d.ambientIdentity)
               .ok(),
          "missing lateral history was inferred");
  auto dynamic = config;
  dynamic.assignments["shock_geometry.orientation_evolution"] = "dynamic-attitude";
  dynamic.physicsFingerprint = ComputePhysicsFingerprint(dynamic.assignments);
  Require(!MaterializeBackgroundEvent(dynamic, d.support, d.numerics, d.ambientIdentity)
               .ok(),
          "research attitude admitted");

  auto tabulated = config;
  tabulated.assignments["shock_geometry.evolution"] = "tabulated-snapshots";
  tabulated.assignments["shock_geometry.snapshots_file"] = "owned-front-v1";
  tabulated.assignments["shock_geometry.snapshot_interpolation"] = "monotone-cubic-c1";
  for (const auto* name : {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                           "lateral_semiaxis_2"}) {
    const auto prefix = "shock_geometry." + std::string(name) + ".";
    tabulated.assignments[prefix + "kinematics"] = "inactive";
    tabulated.assignments[prefix + "after_transition"] = "none";
    for (const auto* field :
         {"initial_value_m", "initial_rate_m_per_s", "final_rate_m_per_s",
          "transition_start_time_s", "transition_duration_s"})
      tabulated.assignments[prefix + field] = "0";
  }
  tabulated.physicsFingerprint = ComputePhysicsFingerprint(tabulated.assignments);
  auto asset = Asset(d.front, true);
  const GeometryAssetReader read = [&](const std::string& name) {
    Require(name == "owned-front-v1", "materializer requested an undeclared asset");
    return SEP::Core::Result<GeometryAsset>::Success(asset);
  };
  const auto table = MaterializeBackgroundEvent(tabulated, d.support, d.numerics,
                                                d.ambientIdentity, read);
  Require(table.ok(), table.status.message);
  const auto middle = table.value.front.At(600);
  Require(
      middle.ok() &&
          middle.value.radialSemiAxisM.value == d.front.components[1].referenceValueM &&
          middle.value.radialSemiAxisM.firstDerivative == 0 &&
          BackgroundEvent::Create(table.value).ok(),
      "verified complete snapshot history changed its explicit constant data");
  asset.expectedSha256 = std::string(64, '0');
  Require(!MaterializeBackgroundEvent(tabulated, d.support, d.numerics,
                                      d.ambientIdentity, read)
               .ok(),
          "wrong expected asset checksum was accepted");
  asset = Asset(d.front, true, "foreign-frame");
  Require(!MaterializeBackgroundEvent(tabulated, d.support, d.numerics,
                                      d.ambientIdentity, read)
               .ok(),
          "foreign snapshot frame was relabeled");
  asset = Asset(d.front, true);
  asset.bytes += "orientation_changed true\n";
  asset.expectedSha256 = ComputeContentChecksum(asset.bytes);
  Require(!MaterializeBackgroundEvent(tabulated, d.support, d.numerics,
                                      d.ambientIdentity, read)
               .ok(),
          "trailing dynamic-attitude asset fields were ignored");
  asset = Asset(d.front, true);
  auto conflicting = tabulated;
  conflicting.assignments["shock_geometry.radial_semiaxis.initial_value_m"] = "1";
  conflicting.physicsFingerprint = ComputePhysicsFingerprint(conflicting.assignments);
  Require(!MaterializeBackgroundEvent(conflicting, d.support, d.numerics,
                                      d.ambientIdentity, read)
               .ok(),
          "table silently overrode a second component-law authority");
}
void S14BG02() {
  const auto d = Attached();
  const auto k = d.front.At(600);
  Require(k.ok(), k.status.message);
  const auto shape = FixedOrientationEllipsoid::FromCenter(d.front.basis, k.value, Rs);
  Require(shape.ok(), shape.status.message);
  const double theta = 1.1, phi = .7;
  const double u[] = {std::cos(theta), std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi)};
  const Vec3 axes[] = {d.front.basis.radial, d.front.basis.firstLateral,
                       d.front.basis.secondLateral};
  const KinematicValue a[] = {k.value.radialSemiAxisM, k.value.firstLateralSemiAxisM,
                              k.value.secondLateralSemiAxisM};
  Vec3 normal{},
      velocity = k.value.centerDistanceM.firstDerivative * d.front.basis.radial;
  for (unsigned i = 0; i < 3; ++i) {
    normal = normal + (u[i] / a[i].value) * axes[i];
    velocity = velocity + (a[i].firstDerivative * u[i]) * axes[i];
  }
  normal = Unit(normal);
  const auto evaluated = shape.value.Evaluate(shape.value.Point(theta, phi));
  Require(evaluated.ok(), evaluated.status.message);
  Require(Norm(evaluated.value.outwardNormal - normal) < 1e-13 &&
              Near(evaluated.value.normalSpeedMPerS, Dot(velocity, normal), 1e-13),
          "independent triaxial gradient/implicit speed mismatch");
  // Independent solid-angle quadrature uses cos(theta), rather than the
  // production midpoint theta Jacobian. No production area routine supplies
  // the reference; all three unequal semiaxes contribute.
  const auto uncut =
      FixedOrientationEllipsoid::FromCenter(d.front.basis, k.value, 0).value;
  double oracle = 0;
  const int n = 256, m = 512;
  for (int i = 0; i < n; ++i) {
    const double mu = -1 + (i + .5) * 2 / n;
    for (int j = 0; j < m; ++j) {
      const double ph = (j + .5) * 2 * Constants::kPi / m;
      oracle += a[0].value * a[1].value * a[2].value *
                std::sqrt(mu * mu / (a[0].value * a[0].value) +
                          (1 - mu * mu) * std::pow(std::cos(ph) / a[1].value, 2) +
                          (1 - mu * mu) * std::pow(std::sin(ph) / a[2].value, 2)) *
                4 * Constants::kPi / (n * m);
    }
  }
  const auto patches = uncut.Tessellate(128, 256);
  Require(patches.ok(), patches.status.message);
  double area = 0;
  for (const auto& p : patches.value)
    area += p.areaM2;
  Require(std::abs(area - oracle) < 2e-4 * oracle,
          "triaxial area differs from independent solid-angle reference");

  // These closed-form sphere caps audit the actual cut-cell rule, including
  // a cap so small that every original parameter-cell midpoint is buried.
  // Promote the exact input doubles before evaluating the independent
  // intersection-circle formula; do not use a production clipping result.
  for (bool tiny : {false, true}) {
    const double radius = .08 * Rs;
    const double center = tiny ? Rs - radius + 1e-6 * Rs : .97 * Rs;
    const EllipsoidKinematics sphereState{
        {center, 0, 0}, {radius, 0, 0}, {radius, 0, 0}, {radius, 0, 0}};
    const auto sphere =
        FixedOrientationEllipsoid::FromCenter(d.front.basis, sphereState, Rs);
    Require(sphere.ok(), sphere.status.message);
    const long double c = center, aSphere = radius, solar = Rs;
    const long double cut =
        ((solar - c) * (solar + c) - aSphere * aSphere) / (2 * c * aSphere);
    const long double exact = 2 * std::acos(-1.L) * aSphere * aSphere * (1 - cut);
    for (int polarCells : {2, 8, 16}) {
      const int azimuthCells = 2 * polarCells + 1;
      const auto clipped =
          TessellateSolarExposed(sphere.value, Rs, polarCells, azimuthCells);
      Require(clipped.ok(), clipped.status.message);
      long double exposedArea = 0;
      std::set<std::uint64_t> ids;
      for (const auto& patch : clipped.value) {
        exposedArea += patch.areaM2;
        Require(patch.areaM2 > 0 && Norm(patch.centerM) >= Rs && patch.physicalId > 0 &&
                    patch.physicalId <=
                        static_cast<std::uint64_t>(polarCells * azimuthCells) &&
                    ids.insert(patch.physicalId).second,
                "cut cell lost its exposed representative/positive area/original ID");
      }
      Require(std::abs(exposedArea - exact) <= 1e-9L * exact,
              "cut-cell area failed the frozen independent sphere-cap tolerance");
      if (tiny) {
        const auto midpoint = sphere.value.Tessellate(polarCells, azimuthCells);
        Require(
            (!midpoint.ok() || midpoint.value.empty()) && !clipped.value.empty(),
            "tiny-cap fixture no longer exercises midpoint-only solar clipping loss");
      }
    }
  }
  const EllipsoidKinematics islandState{
      {.5 * Rs, 0, 0}, {.1 * Rs, 0, 0}, {.87 * Rs, 0, 0}, {.2 * Rs, 0, 0}};
  const auto island =
      FixedOrientationEllipsoid::FromCenter(d.front.basis, islandState, Rs);
  Require(island.ok(), island.status.message);
  const auto midpoint = island.value.Tessellate(8, 6);
  const auto exposed = TessellateSolarExposed(island.value, Rs, 8, 6);
  Require(exposed.ok(), exposed.status.message);
  Require((!midpoint.ok() || midpoint.value.empty()) && !exposed.value.empty(),
          "triaxial exposed islands between all legacy midpoints were dropped");
}
void S14BG03() {
  const auto inside = EvaluateSmoothKinematics(2, 2, 10, 0, 4, 0, 2);
  const auto after = EvaluateSmoothKinematics(6, 6, 10, 0, 4, 0, 2);
  const auto past = EvaluateSmoothKinematics(0, 6, 10, 0, 4, 0, 2);
  Require(inside.ok() && inside.value.value == 10 && after.ok() &&
              after.value.value == 10 && past.ok() && Near(past.value.value, 2),
          "reference value is not honored inside/after a smooth rate transition");
  auto h = Law(1, -2, 2, 4);
  const auto range = h.Range(0, 4);
  Require(range.ok() && Near(range.value.first, -1.75) && Near(range.value.second, 1),
          "smooth between-endpoint minimum was missed");
  h.kind = HistoryKind::CubicHermiteC1;
  h.knots = {{0, 10, -80}, {1, 10, 80}};
  const auto table = h.Range(0, 1);
  Require(table.ok() && Near(table.value.first, -10),
          "Hermite interior negative extremum was missed");
  auto d = Attached();
  d.front.components[2] = h;
  d.support.endS = 1;
  Require(!BackgroundEvent::Create(d).ok(),
          "positive endpoint axes hid an interior invalid axis");
  Require(!h.Range(-1, 1).ok(), "history silently extrapolated");

  // The public asset limit is 100,000 knots. Range must validate the table
  // once and evaluate each segment locally, rather than rescan all knots for
  // every extremum. This full-size exact linear history also audits endpoints.
  h.knots.clear();
  h.knots.reserve(100000);
  for (unsigned i = 0; i < 100000; ++i)
    h.knots.push_back({static_cast<double>(i), 10. + i, 1.});
  const auto large = h.Range(0, 99999);
  Require(large.ok() && large.value.first == 10 && large.value.second == 100009,
          "maximum-size Hermite history lost its exact linear range");
}
void S14BG04() {
  auto d = Attached();
  d.support.plasmaMinimumRadiusM = 1.02 * Rs;
  d.support.plasmaCoverage = PlasmaCoveragePolicy::GeometryOnlyOutsideCoverage;
  const auto event = BackgroundEvent::Create(d);
  Require(event.ok(), event.status.message);
  const auto epoch = event.value->Prepare(0, std::make_shared<UniformEpoch>(d, 0));
  Require(epoch.ok(), epoch.status.message);
  unsigned below = 0, available = 0, qualified = 0;
  for (const auto& p : epoch.value->frontPatches) {
    Require(Norm(p.geometry.centerM) >= Rs, "patch below solar absorption sphere");
    below += p.plasmaSupport == PatchSupport::BelowPlasmaCoverage;
    available += p.plasmaSupport == PatchSupport::PlasmaAvailable;
    qualified += p.sourceRadialSupport == SourceRadialSupport::QualifiedRadius;
  }
  Require(below && available && !qualified,
          "solar/plasma/source masks were conflated at apex=R_in");
  Require(epoch.value->frontRadialExtent.intersectsSolarSurface,
          "topological solar intersection lost");
  auto d2 = d;
  d2.support.plasmaCoverage = PlasmaCoveragePolicy::RequireEveryExposedPatch;
  const auto strict = BackgroundEvent::Create(d2);
  Require(strict.ok(), strict.status.message);
  Require(!strict.value->Prepare(0, std::make_shared<UniformEpoch>(d2, 0)).ok(),
          "missing required guard plasma was hidden");
  auto cap = Attached();
  cap.support.plasmaMaximumRadiusM = 1.05 * Rs - 1;
  const auto capEvent = BackgroundEvent::Create(cap);
  Require(capEvent.ok(), capEvent.status.message);
  Require(!capEvent.value->Prepare(0, std::make_shared<UniformEpoch>(cap, 0)).ok(),
          "an unsupported apex between quadrature centers escaped continuous radial "
          "coverage admission");
  // Sphere-clipping reference is exact; production patch areas converge as
  // midpoint surface quadrature. The independently derived one-cell boundary
  // bound is numerical, not a magnetic/source floor or adjusted tolerance.
  const double center = .97 * Rs, a = .08 * Rs,
               cut = (Rs * Rs - center * center - a * a) / (2 * center * a),
               oracle = 2 * Constants::kPi * a * a * (1 - cut);
  EllipsoidKinematics k{{center, 0, 0}, {a, 0, 0}, {a, 0, 0}, {a, 0, 0}};
  const auto sphere = FixedOrientationEllipsoid::FromCenter(d.front.basis, k, Rs).value;
  for (int n : {32, 64, 128}) {
    const auto p = sphere.Tessellate(n, 2 * n);
    Require(p.ok(), p.status.message);
    double sum = 0;
    for (const auto& q : p.value)
      sum += q.areaM2;
    Require(std::abs(sum - oracle) <= 4 * Constants::kPi * Constants::kPi * a * a / n,
            "solar clipped area exceeds independent midpoint-cell error bound");
  }
}
void S14BG05() {
  const auto d = Attached();
  const auto event = BackgroundEvent::Create(d);
  Require(event.ok(), event.status.message);
  const auto first = event.value->Prepare(0, std::make_shared<UniformEpoch>(d, 0));
  Require(first.ok(), first.status.message);
  const auto oldShock = first.value->shocks;
  auto changed = d;
  changed.front.components[1].finalRateMPerS += 1;
  const auto different = BackgroundEvent::Create(changed);
  Require(different.ok(), different.status.message);
  Require(different.value->Definition().physicsIdentity !=
              event.value->Definition().physicsIdentity,
          "changed physical geometry reused a caller's stale event identity");
  auto changingOwner = std::make_shared<UniformEpoch>(d, 600, 2);
  changingOwner->mutateIdentity = true;
  Require(!event.value->Prepare(600, changingOwner).ok() &&
              event.value->Current() == first.value,
          "ambient identity changed during preparation and published a mixed epoch");
  auto bad = std::make_shared<UniformEpoch>(d, 600, 2);
  bad->throwLate = true;
  Require(!event.value->Prepare(600, bad).ok() &&
              event.value->Current() == first.value &&
              event.value->Current()->shocks == oldShock,
          "failed epoch partially published geometry/background/shock");
  const auto next =
      event.value->Prepare(600, std::make_shared<UniformEpoch>(d, 600, 2));
  Require(next.ok(), next.status.message);
  Require(next.value->generation == 2 && first.value->generation == 1 &&
              oldShock->generation == 1 && first.value->timeS == 0,
          "owning complete epoch lifetime or generation changed");
  Require(
      !event.value->Prepare(1201, std::make_shared<UniformEpoch>(d, 1201, 3)).ok() &&
          event.value->Current() == next.value,
      "out-of-history epoch changed prior publication");
}
void S14BG06() {
  auto d = Attached();
  const auto event = BackgroundEvent::Create(d);
  Require(event.ok(), event.status.message);
  const auto early = event.value->Prepare(0, std::make_shared<UniformEpoch>(d, 0));
  Require(early.ok(), early.status.message);
  Require(early.value->shocks->measures.fastAreaM2 == 0 &&
              early.value->shocks->measures.sourceActiveAreaM2 == 0,
          "sub-fast geometry manufactured an initial shock/source");
  const auto later =
      event.value->Prepare(1200, std::make_shared<UniformEpoch>(d, 1200, 2));
  Require(later.ok(), later.status.message);
  Require(later.value->shocks->measures.fastAreaM2 > 0 &&
              later.value->shocks->measures.sourceActiveAreaM2 == 0,
          "delayed physical shock classification or background-only source suppression "
          "failed");
  for (const auto& patch : later.value->shocks->patches) {
    Require(!patch.sourceActive && !patch.sourceEligibleBeforeClearance &&
                patch.physicalNumberRatePerS == 0,
            "background assembly required/created a particle source");
    if (patch.fast)
      Require(patch.jump.compressionRatio > 1 && patch.jump.entropyLogIncrement > 0 &&
                  patch.jump.residuals.maximum <= 1e-9,
              "assembled RH jump violates unchanged canonical gate");
    else
      Require(patch.jump.compressionRatio == 0,
              "sub-fast patch has invented compression");
  }
  // Preserve the original very-high-beta Cartesian counterexample as an
  // explicit fail-closed test. Its tiny tangential magnetic stress cannot be
  // represented to the unchanged relative RH gate after large Cartesian
  // normal velocities are reconstructed. No field floor or relaxed flux
  // normalization is introduced to admit that unsupported numerical state.
  auto unresolved = std::make_shared<UniformEpoch>(d, 1200, 3);
  unresolved->magneticFieldT = 1e-12;
  const auto rejected = event.value->Prepare(1200, unresolved);
  Require(!rejected.ok() && event.value->Current() == later.value,
          "unresolved high-beta RH state was admitted or changed the prior epoch");
  // An exactly axis-aligned parallel weak limit has independently
  // representable tangential invariants. Its hydrodynamic compression is an
  // independent long-double oracle, while the complete assembled surface
  // above retains nondegenerate oblique MHD fields and the original1e-9 gate.
  const double gamma = 5. / 3.;
  MhdPrimitiveState parallel{
      1,
      1 / gamma,
      {0, 0, 0},
      {.5 * std::sqrt(Constants::kVacuumPermeabilityHPerM), 0, 0}};
  for (double excess : {1e-2, 1e-3, 1e-4}) {
    const auto solved = SolveObliqueFastShock(parallel, {1, 0, 0}, 1 + excess, gamma);
    Require(solved.ok(), solved.status.message);
    const long double mach = 1 + static_cast<long double>(excess), g = 5.L / 3.L;
    const double oracle =
        static_cast<double>((g + 1) * mach * mach / ((g - 1) * mach * mach + 2));
    Require(std::abs(solved.value.compressionRatio - oracle) < 2e-9 &&
                solved.value.residuals.maximum <= 1e-9,
            "resolved axis-aligned weak-shock limit violates independent "
            "compression/RH oracle");
  }
}
void S14BG07() {
  const auto d = Attached();
  const auto event = BackgroundEvent::Create(d);
  Require(event.ok(), event.status.message);
  auto ambient = std::make_shared<UniformEpoch>(d, 0);
  ambient->badPressure = true;
  Require(!event.value->Prepare(0, ambient).ok() && !event.value->Current(),
          "inconsistent composition/pressure accepted");
  ambient = std::make_shared<UniformEpoch>(d, 0);
  ambient->id.frame = "foreign-frame";
  Require(!event.value->Prepare(0, ambient).ok(), "foreign frame accepted");
  ambient = std::make_shared<UniformEpoch>(d, 0);
  ambient->id.gammaAdiabatic = 1.4;
  Require(!event.value->Prepare(0, ambient).ok(),
          "different background and shock EOS admitted");
  ambient = std::make_shared<UniformEpoch>(d, 0);
  ambient->oneSided = true;
  Require(!event.value->Prepare(0, ambient).ok(),
          "one-sided interface state used without an explicit side");
  const auto side = event.value->Prepare(0, ambient, QuerySide::Plus);
  Require(side.ok(), side.status.message);
  for (const auto& p : side.value->frontPatches)
    Require(p.ambient.oneSided && p.ambient.side == QuerySide::Plus,
            "one-sided state/derivative identity lost");
}
void S14BG08() {
  auto d = Attached();
  d.pistonEnabled = true;
  d.piston = d.front;
  // A genuinely independent, slowly growing piston, not a hidden standoff
  // subtraction or an inferred contact/material-sheath closure.
  d.piston.components = {{Law(.97 * Rs, 0, 0), Law(.04 * Rs, 1000, 10000),
                          Law(.04 * Rs, 1000, 10000), Law(.04 * Rs, 1000, 10000)}};
  d.minimumPistonSeparationM = .005 * Rs;
  const auto event = BackgroundEvent::Create(d);
  Require(event.ok(), event.status.message);
  const auto epoch = event.value->Prepare(600, std::make_shared<UniformEpoch>(d, 600));
  Require(epoch.ok(), epoch.status.message);
  Require(epoch.value->pistonEnabled && !epoch.value->pistonPatches.empty(),
          "independent exposed piston not assembled");
  for (const auto& p : epoch.value->pistonPatches) {
    const auto inside = epoch.value->front.Evaluate(p.centerM);
    Require(inside.ok() && inside.value.implicitValue < 0,
            "independent piston escaped the front");
  }
  auto bad = d;
  bad.piston.components[1].referenceValueM = 2 * Rs;
  Require(!BackgroundEvent::Create(bad).ok(),
          "invalid independent piston history accepted");

  // Exercise the public configuration/asset boundary, in addition to the
  // typed-history certificate above. The piston carries its own four laws
  // and does not acquire a center, axes or rates from the front.
  auto config = Configuration();
  config.assignments["cme_piston.enabled"] = "true";
  config.assignments["cme_piston.geometry_source"] = "independent-record";
  config.assignments["cme_piston.geometry_record"] = "piston_geometry";
  config.assignments["piston_geometry.model"] = "triaxial-ellipsoid-versioned-file";
  config.assignments["piston_geometry.record_file"] = "owned-piston-v1";
  config.assignments["piston_geometry.minimum_front_separation_m"] = "3478500";
  config.physicsFingerprint = ComputePhysicsFingerprint(config.assignments);
  auto asset = Asset(d.piston, false);
  const GeometryAssetReader read = [&](const std::string& name) {
    Require(name == "owned-piston-v1", "unexpected piston asset request");
    return SEP::Core::Result<GeometryAsset>::Success(asset);
  };
  const auto loaded = MaterializeBackgroundEvent(config, d.support, d.numerics,
                                                 d.ambientIdentity, read);
  Require(loaded.ok(), loaded.status.message);
  Require(loaded.value.pistonEnabled &&
              loaded.value.piston.components[1].finalRateMPerS == 10000 &&
              loaded.value.front.components[1].finalRateMPerS == 600000 &&
              BackgroundEvent::Create(loaded.value).ok(),
          "independent piston asset was replaced with an inferred front offset");
  auto unbound = config;
  unbound.assignments["cme_piston.geometry_record"] = "none";
  unbound.physicsFingerprint = ComputePhysicsFingerprint(unbound.assignments);
  Require(!MaterializeBackgroundEvent(unbound, d.support, d.numerics,
                                     d.ambientIdentity, read).ok(),
          "enabled piston accepted an explicitly unbound geometry_record=none");
  asset.bytes += "unrecognized-contact-state 1\n";
  asset.expectedSha256 = ComputeContentChecksum(asset.bytes);
  Require(!MaterializeBackgroundEvent(config, d.support, d.numerics, d.ambientIdentity,
                                      read)
               .ok(),
          "unrecognized piston/contact record was silently accepted");
}
void S14BG10() {
  auto d = Attached();
  Require(BackgroundEvent::Create(d).ok(), "documented1200s attached fixture rejected");
  d.front.components[0] = Law(.97 * Rs, 20000, 600000, 120);
  d.front.components[1] = Law(.08 * Rs, 5000, 100000, 120);
  Require(!BackgroundEvent::Create(d).ok(),
          "detachment stress history admitted under attached policy");
  d.support.attachment = AttachmentPolicy::AllowDetachedFront;
  const auto model = BackgroundEvent::Create(d);
  Require(model.ok(), model.status.message);
  const auto before = d.front.At(211.253), after = d.front.At(211.255);
  Require(before.ok() && after.ok(), "stress history unavailable");
  const auto a = EvaluateEllipsoidRadialExtent(before.value, Rs),
             b = EvaluateEllipsoidRadialExtent(after.value, Rs);
  Require(a.ok() && b.ok() && a.value.intersectsSolarSurface &&
              !b.value.intersectsSolarSurface,
          "exact spatial topology missed documented211.254s detachment");
  const double independentTime =
      120 +
      (Rs - (.97 * Rs - .08 * Rs + 120 * ((20000 + 600000) - (5000 + 100000)) / 2)) /
          500000;
  Require(Near(independentTime, 211.254, 1e-12),
          "independent constant-final-rate detachment oracle changed");
  const auto epoch = model.value->Prepare(212, std::make_shared<UniformEpoch>(d, 212));
  Require(epoch.ok() && !epoch.value->frontRadialExtent.intersectsSolarSurface,
          "declared detached-front branch not retained");
}
}  // namespace
// namespace
void RegisterBackgroundEvent(Registry* tests) {
  (*tests)["S14BG01"] = S14BG01;
  (*tests)["S14BG02"] = S14BG02;
  (*tests)["S14BG03"] = S14BG03;
  (*tests)["S14BG04"] = S14BG04;
  (*tests)["S14BG05"] = S14BG05;
  (*tests)["S14BG06"] = S14BG06;
  (*tests)["S14BG07"] = S14BG07;
  (*tests)["S14BG08"] = S14BG08;
  (*tests)["S14BG10"] = S14BG10;
}
}  // namespace SCCMTest
// namespace SCCMTest
