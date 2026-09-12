#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

// Verify the full foreign-state PST03 diagnostic while also retaining PST02's
// model-instance precedence.  "Expected" is the receiver's current
// configuration; "supplied" is the snapshot carried by the prepared state.
void expect_cross_configuration_rejection(
    swcme_test::Context& context, const swcme::ModelStatus& status,
    swcme::ConfigurationDigest expected,
    swcme::ConfigurationDigest supplied, const std::string& label) {
  context.expect_true(status.code==swcme::StatusCode::StateModelMismatch,
                      label+" preserves cross-model precedence");
  context.expect_true(status.has_configuration_digests,
                      label+" carries configuration diagnostics");
  context.expect_true(status.expected_configuration_digest==expected,
                      label+" identifies receiver configuration");
  context.expect_true(status.supplied_configuration_digest==supplied,
                      label+" identifies prepared configuration");
  const std::string summary=status.summary();
  context.expect_true(
      summary.find("expected_configuration_digest=0x")!=std::string::npos &&
      summary.find("supplied_configuration_digest=0x")!=std::string::npos,
      label+" summary prints stable hexadecimal digests");
}

}  // namespace

// PST03: configuration ownership covers every configuration family rather
// than special-casing gamma.  Each foreign model below changes exactly one
// field, followed by a multi-field case and an equal default-versus-explicit
// case.  A mutable 1-D model additionally proves that the dedicated same-owner
// status rejects stale state before output modification.
void test_pst03(swcme_test::Context& context) {
  std::cout << "PST03 cross-configuration prepared-state rejection\n";

  swcme3d::Params baseline;
  baseline.shape=swcme3d::ShockShape::Sphere;
  baseline.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  baseline.r0_Rs=40.0;
  baseline.V0_sh_kms=1200.0;
  baseline.V_sw_kms=400.0;
  baseline.cme_dir[0]=0.73;
  baseline.cme_dir[1]=-0.41;
  baseline.cme_dir[2]=0.547;
  baseline.region_mode=swcme::regions::Mode::ShockOnly;
  baseline.shock_acceleration_mode=swcme::acceleration::Mode::Source;

  const swcme::ConfigurationDigest baseline_digest=
      swcme3d::configuration_digest(baseline);
  // This reviewed golden value detects accidental changes to serialization
  // order, convention tags, or floating-point encoding across compiler runs.
  // Intentional configuration-schema changes must update both the schema tag
  // and this value as an explicit validation-contract revision.
  context.expect_true(baseline_digest==13445060466342315912ULL,
                      "baseline configuration digest is reproducible");
  swcme3d::Model producer(baseline);
  const swcme3d::StepState state=producer.prepare_step(0.0);
  context.expect_true(state.configuration_digest==baseline_digest,
                      "prepared 3-D state records resolved configuration");

  // Independent construction and explicit assignment of default-valued
  // fields must not affect the digest.  The state is still rejected because
  // it belongs to another model instance, and equal digests explain why.
  swcme3d::Params explicit_defaults=baseline;
  explicit_defaults.axis_ratio_y=swcme::defaults::ELLIPSOID_AXIS_RATIO_Y;
  explicit_defaults.axis_ratio_z=swcme::defaults::ELLIPSOID_AXIS_RATIO_Z;
  explicit_defaults.T_K=swcme::defaults::T_K;
  explicit_defaults.gamma_ad=swcme::defaults::GAMMA_AD;
  explicit_defaults.data_time_s=std::vector<double>{};
  explicit_defaults.data_radius_Rs=std::vector<double>{};
  context.expect_true(swcme3d::configuration_digest(explicit_defaults)==
                          baseline_digest,
                      "default and explicit-equal configurations hash equally");
  swcme3d::Model equal_receiver(explicit_defaults);
  double x=swcme::constants::AU_M,y=0.0,z=0.0;
  double density=701.0,vx=702.0,vy=703.0,vz=704.0;
  swcme::ModelStatus status=equal_receiver.evaluate_cartesian_fast_checked(
      state,&x,&y,&z,&density,&vx,&vy,&vz,1);
  expect_cross_configuration_rejection(
      context,status,baseline_digest,baseline_digest,
      "explicit-equal foreign model");
  context.expect_true(density==701.0 && vx==702.0 && vy==703.0 && vz==704.0,
                      "equal foreign rejection preserves outputs");

  // Each variant changes one public field representing a validation-plan
  // configuration family. Global Parker polarity remains a compile-time
  // convention; the closure and Parker source radius are now explicit fields
  // and therefore receive their own one-field provenance cases.
  std::vector<std::pair<std::string,swcme3d::Params>> variants;
  swcme3d::Params gamma=baseline; gamma.gamma_ad=1.55;
  variants.emplace_back("gamma",gamma);
  swcme3d::Params geometry=baseline; geometry.axis_ratio_y+=0.2;
  variants.emplace_back("geometry",geometry);
  swcme3d::Params parker_axis=baseline; parker_axis.solar_rotation_axis[1]=0.25;
  variants.emplace_back("Parker axis",parker_axis);
  swcme3d::Params kinematics=baseline;
  kinematics.kinematics_mode=swcme::kinematics::Mode::DBM;
  variants.emplace_back("kinematics",kinematics);
  swcme3d::Params closure=baseline; closure.T_K+=25000.0;
  variants.emplace_back("thermal closure input",closure);
  swcme3d::Params closure_mode=baseline;
  closure_mode.thermodynamic_closure=
      swcme::solarwind::ThermodynamicClosure::MultiSpecies;
  variants.emplace_back("thermodynamic closure mode",closure_mode);
  swcme3d::Params alpha_abundance=baseline;
  alpha_abundance.alpha_to_proton_ratio=0.05;
  variants.emplace_back("alpha abundance",alpha_abundance);
  swcme3d::Params polarity_normalization=baseline;
  polarity_normalization.B1AU_nT+=1.0;
  variants.emplace_back("Parker polarity normalization",polarity_normalization);
  swcme3d::Params parker_source=baseline;
  parker_source.parker_source_radius_Rs=0.75;
  variants.emplace_back("Parker source radius",parker_source);
  swcme3d::Params region=baseline;
  region.region_mode=swcme::regions::Mode::FullICME;
  variants.emplace_back("region mode",region);

  for (const auto& variant : variants) {
    const swcme::ConfigurationDigest receiver_digest=
        swcme3d::configuration_digest(variant.second);
    context.expect_true(receiver_digest!=baseline_digest,
                        variant.first+" changes the digest");
    swcme3d::Model receiver(variant.second);
    density=711.0; vx=712.0; vy=713.0; vz=714.0;
    status=receiver.evaluate_cartesian_fast_checked(
        state,&x,&y,&z,&density,&vx,&vy,&vz,1);
    expect_cross_configuration_rejection(
        context,status,receiver_digest,baseline_digest,variant.first);
    context.expect_true(
        density==711.0 && vx==712.0 && vy==713.0 && vz==714.0,
        variant.first+" rejection occurs before output modification");
  }

  // Multiple simultaneous changes must use the same general comparison and
  // status path; parameter order must not select a special-case branch.
  swcme3d::Params multiple=baseline;
  multiple.gamma_ad=1.5;
  multiple.half_width_rad*=0.8;
  multiple.solar_rotation_rate_rad_s*=0.9;
  multiple.region_mode=swcme::regions::Mode::FullICME;
  multiple.shock_acceleration_mode=
      swcme::acceleration::Mode::ResolvedCompression;
  const swcme::ConfigurationDigest multiple_digest=
      swcme3d::configuration_digest(multiple);
  swcme3d::Model multiple_receiver(multiple);
  density=721.0; vx=722.0; vy=723.0; vz=724.0;
  status=multiple_receiver.evaluate_cartesian_fast_checked(
      state,&x,&y,&z,&density,&vx,&vy,&vz,1);
  expect_cross_configuration_rejection(
      context,status,multiple_digest,baseline_digest,"multi-field mismatch");
  context.expect_true(density==721.0 && vx==722.0 && vy==723.0 && vz==724.0,
                      "multi-field rejection preserves outputs");

  // PST01 now freezes the public model configuration after preparation, so no
  // supported API can create a same-owner stale state.  Corrupt only a copied
  // digest tag to keep the defensive STATE_CONFIGURATION_MISMATCH path under
  // direct test without weakening the public immutability contract.
  swcme1d::Params one_params;
  one_params.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  one_params.r0_Rs=40.0;
  one_params.V0_sh_kms=1200.0;
  one_params.V_sw_kms=400.0;
  one_params.region_mode=swcme::regions::Mode::ShockOnly;
  one_params.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme1d::Model one_model(one_params);
  swcme1d::StepState one_state=one_model.prepare_step(0.0);
  const swcme::ConfigurationDigest prepared_one_digest=
      one_state.configuration_digest;
  const swcme::ConfigurationDigest current_one_digest=
      swcme1d::configuration_digest(one_model.GetParams());
  one_state.configuration_digest^=1ULL;
  double radius=swcme::constants::AU_M;
  double one_density=731.0,one_velocity=732.0;
  status=one_model.evaluate_radii_fast_checked(
      one_state,&radius,&one_density,&one_velocity,1);
  context.expect_true(
      status.code==swcme::StatusCode::StateConfigurationMismatch,
      "same model rejects a mismatched state configuration tag");
  context.expect_true(status.has_configuration_digests &&
                          status.expected_configuration_digest==current_one_digest &&
                          status.supplied_configuration_digest==
                              (prepared_one_digest^1ULL),
                      "same-model mismatch carries current and prepared digests");
  context.expect_true(one_density==731.0 && one_velocity==732.0,
                      "same-model stale-state rejection preserves outputs");
}
