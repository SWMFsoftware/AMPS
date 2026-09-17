#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace {

// Check the complete PST02 diagnostic contract, not merely the enum value.
// Both process-unique identities must survive status propagation so an AMPS
// log can identify the receiving model and the model that prepared the cache.
void expect_mismatch(swcme_test::Context& context,
                     const swcme::ModelStatus& status,
                     swcme::ModelIdentity expected,
                     swcme::ModelIdentity supplied,
                     const std::string& label) {
  context.expect_true(status.code==swcme::StatusCode::StateModelMismatch,
                      label+" reports STATE_MODEL_MISMATCH");
  context.expect_true(status.has_model_identities,
                      label+" carries ownership diagnostics");
  context.expect_true(status.expected_model_identity==expected,
                      label+" identifies receiving model");
  context.expect_true(status.supplied_model_identity==supplied,
                      label+" identifies prepared-state owner");
  const std::string summary=status.summary();
  context.expect_true(summary.find("expected_model_identity=")!=std::string::npos &&
                      summary.find("supplied_model_identity=")!=std::string::npos,
                      label+" summary prints both identities");
}

// Read a small sentinel file exactly.  PST02 uses this to prove that checked
// writers reject ownership before fopen("w") can truncate prior output.
std::string read_file(const char* path) {
  std::ifstream input(path,std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(input),
                     std::istreambuf_iterator<char>());
}

// Legacy geometry APIs predate ModelStatus returns.  They must still reject a
// foreign state explicitly through the shared exception policy instead of
// returning a plausible radius, mesh, or disconnected connectivity record.
template <typename Function>
void expect_mismatch_exception(swcme_test::Context& context,
                               Function&& function,
                               const std::string& label) {
  bool caught=false;
  std::string message;
  try {
    function();
  } catch (const std::runtime_error& error) {
    caught=true;
    message=error.what();
  }
  context.expect_true(caught,label+" throws on foreign state");
  context.expect_true(message.find("STATE_MODEL_MISMATCH")!=std::string::npos,
                      label+" exception preserves mismatch status");
}

}  // namespace

// PST02: a prepared state is owned by one exact model instance.  The test uses
// both identical and deliberately different configurations, exercises direct
// 1-D/3-D checked physics, SEP-facing source/background adapters, legacy
// geometry wrappers, and checked output.  Every output starts with a sentinel
// so rejection ordering is verified rather than inferred from the status.
void test_pst02(swcme_test::Context& context) {
  std::cout << "PST02 cross-model prepared-state rejection\n";

  // Independently constructed equal models must still have different owner
  // identities; equality of physical parameters is not proof of provenance.
  swcme1d::Params p1;
  p1.r0_Rs=20.0;
  p1.V0_sh_kms=1200.0;
  p1.V_sw_kms=400.0;
  p1.region_mode=swcme::regions::Mode::ShockOnly;
  p1.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme1d::Model one_a(p1),one_b(p1);
  const swcme1d::StepState one_step=one_a.prepare_step(0.0);
  context.expect_true(one_a.model_identity()!=one_b.model_identity(),
                      "equal 1-D model instances have distinct identities");
  context.expect_true(one_step.owner_model_identity==one_a.model_identity(),
                      "1-D prepared state records its owner");

  const double radius=swcme::constants::AU_M;
  double density=101.0,velocity=102.0;
  swcme::ModelStatus status=one_b.evaluate_radii_fast_checked(
      one_step,&radius,&density,&velocity,1);
  expect_mismatch(context,status,one_b.model_identity(),one_a.model_identity(),
                  "1-D fast evaluator");
  context.expect_true(density==101.0 && velocity==102.0,
                      "1-D fast outputs remain unchanged");

  double Br=103.0,Bphi=104.0,Bmag=105.0,divV=106.0;
  status=one_b.evaluate_radii_with_B_div_checked(
      one_step,&radius,&density,&velocity,&Br,&Bphi,&Bmag,&divV,1);
  expect_mismatch(context,status,one_b.model_identity(),one_a.model_identity(),
                  "1-D full evaluator");
  context.expect_true(density==101.0 && velocity==102.0 && Br==103.0 &&
                      Bphi==104.0 && Bmag==105.0 && divV==106.0,
                      "1-D full outputs remain unchanged");

  swcme::acceleration::ShockAccelerationState one_source;
  one_source.time_s=107.0;
  one_source.radius_m=108.0;
  status=one_b.shock_acceleration_state_checked(one_step,one_source);
  expect_mismatch(context,status,one_b.model_identity(),one_a.model_identity(),
                  "1-D shock source");
  context.expect_true(one_source.time_s==107.0 && one_source.radius_m==108.0,
                      "1-D source output remains unchanged");

  // An unprepared state carries the reserved identity zero and must be
  // rejected by the same path rather than treated as an all-zero physics cache.
  const swcme1d::StepState unprepared_one;
  density=111.0; velocity=112.0;
  status=one_b.evaluate_radii_fast_checked(
      unprepared_one,&radius,&density,&velocity,1);
  expect_mismatch(context,status,one_b.model_identity(),0,
                  "default-constructed 1-D state");
  context.expect_true(density==111.0 && velocity==112.0,
                      "unprepared-state rejection preserves outputs");

  // Copy construction creates a new logical owner, so copying a model cannot
  // accidentally copy permission to consume the source model's old states.
  swcme1d::Model one_copy(one_a);
  density=113.0; velocity=114.0;
  status=one_copy.evaluate_radii_fast_checked(
      one_step,&radius,&density,&velocity,1);
  expect_mismatch(context,status,one_copy.model_identity(),one_a.model_identity(),
                  "copied 1-D model evaluator");

  // A different receiver configuration exercises the same ownership branch;
  // the implementation must not fall back to comparing parameter hashes.
  swcme1d::Params p1_different=p1;
  p1_different.V_sw_kms+=75.0;
  swcme1d::Model one_c(p1_different);
  density=109.0; velocity=110.0;
  status=one_c.evaluate_radii_fast_checked(
      one_step,&radius,&density,&velocity,1);
  expect_mismatch(context,status,one_c.model_identity(),one_a.model_identity(),
                  "unequal 1-D model evaluator");
  context.expect_true(density==109.0 && velocity==110.0,
                      "unequal-model rejection preserves outputs");

  // Verify checked output is transactional with a real pre-existing file.
  const char* one_path="output/PST02_preserve_1d.dat";
  const std::string one_bytes="PST02-ONE-D-SENTINEL\n";
  { std::ofstream output(one_path,std::ios::binary); output << one_bytes; }
  status=one_b.write_tecplot_radial_profile_checked(
      one_step,&radius,&density,&velocity,&Br,&Bphi,&Bmag,&divV,1,one_path);
  expect_mismatch(context,status,one_b.model_identity(),one_a.model_identity(),
                  "1-D checked writer");
  context.expect_true(read_file(one_path)==one_bytes,
                      "1-D rejected writer preserves destination bytes");
  std::remove(one_path);

  // Exercise the 3-D model's direct background, shock, source, divergence,
  // geometry, mesh, connectivity, and output ownership boundaries.
  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  // Use the same well-conditioned ballistic surface fixture as the mesh
  // validation suite so this ownership test cannot fail for an unrelated
  // difficult oblique Rankine-Hugoniot branch on the rear hemisphere.
  p3.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  p3.r0_Rs=40.0;
  p3.V0_sh_kms=1200.0;
  p3.V_sw_kms=400.0;
  p3.cme_dir[0]=0.73;
  p3.cme_dir[1]=-0.41;
  p3.cme_dir[2]=0.547;
  p3.region_mode=swcme::regions::Mode::ShockOnly;
  p3.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  swcme3d::Model three_a(p3),three_b(p3);
  const swcme3d::StepState three_step=three_a.prepare_step(0.0);
  context.expect_true(three_a.model_identity()!=three_b.model_identity(),
                      "equal 3-D model instances have distinct identities");
  context.expect_true(three_step.owner_model_identity==three_a.model_identity(),
                      "3-D prepared state records its owner");

  const double x=radius,y=0.0,z=0.0;
  double n3=201.0,vx=202.0,vy=203.0,vz=204.0;
  status=three_b.evaluate_cartesian_fast_checked(
      three_step,&x,&y,&z,&n3,&vx,&vy,&vz,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D fast evaluator");
  context.expect_true(n3==201.0 && vx==202.0 && vy==203.0 && vz==204.0,
                      "3-D fast outputs remain unchanged");

  double bx=205.0,by=206.0,bz=207.0;
  status=three_b.evaluate_cartesian_with_B_checked(
      three_step,&x,&y,&z,&n3,&vx,&vy,&vz,&bx,&by,&bz,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D magnetic evaluator");
  context.expect_true(bx==205.0 && by==206.0 && bz==207.0,
                      "3-D magnetic outputs remain unchanged");

  divV=208.0;
  status=three_b.evaluate_cartesian_with_B_div_checked(
      three_step,&x,&y,&z,&n3,&vx,&vy,&vz,&bx,&by,&bz,&divV,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D combined field/divergence evaluator");
  context.expect_true(n3==201.0 && bx==205.0 && divV==208.0,
                      "3-D combined outputs remain unchanged");

  divV=209.0;
  status=three_b.compute_divV_checked(three_step,&x,&y,&z,&divV,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D divergence evaluator");
  context.expect_true(divV==209.0,"3-D divergence output remains unchanged");

  divV=210.0;
  status=three_b.compute_divV_cartesian_checked(
      three_step,&x,&y,&z,&divV,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D Cartesian divergence evaluator");
  context.expect_true(divV==210.0,
                      "3-D Cartesian divergence output remains unchanged");

  divV=211.0;
  status=three_b.compute_divV_radial_checked(
      three_step,&x,&y,&z,&divV,1);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D radial compatibility evaluator");
  context.expect_true(divV==211.0,
                      "3-D radial divergence output remains unchanged");

  const double direction[3]={1.0,0.0,0.0};
  swcme3d::LocalShockState shock;
  shock.surface_exists=true;
  shock.Rdir_m=212.0;
  status=three_b.shock_state_direction_checked(three_step,direction,shock);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D shock evaluator");
  context.expect_true(shock.surface_exists && shock.Rdir_m==212.0,
                      "3-D shock output remains unchanged");

  swcme::acceleration::ShockAccelerationState three_source;
  three_source.time_s=213.0;
  status=three_b.shock_acceleration_state_checked(
      three_step,direction,three_source);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D acceleration source");
  context.expect_true(three_source.time_s==213.0,
                      "3-D acceleration output remains unchanged");

  double Rdir=214.0,normal[3]={215.0,216.0,217.0};
  expect_mismatch_exception(context,[&] {
    (void)three_b.shape_radius_normal(
        three_step,1.0,0.0,0.0,Rdir,normal);
  },"3-D geometry wrapper");
  context.expect_true(Rdir==214.0 && normal[0]==215.0 &&
                      normal[1]==216.0 && normal[2]==217.0,
                      "3-D geometry outputs remain unchanged");

  expect_mismatch_exception(context,[&] {
    (void)three_b.build_shock_mesh(three_step,4,8);
  },"3-D mesh builder");
  const double observer[3]={radius,0.0,0.0};
  expect_mismatch_exception(context,[&] {
    (void)three_b.observer_connectivity(three_step,observer);
  },"3-D connectivity solver");

  const swcme3d::ShockMesh mesh=three_a.build_shock_mesh(three_step,4,8);
  swcme3d::TriMetrics metrics;
  three_a.compute_triangle_metrics(mesh,metrics);
  const swcme3d::BoxSpec box=three_a.default_apex_box(three_step,0.02,2);
  const char* three_path="output/PST02_preserve_3d.dat";
  const std::string three_bytes="PST02-THREE-D-SENTINEL\n";
  { std::ofstream output(three_path,std::ios::binary); output << three_bytes; }
  status=three_b.write_tecplot_dataset_bundle_checked(
      mesh,metrics,three_step,box,three_path);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D checked writer");
  context.expect_true(read_file(three_path)==three_bytes,
                      "3-D rejected writer preserves destination bytes");

  status=three_b.write_box_face_minX_tecplot_structured_checked(
      three_step,box,three_path);
  expect_mismatch(context,status,three_b.model_identity(),three_a.model_identity(),
                  "3-D checked face writer");
  context.expect_true(read_file(three_path)==three_bytes,
                      "3-D rejected face writer preserves destination bytes");
  std::remove(three_path);

  // The public AMPS-facing adapters must preserve destination records too.
  // These calls also prove that state ownership survives through typedefs and
  // is not lost at the interface boundary.
  swcme::sep::Interface1D interface1_a(p1),interface1_b(p1);
  const auto interface1_step=interface1_a.prepare(0.0);
  swcme::sep::BackgroundState background1;
  background1.density_m3=301.0;
  status=interface1_b.evaluate_background(interface1_step,radius,background1);
  expect_mismatch(context,status,interface1_b.model().model_identity(),
                  interface1_a.model().model_identity(),"SEP 1-D background");
  context.expect_true(background1.density_m3==301.0,
                      "SEP 1-D background remains unchanged");
  swcme::sep::SEPSourceState sep1;
  sep1.time_s=302.0;
  status=interface1_b.source_at_shock(interface1_step,sep1);
  expect_mismatch(context,status,interface1_b.model().model_identity(),
                  interface1_a.model().model_identity(),"SEP 1-D source");
  context.expect_true(sep1.time_s==302.0,
                      "SEP 1-D source remains unchanged");

  swcme::sep::Interface3D interface3_a(p3),interface3_b(p3);
  const auto interface3_step=interface3_a.prepare(0.0);
  swcme::sep::BackgroundState background3;
  background3.density_m3=303.0;
  status=interface3_b.evaluate_background(
      interface3_step,{{radius,0.0,0.0}},background3);
  expect_mismatch(context,status,interface3_b.model().model_identity(),
                  interface3_a.model().model_identity(),"SEP 3-D background");
  context.expect_true(background3.density_m3==303.0,
                      "SEP 3-D background remains unchanged");

  swcme::sep::SEPSourceState sep3;
  sep3.time_s=304.0;
  status=interface3_b.source_at_direction(
      interface3_step,{{1.0,0.0,0.0}},sep3);
  expect_mismatch(context,status,interface3_b.model().model_identity(),
                  interface3_a.model().model_identity(),"SEP 3-D source");
  context.expect_true(sep3.time_s==304.0,
                      "SEP 3-D source remains unchanged");

  swcme::sep::SourceSurface surface;
  surface.time_s=305.0;
  surface.patch_count=306;
  status=interface3_b.build_shock_surface_source(
      interface3_step,4,8,surface);
  expect_mismatch(context,status,interface3_b.model().model_identity(),
                  interface3_a.model().model_identity(),"SEP 3-D surface source");
  context.expect_true(surface.time_s==305.0 && surface.patch_count==306,
                      "SEP surface output remains unchanged");

  swcme3d::ConnectivityState connectivity;
  connectivity.observer_radius_m=307.0;
  sep3.time_s=308.0;
  status=interface3_b.source_at_observer_cobpoint(
      interface3_step,{{radius,0.0,0.0}},sep3,&connectivity);
  expect_mismatch(context,status,interface3_b.model().model_identity(),
                  interface3_a.model().model_identity(),"SEP cobpoint source");
  context.expect_true(sep3.time_s==308.0 &&
                      connectivity.observer_radius_m==307.0,
                      "SEP cobpoint outputs remain unchanged");
}
