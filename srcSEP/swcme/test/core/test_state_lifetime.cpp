#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <array>
#include <future>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

// Capture only public production outputs and compare them exactly.  PST05 is
// about ownership/lifetime behavior, not numerical tolerance, so a successful
// relocation must preserve the same bits produced before the move.
struct Probe1D {
  swcme::StatusCode code=swcme::StatusCode::NonFiniteResult;
  double density=0.0;
  double velocity=0.0;
  double Br=0.0;
  double Bphi=0.0;
  double Bmag=0.0;
  double divV=0.0;
};

Probe1D probe(const swcme1d::Model& model,
              const swcme1d::StepState& step) {
  Probe1D out;
  const double radius=swcme::constants::AU_M;
  out.code=model.evaluate_radii_with_B_div_checked(
      step,&radius,&out.density,&out.velocity,&out.Br,&out.Bphi,&out.Bmag,
      &out.divV,1).code;
  return out;
}

bool identical(const Probe1D& left,const Probe1D& right) {
  return left.code==right.code && left.density==right.density &&
         left.velocity==right.velocity && left.Br==right.Br &&
         left.Bphi==right.Bphi && left.Bmag==right.Bmag &&
         left.divV==right.divV;
}

struct Probe3D {
  swcme::StatusCode code=swcme::StatusCode::NonFiniteResult;
  double density=0.0;
  std::array<double,3> velocity{{0.0,0.0,0.0}};
  std::array<double,3> magnetic{{0.0,0.0,0.0}};
  double divV=0.0;
};

Probe3D probe(const swcme3d::Model& model,
              const swcme3d::StepState& step) {
  Probe3D out;
  const double x=swcme::constants::AU_M,y=0.0,z=0.0;
  out.code=model.evaluate_cartesian_with_B_div_checked(
      step,&x,&y,&z,&out.density,&out.velocity[0],&out.velocity[1],
      &out.velocity[2],&out.magnetic[0],&out.magnetic[1],&out.magnetic[2],
      &out.divV,1).code;
  return out;
}

bool identical(const Probe3D& left,const Probe3D& right) {
  return left.code==right.code && left.density==right.density &&
         left.velocity==right.velocity && left.magnetic==right.magnetic &&
         left.divV==right.divV;
}

// Check the deterministic orphan/moved-from rejection contract, including the
// diagnostic IDs required to distinguish lifetime misuse from physics failure.
void expect_owner_mismatch(swcme_test::Context& context,
                           const swcme::ModelStatus& status,
                           swcme::ModelIdentity expected,
                           swcme::ModelIdentity supplied,
                           const std::string& label) {
  context.expect_true(status.code==swcme::StatusCode::StateModelMismatch,
                      label+" reports STATE_MODEL_MISMATCH");
  context.expect_true(status.has_model_identities,
                      label+" carries ownership diagnostics");
  context.expect_true(status.expected_model_identity==expected &&
                          status.supplied_model_identity==supplied,
                      label+" reports live and orphan owner identities");
}

// Compare the transport fields most likely to expose a dangling Params/cache
// reference.  The complete record is already covered field-by-field by PST07;
// PST05 requires those values to survive adapter relocation and async use.
bool identical(const swcme::sep::BackgroundState& left,
               const swcme::sep::BackgroundState& right) {
  return left.status.code==right.status.code &&
         left.owner_model_identity==right.owner_model_identity &&
         left.configuration_digest==right.configuration_digest &&
         left.position_m==right.position_m && left.density_m3==right.density_m3 &&
         left.pressure_Pa==right.pressure_Pa &&
         left.velocity_m_s==right.velocity_m_s &&
         left.magnetic_T==right.magnetic_T &&
         left.magnetic_magnitude_T==right.magnetic_magnitude_T &&
         left.div_velocity_s_inv==right.div_velocity_s_inv &&
         left.focusing_length_m==right.focusing_length_m;
}

}  // namespace

// PST05 defines prepared-state lifetime explicitly.  StepState is a
// self-contained value after owner destruction, but evaluation permission is
// carried by one live logical Model.  Moves relocate that owner; copies and new
// equal models do not.  The matrix below exercises direct and AMPS-facing APIs,
// destruction, both move operations, vector relocation, and asynchronous use.
void test_pst05(swcme_test::Context& context) {
  std::cout << "PST05 prepared-state lifetime contract\n";

  // These compile-time requirements guard the value/relocation properties used
  // by the contract.  The explicit StepState field audit separately establishes
  // that those values contain no back-pointer to Model-owned storage.
  context.expect_true(std::is_copy_constructible<swcme1d::StepState>::value &&
                          std::is_move_constructible<swcme1d::StepState>::value,
                      "1-D prepared state is a copyable and movable value");
  context.expect_true(std::is_copy_constructible<swcme3d::StepState>::value &&
                          std::is_move_constructible<swcme3d::StepState>::value,
                      "3-D prepared state is a copyable and movable value");
  context.expect_true(std::is_nothrow_move_constructible<swcme1d::Model>::value &&
                          std::is_nothrow_move_constructible<swcme3d::Model>::value,
                      "model move construction is non-throwing for relocation");

  // Move construction transfers the prepared owner identity in 1-D.  The old
  // source gets a fresh identity and cannot accidentally remain a co-owner.
  swcme1d::Params one_params;
  swcme1d::Model one_source(one_params);
  const swcme1d::StepState one_step=one_source.prepare_step(3600.0);
  const Probe1D one_before=probe(one_source,one_step);
  const swcme::ModelIdentity one_owner=one_source.model_identity();
  swcme1d::Model one_moved(std::move(one_source));
  context.expect_true(one_moved.model_identity()==one_owner &&
                          one_step.owner_model_identity==one_owner,
                      "1-D move construction transfers logical ownership");
  context.expect_true(one_source.model_identity()!=one_owner,
                      "1-D moved-from object receives a fresh identity");
  context.expect_true(identical(one_before,probe(one_moved,one_step)),
                      "1-D state evaluates identically after move construction");
  double one_n=701.0,one_v=702.0;
  const double one_r=swcme::constants::AU_M;
  swcme::ModelStatus status=one_source.evaluate_radii_fast_checked(
      one_step,&one_r,&one_n,&one_v,1);
  expect_owner_mismatch(context,status,one_source.model_identity(),one_owner,
                        "1-D moved-from model");
  context.expect_true(one_n==701.0 && one_v==702.0,
                      "moved-from rejection preserves 1-D outputs");

  // Move assignment is allowed only into an unprepared destination.  It
  // transfers both identity and frozen phase so the original state follows the
  // source's logical lifetime into the destination.
  swcme3d::Params three_params;
  three_params.shape=swcme3d::ShockShape::Sphere;
  swcme3d::Model three_source(three_params);
  const swcme3d::StepState three_step=three_source.prepare_step(3600.0);
  const Probe3D three_before=probe(three_source,three_step);
  const swcme::ModelIdentity three_owner=three_source.model_identity();
  swcme3d::Model three_target(three_params);
  three_target=std::move(three_source);
  context.expect_true(three_target.model_identity()==three_owner &&
                          three_target.configuration_locked(),
                      "3-D move assignment transfers owner and prepared phase");
  context.expect_true(identical(three_before,probe(three_target,three_step)),
                      "3-D state evaluates identically after move assignment");
  double d=711.0,vx=712.0,vy=713.0,vz=714.0;
  double bx=715.0,by=716.0,bz=717.0,div=718.0;
  const double x=one_r,y=0.0,z=0.0;
  status=three_source.evaluate_cartesian_with_B_div_checked(
      three_step,&x,&y,&z,&d,&vx,&vy,&vz,&bx,&by,&bz,&div,1);
  expect_owner_mismatch(context,status,three_source.model_identity(),three_owner,
                        "3-D moved-from model");
  context.expect_true(d==711.0 && vx==712.0 && vy==713.0 && vz==714.0 &&
                          bx==715.0 && by==716.0 && bz==717.0 && div==718.0,
                      "moved-from rejection preserves 3-D outputs");

  // A frozen destination must not lose its existing states or partially steal
  // the source when move assignment is attempted.
  swcme1d::Model frozen_destination(one_params);
  const swcme1d::StepState frozen_step=frozen_destination.prepare_step(1800.0);
  swcme1d::Model assignment_source(one_params);
  const swcme1d::StepState assignment_step=assignment_source.prepare_step(2400.0);
  const swcme::ModelIdentity frozen_id=frozen_destination.model_identity();
  const swcme::ModelIdentity assignment_id=assignment_source.model_identity();
  bool rejected=false;
  try {
    frozen_destination=std::move(assignment_source);
  } catch (const std::logic_error&) {
    rejected=true;
  }
  context.expect_true(rejected && frozen_destination.model_identity()==frozen_id &&
                          assignment_source.model_identity()==assignment_id,
                      "move assignment into prepared model is transactional");
  context.expect_true(probe(frozen_destination,frozen_step).code==
                          swcme::StatusCode::Ok &&
                          probe(assignment_source,assignment_step).code==
                          swcme::StatusCode::Ok,
                      "failed move leaves both owners and states usable");

  // Destroying an owner cannot dangle the value-owned state.  Its diagnostic
  // fields and private seal remain readable/copyable, while an equal replacement
  // is a different owner and rejects evaluation before touching output.
  swcme1d::StepState orphan_one;
  swcme::ModelIdentity destroyed_one_id=0;
  swcme::ConfigurationDigest orphan_one_seal=0;
  {
    swcme1d::Model temporary(one_params);
    orphan_one=temporary.prepare_step(1200.0);
    destroyed_one_id=temporary.model_identity();
    orphan_one_seal=orphan_one.integrity_digest();
  }
  const swcme1d::StepState orphan_one_copy=orphan_one;
  context.expect_true(orphan_one.owner_model_identity==destroyed_one_id &&
                          orphan_one_copy.integrity_digest()==orphan_one_seal &&
                          swcme1d::prepared_state_integrity(orphan_one_copy)==
                              orphan_one_seal,
                      "1-D state remains a valid inert snapshot after destruction");
  swcme1d::Model replacement_one(one_params);
  one_n=721.0; one_v=722.0;
  status=replacement_one.evaluate_radii_fast_checked(
      orphan_one_copy,&one_r,&one_n,&one_v,1);
  expect_owner_mismatch(context,status,replacement_one.model_identity(),
                        destroyed_one_id,"destroyed 1-D owner replacement");
  context.expect_true(one_n==721.0 && one_v==722.0,
                      "orphan rejection preserves direct outputs");

  // The same destruction behavior is required at the AMPS-facing boundary;
  // retaining PreparedStep is safe, but constructing an equal adapter does not
  // resurrect its former embedded Model identity.
  swcme3d::StepState orphan_adapter_step;
  swcme::ModelIdentity destroyed_adapter_id=0;
  {
    swcme::sep::Interface3D temporary(three_params);
    orphan_adapter_step=temporary.prepare(1200.0);
    destroyed_adapter_id=temporary.model().model_identity();
  }
  swcme::sep::Interface3D replacement_adapter(three_params);
  swcme::sep::BackgroundState preserved_background;
  preserved_background.density_m3=731.0;
  status=replacement_adapter.evaluate_background(
      orphan_adapter_step,{{one_r,0.0,0.0}},preserved_background);
  expect_owner_mismatch(context,status,
                        replacement_adapter.model().model_identity(),
                        destroyed_adapter_id,"destroyed AMPS owner replacement");
  context.expect_true(preserved_background.density_m3==731.0,
                      "orphan rejection preserves AMPS output");

  // A noexcept move constructor makes standard-container growth a supported
  // relocation mechanism.  Test both a direct model and the AMPS adapter so a
  // future special-member change cannot silently select copy construction.
  std::vector<swcme1d::Model> one_models;
  one_models.reserve(1);
  one_models.emplace_back(one_params);
  const swcme1d::StepState vector_one_step=one_models[0].prepare_step(900.0);
  const Probe1D vector_one_before=probe(one_models[0],vector_one_step);
  const swcme::ModelIdentity vector_one_id=one_models[0].model_identity();
  one_models.emplace_back(one_params);  // forces growth beyond reserve(1)
  context.expect_true(one_models[0].model_identity()==vector_one_id &&
                          identical(vector_one_before,
                                    probe(one_models[0],vector_one_step)),
                      "std::vector relocation preserves direct model lifetime");

  std::vector<swcme::sep::Interface3D> adapters;
  adapters.reserve(1);
  adapters.emplace_back(three_params);
  const swcme3d::StepState vector_adapter_step=adapters[0].prepare(900.0);
  swcme::sep::BackgroundState adapter_before;
  status=adapters[0].evaluate_background(
      vector_adapter_step,{{one_r,0.0,0.0}},adapter_before);
  const swcme::ModelIdentity vector_adapter_id=
      adapters[0].model().model_identity();
  adapters.emplace_back(three_params);  // relocates the prepared first adapter
  swcme::sep::BackgroundState adapter_after;
  const swcme::ModelStatus adapter_status=adapters[0].evaluate_background(
      vector_adapter_step,{{one_r,0.0,0.0}},adapter_after);
  context.expect_true(status.ok() && adapter_status.ok() &&
                          adapters[0].model().model_identity()==vector_adapter_id &&
                          identical(adapter_before,adapter_after),
                      "std::vector relocation preserves AMPS adapter lifetime");

  // Prepared-state copies are safe to hand to asynchronous work while the
  // relocated owner remains alive.  Capturing the states by value also checks
  // that no hidden address relationship is required by either public path.
  auto direct_future=std::async(
      std::launch::async,[&one_moved,one_step] { return probe(one_moved,one_step); });
  auto adapter_future=std::async(
      std::launch::async,[&adapters,vector_adapter_step,one_r] {
        swcme::sep::BackgroundState out;
        adapters[0].evaluate_background(
            vector_adapter_step,{{one_r,0.0,0.0}},out);
        return out;
      });
  context.expect_true(identical(one_before,direct_future.get()),
                      "async direct evaluation survives owner relocation");
  context.expect_true(identical(adapter_before,adapter_future.get()),
                      "async AMPS evaluation survives adapter relocation");
}
