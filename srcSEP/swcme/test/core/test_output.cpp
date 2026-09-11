#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_output.hpp>

#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

// FaultSink is a deterministic substitute for a stdio stream.  Each switch
// represents one lifecycle failure that is difficult to trigger portably on a
// real filesystem, while fail_after_bytes can split any formatted record at
// an exact byte.  Counters also prove that cleanup still closes the handle.
struct FaultSink {
  bool fail_open=false;
  std::size_t fail_after_bytes=std::numeric_limits<std::size_t>::max();
  bool fail_flush=false;
  bool report_stream_error=false;
  bool fail_close=false;
  bool opened=false;
  bool flushed=false;
  bool closed=false;
  std::size_t accepted_bytes=0;
};

// These callbacks deliberately use only the public FileOperations ABI.  The
// model writers therefore see exactly the same call sequence as production
// stdio without test-only branches in CheckedTextFile or the physics models.
void* fault_open(void* user_data,const char*,const char*) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  if (sink.fail_open) return nullptr;
  sink.opened=true;
  return &sink;
}

std::size_t fault_write(void* user_data,void*,const char*,std::size_t count) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  if (sink.accepted_bytes>=sink.fail_after_bytes) return 0;
  const std::size_t capacity=sink.fail_after_bytes-sink.accepted_bytes;
  const std::size_t accepted=count<capacity ? count : capacity;
  sink.accepted_bytes+=accepted;
  return accepted;
}

int fault_flush(void* user_data,void*) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  sink.flushed=true;
  return sink.fail_flush ? -1 : 0;
}

int fault_error(void* user_data,void*) {
  const FaultSink& sink=*static_cast<const FaultSink*>(user_data);
  return sink.report_stream_error ? 1 : 0;
}

int fault_close(void* user_data,void*) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  sink.closed=true;
  return sink.fail_close ? -1 : 0;
}

swcme::output::FileOperations operations_for(FaultSink& sink) {
  // Returning the table by value is intentional: CheckedTextFile copies the
  // callback pointers and user-data address, so no global mutable test state
  // is shared across cases or future parallel test execution.
  return {&sink,fault_open,fault_write,fault_flush,fault_error,fault_close};
}

struct OneDimensionalFixture {
  swcme1d::Model model;
  swcme1d::StepState step;
  double radius=swcme::constants::AU_M;
  double density=0.0;
  double velocity=0.0;
  double radial_field=0.0;
  double azimuthal_field=0.0;
  double field_magnitude=0.0;
  double divergence=0.0;

  OneDimensionalFixture() : model(make_params()),step(model.prepare_step(0.0)) {
    model.evaluate_radii_with_B_div(
        step,&radius,&density,&velocity,&radial_field,&azimuthal_field,
        &field_magnitude,&divergence,1);
  }

  static swcme1d::Params make_params() {
    swcme1d::Params params;
    params.r0_Rs=20.0;
    params.V0_sh_kms=1200.0;
    params.V_sw_kms=400.0;
    params.region_mode=swcme::regions::Mode::ShockOnly;
    params.shock_acceleration_mode=swcme::acceleration::Mode::Source;
    return params;
  }

  swcme::ModelStatus write(const swcme::output::FileOperations& operations) {
    return model.write_tecplot_radial_profile_checked(
        step,&radius,&density,&velocity,&radial_field,&azimuthal_field,
        &field_magnitude,&divergence,1,"unused-injected-output.dat",0.0,
        &operations);
  }
};

struct ThreeDimensionalFixture {
  swcme3d::Model model;
  swcme3d::StepState step;
  swcme3d::ShockMesh mesh;
  swcme3d::TriMetrics metrics;
  swcme3d::BoxSpec box;

  ThreeDimensionalFixture()
      : model(make_params()),step(model.prepare_step(0.0)),
        mesh(model.build_shock_mesh(step,4,8)),
        box(model.default_apex_box(step,0.02,2)) {
    model.compute_triangle_metrics(mesh,metrics);
  }

  static swcme3d::Params make_params() {
    // A spherical ballistic fixture is deliberately well conditioned over the
    // complete mesh, keeping OUT02 isolated from nonlinear shock-solver tests.
    swcme3d::Params params;
    params.shape=swcme3d::ShockShape::Sphere;
    params.kinematics_mode=swcme::kinematics::Mode::Ballistic;
    params.r0_Rs=40.0;
    params.V0_sh_kms=1200.0;
    params.V_sw_kms=400.0;
    params.cme_dir[0]=0.73;
    params.cme_dir[1]=-0.41;
    params.cme_dir[2]=0.547;
    params.region_mode=swcme::regions::Mode::ShockOnly;
    params.shock_acceleration_mode=swcme::acceleration::Mode::Source;
    return params;
  }

  swcme::ModelStatus write_face(
      const swcme::output::FileOperations& operations) {
    return model.write_box_face_minX_tecplot_structured_checked(
        step,box,"unused-injected-face.dat",&operations);
  }
};

void expect_write_failure(swcme_test::Context& context,
                          const swcme::ModelStatus& status,
                          std::size_t expected_offset,
                          const std::string& label) {
  context.expect_true(status.code==swcme::StatusCode::FileWriteFailure,
                      label+" reports FILE_WRITE_FAILURE");
  context.expect_true(status.has_io_byte_offset,
                      label+" carries an output byte offset");
  context.expect_true(status.io_byte_offset==expected_offset,
                      label+" reports the first unaccepted byte");
  context.expect_true(status.summary().find("io_byte_offset=")!=
                          std::string::npos,
                      label+" summary prints the byte offset");
}

template <typename Writer>
void expect_dev_full_failure(swcme_test::Context& context,Writer&& writer,
                             const std::string& label) {
  const swcme::ModelStatus status=writer();
  context.expect_true(status.code==swcme::StatusCode::FileWriteFailure,
                      label+" detects /dev/full write failure");
  context.expect_true(status.has_io_byte_offset,
                      label+" reports /dev/full byte context");
}

}  // namespace

// OUT02: every output layer must distinguish open failure from data-loss
// failure and must inspect delayed stdio failures at flush, stream-error, and
// close.  The injected matrix checks exact diagnostics; /dev/full then proves
// the production stdio backend propagates a real operating-system failure for
// every public Tecplot product.
void test_out02(swcme_test::Context& context) {
  std::cout << "OUT02 write failure detection and propagation\n";
  OneDimensionalFixture one;
  ThreeDimensionalFixture three;

  // A complete injected run guards against false-positive failure handling in
  // both the inline 1-D writer and the separately compiled 3-D writer.
  FaultSink one_success;
  swcme::output::FileOperations one_success_ops=operations_for(one_success);
  swcme::ModelStatus status=one.write(one_success_ops);
  context.expect_true(status.ok() && one_success.opened && one_success.flushed &&
                          one_success.closed && one_success.accepted_bytes>0,
                      "1-D injected writer completes the full lifecycle");

  FaultSink three_success;
  swcme::output::FileOperations three_success_ops=operations_for(three_success);
  status=three.write_face(three_success_ops);
  context.expect_true(status.ok() && three_success.opened &&
                          three_success.flushed && three_success.closed &&
                          three_success.accepted_bytes>0,
                      "3-D injected writer completes the full lifecycle");

  FaultSink open_failure;
  open_failure.fail_open=true;
  swcme::output::FileOperations open_failure_ops=operations_for(open_failure);
  status=one.write(open_failure_ops);
  context.expect_true(status.code==swcme::StatusCode::FileOpenFailure,
                      "open rejection reports FILE_OPEN_FAILURE");
  context.expect_true(!status.has_io_byte_offset && !open_failure.closed,
                      "open rejection is not mislabeled as a write failure");

  FaultSink first_write;
  first_write.fail_after_bytes=0;
  swcme::output::FileOperations first_write_ops=operations_for(first_write);
  status=one.write(first_write_ops);
  expect_write_failure(context,status,0,"first formatted write");
  context.expect_true(std::string(status.context).find("title")!=std::string::npos,
                      "first-write diagnostic identifies the title phase");
  context.expect_true(first_write.closed,
                      "first-write failure still closes the output handle");

  FaultSink partial_write;
  partial_write.fail_after_bytes=57;
  partial_write.fail_close=true;
  swcme::output::FileOperations partial_write_ops=operations_for(partial_write);
  status=one.write(partial_write_ops);
  expect_write_failure(context,status,57,"selected-byte partial write");
  context.expect_true(std::string(status.context).find("variables")!=
                          std::string::npos,
                      "partial-write diagnostic preserves record context");
  context.expect_true(partial_write.closed,
                      "close is attempted after a partial write");

  FaultSink row_write;
  row_write.fail_after_bytes=250;
  swcme::output::FileOperations row_write_ops=operations_for(row_write);
  status=one.write(row_write_ops);
  expect_write_failure(context,status,250,"1-D data-row partial write");
  context.expect_true(status.sample_index==0 &&
                          std::string(status.context).find("row")!=
                              std::string::npos,
                      "row failure records its row index and phase");

  FaultSink three_partial;
  three_partial.fail_after_bytes=41;
  swcme::output::FileOperations three_partial_ops=operations_for(three_partial);
  status=three.write_face(three_partial_ops);
  expect_write_failure(context,status,41,"compiled 3-D partial write");
  context.expect_true(three_partial.closed,
                      "compiled 3-D failure closes the output handle");

  FaultSink flush_failure;
  flush_failure.fail_flush=true;
  swcme::output::FileOperations flush_failure_ops=operations_for(flush_failure);
  status=one.write(flush_failure_ops);
  expect_write_failure(context,status,flush_failure.accepted_bytes,
                       "flush failure");
  context.expect_true(std::string(status.context).find("flush")!=std::string::npos,
                      "flush diagnostic identifies lifecycle phase");
  context.expect_true(flush_failure.closed,
                      "flush failure still closes the output handle");

  FaultSink stream_failure;
  stream_failure.report_stream_error=true;
  swcme::output::FileOperations stream_failure_ops=operations_for(stream_failure);
  status=one.write(stream_failure_ops);
  expect_write_failure(context,status,stream_failure.accepted_bytes,
                       "persistent stream error");
  context.expect_true(std::string(status.context).find("stream error")!=
                          std::string::npos,
                      "stream diagnostic identifies lifecycle phase");
  context.expect_true(stream_failure.closed,
                      "stream error still closes the output handle");

  FaultSink close_failure;
  close_failure.fail_close=true;
  swcme::output::FileOperations close_failure_ops=operations_for(close_failure);
  status=one.write(close_failure_ops);
  expect_write_failure(context,status,close_failure.accepted_bytes,
                       "close failure");
  context.expect_true(std::string(status.context).find("close")!=std::string::npos,
                      "close diagnostic identifies lifecycle phase");
  context.expect_true(close_failure.closed,
                      "failed close was attempted exactly at cleanup");

  // /dev/full is the operating-system integration check for buffered stdio:
  // writes may appear successful until fflush, so a passing result here proves
  // production code checks more than fopen and formatted-write return values.
  std::FILE* full_probe=std::fopen("/dev/full","w");
  if (!full_probe) {
    context.record_skip();
  } else {
    (void)std::fclose(full_probe);
    expect_dev_full_failure(context,[&] {
      return one.model.write_tecplot_radial_profile_checked(
          one.step,&one.radius,&one.density,&one.velocity,&one.radial_field,
          &one.azimuthal_field,&one.field_magnitude,&one.divergence,1,
          "/dev/full",0.0);
    },"1-D radial profile");
    expect_dev_full_failure(context,[&] {
      return one.model.write_tecplot_shock_vs_time_checked(10.0,2,"/dev/full");
    },"1-D shock history");
    expect_dev_full_failure(context,[&] {
      return three.model.write_shock_surface_center_metrics_tecplot_checked(
          three.mesh,three.metrics,"/dev/full");
    },"3-D surface");
    expect_dev_full_failure(context,[&] {
      return three.model.write_tecplot_dataset_bundle_checked(
          three.mesh,three.metrics,three.step,three.box,"/dev/full");
    },"3-D dataset bundle");
    expect_dev_full_failure(context,[&] {
      return three.model.write_box_face_minX_tecplot_structured_checked(
          three.step,three.box,"/dev/full");
    },"3-D box face");

    // The source-compatible bool APIs intentionally discard detailed status,
    // but they must still delegate to the checked implementation and return
    // false rather than reproducing the historical false-success behavior.
    context.expect_true(!one.model.write_tecplot_radial_profile(
        one.step,&one.radius,&one.density,&one.velocity,&one.radial_field,
        &one.azimuthal_field,&one.field_magnitude,&one.divergence,1,
        "/dev/full",0.0),"legacy 1-D radial writer rejects /dev/full");
    context.expect_true(!one.model.write_tecplot_shock_vs_time(
        10.0,2,"/dev/full"),"legacy 1-D history writer rejects /dev/full");
    context.expect_true(
        !three.model.write_shock_surface_center_metrics_tecplot(
            three.mesh,three.metrics,"/dev/full"),
        "legacy 3-D surface writer rejects /dev/full");
    context.expect_true(!three.model.write_tecplot_dataset_bundle(
        three.mesh,three.metrics,three.step,three.box,"/dev/full"),
        "legacy 3-D bundle writer rejects /dev/full");
    context.expect_true(!three.model.write_box_face_minX_tecplot_structured(
        three.step,three.box,"/dev/full"),
        "legacy 3-D face writer rejects /dev/full");
  }
}
