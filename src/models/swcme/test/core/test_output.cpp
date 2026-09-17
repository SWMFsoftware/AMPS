#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_output.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
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
  unsigned open_failures_remaining=0;
  std::size_t fail_after_bytes=std::numeric_limits<std::size_t>::max();
  bool fail_flush=false;
  bool report_stream_error=false;
  bool fail_close=false;
  bool fail_commit=false;
  bool fail_remove=false;
  bool opened=false;
  unsigned open_attempts=0;
  bool flushed=false;
  bool closed=false;
  bool commit_called=false;
  bool remove_called=false;
  std::size_t accepted_bytes=0;
  std::string opened_path;
  std::string opened_mode;
  std::string committed_destination;
  std::string temporary_bytes;
  std::string destination_bytes="PREVIOUS-COMPLETE-OUTPUT\n";
};

// These callbacks deliberately use only the public FileOperations ABI.  The
// model writers therefore see exactly the same call sequence as production
// stdio without test-only branches in CheckedTextFile or the physics models.
void* fault_open(void* user_data,const char* path,const char* mode) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  ++sink.open_attempts;
  if (sink.fail_open) return nullptr;
  if (sink.open_failures_remaining>0) {
    --sink.open_failures_remaining;
    return nullptr;
  }
  sink.opened=true;
  sink.opened_path=path ? path : "";
  sink.opened_mode=mode ? mode : "";
  sink.temporary_bytes.clear();
  return &sink;
}

std::size_t fault_write(void* user_data,void*,const char* bytes,
                        std::size_t count) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  if (sink.accepted_bytes>=sink.fail_after_bytes) return 0;
  const std::size_t capacity=sink.fail_after_bytes-sink.accepted_bytes;
  const std::size_t accepted=count<capacity ? count : capacity;
  // Retain exactly the accepted prefix so OUT03 can verify that commit sees a
  // complete staged product and that failure cleanup never modifies the
  // simulated pre-existing destination.
  sink.temporary_bytes.append(bytes,accepted);
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

int fault_commit(void* user_data,const char*,const char* destination_path) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  sink.commit_called=true;
  sink.committed_destination=destination_path ? destination_path : "";
  if (sink.fail_commit) return -1;
  sink.destination_bytes=sink.temporary_bytes;
  sink.temporary_bytes.clear();
  return 0;
}

int fault_remove(void* user_data,const char*) {
  FaultSink& sink=*static_cast<FaultSink*>(user_data);
  sink.remove_called=true;
  if (sink.fail_remove) return -1;
  sink.temporary_bytes.clear();
  return 0;
}

swcme::output::FileOperations operations_for(FaultSink& sink) {
  // Returning the table by value is intentional: CheckedTextFile copies the
  // callback pointers and user-data address, so no global mutable test state
  // is shared across cases or future parallel test execution.
  return {&sink,fault_open,fault_write,fault_flush,fault_error,fault_close,
          fault_commit,fault_remove};
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

std::string read_file(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(input),
                     std::istreambuf_iterator<char>());
}

void write_file(const std::filesystem::path& path,const std::string& bytes) {
  std::ofstream output(path,std::ios::binary|std::ios::trunc);
  output << bytes;
}

std::size_t count_staging_files(const std::filesystem::path& destination) {
  // Transaction names are siblings formed from the destination filename, so
  // a directory scan can detect leaked staging files without depending on the
  // process-global sequence number used to make each name unique.
  const std::filesystem::path directory=destination.parent_path().empty()
      ? std::filesystem::path(".") : destination.parent_path();
  const std::string prefix=destination.filename().string()+".swcme-tmp-";
  std::size_t count=0;
  for (const std::filesystem::directory_entry& entry :
       std::filesystem::directory_iterator(directory)) {
    if (entry.path().filename().string().rfind(prefix,0)==0) ++count;
  }
  return count;
}

}  // namespace

// OUT02: every output layer must distinguish open failure from data-loss
// failure and must inspect delayed stdio failures at flush, stream-error, and
// close.  The injected matrix checks exact diagnostics; /dev/full then proves
// the production stdio backend propagates a real operating-system failure at
// the shared byte-stream layer used beneath every public Tecplot product.
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

  // /dev/full is tested directly at the shared stream layer.  OUT03 model
  // writers stage regular files rather than opening their public destination,
  // so attempting to route them through a device would test commit policy
  // instead of the buffered-write failure that OUT02 is designed to isolate.
  std::FILE* full_probe=std::fopen("/dev/full","w");
  if (!full_probe) {
    context.record_skip();
  } else {
    (void)std::fclose(full_probe);
    swcme::output::CheckedTextFile full_output(
        swcme::output::stdio_file_operations());
    context.expect_true(full_output.open("/dev/full"),
                        "/dev/full opens for the OUT02 stream probe");
    (void)full_output.print("OUT02 /dev/full record",0,
                            "buffered output probe %d\n",42);
    const swcme::ModelStatus full_status=full_output.finish(
        "OUT02 /dev/full flush","OUT02 /dev/full stream error",
        "OUT02 /dev/full close");
    context.expect_true(
        full_status.code==swcme::StatusCode::FileWriteFailure &&
            full_status.has_io_byte_offset,
        "stdio stream detects delayed /dev/full failure");
  }
}

// OUT03: public products are invisible until every byte and lifecycle check
// succeeds.  The deterministic backend proves commit ordering and preservation
// without timing assumptions; real-filesystem checks then verify replacement
// and cleanup behavior with the production backend.
void test_out03(swcme_test::Context& context) {
  std::cout << "OUT03 transactional output commit\n";
  OneDimensionalFixture one;
  ThreeDimensionalFixture three;
  const std::string sentinel="PREVIOUS-COMPLETE-OUTPUT\n";

  FaultSink success;
  success.destination_bytes=sentinel;
  swcme::output::FileOperations success_ops=operations_for(success);
  swcme::ModelStatus status=one.write(success_ops);
  context.expect_true(status.ok() && success.commit_called,
                      "complete staged output is committed exactly once");
  context.expect_true(success.opened_mode=="wx" &&
                          success.opened_path.find(
                              "unused-injected-output.dat.swcme-tmp-")==0,
                      "writer exclusively opens a same-directory staging name");
  context.expect_true(success.committed_destination==
                          "unused-injected-output.dat" &&
                          success.destination_bytes.find(
                              "TITLE=\"1D SW+CME radial profile\"")==0,
                      "commit replaces the destination with complete output");
  context.expect_true(!success.remove_called && success.temporary_bytes.empty(),
                      "successful commit needs no failure cleanup");

  FaultSink collision;
  collision.destination_bytes=sentinel;
  collision.open_failures_remaining=1;
  swcme::output::FileOperations collision_ops=operations_for(collision);
  status=one.write(collision_ops);
  context.expect_true(status.ok() && collision.open_attempts==2 &&
                          collision.commit_called,
                      "exclusive staging retries after one name collision");

  // Exercise every remaining public product through the injected transaction
  // boundary.  The detailed ordering checks above need not be duplicated, but
  // each writer must demonstrably request commit after producing nonempty data.
  FaultSink history_success;
  swcme::output::FileOperations history_success_ops=
      operations_for(history_success);
  status=one.model.write_tecplot_shock_vs_time_checked(
      10.0,2,"unused-injected-history.dat",&history_success_ops);
  context.expect_true(status.ok() && history_success.commit_called &&
                          history_success.destination_bytes.find(
                              "TITLE=\"Shock kinematics vs time\"")==0,
                      "1-D shock history commits a complete transaction");

  FaultSink surface_success;
  swcme::output::FileOperations surface_success_ops=
      operations_for(surface_success);
  status=three.model.write_shock_surface_center_metrics_tecplot_checked(
      three.mesh,three.metrics,"unused-injected-surface.dat",
      &surface_success_ops);
  context.expect_true(status.ok() && surface_success.commit_called &&
                          surface_success.destination_bytes.find(
                              "TITLE=\"Shock surface")==0,
                      "3-D surface commits a complete transaction");

  FaultSink bundle_success;
  swcme::output::FileOperations bundle_success_ops=
      operations_for(bundle_success);
  status=three.model.write_tecplot_dataset_bundle_checked(
      three.mesh,three.metrics,three.step,three.box,
      "unused-injected-bundle.dat",&bundle_success_ops);
  context.expect_true(status.ok() && bundle_success.commit_called &&
                          bundle_success.destination_bytes.find(
                              "TITLE = \"SW+CME dataset\"")==0,
                      "3-D dataset bundle commits a complete transaction");

  FaultSink face_success;
  swcme::output::FileOperations face_success_ops=operations_for(face_success);
  status=three.write_face(face_success_ops);
  context.expect_true(status.ok() && face_success.commit_called &&
                          face_success.destination_bytes.find(
                              "TITLE = \"Box face (minX)\"")==0,
                      "3-D standalone face commits a complete transaction");

  FaultSink write_failure;
  write_failure.destination_bytes=sentinel;
  write_failure.fail_after_bytes=57;
  swcme::output::FileOperations write_failure_ops=operations_for(write_failure);
  status=one.write(write_failure_ops);
  expect_write_failure(context,status,57,"transactional partial write");
  context.expect_true(!write_failure.commit_called &&
                          write_failure.remove_called &&
                          write_failure.destination_bytes==sentinel &&
                          write_failure.temporary_bytes.empty(),
                      "partial output is removed without changing destination");

  FaultSink close_failure;
  close_failure.destination_bytes=sentinel;
  close_failure.fail_close=true;
  swcme::output::FileOperations close_failure_ops=operations_for(close_failure);
  status=one.write(close_failure_ops);
  expect_write_failure(context,status,close_failure.accepted_bytes,
                       "transactional close failure");
  context.expect_true(!close_failure.commit_called &&
                          close_failure.remove_called &&
                          close_failure.destination_bytes==sentinel,
                      "close failure cannot publish staged output");

  FaultSink commit_failure;
  commit_failure.destination_bytes=sentinel;
  commit_failure.fail_commit=true;
  swcme::output::FileOperations commit_failure_ops=operations_for(commit_failure);
  status=three.write_face(commit_failure_ops);
  context.expect_true(status.code==swcme::StatusCode::FileCommitFailure &&
                          status.has_io_byte_offset &&
                          status.io_byte_offset==commit_failure.accepted_bytes,
                      "failed atomic replacement reports FILE_COMMIT_FAILURE");
  context.expect_true(std::string(status.context).find("commit")!=
                          std::string::npos &&
                          status.summary().find("FILE_COMMIT_FAILURE")!=
                              std::string::npos,
                      "commit diagnostic identifies phase and staged size");
  context.expect_true(commit_failure.commit_called &&
                          commit_failure.remove_called &&
                          commit_failure.destination_bytes==sentinel &&
                          commit_failure.temporary_bytes.empty(),
                      "failed commit preserves destination and removes staging");

  FaultSink open_failure;
  open_failure.destination_bytes=sentinel;
  open_failure.fail_open=true;
  swcme::output::FileOperations open_failure_ops=operations_for(open_failure);
  status=one.write(open_failure_ops);
  context.expect_true(status.code==swcme::StatusCode::FileOpenFailure &&
                          !open_failure.commit_called &&
                          !open_failure.remove_called &&
                          open_failure.destination_bytes==sentinel,
                      "staging-open failure leaves destination untouched");

  FaultSink cleanup_failure;
  cleanup_failure.destination_bytes=sentinel;
  cleanup_failure.fail_after_bytes=57;
  cleanup_failure.fail_remove=true;
  swcme::output::FileOperations cleanup_failure_ops=operations_for(cleanup_failure);
  status=one.write(cleanup_failure_ops);
  context.expect_true(status.code==swcme::StatusCode::FileWriteFailure &&
                          cleanup_failure.remove_called &&
                          !cleanup_failure.commit_called &&
                          cleanup_failure.destination_bytes==sentinel,
                      "cleanup failure cannot mask the original write failure");

  const std::filesystem::path file_path="output/OUT03_transaction.dat";
  const std::filesystem::path directory_path=
      "output/OUT03_nonregular_destination";
  std::error_code cleanup_error;
  std::filesystem::remove(file_path,cleanup_error);
  std::filesystem::remove(directory_path,cleanup_error);
  write_file(file_path,sentinel);

  status=one.model.write_tecplot_radial_profile_checked(
      one.step,&one.radius,&one.density,&one.velocity,&one.radial_field,
      &one.azimuthal_field,&one.field_magnitude,&one.divergence,1,
      file_path.string().c_str(),0.0);
  const std::string committed_bytes=read_file(file_path);
  context.expect_true(status.ok() && committed_bytes!=sentinel &&
                          committed_bytes.find(
                              "TITLE=\"1D SW+CME radial profile\"")==0,
                      "production rename replaces an existing regular file");
  context.expect_true(count_staging_files(file_path)==0,
                      "successful production commit leaves no staging file");

  const std::filesystem::path new_file_path="output/OUT03_new_product.dat";
  std::filesystem::remove(new_file_path,cleanup_error);
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,three.box,new_file_path.string().c_str());
  context.expect_true(status.ok() && std::filesystem::is_regular_file(
                          new_file_path) &&
                          read_file(new_file_path).find(
                              "TITLE = \"Box face (minX)\"")==0,
                      "production commit installs a new destination file");
  context.expect_true(count_staging_files(new_file_path)==0,
                      "new-file commit leaves no staging file");

  // Legacy boolean APIs cannot expose FILE_COMMIT_FAILURE details, but they
  // must retain the same publish-on-success transaction as checked callers.
  const std::filesystem::path legacy_path="output/OUT03_legacy_product.dat";
  write_file(legacy_path,sentinel);
  const bool legacy_ok=one.model.write_tecplot_radial_profile(
      one.step,&one.radius,&one.density,&one.velocity,&one.radial_field,
      &one.azimuthal_field,&one.field_magnitude,&one.divergence,1,
      legacy_path.string().c_str(),0.0);
  context.expect_true(legacy_ok && read_file(legacy_path)!=sentinel &&
                          count_staging_files(legacy_path)==0,
                      "legacy writer commits one complete replacement");

  std::filesystem::create_directory(directory_path);
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,three.box,directory_path.string().c_str());
  context.expect_true(status.code==swcme::StatusCode::FileCommitFailure &&
                          std::filesystem::is_directory(directory_path),
                      "production commit refuses to replace a nonregular target");
  context.expect_true(count_staging_files(directory_path)==0,
                      "rejected production commit removes its staging file");
  context.expect_true(
      !three.model.write_box_face_minX_tecplot_structured(
          three.step,three.box,directory_path.string().c_str()) &&
          std::filesystem::is_directory(directory_path) &&
          count_staging_files(directory_path)==0,
      "legacy writer returns false without replacing nonregular target");

  std::filesystem::remove(file_path,cleanup_error);
  std::filesystem::remove(new_file_path,cleanup_error);
  std::filesystem::remove(legacy_path,cleanup_error);
  std::filesystem::remove(directory_path,cleanup_error);
}

// OUT05: every point that a model writer will emit must be finite and inside
// the shared solar-wind domain before the output backend is touched.  The
// injected sink makes "no filesystem side effect" observable as zero open,
// commit, and remove callbacks, while one production-file check verifies the
// same byte-preservation and no-staging guarantee on the real backend.
void test_out05(swcme_test::Context& context) {
  std::cout << "OUT05 model-domain output preflight\n";
  OneDimensionalFixture one;
  ThreeDimensionalFixture three;
  const double minimum=swcme::solarwind::MIN_RADIUS_M;
  const std::string sentinel="PREVIOUS-COMPLETE-OUTPUT\n";

  // Put the invalid radius in the middle of an otherwise valid profile.  This
  // proves the implementation scans the complete request rather than checking
  // only its first point before opening the transaction.
  const double radii[]={one.radius,0.5*minimum,one.radius};
  const double density[]={one.density,one.density,one.density};
  const double velocity[]={one.velocity,one.velocity,one.velocity};
  const double radial_field[]={one.radial_field,one.radial_field,
                               one.radial_field};
  const double azimuthal_field[]={one.azimuthal_field,one.azimuthal_field,
                                  one.azimuthal_field};
  const double magnitude[]={one.field_magnitude,one.field_magnitude,
                            one.field_magnitude};
  const double divergence[]={one.divergence,one.divergence,one.divergence};
  FaultSink radial_domain;
  radial_domain.destination_bytes=sentinel;
  swcme::output::FileOperations radial_domain_ops=operations_for(radial_domain);
  swcme::ModelStatus status=one.model.write_tecplot_radial_profile_checked(
      one.step,radii,density,velocity,radial_field,azimuthal_field,magnitude,
      divergence,3,"unused-out05-radial.dat",0.0,&radial_domain_ops);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain &&
                          status.sample_index==1 &&
                          status.has_offending_value &&
                          status.offending_value==radii[1],
                      "1-D profile identifies the first out-of-domain row");
  context.expect_true(radial_domain.open_attempts==0 &&
                          !radial_domain.commit_called &&
                          !radial_domain.remove_called &&
                          radial_domain.destination_bytes==sentinel,
                      "1-D domain rejection occurs before any output callback");

  // Non-finite caller coordinates and precomputed fields are distinct failure
  // classes.  Both must retain their row and fail before an injected open.
  double nonfinite_radii[]={one.radius,one.radius,
                            std::numeric_limits<double>::quiet_NaN()};
  FaultSink nonfinite_coordinate;
  swcme::output::FileOperations nonfinite_coordinate_ops=
      operations_for(nonfinite_coordinate);
  status=one.model.write_tecplot_radial_profile_checked(
      one.step,nonfinite_radii,density,velocity,radial_field,azimuthal_field,
      magnitude,divergence,3,"unused-out05-coordinate.dat",0.0,
      &nonfinite_coordinate_ops);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput &&
                          status.sample_index==2 &&
                          nonfinite_coordinate.open_attempts==0,
                      "1-D non-finite radius is rejected as input before open");

  double bad_density[]={one.density,
                        std::numeric_limits<double>::infinity(),one.density};
  FaultSink nonfinite_field;
  swcme::output::FileOperations nonfinite_field_ops=operations_for(nonfinite_field);
  status=one.model.write_tecplot_radial_profile_checked(
      one.step,nonfinite_radii,bad_density,velocity,radial_field,
      azimuthal_field,magnitude,divergence,2,"unused-out05-field.dat",0.0,
      &nonfinite_field_ops);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteResult &&
                          status.sample_index==1 &&
                          nonfinite_field.open_attempts==0,
                      "1-D non-finite field is rejected as a result before open");

  // A data-driven history can fail only at a late requested time.  Preparing
  // the entire immutable series first prevents the earlier valid samples from
  // causing even a temporary output artifact.
  swcme1d::Params history_params=OneDimensionalFixture::make_params();
  history_params.kinematics_mode=swcme::kinematics::Mode::DataDriven;
  history_params.data_time_s={0.0,5.0};
  history_params.data_radius_Rs={40.0,41.0};
  swcme1d::Model history_model(history_params);
  FaultSink history_domain;
  swcme::output::FileOperations history_domain_ops=operations_for(history_domain);
  status=history_model.write_tecplot_shock_vs_time_checked(
      10.0,3,"unused-out05-history.dat",&history_domain_ops);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteResult &&
                          status.sample_index==2 &&
                          status.has_offending_value &&
                          status.offending_value==10.0 &&
                          history_domain.open_attempts==0,
                      "late invalid history time is identified before open");

  // Surface output does not evaluate the ambient background, so it needs an
  // explicit node-domain pass.  Corrupting a non-first vertex confirms the
  // pass records the mesh-node index and precedes metric/output work.
  swcme3d::ShockMesh bad_mesh=three.mesh;
  bad_mesh.x[1]=0.5*minimum;
  bad_mesh.y[1]=0.0;
  bad_mesh.z[1]=0.0;
  swcme3d::TriMetrics bad_mesh_metrics;
  three.model.compute_triangle_metrics(bad_mesh,bad_mesh_metrics);
  FaultSink surface_domain;
  swcme::output::FileOperations surface_domain_ops=operations_for(surface_domain);
  status=three.model.write_shock_surface_center_metrics_tecplot_checked(
      bad_mesh,bad_mesh_metrics,"unused-out05-surface.dat",&surface_domain_ops);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain &&
                          status.sample_index==1 &&
                          surface_domain.open_attempts==0,
                      "3-D surface rejects an invalid late vertex before open");

  // Endpoints at twice the lower bound are valid, but the middle point of this
  // three-point line crosses the origin.  OUT05 must report flattened volume
  // row one before the bundle opens or begins evaluating physics.
  swcme3d::BoxSpec crossing_box;
  crossing_box.cx=0.0;
  crossing_box.cy=0.0;
  crossing_box.cz=0.0;
  crossing_box.hx=2.0*minimum;
  crossing_box.hy=0.0;
  crossing_box.hz=0.0;
  crossing_box.Ni=3;
  crossing_box.Nj=2;
  crossing_box.Nk=2;
  FaultSink bundle_domain;
  swcme::output::FileOperations bundle_domain_ops=operations_for(bundle_domain);
  status=three.model.write_tecplot_dataset_bundle_checked(
      three.mesh,three.metrics,three.step,crossing_box,
      "unused-out05-bundle.dat",&bundle_domain_ops);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain &&
                          status.sample_index==1 &&
                          std::string(status.context).find("volume")!=
                              std::string::npos &&
                          bundle_domain.open_attempts==0,
                      "3-D bundle preflights every structured volume point");

  swcme3d::BoxSpec face_box=three.box;
  face_box.cx=0.5*minimum;
  face_box.cy=0.0;
  face_box.cz=0.0;
  face_box.hx=0.0;
  face_box.hy=0.0;
  face_box.hz=0.0;
  FaultSink face_domain;
  swcme::output::FileOperations face_domain_ops=operations_for(face_domain);
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,face_box,"unused-out05-face.dat",&face_domain_ops);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain &&
                          status.sample_index==0 &&
                          face_domain.open_attempts==0,
                      "3-D face rejects its first invalid generated point");
  context.expect_true(
      !three.model.write_box_face_minX_tecplot_structured(
          three.step,face_box,"unused-out05-legacy-face.dat"),
      "legacy face writer delegates to the same domain preflight");

  // The lower boundary is inclusive in the public evaluator contract.  A
  // successful exact-boundary write guards against accidentally tightening
  // OUT05 to radius > MIN_RADIUS_M.
  swcme3d::BoxSpec boundary_box=face_box;
  boundary_box.cx=minimum;
  FaultSink boundary_success;
  swcme::output::FileOperations boundary_success_ops=
      operations_for(boundary_success);
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,boundary_box,"unused-out05-boundary.dat",
      &boundary_success_ops);
  context.expect_true(status.ok() && boundary_success.open_attempts==1 &&
                          boundary_success.commit_called,
                      "exact domain-boundary face passes and commits");

  // Finally exercise the production backend against an existing destination.
  // A true preflight leaves both its bytes and sibling namespace untouched.
  const std::filesystem::path preserved_path="output/OUT05_preserved.dat";
  std::error_code cleanup_error;
  std::filesystem::remove(preserved_path,cleanup_error);
  write_file(preserved_path,sentinel);
  status=one.model.write_tecplot_radial_profile_checked(
      one.step,radii,density,velocity,radial_field,azimuthal_field,magnitude,
      divergence,3,preserved_path.string().c_str(),0.0);
  context.expect_true(status.code==swcme::StatusCode::OutsideModelDomain &&
                          read_file(preserved_path)==sentinel &&
                          count_staging_files(preserved_path)==0,
                      "production preflight preserves destination and creates no staging");
  std::filesystem::remove(preserved_path,cleanup_error);
}

// OUT04: BoxSpec is a structural input contract, distinct from OUT05's scan of
// otherwise valid generated points.  Every malformed-box case uses FaultSink
// to prove validation completes before opening, and the matrix deliberately
// covers both bundle and standalone-face entry points plus the box factory.
void test_out04(swcme_test::Context& context) {
  std::cout << "OUT04 box specification validation\n";
  ThreeDimensionalFixture three;
  const double minimum=swcme::solarwind::MIN_RADIUS_M;
  const std::string sentinel="PREVIOUS-COMPLETE-OUTPUT\n";

  // Run a standalone-face rejection while retaining its complete observable
  // backend state.  Returning the status lets each caller verify the precise
  // classification without duplicating the no-open/no-commit assertions.
  const auto reject_face=[&](const swcme3d::BoxSpec& box,
                             FaultSink& sink) {
    sink.destination_bytes=sentinel;
    swcme::output::FileOperations operations=operations_for(sink);
    const swcme::ModelStatus status=
        three.model.write_box_face_minX_tecplot_structured_checked(
            three.step,box,"unused-out04-face.dat",&operations);
    context.expect_true(sink.open_attempts==0 && !sink.commit_called &&
                            !sink.remove_called &&
                            sink.destination_bytes==sentinel,
                        "invalid BoxSpec causes no output callback");
    return status;
  };

  swcme3d::BoxSpec invalid=three.box;
  invalid.cx=std::numeric_limits<double>::quiet_NaN();
  FaultSink nonfinite_center;
  swcme::ModelStatus status=reject_face(invalid,nonfinite_center);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput &&
                          status.has_offending_value &&
                          std::isnan(status.offending_value),
                      "non-finite box center is rejected as input");

  invalid=three.box;
  invalid.hy=std::numeric_limits<double>::infinity();
  FaultSink nonfinite_extent;
  status=reject_face(invalid,nonfinite_extent);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput &&
                          status.has_offending_value &&
                          std::isinf(status.offending_value),
                      "non-finite box half extent is rejected as input");

  // Negative half sizes invert grid orientation and do not describe the
  // documented center-plus-half-size box.  Check the bundle specifically so
  // both BoxSpec-consuming checked APIs are held to the shared validator.
  invalid=three.box;
  invalid.hx=-1.0;
  FaultSink negative_extent;
  negative_extent.destination_bytes=sentinel;
  swcme::output::FileOperations negative_extent_ops=
      operations_for(negative_extent);
  status=three.model.write_tecplot_dataset_bundle_checked(
      three.mesh,three.metrics,three.step,invalid,"unused-out04-bundle.dat",
      &negative_extent_ops);
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          status.has_offending_value &&
                          status.offending_value==-1.0 &&
                          negative_extent.open_attempts==0 &&
                          !negative_extent.commit_called,
                      "bundle rejects a negative half extent before open");

  // Although a face does not use Ni to choose its Y/Z resolution, Ni remains
  // part of the BoxSpec invariant.  This case guards against the former API
  // discrepancy where the face accepted a box that the bundle could reject.
  invalid=three.box;
  invalid.Ni=1;
  FaultSink small_ni;
  status=reject_face(invalid,small_ni);
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          status.has_offending_value &&
                          status.offending_value==1.0,
                      "standalone face enforces Ni >= 2");

  invalid=three.box;
  invalid.Nj=0;
  FaultSink small_nj;
  status=reject_face(invalid,small_nj);
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          status.offending_value==0.0,
                      "standalone face enforces Nj >= 2");

  invalid=three.box;
  invalid.Nk=-4;
  FaultSink small_nk;
  status=reject_face(invalid,small_nk);
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          status.offending_value==-4.0,
                      "standalone face enforces Nk >= 2");

  // These two cases exercise different floating-point hazards: overflowing a
  // center +/- extent bound, and overflowing 2*h even though both symmetric
  // bounds remain finite.  Both are malformed specifications, not point-domain
  // failures, and must be caught before OUT05 begins its grid traversal.
  invalid=three.box;
  invalid.cx=std::numeric_limits<double>::max();
  invalid.hx=std::numeric_limits<double>::max();
  FaultSink bound_overflow;
  status=reject_face(invalid,bound_overflow);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput &&
                          status.has_offending_value,
                      "overflowed box bound is rejected before sampling");

  invalid=three.box;
  invalid.cx=0.0;
  invalid.hx=0.75*std::numeric_limits<double>::max();
  FaultSink span_overflow;
  status=reject_face(invalid,span_overflow);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteInput &&
                          status.has_offending_value,
                      "overflowed doubled span is rejected before sampling");

  // INT_MAX cubed exceeds size_t on the supported 64-bit platform.  The test
  // would be infeasible if validation entered the grid loops, so a prompt
  // INVALID_CONFIGURATION result also proves multiplication is overflow-safe.
  invalid=three.box;
  invalid.Ni=std::numeric_limits<int>::max();
  invalid.Nj=std::numeric_limits<int>::max();
  invalid.Nk=std::numeric_limits<int>::max();
  FaultSink count_overflow;
  status=reject_face(invalid,count_overflow);
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          !status.has_offending_value,
                      "unrepresentable grid point count is rejected");

  // Zero half extents remain intentional and useful for repeated-point
  // diagnostics; the structural rule is nonnegative, while OUT05 still checks
  // the resulting location.  Place the collapsed box safely above the floor.
  swcme3d::BoxSpec collapsed=three.box;
  collapsed.cx=2.0*minimum;
  collapsed.cy=0.0;
  collapsed.cz=0.0;
  collapsed.hx=collapsed.hy=collapsed.hz=0.0;
  FaultSink collapsed_success;
  swcme::output::FileOperations collapsed_success_ops=
      operations_for(collapsed_success);
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,collapsed,"unused-out04-collapsed.dat",
      &collapsed_success_ops);
  context.expect_true(status.ok() && collapsed_success.open_attempts==1 &&
                          collapsed_success.commit_called,
                      "nonnegative zero extents remain valid");

  // Legacy wrappers cannot expose ModelStatus, but must delegate to the same
  // validator and report false without creating the requested destination.
  invalid=three.box;
  invalid.hz=-0.5;
  const std::filesystem::path legacy_path="output/OUT04_legacy_invalid.dat";
  std::error_code cleanup_error;
  std::filesystem::remove(legacy_path,cleanup_error);
  context.expect_true(
      !three.model.write_box_face_minX_tecplot_structured(
          three.step,invalid,legacy_path.string().c_str()) &&
          !std::filesystem::exists(legacy_path),
      "legacy face writer rejects malformed BoxSpec without a file");

  // The convenience factory is part of the same contract: it must reject bad
  // source units or resolution rather than returning a box that fails later.
  const auto factory_rejects=[&](double half_au,int resolution) {
    try {
      (void)three.model.default_apex_box(three.step,half_au,resolution);
    } catch (const std::exception&) {
      return true;
    }
    return false;
  };
  context.expect_true(factory_rejects(-0.01,4),
                      "default box factory rejects negative half size");
  context.expect_true(
      factory_rejects(std::numeric_limits<double>::quiet_NaN(),4),
      "default box factory rejects non-finite half size");
  context.expect_true(factory_rejects(0.01,1),
                      "default box factory rejects resolution below two");
  const swcme3d::BoxSpec factory_box=
      three.model.default_apex_box(three.step,0.01,3);
  context.expect_true(factory_box.Ni==3 && factory_box.Nj==3 &&
                          factory_box.Nk==3 && factory_box.hx>0.0,
                      "default box factory still returns a valid specification");

  // A real pre-existing product must remain untouched by structural rejection,
  // complementing FaultSink's callback proof with production namespace checks.
  const std::filesystem::path preserved_path="output/OUT04_preserved.dat";
  write_file(preserved_path,sentinel);
  invalid=three.box;
  invalid.Nk=1;
  status=three.model.write_box_face_minX_tecplot_structured_checked(
      three.step,invalid,preserved_path.string().c_str());
  context.expect_true(status.code==swcme::StatusCode::InvalidConfiguration &&
                          read_file(preserved_path)==sentinel &&
                          count_staging_files(preserved_path)==0,
                      "production BoxSpec rejection preserves destination and staging namespace");
  std::filesystem::remove(preserved_path,cleanup_error);
}

// OUT06: surface output must accept only one self-consistent mesh-plus-metrics
// record.  FaultSink proves every structural, physical, topology, and stale-
// metric rejection occurs before staging; positive cases verify both supplied
// canonical metrics and the documented all-empty auto-compute request.
void test_out06(swcme_test::Context& context) {
  std::cout << "OUT06 mesh output validation\n";
  ThreeDimensionalFixture three;
  const std::string sentinel="PREVIOUS-COMPLETE-OUTPUT\n";

  // Establish that the stricter checks accept the ordinary production record
  // before fault cases mutate its arrays.  This also exercises comparison of
  // complete caller-supplied metrics with their canonical derivation.
  FaultSink supplied_metrics_success;
  swcme::output::FileOperations supplied_metrics_success_ops=
      operations_for(supplied_metrics_success);
  swcme::ModelStatus status=
      three.model.write_shock_surface_center_metrics_tecplot_checked(
          three.mesh,three.metrics,"unused-out06-supplied-metrics.dat",
          &supplied_metrics_success_ops);
  context.expect_true(status.ok() &&
                          supplied_metrics_success.open_attempts==1 &&
                          supplied_metrics_success.commit_called,
                      "canonical supplied mesh and metrics commit successfully");

  const auto reject_surface=[&](const swcme3d::ShockMesh& mesh,
                                const swcme3d::TriMetrics& metrics,
                                FaultSink& sink) {
    sink.destination_bytes=sentinel;
    swcme::output::FileOperations operations=operations_for(sink);
    const swcme::ModelStatus status=
        three.model.write_shock_surface_center_metrics_tecplot_checked(
            mesh,metrics,"unused-out06-surface.dat",&operations);
    context.expect_true(sink.open_attempts==0 && !sink.commit_called &&
                            !sink.remove_called &&
                            sink.destination_bytes==sentinel,
                        "invalid mesh record causes no output callback");
    return status;
  };

  // Parallel nodal arrays and connectivity arrays are part of one record; a
  // missing element must be INVALID_MESH rather than being discovered as an
  // out-of-bounds access during serialization.
  swcme3d::ShockMesh invalid_mesh=three.mesh;
  invalid_mesh.n_hat_z.pop_back();
  FaultSink nodal_size;
  status=reject_surface(invalid_mesh,three.metrics,nodal_size);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh,
                      "mismatched nodal array size is INVALID_MESH");

  invalid_mesh=three.mesh;
  invalid_mesh.tri_k.pop_back();
  FaultSink connectivity_size;
  status=reject_surface(invalid_mesh,three.metrics,connectivity_size);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh,
                      "mismatched connectivity array size is INVALID_MESH");

  swcme3d::ShockMesh empty_mesh;
  swcme3d::TriMetrics empty_metrics;
  FaultSink empty_record;
  status=reject_surface(empty_mesh,empty_metrics,empty_record);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh,
                      "empty surface is rejected before output");

  // Finiteness alone does not make a nodal state valid.  Normals must remain
  // unit length, compression cannot fall below one, and outward shock-normal
  // speed cannot be negative.  Each diagnostic retains the bad node index.
  invalid_mesh=three.mesh;
  invalid_mesh.n_hat_x[2]*=2.0;
  invalid_mesh.n_hat_y[2]*=2.0;
  invalid_mesh.n_hat_z[2]*=2.0;
  FaultSink bad_nodal_normal;
  status=reject_surface(invalid_mesh,three.metrics,bad_nodal_normal);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==2 &&
                          status.has_offending_value,
                      "non-unit nodal normal identifies its node");

  invalid_mesh=three.mesh;
  invalid_mesh.rc[3]=0.75;
  FaultSink bad_compression;
  status=reject_surface(invalid_mesh,three.metrics,bad_compression);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==3 &&
                          status.offending_value==0.75,
                      "sub-unity compression identifies its node");

  invalid_mesh=three.mesh;
  invalid_mesh.Vsh_n[4]=-1.0;
  FaultSink bad_speed;
  status=reject_surface(invalid_mesh,three.metrics,bad_speed);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==4 &&
                          status.offending_value==-1.0,
                      "negative normal speed identifies its node");

  invalid_mesh=three.mesh;
  invalid_mesh.rc[1]=std::numeric_limits<double>::quiet_NaN();
  FaultSink nonfinite_node;
  status=reject_surface(invalid_mesh,three.metrics,nonfinite_node);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteResult &&
                          status.sample_index==1 &&
                          std::isnan(status.offending_value),
                      "non-finite nodal physics is NONFINITE_RESULT");

  // Validate the exact one-based integers that Tecplot will receive.  Use a
  // non-first triangle so sample_index proves the complete connectivity list
  // is scanned rather than only its leading record.
  invalid_mesh=three.mesh;
  invalid_mesh.tri_i[1]=0;
  FaultSink zero_connectivity;
  status=reject_surface(invalid_mesh,three.metrics,zero_connectivity);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==1,
                      "zero one-based connectivity identifies its triangle");

  invalid_mesh=three.mesh;
  invalid_mesh.tri_j[2]=static_cast<int>(invalid_mesh.x.size()+1);
  FaultSink range_connectivity;
  status=reject_surface(invalid_mesh,three.metrics,range_connectivity);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==2,
                      "out-of-range connectivity identifies its triangle");

  invalid_mesh=three.mesh;
  invalid_mesh.tri_k[3]=invalid_mesh.tri_i[3];
  FaultSink repeated_connectivity;
  status=reject_surface(invalid_mesh,three.metrics,repeated_connectivity);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==3,
                      "repeated triangle node is rejected");

  // Unique indices can still form an unusable geometric cell.  Coincident
  // coordinates exercise the production metric builder's scale-relative area
  // and winding checks through the checked writer status boundary.
  invalid_mesh=three.mesh;
  invalid_mesh.x[1]=invalid_mesh.x[0];
  invalid_mesh.y[1]=invalid_mesh.y[0];
  invalid_mesh.z[1]=invalid_mesh.z[0];
  FaultSink degenerate_geometry;
  status=reject_surface(invalid_mesh,empty_metrics,degenerate_geometry);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh,
                      "geometrically degenerate triangle is rejected");

  // TriMetrics is either completely absent (requesting canonical computation)
  // or complete.  Partial records must not silently discard caller data.
  swcme3d::TriMetrics partial_metrics;
  partial_metrics.area=three.metrics.area;
  FaultSink partial_record;
  status=reject_surface(three.mesh,partial_metrics,partial_record);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh,
                      "partially populated metrics are rejected");

  swcme3d::TriMetrics corrupt_metrics=three.metrics;
  corrupt_metrics.area[1]=std::numeric_limits<double>::infinity();
  FaultSink nonfinite_metric;
  status=reject_surface(three.mesh,corrupt_metrics,nonfinite_metric);
  context.expect_true(status.code==swcme::StatusCode::NonFiniteResult &&
                          status.sample_index==1 &&
                          std::isinf(status.offending_value),
                      "non-finite metric identifies its triangle");

  corrupt_metrics=three.metrics;
  corrupt_metrics.area[2]*=1.01;
  FaultSink stale_metric;
  status=reject_surface(three.mesh,corrupt_metrics,stale_metric);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==2 &&
                          status.has_offending_value,
                      "finite stale metric is rejected against canonical geometry");

  // Changing a valid nodal physical value without regenerating its cell means
  // is another stale-record path that array-size and finiteness checks miss.
  invalid_mesh=three.mesh;
  invalid_mesh.rc[0]+=0.1;
  FaultSink stale_nodal_metric;
  status=reject_surface(invalid_mesh,three.metrics,stale_nodal_metric);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.has_offending_value,
                      "metrics stale after nodal-state change are rejected");

  // A completely empty metrics record is the sole auto-compute request.  It
  // must produce and commit the same valid surface rather than being confused
  // with a malformed partial record.
  FaultSink auto_metrics;
  swcme::output::FileOperations auto_metrics_ops=operations_for(auto_metrics);
  status=three.model.write_shock_surface_center_metrics_tecplot_checked(
      three.mesh,empty_metrics,"unused-out06-auto-metrics.dat",
      &auto_metrics_ops);
  context.expect_true(status.ok() && auto_metrics.open_attempts==1 &&
                          auto_metrics.commit_called &&
                          auto_metrics.destination_bytes.find(
                              "TITLE=\"Shock surface")==0,
                      "all-empty metrics request canonical computation and commits");

  // Hold the bundle to the same topology validator.  Its valid BoxSpec and
  // state cannot mask a bad surface record or allow the file to open.
  invalid_mesh=three.mesh;
  invalid_mesh.tri_i[1]=0;
  FaultSink bundle_mesh;
  swcme::output::FileOperations bundle_mesh_ops=operations_for(bundle_mesh);
  status=three.model.write_tecplot_dataset_bundle_checked(
      invalid_mesh,three.metrics,three.step,three.box,
      "unused-out06-bundle.dat",&bundle_mesh_ops);
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          status.sample_index==1 &&
                          bundle_mesh.open_attempts==0,
                      "bundle shares complete mesh validation before open");

  // Legacy behavior maps the checked status to false and must not leave a new
  // path.  A separate production case verifies preservation of existing bytes
  // and absence of a transaction sibling.
  const std::filesystem::path legacy_path="output/OUT06_legacy_invalid.dat";
  std::error_code cleanup_error;
  std::filesystem::remove(legacy_path,cleanup_error);
  context.expect_true(
      !three.model.write_shock_surface_center_metrics_tecplot(
          invalid_mesh,three.metrics,legacy_path.string().c_str()) &&
          !std::filesystem::exists(legacy_path),
      "legacy surface writer rejects invalid mesh without a file");

  const std::filesystem::path preserved_path="output/OUT06_preserved.dat";
  write_file(preserved_path,sentinel);
  status=three.model.write_shock_surface_center_metrics_tecplot_checked(
      invalid_mesh,three.metrics,preserved_path.string().c_str());
  context.expect_true(status.code==swcme::StatusCode::InvalidMesh &&
                          read_file(preserved_path)==sentinel &&
                          count_staging_files(preserved_path)==0,
                      "production mesh rejection preserves destination and staging namespace");
  std::filesystem::remove(preserved_path,cleanup_error);
}
