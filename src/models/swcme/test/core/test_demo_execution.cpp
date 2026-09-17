#include "test_framework.hpp"
#include "tecplot_parser.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <locale>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

namespace {

struct CsvTable {
  std::vector<std::string> header;
  std::vector<std::vector<double>> rows;
};

struct DemoRun {
  bool exited_normally=false;
  int exit_code=-1;
  std::string stdout_text;
  std::string stderr_text;
};

std::string read_text(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(input),
                     std::istreambuf_iterator<char>());
}

std::vector<std::string> split_commas(const std::string& line) {
  // Demo CSV fields contain no quoting or embedded commas, so a deliberately
  // small splitter is preferable to a permissive general CSV library: empty
  // fields and wrong widths remain visible to OUT07.
  std::vector<std::string> fields;
  std::size_t begin=0;
  for (;;) {
    const std::size_t comma=line.find(',',begin);
    fields.push_back(line.substr(
        begin,comma==std::string::npos ? std::string::npos : comma-begin));
    if (comma==std::string::npos) break;
    begin=comma+1;
  }
  return fields;
}

bool parse_csv(const std::filesystem::path& path,
               const std::vector<std::string>& expected_header,
               CsvTable& table,std::string& error) {
  const std::string text=read_text(path);
  if (text.empty()) {
    error="empty or unreadable file";
    return false;
  }
  if (text.back()!='\n') {
    error="missing final newline";
    return false;
  }

  std::istringstream input(text);
  input.imbue(std::locale::classic());
  std::string line;
  if (!std::getline(input,line)) {
    error="missing header";
    return false;
  }
  if (!line.empty() && line.back()=='\r') line.pop_back();
  table.header=split_commas(line);
  if (table.header!=expected_header) {
    error="header names, units, or order differ from the demo contract";
    return false;
  }

  std::size_t row_index=0;
  while (std::getline(input,line)) {
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line.empty()) {
      error="empty record at data row "+std::to_string(row_index);
      return false;
    }
    const std::vector<std::string> fields=split_commas(line);
    if (fields.size()!=expected_header.size()) {
      error="wrong field count at data row "+std::to_string(row_index);
      return false;
    }
    std::vector<double> values;
    values.reserve(fields.size());
    for (const std::string& field : fields) {
      std::istringstream number(field);
      number.imbue(std::locale::classic());
      double value=0.0;
      char trailing='\0';
      if (!(number>>value) || (number>>trailing) || !std::isfinite(value)) {
        error="invalid numeric field at data row "+std::to_string(row_index);
        return false;
      }
      values.push_back(value);
    }
    table.rows.push_back(values);
    ++row_index;
  }
  if (table.rows.empty()) {
    error="no data rows";
    return false;
  }
  return true;
}

bool has_time_axis(const CsvTable& table,std::size_t expected_rows,
                   double expected_step,double expected_end) {
  if (table.rows.size()!=expected_rows || table.rows.front().empty())
    return false;
  for (std::size_t row=0; row<table.rows.size(); ++row) {
    const double expected=expected_step*static_cast<double>(row);
    if (std::abs(table.rows[row][0]-expected)>1.0e-8) return false;
  }
  return std::abs(table.rows.back()[0]-expected_end)<=1.0e-8;
}

std::string shell_quote(const std::string& value) {
  // Executable/work paths are generated locally, but robust POSIX quoting
  // prevents whitespace or apostrophes in a checkout path from changing the
  // command.  No user-controlled text is passed to the shell.
  std::string quoted="'";
  for (char character : value) {
    if (character=='\'') quoted+="'\"'\"'";
    else quoted.push_back(character);
  }
  quoted.push_back('\'');
  return quoted;
}

DemoRun run_demo(const std::filesystem::path& executable,
                 const std::filesystem::path& directory) {
  // Each program receives a new empty working directory, and stdout/stderr are
  // captured beside its products.  LC_ALL=C fixes decimal formatting for the
  // independent text parsers without altering model parameters.
  const std::string command="cd "+shell_quote(directory.string())+
      " && LC_ALL=C "+shell_quote(executable.string())+
      " > stdout.log 2> stderr.log";
  const int status=std::system(command.c_str());
  DemoRun run;
  if (status!=-1 && WIFEXITED(status)) {
    run.exited_normally=true;
    run.exit_code=WEXITSTATUS(status);
  }
  run.stdout_text=read_text(directory/"stdout.log");
  run.stderr_text=read_text(directory/"stderr.log");
  return run;
}

bool exact_manifest(const std::filesystem::path& directory,
                    const std::set<std::string>& expected,
                    std::string& diagnostic) {
  std::set<std::string> actual;
  std::error_code error;
  for (std::filesystem::directory_iterator entries(directory,error),end;
       !error && entries!=end; entries.increment(error)) {
    if (!entries->is_regular_file()) {
      diagnostic="non-regular artifact "+entries->path().filename().string();
      return false;
    }
    actual.insert(entries->path().filename().string());
  }
  if (error) {
    diagnostic="directory scan failed: "+error.message();
    return false;
  }
  if (actual!=expected) {
    diagnostic="actual artifact manifest differs from declaration";
    return false;
  }
  for (const std::string& name : actual)
    if (name.find(".swcme-tmp-")!=std::string::npos) {
      diagnostic="transaction staging artifact leaked";
      return false;
    }
  return true;
}

const std::vector<std::string> kProfileVariables={
    "r[m]","R[AU]","rSun[R_s]","n[m^-3]","V[m/s]","Br[T]",
    "Bphi[T]","Bmag[T]","divV[s^-1]","rc","R_sh[m]","R_LE[m]",
    "R_TE[m]"};

const std::vector<std::string> kPointVariables={
    "X[m]","Y[m]","Z[m]","n[m^-3]","Vx[m/s]","Vy[m/s]","Vz[m/s]",
    "Bx[T]","By[T]","Bz[T]","rc[-]","Vsh_n[m/s]"};

const std::vector<std::string> kTimeHeader={
    "t_s","n_m3","Vx_ms","Vy_ms","Vz_ms","V_mag_ms"};
const std::vector<std::string> kShockHeader={
    "t_s","R_sh_AU","V_sh_km_s","rc","R_LE_AU","R_TE_AU"};
const std::vector<std::string> kStrengthHeader={
    "t_s","has_shock","rc_apex","M_fast","U1n_ms","Vsh_n_ms",
    "theta_Bn_rad","B2overB1","field_rotation_rad","mass_residual",
    "normal_B_residual","electric_residual","momentum_residual",
    "energy_residual"};

bool parse_bundle_file(const std::filesystem::path& path,
                       swcme_test::tecplot::ParsedDocument& document,
                       std::string& error) {
  return swcme_test::tecplot::parse_bundle_text(
      read_text(path),document,error);
}

bool valid_demo_bundle(const swcme_test::tecplot::ParsedDocument& document) {
  // Both 3-D examples request nTheta=24, nPhi=48 and a 12^3 box.  A finite SSE
  // cap therefore has one apex plus 24 unique rings and 48 + 23*96 triangles.
  // Computing these counts here, rather than reading them from production mesh
  // metadata, verifies that the example and its documentation remain aligned.
  constexpr std::size_t nodes=1+24*48;
  constexpr std::size_t elements=48+23*2*48;
  return document.zones.size()==4 && document.zones[0].nodes==nodes &&
      document.zones[0].elements==elements &&
      document.zones[1].nodes==nodes &&
      document.zones[1].elements==elements &&
      document.zones[2].ni==12 && document.zones[2].nj==12 &&
      document.zones[2].nk==12 && document.zones[2].rows.size()==12*12*12 &&
      document.zones[3].ni==12 && document.zones[3].nj==12 &&
      document.zones[3].rows.size()==12*12;
}

}  // namespace

// OUT07 is the executable-document-documentation gate.  It builds on OUT01's
// independent parser but drives the actual user-facing binaries in isolated
// directories, so stale comments, ignored return values, invalid default
// geometry, missing products, and unparseable output all become test failures.
void test_out07(swcme_test::Context& context) {
  std::cout << "OUT07 demonstration program execution\n";
  const int failures_before=context.failures();
  const std::filesystem::path test_directory=std::filesystem::current_path();
  // Include both the process identifier and a high-resolution start stamp in
  // the run root.  Some container runtimes reuse namespace-local process IDs
  // across separate command sessions, while the stamp remains distinct; the
  // pair therefore lets concurrent campaigns run without deleting or
  // overwriting one another's products.  A retained directory from an earlier
  // failed campaign likewise cannot contaminate the exact-manifest checks.
  const auto start_stamp=std::chrono::high_resolution_clock::now()
                             .time_since_epoch().count();
  const std::filesystem::path run_root=
      test_directory/"output"/
      ("OUT07_demo_runs_"+std::to_string(static_cast<long long>(::getpid()))+
       "_"+std::to_string(start_stamp));

  std::error_code filesystem_error;
  std::filesystem::remove_all(run_root,filesystem_error);
  context.expect_true(!filesystem_error,
                      "OUT07 removes any prior isolated run directory");
  filesystem_error.clear();
  std::filesystem::create_directories(run_root/"demo1d",filesystem_error);
  std::filesystem::create_directories(run_root/"demo3d_1",filesystem_error);
  std::filesystem::create_directories(run_root/"demo3d_2",filesystem_error);
  context.expect_true(!filesystem_error,
                      "OUT07 creates three empty demonstration directories");

  const DemoRun one=run_demo(
      std::filesystem::absolute(test_directory/"output"/"demo1d"),
      run_root/"demo1d");
  const DemoRun three_basic=run_demo(
      std::filesystem::absolute(test_directory/"output"/"demo3d_1"),
      run_root/"demo3d_1");
  const DemoRun three_extended=run_demo(
      std::filesystem::absolute(test_directory/"output"/"demo3d_2"),
      run_root/"demo3d_2");
  context.expect_true(one.exited_normally && one.exit_code==0 &&
                          one.stderr_text.empty() &&
                          one.stdout_text.find("profile_1d.dat")!=std::string::npos,
                      "demo1d exits zero, reports its product, and has empty stderr");
  context.expect_true(three_basic.exited_normally &&
                          three_basic.exit_code==0 &&
                          three_basic.stderr_text.empty() &&
                          three_basic.stdout_text.find(
                              "sse_apex_bundle_tecplot.dat")!=std::string::npos,
                      "demo3d_1 exits zero and reports its supported-domain bundle");
  context.expect_true(three_extended.exited_normally &&
                          three_extended.exit_code==0 &&
                          three_extended.stderr_text.empty() &&
                          three_extended.stdout_text.find(
                              "max_RH_residual")!=std::string::npos,
                      "demo3d_2 exits zero and reports ideal-MHD shock diagnostics");

  std::string manifest_error;
  context.expect_true(exact_manifest(
      run_root/"demo1d",{"profile_1d.dat","stdout.log","stderr.log"},
      manifest_error),"demo1d produces exactly its declared files: "+manifest_error);
  manifest_error.clear();
  context.expect_true(exact_manifest(
      run_root/"demo3d_1",{"ts_cone.csv","shock_cone.csv",
                           "sse_apex_bundle_tecplot.dat","stdout.log","stderr.log"},
      manifest_error),"demo3d_1 produces exactly its declared files: "+manifest_error);
  manifest_error.clear();
  context.expect_true(exact_manifest(
      run_root/"demo3d_2",{"strength_summary.csv","ts_cone.csv",
                           "shock_cone.csv","sse_apex_bundle_tecplot.dat",
                           "predefined_points_tecplot.dat",
                           "surface_random_samples_tecplot.dat",
                           "stdout.log","stderr.log"},manifest_error),
      "demo3d_2 produces exactly its declared files: "+manifest_error);

  using swcme_test::tecplot::ParsedDocument;
  ParsedDocument profile;
  std::string parse_error;
  const bool profile_valid=swcme_test::tecplot::parse_point_text(
      read_text(run_root/"demo1d"/"profile_1d.dat"),
      "1D SW+CME radial profile",kProfileVariables,"radial","I",profile,
      parse_error);
  context.expect_true(profile_valid && profile.zones.size()==1 &&
                          profile.zones[0].ni==1200 &&
                          profile.zones[0].rows.size()==1200,
                      "demo1d profile parses with 1200 finite 13-column rows: "+
                          parse_error);

  ParsedDocument basic_bundle;
  parse_error.clear();
  const bool basic_bundle_valid=parse_bundle_file(
      run_root/"demo3d_1"/"sse_apex_bundle_tecplot.dat",basic_bundle,
      parse_error);
  context.expect_true(basic_bundle_valid && valid_demo_bundle(basic_bundle),
                      "demo3d_1 bundle has documented zones and dimensions: "+
                          parse_error);

  ParsedDocument extended_bundle;
  parse_error.clear();
  const bool extended_bundle_valid=parse_bundle_file(
      run_root/"demo3d_2"/"sse_apex_bundle_tecplot.dat",extended_bundle,
      parse_error);
  context.expect_true(extended_bundle_valid && valid_demo_bundle(extended_bundle),
                      "demo3d_2 bundle has documented zones and dimensions: "+
                          parse_error);

  ParsedDocument predefined;
  parse_error.clear();
  const bool predefined_valid=swcme_test::tecplot::parse_point_text(
      read_text(run_root/"demo3d_2"/"predefined_points_tecplot.dat"),
      "Point cloud (plasma + B + shock diagnostics)",kPointVariables,"points",
      "N",predefined,parse_error);
  context.expect_true(predefined_valid && predefined.zones[0].nodes==6 &&
                          predefined.zones[0].rows.size()==6,
                      "predefined point cloud parses as six finite rows: "+
                          parse_error);

  ParsedDocument random_samples;
  parse_error.clear();
  const bool samples_valid=swcme_test::tecplot::parse_point_text(
      read_text(run_root/"demo3d_2"/"surface_random_samples_tecplot.dat"),
      "Point cloud (plasma + B + shock diagnostics)",kPointVariables,"points",
      "N",random_samples,parse_error);
  constexpr std::size_t expected_triangles=48+23*2*48;
  context.expect_true(samples_valid &&
                          random_samples.zones[0].nodes==expected_triangles*10 &&
                          random_samples.zones[0].rows.size()==expected_triangles*10,
                      "surface point cloud contains ten finite samples per cell: "+
                          parse_error);

  CsvTable basic_time,basic_shock,extended_time,extended_shock,strength;
  std::string csv_error;
  const bool basic_time_valid=parse_csv(
      run_root/"demo3d_1"/"ts_cone.csv",kTimeHeader,basic_time,csv_error);
  context.expect_true(basic_time_valid &&
                          has_time_axis(basic_time,865,300.0,72.0*3600.0),
                      "demo3d_1 plasma CSV has a complete five-minute axis: "+
                          csv_error);
  csv_error.clear();
  const bool basic_shock_valid=parse_csv(
      run_root/"demo3d_1"/"shock_cone.csv",kShockHeader,basic_shock,csv_error);
  context.expect_true(basic_shock_valid &&
                          has_time_axis(basic_shock,865,300.0,72.0*3600.0),
                      "demo3d_1 shock CSV has a complete five-minute axis: "+
                          csv_error);
  csv_error.clear();
  const bool extended_time_valid=parse_csv(
      run_root/"demo3d_2"/"ts_cone.csv",kTimeHeader,extended_time,csv_error);
  context.expect_true(extended_time_valid &&
                          has_time_axis(extended_time,865,300.0,72.0*3600.0),
                      "demo3d_2 plasma CSV has a complete five-minute axis: "+
                          csv_error);
  csv_error.clear();
  const bool extended_shock_valid=parse_csv(
      run_root/"demo3d_2"/"shock_cone.csv",kShockHeader,extended_shock,csv_error);
  context.expect_true(extended_shock_valid &&
                          has_time_axis(extended_shock,865,300.0,72.0*3600.0),
                      "demo3d_2 shock CSV has a complete five-minute axis: "+
                          csv_error);
  csv_error.clear();
  const bool strength_valid=parse_csv(
      run_root/"demo3d_2"/"strength_summary.csv",kStrengthHeader,strength,
      csv_error);
  context.expect_true(strength_valid &&
                          has_time_axis(strength,865,300.0,72.0*3600.0),
                      "ideal-MHD strength CSV has finite states and residuals: "+
                          csv_error);

  // Passing runs leave no bulky demo products in the source tree.  A failed
  // run is retained under test/output/OUT07_demo_runs_<run-id> so its captured
  // logs and partial manifest remain available for diagnosis without blocking
  // another concurrent campaign.
  if (context.failures()==failures_before) {
    filesystem_error.clear();
    std::filesystem::remove_all(run_root,filesystem_error);
    context.expect_true(!filesystem_error,
                        "OUT07 removes successful temporary run artifacts");
  }
}
