// Portable regression for the real srcSEP3D output callbacks. The runner
// extracts their unchanged function bodies from main_lib.cpp into an include.
// Only AMPS byte-buffer/channel services are test doubles; the production
// interpolation and presentation arithmetic is compiled from sampling.cpp.
// This verifies callback behavior, not MPI transport or AMPS cut-cell geometry.
#include "../../output/sampling.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace PIC {
int ThisThread = 0;
long int LastSampleLength = 0;
namespace Mesh {
int DatumParticleNumber = 0;
class cDataCenterNode {
 public:
  std::vector<char> bytes = std::vector<char>(512, char(0x5a));
  double position[3] = {2.0, 0.0, 0.0};
  double sample = 0.0;
  char* GetAssociatedDataBufferPointer() { return bytes.data(); }
  void GetX(double* x) { std::copy(position, position + 3, x); }
  double GetDatumAverage(int, int) { return sample; }
};
}  // namespace Mesh
}  // namespace PIC

struct CMPI_channel {
  int ThisThread = 0;
  std::vector<double> sent;
  void send(double* data, int n) { sent.assign(data, data + n); }
  void recv(double* data, int n, int) {
    if (sent.size() != static_cast<std::size_t>(n))
      throw std::runtime_error("owner/root row lengths differ");
    std::copy(sent.begin(), sent.end(), data);
  }
};

namespace {
namespace C = SEP3D::Core;
namespace R = SEP3D::RuntimeModel;
namespace O = SEP3D::Output;
using Node = PIC::Mesh::cDataCenterNode;

struct ConfigurationFixture {
  R::StorageLayout layout;
  R::RunConfiguration3DOptions settings;
  const R::StorageLayout& storage_layout() const { return layout; }
  const R::RunConfiguration3DOptions& options() const { return settings; }
} configuration;
ConfigurationFixture& Configuration() { return configuration; }
int gStaticCellDataOffset = 3;  // deliberately unaligned, as permitted by AMPS
struct Rejected : std::runtime_error {
  C::StatusCode code;
  explicit Rejected(const C::Status& status)
      : std::runtime_error(status.message), code(status.code) {}
};
[[noreturn]] void StopWithStatus(const char*, const C::Status& status) {
  throw Rejected(status);
}
void LoadBytes(Node* node, std::size_t offset, void* target, std::size_t bytes) {
  std::memcpy(target, node->GetAssociatedDataBufferPointer() +
              gStaticCellDataOffset + offset, bytes);
}

// Provider/event physics remains in a native helper in main_lib.cpp.  This
// deterministic double lets the portable probe verify the callback's optional
// row sizing, column order, and owner/root transfer without replacing or
// reimplementing the production mass-density calculation.
bool reducedProductionColumns = false;
bool ReducedProductionColumnsEnabled() { return reducedProductionColumns; }
void AppendReducedProductionColumns(
    Node*, std::vector<double>* values, std::size_t* cursor) {
  if (values == nullptr || cursor == nullptr || *cursor + 4 > values->size())
    throw std::runtime_error(
        "reduced-column test double received an incomplete row");
  for (double value : {101.0, 102.0, 103.0, 104.0})
    (*values)[(*cursor)++] = value;
}

// Generated from main_lib.cpp by _check_output_boundary; do not maintain a
// second implementation of either native callback in this probe.
#include "native_output_callbacks.inc"

void Require(bool good, const char* reason) {
  if (!good) throw std::runtime_error(reason);
}
void Layout(bool gradients) {
  auto& l = configuration.layout;
  l = R::StorageLayout();
  std::size_t next = 0;
  auto reserve = [&](std::size_t& offset, std::size_t count) {
    offset = next;
    next += count * sizeof(double);
  };
  reserve(l.magneticFieldOffset, 3); reserve(l.bulkVelocityOffset, 3);
  reserve(l.numberDensityOffset, 1); reserve(l.velocityDivergenceOffset, 1);
  reserve(l.temperatureOffset, 1); reserve(l.pressureOffset, 1);
  reserve(l.alfvenSpeedOffset, 1); reserve(l.divBhatOffset, 1);
  reserve(l.focusingLengthOffset, 1); reserve(l.curvatureOffset, 3);
  reserve(l.fieldAlignedStrainOffset, 1);
  if (gradients) {
    reserve(l.magneticGradientOffset, 9); reserve(l.velocityGradientOffset, 9);
  }
  reserve(l.waveEnergyOffset, 2);
  l.cellAssociatedBytes = next;
  configuration.settings.innerRadiusM = 1.0;
  configuration.settings.outerRadiusM = 10.0;
  configuration.settings.coordinateOriginM = C::Vec3();
}
void Fill(Node& node, double value) {
  const auto n = Configuration().layout.cellAssociatedBytes / sizeof(double);
  std::vector<double> values(n, value);
  std::memcpy(node.bytes.data() + gStaticCellDataOffset, values.data(),
              n * sizeof(double));
}
std::vector<double> State(Node& node) {
  std::vector<double> values(Configuration().layout.cellAssociatedBytes / sizeof(double));
  LoadBytes(&node, 0, values.data(), values.size() * sizeof(double));
  return values;
}
template <class F> void Reject(F action, C::StatusCode code) {
  try { action(); } catch (const Rejected& error) {
    Require(error.code == code, "malformed input has the wrong diagnostic");
    return;
  }
  throw std::runtime_error("malformed interpolation input was accepted");
}
std::vector<double> Print(Node& node, CMPI_channel* pipe = nullptr, int owner = 0) {
  FILE* file = std::tmpfile();
  Require(file != nullptr, "tmpfile unavailable");
  PrintInitializationCellData(file, 0, pipe, owner, &node);
  std::rewind(file);
  std::vector<double> values;
  double value;
  while (std::fscanf(file, "%lf", &value) == 1) values.push_back(value);
  std::fclose(file);
  return values;
}

void EmptyInterpolation() {
  for (bool gradients : {false, true}) {
    Layout(gradients);
    Node destination;
    InterpolateInitializationCellData(nullptr, nullptr, 0, &destination);
    const auto values = State(destination);
    Require(std::all_of(values.begin(), values.end(),
                        [](double x) { return x == 0.0; }),
            "no-donor interpolation retained stale application bytes");
    Require(destination.bytes[0] == char(0x5a) &&
            destination.bytes[gStaticCellDataOffset - 1] == char(0x5a) &&
            destination.bytes[gStaticCellDataOffset + Configuration().layout.cellAssociatedBytes] == char(0x5a),
            "interpolation overwrote an AMPS-owned slice");
    // Unused arrays must not be read even when they contain invalid values.
    double unused = std::numeric_limits<double>::quiet_NaN();
    InterpolateInitializationCellData(nullptr, &unused, 0, &destination);
    Reject([&] { InterpolateInitializationCellData(nullptr, nullptr, -1, &destination); }, C::StatusCode::InvalidInput);
    Reject([&] { InterpolateInitializationCellData(nullptr, nullptr, 1, &destination); }, C::StatusCode::InvalidInput);
    Node first, second;
    Fill(first, 2.0); Fill(second, 6.0);
    Node* sources[] = {&first, &second};
    double weights[] = {0.25, 0.75};
    InterpolateInitializationCellData(sources, weights, 2, &destination);
    for (double x : State(destination)) Require(x == 5.0, "populated interpolation changed");
    sources[0] = nullptr;
    Reject([&] { InterpolateInitializationCellData(sources, weights, 2, &destination); }, C::StatusCode::InvalidInput);
    const int offset = gStaticCellDataOffset;
    gStaticCellDataOffset = -1;
    Reject([&] { InterpolateInitializationCellData(nullptr, nullptr, 0, &destination); }, C::StatusCode::LayoutMismatch);
    gStaticCellDataOffset = offset;
  }
  double value = 42.0;
  Require(!O::InterpolateStaticCenterState(nullptr, nullptr, 0, 0, &value).ok(), "zero-size output layout accepted");
  Require(!O::InterpolateStaticCenterState(nullptr, nullptr, 0, 1, nullptr).ok(), "null output buffer accepted");
}

void ExcludedRows() {
  for (bool gradients : {false, true}) {
    Layout(gradients);
    const std::size_t rowSize = gradients ? 44 : 26;
    Node node;
    InterpolateInitializationCellData(nullptr, nullptr, 0, &node);
    PIC::LastSampleLength = 8;
    auto row = Print(node);
    Require(row.size() == rowSize, "empty row violated the variable-count contract");
    for (double x : row) Require(std::isfinite(x), "empty output row is non-finite");
    Require(row[rowSize - 3] == 0.0 && row[rowSize - 2] == 1.0 && row.back() == 0.0,
            "in-shell excluded row or sampling availability is wrong");
    Fill(node, 2.0);
    row = Print(node);
    Require(row[rowSize - 3] == 1.0 && row.back() == 0.0,
            "a populated background with no particles was invalidated");
    node.sample = 3.0;
    row = Print(node);
    Require(row.back() == 1.0, "occupied particle window was lost");
    node.position[0] = 0.5;  // solar interior
    row = Print(node);
    Require(row[rowSize - 3] == 0.0, "solar-interior background became valid");
    for (std::size_t i = 0; i < rowSize - 3; ++i)
      Require(row[i] == 0.0, "excluded background has nonzero placeholders");

    // Exercise the actual callback's owner-send/root-receive branches with
    // a channel double. The nonroot FILE* is normally null in AMPS output.
    CMPI_channel channel;
    channel.ThisThread = PIC::ThisThread = 6;
    PrintInitializationCellData(nullptr, 0, &channel, 6, &node);
    Require(channel.sent == row, "nonroot callback changed the row");
    channel.ThisThread = PIC::ThisThread = 0;
    Node rootTemporary;
    Require(Print(rootTemporary, &channel, 6) == row, "root did not receive the owner row");

    // The reduced production mode appends mass density, time, ambient
    // generation, and front generation before the three validity flags.  Its
    // optional columns must not shift or drop the generic flags, and the same
    // row must still traverse the owner/root channel unchanged.
    node.position[0] = 2.0;
    node.sample = 0.0;
    Fill(node, 2.0);
    reducedProductionColumns = true;
    const auto reducedRow = Print(node);
    Require(reducedRow.size() == rowSize + 4,
            "reduced metadata violated the variable-count contract");
    const std::size_t reducedStart = rowSize - 3;
    for (std::size_t index = 0; index < 4; ++index)
      Require(reducedRow[reducedStart + index] == 101.0 + index,
              "reduced metadata column order changed");
    Require(reducedRow[reducedStart + 4] == 1.0 &&
            reducedRow[reducedStart + 5] == 1.0 &&
            reducedRow[reducedStart + 6] == 0.0,
            "reduced metadata shifted the availability flags");
    reducedProductionColumns = false;
  }
}
}  // namespace

int main(int argc, char** argv) {
  try {
    const std::string id = argc > 1 ? argv[1] : "";
    if (id == "OUT3D01") EmptyInterpolation();
    else if (id == "OUT3D02") ExcludedRows();
    else throw std::runtime_error("expected OUT3D01 or OUT3D02");
    std::cout << id << " PASS: actual output callbacks accept excluded vertices; portable buffer/channel probe\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "output_boundary_probe FAIL: " << error.what() << '\n';
    return 1;
  }
}
