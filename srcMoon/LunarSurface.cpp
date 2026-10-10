#include "LunarSurface.h"

#include "pic.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

namespace Moon {
namespace Surface {
namespace {

namespace fs = std::filesystem;

constexpr double kPi = 3.141592653589793238462643383279502884;
// Bound memory and generation time for accidental metre/sub-metre requests.
// The limit is checked before a refinement can exceed it.
constexpr std::size_t kMaximumFaces = 6u * 1024u * 1024u;
// Initialization-only state consumed by main_lib.cpp after the collective
// generation/load step.  The surface remains static for the lifetime of AMR.
bool gRealisticSurfaceActive = false;
double gMinimumSurfaceRadiusM = 0.0;

// PDS scalar strings may be quoted and numeric fields may carry <UNIT> text.
// This helper removes only surrounding whitespace/quotes; it is not a general
// PDS parser and is intentionally paired with strict required-key checks.
std::string Trim(const std::string& text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n\"'");
  if (first == std::string::npos) return std::string();
  const std::size_t last = text.find_last_not_of(" \t\r\n\"'");
  return text.substr(first, last - first + 1);
}

bool LabelValue(const std::string& label, const std::string& key,
                std::string* value) {
  // LDEM_4 uses one KEY = VALUE assignment per line.  Strip the supported PDS
  // block-comment suffix and match complete keys to avoid confusing similarly
  // named metadata.  Unsupported label layouts fail as missing metadata.
  std::istringstream lines(label);
  std::string line;
  while (std::getline(lines, line)) {
    const std::size_t comment = line.find("/*");
    if (comment != std::string::npos) line.erase(comment);
    const std::size_t equal = line.find('=');
    if (equal == std::string::npos) continue;
    if (Trim(line.substr(0, equal)) == key) {
      *value = Trim(line.substr(equal + 1));
      const std::size_t unit = value->find('<');
      if (unit != std::string::npos) *value = Trim(value->substr(0, unit));
      return true;
    }
  }
  return false;
}

bool ParseUnsigned(const std::string& text, std::size_t* value) {
  char* end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || end == text.c_str() || *end != '\0') return false;
  *value = static_cast<std::size_t>(parsed);
  return true;
}

// Numeric label conversion is locale-independent under the process C locale
// used by AMPS and rejects non-finite values and trailing characters.
bool ParseDouble(const std::string& text, double* value) {
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

bool ReadString(const std::string& label, const char* key, std::string* value,
                std::string* error) {
  if (!LabelValue(label, key, value)) {
    *error = std::string("LOLA label is missing ") + key;
    return false;
  }
  return true;
}

bool ReadSize(const std::string& label, const char* key, std::size_t* value,
              std::string* error) {
  std::string text;
  if (!LabelValue(label, key, &text) || !ParseUnsigned(text, value)) {
    *error = std::string("LOLA label has missing or invalid ") + key;
    return false;
  }
  return true;
}

bool ReadNumber(const std::string& label, const char* key, double* value,
                std::string* error) {
  std::string text;
  if (!LabelValue(label, key, &text) || !ParseDouble(text, value)) {
    *error = std::string("LOLA label has missing or invalid ") + key;
    return false;
  }
  return true;
}

double Norm(const std::array<double, 3>& x) {
  return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

std::array<double, 3> Unit(std::array<double, 3> x) {
  // Callers supply non-antipodal icosphere vertices/midpoints; a zero vector
  // would indicate corrupt topology and is impossible for the fixed seed.
  const double norm = Norm(x);
  for (double& component : x) component /= norm;
  return x;
}

double GreatCircleEdge(const std::array<double, 3>& a,
                       const std::array<double, 3>& b, double radius) {
  // Clamp the normalized dot product to absorb round-off before acos().
  double dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  dot = std::max(-1.0, std::min(1.0, dot / (Norm(a) * Norm(b))));
  return radius * std::acos(dot);
}

std::size_t Midpoint(std::size_t a, std::size_t b,
                     std::vector<std::array<double, 3>>* vertices,
                     std::map<std::pair<std::size_t, std::size_t>,
                              std::size_t>* cache) {
  // Canonicalize the undirected edge and cache its midpoint.  Without the
  // cache, adjacent triangles would create coincident but disconnected nodes,
  // breaking watertightness and ray-intersection parity.
  if (a > b) std::swap(a, b);
  const std::pair<std::size_t, std::size_t> key(a, b);
  const auto found = cache->find(key);
  if (found != cache->end()) return found->second;
  const auto& x = (*vertices)[a];
  const auto& y = (*vertices)[b];
  const std::size_t index = vertices->size();
  vertices->push_back(Unit({x[0] + y[0], x[1] + y[1], x[2] + y[2]}));
  cache->emplace(key, index);
  return index;
}

void OrientOutward(Triangulation* mesh) {
  // For a body enclosing the origin, an outward face satisfies
  // cross(b-a,c-a) dot (a+b+c) > 0.  Radial DEM displacement preserves this
  // star-shaped convention, but rechecking after displacement is inexpensive.
  for (auto& face : mesh->faces) {
    const auto& a = mesh->verticesM[face[0]];
    const auto& b = mesh->verticesM[face[1]];
    const auto& c = mesh->verticesM[face[2]];
    const double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const double normal[3] = {
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0]};
    const double center[3] = {
        a[0] + b[0] + c[0], a[1] + b[1] + c[1],
        a[2] + b[2] + c[2]};
    if (normal[0] * center[0] + normal[1] * center[1] +
        normal[2] * center[2] < 0.0) {
      std::swap(face[1], face[2]);
    }
  }
}

bool EnsureParent(const std::string& path, std::string* error) {
  std::error_code ec;
  const fs::path parent = fs::path(path).parent_path();
  if (!parent.empty()) fs::create_directories(parent, ec);
  if (ec) {
    *error = "cannot create surface output directory '" + parent.string() +
        "': " + ec.message();
    return false;
  }
  return true;
}

template <class Writer>
bool AtomicWrite(const std::string& path, Writer writer, std::string* error) {
  // A reader must never observe a half-written mesh.  Rank zero writes a
  // sibling temporary file, closes/checks it, then renames it into place
  // before the MPI barrier in PrepareProductionSurface().
  if (!EnsureParent(path, error)) return false;
  const fs::path target(path);
  const fs::path temporary = target.string() + ".tmp";
  std::ofstream output(temporary, std::ios::out | std::ios::trunc);
  if (!output) {
    *error = "cannot open surface output '" + temporary.string() + "'";
    return false;
  }
  writer(output);
  output.close();
  if (!output) {
    *error = "I/O failure while writing surface output '" +
        temporary.string() + "'";
    return false;
  }
  std::error_code ec;
  fs::rename(temporary, target, ec);
  if (ec) {
    fs::remove(target, ec);
    ec.clear();
    fs::rename(temporary, target, ec);
  }
  if (ec) {
    *error = "cannot install surface output '" + target.string() + "': " +
        ec.message();
    return false;
  }
  return true;
}

bool InitializeLoadedNormalsFromWinding(std::string* error) {
  // The generated mesh has a known outward winding.  AMPS's generic
  // InitExternalNormalVector() determines orientation by ray tracing and
  // requires PIC::Mesh::mesh->xGlobalMin/xGlobalMax, which are not available
  // until after the surface has already been loaded and the AMR mesh is
  // initialized.  Preserve the deterministic generator orientation here and
  // populate the node-averaged normals needed by surface utilities.
  using PIC::Mesh::IrregularSurface::BoundaryTriangleFaces;
  using PIC::Mesh::IrregularSurface::BoundaryTriangleNodes;
  using PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;
  using PIC::Mesh::IrregularSurface::nBoundaryTriangleNodes;

  for (int node = 0; node < nBoundaryTriangleNodes; ++node) {
    for (int dimension = 0; dimension < 3; ++dimension)
      BoundaryTriangleNodes[node].BallAveragedExternalNormal[dimension] = 0.0;
  }

  for (int face = 0; face < nBoundaryTriangleFaces; ++face) {
    double center[3];
    BoundaryTriangleFaces[face].GetCenterPosition(center);
    const double* normal = BoundaryTriangleFaces[face].ExternalNormal;
    const double outward = center[0] * normal[0] + center[1] * normal[1] +
        center[2] * normal[2];
    if (!(outward > 0.0) || !std::isfinite(outward)) {
      *error = "loaded LOLA surface contains a non-outward face normal";
      return false;
    }
    for (int corner = 0; corner < 3; ++corner) {
      for (int dimension = 0; dimension < 3; ++dimension) {
        BoundaryTriangleFaces[face].node[corner]
            ->BallAveragedExternalNormal[dimension] += normal[dimension];
      }
    }
  }

  for (int node = 0; node < nBoundaryTriangleNodes; ++node) {
    double lengthSquared = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double component =
          BoundaryTriangleNodes[node].BallAveragedExternalNormal[dimension];
      lengthSquared += component * component;
    }
    if (!(lengthSquared > 0.0) || !std::isfinite(lengthSquared)) {
      *error = "loaded LOLA surface contains an undefined vertex normal";
      return false;
    }
    const double inverseLength = 1.0 / std::sqrt(lengthSquared);
    for (int dimension = 0; dimension < 3; ++dimension) {
      BoundaryTriangleNodes[node].BallAveragedExternalNormal[dimension] *=
          inverseLength;
    }
  }
  return true;
}

}  // namespace

bool LolaDem::Load(const Runtime::Configuration& configuration,
                   std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  std::ifstream labelStream(configuration.lolaLabelFile);
  if (!labelStream) {
    *error = "cannot open LOLA label '" + configuration.lolaLabelFile + "'";
    return false;
  }
  const std::string label((std::istreambuf_iterator<char>(labelStream)),
                          std::istreambuf_iterator<char>());

  // Decode into temporaries and commit at the end so a failed reload cannot
  // leave metadata and pixels from different products in the same object.
  LolaMetadata candidate;
  std::size_t sampleBits = 0;
  if (!ReadString(label, "PRODUCT_ID", &candidate.productId, error) ||
      !ReadString(label, "DATA_SET_ID", &candidate.dataSetId, error) ||
      !ReadString(label, "MAP_PROJECTION_TYPE", &candidate.projection, error) ||
      !ReadString(label, "COORDINATE_SYSTEM_TYPE",
                  &candidate.coordinateSystemType, error) ||
      !ReadString(label, "COORDINATE_SYSTEM_NAME",
                  &candidate.coordinateSystemName, error) ||
      !ReadString(label, "POSITIVE_LONGITUDE_DIRECTION",
                  &candidate.positiveLongitudeDirection, error) ||
      !ReadString(label, "SAMPLE_TYPE", &candidate.sampleType, error) ||
      !ReadString(label, "UNIT", &candidate.unit, error) ||
      !ReadSize(label, "LINES", &candidate.lines, error) ||
      !ReadSize(label, "LINE_SAMPLES", &candidate.lineSamples, error) ||
      !ReadSize(label, "SAMPLE_BITS", &sampleBits, error) ||
      !ReadNumber(label, "SCALING_FACTOR", &candidate.scalingFactor, error) ||
      !ReadNumber(label, "OFFSET", &candidate.referenceRadiusM, error) ||
      !ReadNumber(label, "MAP_RESOLUTION",
                  &candidate.mapResolutionPixelsPerDegree, error) ||
      !ReadNumber(label, "MINIMUM_LATITUDE",
                  &candidate.minimumLatitudeDeg, error) ||
      !ReadNumber(label, "MAXIMUM_LATITUDE",
                  &candidate.maximumLatitudeDeg, error) ||
      !ReadNumber(label, "WESTERNMOST_LONGITUDE",
                  &candidate.westernmostLongitudeDeg, error) ||
      !ReadNumber(label, "EASTERNMOST_LONGITUDE",
                  &candidate.easternmostLongitudeDeg, error)) {
    return false;
  }
  candidate.sampleBits = static_cast<unsigned int>(sampleBits);

  // The native label is authoritative.  These checks prevent a different PDS
  // raster from being silently interpreted with LDEM_4 byte order, units,
  // coordinate direction, or projection assumptions.
  if (candidate.productId != configuration.lolaProductId) {
    *error = "LOLA PRODUCT_ID '" + candidate.productId +
        "' does not match configured lola_product_id '" +
        configuration.lolaProductId + "'";
    return false;
  }
  if (candidate.projection != "SIMPLE CYLINDRICAL" ||
      candidate.coordinateSystemType != "BODY-FIXED ROTATING" ||
      candidate.positiveLongitudeDirection != "EAST" ||
      candidate.sampleType != "LSB_INTEGER" || candidate.sampleBits != 16 ||
      candidate.unit != "METER" || candidate.lines < 2 ||
      candidate.lineSamples < 2 || candidate.scalingFactor <= 0.0 ||
      candidate.referenceRadiusM <= 0.0 ||
      candidate.mapResolutionPixelsPerDegree <= 0.0) {
    *error = "unsupported LOLA label: require SIMPLE CYLINDRICAL, "
        "BODY-FIXED ROTATING, EAST longitude, 16-bit LSB_INTEGER metres, "
        "and positive scale/resolution/radius";
    return false;
  }

  // LSB_INTEGER/SAMPLE_BITS=16 has exactly two bytes per pixel and no attached
  // label bytes in this detached IMG.  Exact size rejects truncation, an
  // accidentally selected higher-resolution product, and record padding.
  const std::uintmax_t expectedBytes =
      static_cast<std::uintmax_t>(candidate.lines) *
      static_cast<std::uintmax_t>(candidate.lineSamples) * 2u;
  std::error_code ec;
  const std::uintmax_t actualBytes =
      fs::file_size(configuration.lolaImageFile, ec);
  if (ec || actualBytes != expectedBytes) {
    std::ostringstream message;
    message << "LOLA IMG byte count mismatch: expected " << expectedBytes
            << ", found " << (ec ? 0 : actualBytes);
    *error = message.str();
    return false;
  }

  std::ifstream image(configuration.lolaImageFile, std::ios::binary);
  if (!image) {
    *error = "cannot open LOLA image '" + configuration.lolaImageFile + "'";
    return false;
  }
  std::vector<unsigned char> bytes(static_cast<std::size_t>(actualBytes));
  image.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
  if (!image) {
    *error = "I/O failure while reading LOLA image";
    return false;
  }
  // Assemble little-endian values explicitly so decoding is independent of
  // host endianness and implementation-defined pointer reinterpretation.
  std::vector<std::int16_t> decoded(expectedBytes / 2u);
  for (std::size_t index = 0; index < decoded.size(); ++index) {
    const std::uint16_t raw = static_cast<std::uint16_t>(bytes[2 * index]) |
        (static_cast<std::uint16_t>(bytes[2 * index + 1]) << 8u);
    decoded[index] = raw < 0x8000u ? static_cast<std::int16_t>(raw) :
        static_cast<std::int16_t>(static_cast<std::int32_t>(raw) - 0x10000);
  }
  metadata_ = std::move(candidate);
  elevationDn_ = std::move(decoded);
  return true;
}

double LolaDem::ElevationM(double longitudeDeg, double latitudeDeg) const {
  // The LDEM_4 longitude interval spans one full periodic revolution.  Reduce
  // any caller longitude to that native interval before computing pixel
  // coordinates, so interpolation across 0/360 uses adjacent raster columns.
  const double longitudeWidth = metadata_.easternmostLongitudeDeg -
      metadata_.westernmostLongitudeDeg;
  longitudeDeg = std::fmod(longitudeDeg - metadata_.westernmostLongitudeDeg,
                           longitudeWidth);
  if (longitudeDeg < 0.0) longitudeDeg += longitudeWidth;
  longitudeDeg += metadata_.westernmostLongitudeDeg;
  latitudeDeg = std::max(metadata_.minimumLatitudeDeg,
      std::min(metadata_.maximumLatitudeDeg, latitudeDeg));

  // The label bounds describe pixel edges.  Subtracting 0.5 maps the centre of
  // the first pixel to integer index zero.  Rows run north-to-south in the IMG.
  const double longitudeStep = longitudeWidth / metadata_.lineSamples;
  const double latitudeStep = (metadata_.maximumLatitudeDeg -
      metadata_.minimumLatitudeDeg) / metadata_.lines;
  const double column = (longitudeDeg - metadata_.westernmostLongitudeDeg) /
      longitudeStep - 0.5;
  const double row = (metadata_.maximumLatitudeDeg - latitudeDeg) /
      latitudeStep - 0.5;
  const long column0Raw = static_cast<long>(std::floor(column));
  const double columnWeight = column - std::floor(column);
  // There is no row beyond a pole.  Clamp to the centre of the polar row,
  // whereas columns wrap because longitude is periodic.
  const double rowClamped = std::max(0.0,
      std::min(static_cast<double>(metadata_.lines - 1), row));
  const long row0 = static_cast<long>(std::floor(rowClamped));
  const long row1 = std::min<long>(row0 + 1, metadata_.lines - 1);
  const double rowWeight = rowClamped - row0;
  const auto wrap = [this](long columnIndex) {
    const long width = static_cast<long>(metadata_.lineSamples);
    columnIndex %= width;
    if (columnIndex < 0) columnIndex += width;
    return columnIndex;
  };
  const long column0 = wrap(column0Raw);
  const long column1 = wrap(column0Raw + 1);
  const auto value = [this](long r, long c) {
    return metadata_.scalingFactor * elevationDn_[
        static_cast<std::size_t>(r) * metadata_.lineSamples +
        static_cast<std::size_t>(c)];
  };
  const double north = (1.0 - columnWeight) * value(row0, column0) +
      columnWeight * value(row0, column1);
  const double south = (1.0 - columnWeight) * value(row1, column0) +
      columnWeight * value(row1, column1);
  return (1.0 - rowWeight) * north + rowWeight * south;
}

bool BuildIcosphere(double referenceRadiusM, double requestedResolutionM,
                    Triangulation* result, std::string* error) {
  if (result == nullptr || error == nullptr || referenceRadiusM <= 0.0 ||
      requestedResolutionM <= 0.0 || !std::isfinite(referenceRadiusM) ||
      !std::isfinite(requestedResolutionM)) {
    if (error != nullptr) *error = "invalid icosphere radius or resolution";
    return false;
  }
  // Seed with a unit regular icosahedron.  Its near-uniform faces avoid the
  // vanishing east-west edges and polar singularity of a latitude/longitude
  // tessellation.  Connectivity below is the standard closed 20-face shell.
  const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  Triangulation mesh;
  mesh.verticesM = {
      Unit({-1, phi, 0}), Unit({1, phi, 0}), Unit({-1, -phi, 0}),
      Unit({1, -phi, 0}), Unit({0, -1, phi}), Unit({0, 1, phi}),
      Unit({0, -1, -phi}), Unit({0, 1, -phi}), Unit({phi, 0, -1}),
      Unit({phi, 0, 1}), Unit({-phi, 0, -1}), Unit({-phi, 0, 1})};
  mesh.faces = {
      {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
      {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
      {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
      {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};

  // Resolution is defined on the reference sphere before topography.  That
  // makes it deterministic and independent of local positive/negative relief.
  auto maximumEdge = [&]() {
    double maximum = 0.0;
    for (const auto& face : mesh.faces) {
      maximum = std::max(maximum, GreatCircleEdge(mesh.verticesM[face[0]],
          mesh.verticesM[face[1]], referenceRadiusM));
      maximum = std::max(maximum, GreatCircleEdge(mesh.verticesM[face[1]],
          mesh.verticesM[face[2]], referenceRadiusM));
      maximum = std::max(maximum, GreatCircleEdge(mesh.verticesM[face[2]],
          mesh.verticesM[face[0]], referenceRadiusM));
    }
    return maximum;
  };

  // Each refinement splits every face into four and projects shared edge
  // midpoints back to the unit sphere.  Global levels retain approximately
  // uniform resolution over the whole lunar surface.
  while (maximumEdge() > requestedResolutionM) {
    if (mesh.faces.size() > kMaximumFaces / 4u) {
      *error = "surface_mesh_resolution_m would require more than " +
          std::to_string(kMaximumFaces) + " triangular faces";
      return false;
    }
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> cache;
    std::vector<std::array<std::size_t, 3>> refined;
    refined.reserve(mesh.faces.size() * 4u);
    for (const auto& face : mesh.faces) {
      const std::size_t ab = Midpoint(face[0], face[1], &mesh.verticesM, &cache);
      const std::size_t bc = Midpoint(face[1], face[2], &mesh.verticesM, &cache);
      const std::size_t ca = Midpoint(face[2], face[0], &mesh.verticesM, &cache);
      refined.push_back({face[0], ab, ca});
      refined.push_back({face[1], bc, ab});
      refined.push_back({face[2], ca, bc});
      refined.push_back({ab, bc, ca});
    }
    mesh.faces = std::move(refined);
    ++mesh.subdivisionLevel;
  }
  mesh.maximumEdgeLengthM = maximumEdge();
  // Store final geometry in AMPS SI units only after the angular refinement is
  // complete; all connectivity and directions remain unchanged.
  for (auto& vertex : mesh.verticesM) {
    for (double& component : vertex) component *= referenceRadiusM;
  }
  mesh.minimumRadiusM = referenceRadiusM;
  mesh.maximumRadiusM = referenceRadiusM;
  OrientOutward(&mesh);
  *result = std::move(mesh);
  error->clear();
  return true;
}

bool ApplyLolaTopography(const LolaDem& dem, Triangulation* mesh,
                         std::string* error) {
  if (mesh == nullptr || error == nullptr || mesh->verticesM.empty()) return false;
  double minimum = std::numeric_limits<double>::max();
  double maximum = 0.0;
  // Vertex directions are interpreted as planetocentric latitude and
  // east-positive longitude in the label's MEAN EARTH/POLAR AXIS OF DE421
  // axes.  No IAU_MOON or epoch-dependent transform is applied here.
  for (auto& vertex : mesh->verticesM) {
    const double radius = Norm(vertex);
    if (!(radius > 0.0)) {
      *error = "surface triangulation contains a zero-radius vertex";
      return false;
    }
    const double latitude = std::asin(vertex[2] / radius) * 180.0 / kPi;
    double longitude = std::atan2(vertex[1], vertex[0]) * 180.0 / kPi;
    if (longitude < 0.0) longitude += 360.0;
    const double adjustedRadius = dem.metadata().referenceRadiusM +
        dem.ElevationM(longitude, latitude);
    if (!(adjustedRadius > 0.0) || !std::isfinite(adjustedRadius)) {
      *error = "LOLA interpolation produced an invalid planetary radius";
      return false;
    }
    const double scale = adjustedRadius / radius;
    for (double& component : vertex) component *= scale;
    minimum = std::min(minimum, adjustedRadius);
    maximum = std::max(maximum, adjustedRadius);
  }
  mesh->minimumRadiusM = minimum;
  mesh->maximumRadiusM = maximum;
  OrientOutward(mesh);
  error->clear();
  return true;
}

bool WriteCeaSurface(const Triangulation& mesh, const std::string& path,
                     std::string* error) {
  // CEA long format: counts, Cartesian nodes, then triangle connectivity.
  // Coordinates are metres and connectivity is one-based on disk.
  return AtomicWrite(path, [&mesh](std::ostream& output) {
    output << mesh.verticesM.size() << ' ' << mesh.faces.size() << '\n';
    output << std::setprecision(17);
    for (const auto& vertex : mesh.verticesM)
      output << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << '\n';
    // The AMPS CEA reader subtracts one from each stored node number.
    for (const auto& face : mesh.faces)
      output << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1
             << '\n';
  }, error);
}

bool WriteTecplotSurface(const Triangulation& mesh, const std::string& path,
                         std::string* error) {
  // Keep this representation directly viewable: point-packed Cartesian SI
  // coordinates plus radius, followed by one-based FETRIANGLE connectivity.
  return AtomicWrite(path, [&mesh](std::ostream& output) {
    output << "VARIABLES=\"X [m]\",\"Y [m]\",\"Z [m]\",\"Radius [m]\"\n";
    output << "ZONE N=" << mesh.verticesM.size() << ", E=" << mesh.faces.size()
           << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n";
    output << std::setprecision(17);
    for (const auto& vertex : mesh.verticesM)
      output << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << ' '
             << Norm(vertex) << '\n';
    for (const auto& face : mesh.faces)
      output << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1
             << '\n';
  }, error);
}

bool PrepareProductionSurface(const Runtime::Configuration& configuration,
                              std::string* error) {
  // Sphere mode is a true legacy path: do not generate files or register a
  // triangulation, and leave the historical analytic boundary initialization
  // to main_lib.cpp.
  if (configuration.surfaceGeometry == Runtime::SurfaceGeometry::Sphere) {
    gRealisticSurfaceActive = false;
    gMinimumSurfaceRadiusM = 0.0;
    return true;
  }

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  // The cut-cell triangulation is static in the AMR frame.  Enabling lunar
  // rotation without moving/rebuilding it would mix the label body-fixed frame
  // with LSO and produce physically incorrect impacts, so fail explicitly.
  *error = "LOLA surface mode currently requires orbit calculation off: the "
      "static AMR triangulation is in the LDEM body-fixed frame";
  return false;
#endif

  // Only rank zero reads the large raster and writes shared products.  The
  // success flag and minimum radius are broadcast before any rank invokes the
  // normal AMPS CEA reader, ensuring all ranks see one complete mesh.
  int generated = 1;
  std::string localError;
  double minimumRadius = 0.0;
  if (PIC::ThisThread == 0) {
    LolaDem dem;
    Triangulation mesh;
    generated = dem.Load(configuration, &localError) &&
        BuildIcosphere(dem.metadata().referenceRadiusM,
                       configuration.surfaceMeshResolutionM, &mesh,
                       &localError) &&
        ApplyLolaTopography(dem, &mesh, &localError) &&
        WriteCeaSurface(mesh, configuration.surfaceCeaFile, &localError) &&
        WriteTecplotSurface(mesh, configuration.surfaceTecplotFile,
                            &localError);
    if (generated) {
      minimumRadius = mesh.minimumRadiusM;
      std::cout << "$PREFIX: srcMoon LOLA surface: product="
                << dem.metadata().productId << ", frame=\""
                << dem.metadata().coordinateSystemName << "\", vertices="
                << mesh.verticesM.size() << ", faces=" << mesh.faces.size()
                << ", subdivision=" << mesh.subdivisionLevel
                << ", maximum_edge_m=" << mesh.maximumEdgeLengthM << '\n';
    }
  }
  MPI_Bcast(&generated, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (!generated) {
    if (PIC::ThisThread == 0) *error = localError;
    else *error = "rank zero failed to generate the LOLA surface";
    return false;
  }
  MPI_Bcast(&minimumRadius, 1, MPI_DOUBLE, 0, MPI_GLOBAL_COMMUNICATOR);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  // Use the same generic loader as a pre-generated AMPS surface.  This keeps
  // production ray tracing and cut-cell construction on the standard path.
  PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(
      configuration.surfaceCeaFile.c_str(), 1.0);
  if (!InitializeLoadedNormalsFromWinding(error)) return false;
  gMinimumSurfaceRadiusM = minimumRadius;
  gRealisticSurfaceActive = true;
  error->clear();
  return true;
}

bool RealisticSurfaceActive() { return gRealisticSurfaceActive; }
double MinimumLoadedSurfaceRadiusM() { return gMinimumSurfaceRadiusM; }

}  // namespace Surface
}  // namespace Moon
