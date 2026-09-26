#ifndef _SRC_EARTH_UTIL_SWMF_SNAPSHOT_CONTRACT_H_
#define _SRC_EARTH_UTIL_SWMF_SNAPSHOT_CONTRACT_H_

//======================================================================================
// SWMFSnapshotContract.h -- Roadmap Step 9
//======================================================================================
// Dependency-free interchange contract for one frozen SWMF magnetic-field/plasma-
// velocity state.  The live AMPS/SWMF component and standalone Mode3D replay use this
// exact format, so a cross-path comparison is a comparison of the same B/u samples,
// not a comparison of two independently reconstructed phenomenological fields.
//
// Phase-1 physics scope
// ---------------------
// A file represents one instantaneous/quasi-static GSM snapshot.  B and bulk velocity
// are frozen at cell centres; E is defined (not independently fitted) by ideal MHD,
//
//                         E = -u x B .
//
// The contract does not claim time-dependent characteristics, electric acceleration,
// trapping/loss evolution, or interpolation between coupling times.  Those are later
// roadmap capabilities.  Step 9 only guarantees that every trajectory in one cutoff/
// access/flux batch sees one validated, immutable state.
//
// Why this header has no AMPS/PIC/MPI dependency
// ----------------------------------------------
// File parsing, dimensional validation, deterministic identity, restart reproducibility,
// and live-versus-replay comparison must be testable without an SWMF installation.  The
// production AMR adapter is in 3d/GlobalMagneticField.cpp; this header owns only the
// portable records and strict validation rules.
//======================================================================================

#include "FieldProvider.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace Earth {
namespace SWMFSnapshot {

static const char* const kSchema="sep-in-geospace/swmf-field-snapshot/v1";
static const char* const kFrame="GSM";
static const char* const kPositionUnit="m";
static const char* const kMagneticFieldUnit="T";
static const char* const kPlasmaVelocityUnit="m/s";
static const char* const kElectricFieldConvention="E=-u_cross_B";
static const char* const kMagneticOnlyMode="MAGNETIC_ONLY";
static const char* const kExperimentalIdealMhdMode="EXPERIMENTAL_IDEAL_MHD";

struct Cell {
  long int blockId;
  int i,j,k;
  double position_m[3];
  double magneticField_T[3];
  double plasmaVelocity_m_s[3];

  Cell() : blockId(-1),i(-1),j(-1),k(-1) {
    for (int d=0;d<3;++d) {
      position_m[d]=0.0;
      magneticField_T[d]=0.0;
      plasmaVelocity_m_s[d]=0.0;
    }
  }
};

inline bool CellLess(const Cell& a,const Cell& b) {
  return std::tie(a.blockId,a.i,a.j,a.k)<std::tie(b.blockId,b.i,b.j,b.k);
}

struct Snapshot {
  std::string schema{kSchema};
  std::string snapshotId;
  std::string contentFingerprint;
  // Content-derived identity of the AMR keys, cell centres, dimensions, and box.
  // Keeping it separate from the B/u fingerprint distinguishes a stale mesh from a
  // stale field receive without relying on an implementation-specific pointer value.
  std::string meshRevision;
  std::string epochUTC;
  // Authoritative PT time relative to the configured coupling epoch [s].
  double simulationTime_s{0.0};
  std::string frame{kFrame};
  std::string positionUnit{kPositionUnit};
  std::string magneticFieldUnit{kMagneticFieldUnit};
  std::string plasmaVelocityUnit{kPlasmaVelocityUnit};
  std::string electricFieldConvention{kElectricFieldConvention};
  // Phase 1 is magnetic-only by default.  Ideal-MHD E can be published only under
  // the explicitly experimental mode pending its later validation gates.
  std::string electricFieldMode{kMagneticOnlyMode};
  long int usedLeafBlocks{0};
  int blockCellsX{0},blockCellsY{0},blockCellsZ{0};
  Earth::Field::AxisAlignedDomainSI domain;
  std::vector<Cell> cells;
};

inline void DeriveElectricField(const double velocity_m_s[3],
                                const double magneticField_T[3],
                                double electricField_V_m[3]) {
  // E=-u x B.  Keeping this formula in the portable contract ensures the live gather,
  // exported file, standalone replay, and comparison tests use one sign convention.
  electricField_V_m[0]=-(velocity_m_s[1]*magneticField_T[2]-
                          velocity_m_s[2]*magneticField_T[1]);
  electricField_V_m[1]=-(velocity_m_s[2]*magneticField_T[0]-
                          velocity_m_s[0]*magneticField_T[2]);
  electricField_V_m[2]=-(velocity_m_s[0]*magneticField_T[1]-
                          velocity_m_s[1]*magneticField_T[0]);
}

namespace Detail {

class StableHash64 {
public:
  StableHash64() : value_(14695981039346656037ULL) {}

  void AddByte(unsigned char byte) {
    value_^=static_cast<std::uint64_t>(byte);
    value_*=1099511628211ULL;
  }

  void AddUnsigned(std::uint64_t value) {
    // Always feed least-significant byte first.  The result is independent of host
    // endianness and therefore suitable for restart/cross-layout identity checks.
    for (int shift=0;shift<64;shift+=8)
      AddByte(static_cast<unsigned char>((value>>shift)&0xffU));
  }

  void AddSigned(std::int64_t value) {
    AddUnsigned(static_cast<std::uint64_t>(value));
  }

  void AddString(const std::string& value) {
    AddUnsigned(static_cast<std::uint64_t>(value.size()));
    for (std::string::const_iterator it=value.begin();it!=value.end();++it)
      AddByte(static_cast<unsigned char>(*it));
  }

  void AddDouble(double value) {
    if (!std::isfinite(value))
      throw std::invalid_argument("Cannot hash a non-finite SWMF snapshot value");
    // Canonicalize signed zero: +0 and -0 are the same physical state.
    if (value==0.0) value=0.0;
    static_assert(sizeof(double)==sizeof(std::uint64_t),
                  "SWMF snapshot identity requires 64-bit double precision");
    std::uint64_t bits=0;
    std::memcpy(&bits,&value,sizeof(bits));
    AddUnsigned(bits);
  }

  std::uint64_t Value() const { return value_; }

  std::string Hex(const char* prefix) const {
    std::ostringstream out;
    out << (prefix!=NULL ? prefix : "") << std::hex << std::setw(16)
        << std::setfill('0') << value_;
    return out.str();
  }

private:
  std::uint64_t value_;
};

inline std::string Trim(const std::string& value) {
  std::string::size_type begin=0,end=value.size();
  while (begin<end && (value[begin]==' ' || value[begin]=='\t' ||
                       value[begin]=='\r' || value[begin]=='\n')) ++begin;
  while (end>begin && (value[end-1]==' ' || value[end-1]=='\t' ||
                       value[end-1]=='\r' || value[end-1]=='\n')) --end;
  return value.substr(begin,end-begin);
}

inline std::vector<std::string> SplitCsv(const std::string& line) {
  // Snapshot tokens never contain commas.  Rejecting CSV quoting keeps the format
  // deterministic and makes accidental locale-dependent output impossible.
  std::vector<std::string> result;
  std::string token;
  std::istringstream input(line);
  while (std::getline(input,token,',')) result.push_back(Trim(token));
  if (!line.empty() && line[line.size()-1]==',') result.push_back("");
  return result;
}

inline long int ParseLong(const std::string& token,const std::string& label) {
  errno=0;
  char* end=NULL;
  const long value=std::strtol(token.c_str(),&end,10);
  if (errno!=0 || end==token.c_str() || *end!='\0')
    throw std::invalid_argument("Invalid integer for "+label+": '"+token+"'");
  return value;
}

inline double ParseDouble(const std::string& token,const std::string& label) {
  errno=0;
  char* end=NULL;
  const double value=std::strtod(token.c_str(),&end);
  if (errno!=0 || end==token.c_str() || *end!='\0' || !std::isfinite(value))
    throw std::invalid_argument("Invalid finite number for "+label+": '"+token+"'");
  return value;
}

inline std::string RequireKey(const std::map<std::string,std::string>& values,
                              const std::string& key) {
  const std::map<std::string,std::string>::const_iterator it=values.find(key);
  if (it==values.end() || it->second.empty())
    throw std::invalid_argument("SWMF snapshot is missing metadata key '"+key+"'");
  return it->second;
}

inline std::string JsonEscape(const std::string& value) {
  std::string out;
  for (std::string::const_iterator it=value.begin();it!=value.end();++it) {
    const char c=*it;
    if (c=='\\') out+="\\\\";
    else if (c=='\"') out+="\\\"";
    else if (c=='\n') out+="\\n";
    else if (c=='\r') out+="\\r";
    else if (c=='\t') out+="\\t";
    else out.push_back(c);
  }
  return out;
}

} // namespace Detail

inline std::vector<Cell> CanonicalCells(const Snapshot& snapshot) {
  std::vector<Cell> cells=snapshot.cells;
  std::sort(cells.begin(),cells.end(),CellLess);
  return cells;
}

inline bool UsesExperimentalIdealMhdElectricField(const Snapshot& snapshot) {
  return snapshot.electricFieldMode==kExperimentalIdealMhdMode;
}

inline std::string ComputeMeshRevision(const Snapshot& snapshot) {
  // The mesh identity excludes B, u, epoch, and product options. Repartitioning and
  // traversal order leave it unchanged; topology/domain/centre changes do not.
  Detail::StableHash64 hash;
  hash.AddString("sep-in-geospace/swmf-mesh/v1");
  hash.AddSigned(snapshot.usedLeafBlocks);
  hash.AddSigned(snapshot.blockCellsX);
  hash.AddSigned(snapshot.blockCellsY);
  hash.AddSigned(snapshot.blockCellsZ);
  hash.AddUnsigned(snapshot.domain.enabled ? 1U : 0U);
  for (int d=0;d<3;++d) {
    hash.AddDouble(snapshot.domain.minimum_m[d]);
    hash.AddDouble(snapshot.domain.maximum_m[d]);
  }
  const std::vector<Cell> cells=CanonicalCells(snapshot);
  hash.AddUnsigned(static_cast<std::uint64_t>(cells.size()));
  for (std::vector<Cell>::const_iterator it=cells.begin();it!=cells.end();++it) {
    hash.AddSigned(it->blockId);
    hash.AddSigned(it->i);
    hash.AddSigned(it->j);
    hash.AddSigned(it->k);
    for (int d=0;d<3;++d) hash.AddDouble(it->position_m[d]);
  }
  return hash.Hex("swmf-mesh-v1-");
}

inline std::string ComputeContentFingerprint(const Snapshot& snapshot) {
  Detail::StableHash64 hash;
  hash.AddString(snapshot.schema);
  hash.AddString(snapshot.epochUTC);
  hash.AddDouble(snapshot.simulationTime_s);
  hash.AddString(snapshot.frame);
  hash.AddString(snapshot.positionUnit);
  hash.AddString(snapshot.magneticFieldUnit);
  hash.AddString(snapshot.plasmaVelocityUnit);
  hash.AddString(snapshot.electricFieldConvention);
  hash.AddString(snapshot.electricFieldMode);
  hash.AddString(snapshot.meshRevision);
  hash.AddSigned(snapshot.usedLeafBlocks);
  hash.AddSigned(snapshot.blockCellsX);
  hash.AddSigned(snapshot.blockCellsY);
  hash.AddSigned(snapshot.blockCellsZ);
  hash.AddUnsigned(snapshot.domain.enabled ? 1U : 0U);
  for (int d=0;d<3;++d) {
    hash.AddDouble(snapshot.domain.minimum_m[d]);
    hash.AddDouble(snapshot.domain.maximum_m[d]);
  }

  const std::vector<Cell> cells=CanonicalCells(snapshot);
  hash.AddUnsigned(static_cast<std::uint64_t>(cells.size()));
  for (std::vector<Cell>::const_iterator it=cells.begin();it!=cells.end();++it) {
    hash.AddSigned(it->blockId);
    hash.AddSigned(it->i);
    hash.AddSigned(it->j);
    hash.AddSigned(it->k);
    for (int d=0;d<3;++d) hash.AddDouble(it->position_m[d]);
    for (int d=0;d<3;++d) hash.AddDouble(it->magneticField_T[d]);
    for (int d=0;d<3;++d) hash.AddDouble(it->plasmaVelocity_m_s[d]);
  }
  return hash.Hex("swmf-state-v1-");
}

inline std::string SnapshotIdFromFingerprint(const std::string& epochUTC,
                                             const std::string& fingerprint) {
  return Earth::Field::MakeSnapshotId(
      "PIC::CPLR:SWMF",epochUTC,"content="+fingerprint);
}

inline void FinalizeIdentity(Snapshot& snapshot) {
  snapshot.meshRevision=ComputeMeshRevision(snapshot);
  snapshot.contentFingerprint=ComputeContentFingerprint(snapshot);
  snapshot.snapshotId=SnapshotIdFromFingerprint(
      snapshot.epochUTC,snapshot.contentFingerprint);
}

inline void Validate(const Snapshot& snapshot,bool requireDeclaredIdentity=true) {
  if (snapshot.schema!=kSchema)
    throw std::invalid_argument("Unsupported SWMF snapshot schema '"+
                                snapshot.schema+"'");
  if (snapshot.epochUTC.empty())
    throw std::invalid_argument("SWMF snapshot epoch_utc is empty");
  if (!std::isfinite(snapshot.simulationTime_s) || snapshot.simulationTime_s<0.0)
    throw std::invalid_argument(
        "SWMF snapshot simulation_time_s must be finite and nonnegative");
  if (snapshot.frame!=kFrame || snapshot.positionUnit!=kPositionUnit ||
      snapshot.magneticFieldUnit!=kMagneticFieldUnit ||
      snapshot.plasmaVelocityUnit!=kPlasmaVelocityUnit ||
      snapshot.electricFieldConvention!=kElectricFieldConvention)
    throw std::invalid_argument(
        "SWMF snapshot must use GSM, position m, B tesla, velocity m/s, and E=-u_cross_B");
  if (snapshot.electricFieldMode!=kMagneticOnlyMode &&
      snapshot.electricFieldMode!=kExperimentalIdealMhdMode)
    throw std::invalid_argument(
        "SWMF snapshot electric_field_mode must be MAGNETIC_ONLY or "
        "EXPERIMENTAL_IDEAL_MHD");
  if (!snapshot.domain.enabled)
    throw std::invalid_argument("SWMF snapshot requires an explicit validity domain");
  for (int d=0;d<3;++d) {
    if (!std::isfinite(snapshot.domain.minimum_m[d]) ||
        !std::isfinite(snapshot.domain.maximum_m[d]) ||
        snapshot.domain.minimum_m[d]>=snapshot.domain.maximum_m[d])
      throw std::invalid_argument("SWMF snapshot has an invalid domain bound");
  }
  if (snapshot.usedLeafBlocks<=0 || snapshot.blockCellsX<=0 ||
      snapshot.blockCellsY<=0 || snapshot.blockCellsZ<=0)
    throw std::invalid_argument("SWMF snapshot has invalid AMR topology dimensions");

  const std::uint64_t cellsPerBlock=
      static_cast<std::uint64_t>(snapshot.blockCellsX)*
      static_cast<std::uint64_t>(snapshot.blockCellsY)*
      static_cast<std::uint64_t>(snapshot.blockCellsZ);
  const std::uint64_t expected=
      static_cast<std::uint64_t>(snapshot.usedLeafBlocks)*cellsPerBlock;
  if (expected!=static_cast<std::uint64_t>(snapshot.cells.size()))
    throw std::invalid_argument(
        "SWMF snapshot cell count does not equal blocks*Nx*Ny*Nz");

  std::set<std::tuple<long int,int,int,int> > keys;
  for (std::vector<Cell>::const_iterator it=snapshot.cells.begin();
       it!=snapshot.cells.end();++it) {
    if (it->blockId<0 || it->blockId>=snapshot.usedLeafBlocks ||
        it->i<0 || it->i>=snapshot.blockCellsX ||
        it->j<0 || it->j>=snapshot.blockCellsY ||
        it->k<0 || it->k>=snapshot.blockCellsZ)
      throw std::invalid_argument("SWMF snapshot contains an out-of-range cell index");
    const std::tuple<long int,int,int,int> key(it->blockId,it->i,it->j,it->k);
    if (!keys.insert(key).second)
      throw std::invalid_argument("SWMF snapshot contains a duplicate cell index");
    if (!Earth::Field::FiniteVector3(it->position_m) ||
        !Earth::Field::FiniteVector3(it->magneticField_T) ||
        !Earth::Field::FiniteVector3(it->plasmaVelocity_m_s))
      throw std::invalid_argument("SWMF snapshot contains a non-finite cell value");
    for (int d=0;d<3;++d) {
      const double span=snapshot.domain.maximum_m[d]-snapshot.domain.minimum_m[d];
      const double tolerance=64.0*std::numeric_limits<double>::epsilon()*
          std::max(1.0,std::fabs(span));
      if (it->position_m[d]<snapshot.domain.minimum_m[d]-tolerance ||
          it->position_m[d]>snapshot.domain.maximum_m[d]+tolerance)
        throw std::invalid_argument("SWMF snapshot cell centre is outside its domain");
    }
    double electric[3];
    DeriveElectricField(it->plasmaVelocity_m_s,it->magneticField_T,electric);
    if (!Earth::Field::FiniteVector3(electric))
      throw std::invalid_argument("SWMF snapshot produces non-finite E=-u x B");
  }

  const std::string meshRevision=ComputeMeshRevision(snapshot);
  if (requireDeclaredIdentity && snapshot.meshRevision!=meshRevision)
    throw std::invalid_argument(
        "SWMF snapshot mesh_revision does not match its topology/cell centres");
  const std::string fingerprint=ComputeContentFingerprint(snapshot);
  const std::string snapshotId=SnapshotIdFromFingerprint(snapshot.epochUTC,fingerprint);
  if (requireDeclaredIdentity &&
      (snapshot.contentFingerprint!=fingerprint || snapshot.snapshotId!=snapshotId))
    throw std::invalid_argument(
        "SWMF snapshot content fingerprint or snapshot_id does not match its data");
}

inline void Write(const Snapshot& input,const std::string& fileName) {
  Snapshot snapshot=input;
  FinalizeIdentity(snapshot);
  Validate(snapshot,true);
  snapshot.cells=CanonicalCells(snapshot);

  std::ofstream output(fileName.c_str(),std::ios::out|std::ios::trunc);
  if (!output) throw std::runtime_error("Cannot write SWMF snapshot file: "+fileName);
  output.imbue(std::locale::classic());
  output << "# schema=" << snapshot.schema << "\n"
         << "# snapshot_id=" << snapshot.snapshotId << "\n"
         << "# content_fingerprint=" << snapshot.contentFingerprint << "\n"
         << "# mesh_revision=" << snapshot.meshRevision << "\n"
         << "# epoch_utc=" << snapshot.epochUTC << "\n"
         << std::setprecision(17)
         << "# simulation_time_s=" << snapshot.simulationTime_s << "\n"
         << "# frame=" << snapshot.frame << "\n"
         << "# position_unit=" << snapshot.positionUnit << "\n"
         << "# magnetic_field_unit=" << snapshot.magneticFieldUnit << "\n"
         << "# plasma_velocity_unit=" << snapshot.plasmaVelocityUnit << "\n"
         << "# electric_field_convention=" << snapshot.electricFieldConvention << "\n"
         << "# electric_field_mode=" << snapshot.electricFieldMode << "\n"
         << "# used_leaf_blocks=" << snapshot.usedLeafBlocks << "\n"
         << "# block_cells_x=" << snapshot.blockCellsX << "\n"
         << "# block_cells_y=" << snapshot.blockCellsY << "\n"
         << "# block_cells_z=" << snapshot.blockCellsZ << "\n";
  output << std::setprecision(17)
         << "# domain_min_x_m=" << snapshot.domain.minimum_m[0] << "\n"
         << "# domain_min_y_m=" << snapshot.domain.minimum_m[1] << "\n"
         << "# domain_min_z_m=" << snapshot.domain.minimum_m[2] << "\n"
         << "# domain_max_x_m=" << snapshot.domain.maximum_m[0] << "\n"
         << "# domain_max_y_m=" << snapshot.domain.maximum_m[1] << "\n"
         << "# domain_max_z_m=" << snapshot.domain.maximum_m[2] << "\n"
         << "block_id,i,j,k,x_m,y_m,z_m,bx_T,by_T,bz_T,ux_m_s,uy_m_s,uz_m_s\n";
  for (std::vector<Cell>::const_iterator it=snapshot.cells.begin();
       it!=snapshot.cells.end();++it) {
    output << it->blockId << ',' << it->i << ',' << it->j << ',' << it->k;
    for (int d=0;d<3;++d) output << ',' << it->position_m[d];
    for (int d=0;d<3;++d) output << ',' << it->magneticField_T[d];
    for (int d=0;d<3;++d) output << ',' << it->plasmaVelocity_m_s[d];
    output << '\n';
  }
  output.flush();
  if (!output) throw std::runtime_error("Failed while writing SWMF snapshot: "+fileName);
}

inline Snapshot Read(const std::string& fileName) {
  std::ifstream input(fileName.c_str());
  if (!input) throw std::runtime_error("Cannot open SWMF snapshot file: "+fileName);
  input.imbue(std::locale::classic());

  std::map<std::string,std::string> metadata;
  Snapshot snapshot;
  std::string line;
  bool columnsSeen=false;
  long int lineNo=0;
  const std::string expectedColumns=
      "block_id,i,j,k,x_m,y_m,z_m,bx_T,by_T,bz_T,ux_m_s,uy_m_s,uz_m_s";
  while (std::getline(input,line)) {
    ++lineNo;
    line=Detail::Trim(line);
    if (line.empty()) continue;
    if (!columnsSeen && line[0]=='#') {
      const std::string item=Detail::Trim(line.substr(1));
      const std::string::size_type equal=item.find('=');
      if (equal==std::string::npos)
        throw std::invalid_argument("Malformed SWMF metadata line "+
                                    std::to_string(lineNo));
      const std::string key=Detail::Trim(item.substr(0,equal));
      const std::string value=Detail::Trim(item.substr(equal+1));
      if (key.empty() || value.empty() || !metadata.insert(std::make_pair(key,value)).second)
        throw std::invalid_argument("Duplicate/empty SWMF metadata at line "+
                                    std::to_string(lineNo));
      continue;
    }
    if (!columnsSeen) {
      if (line!=expectedColumns)
        throw std::invalid_argument("Unexpected SWMF snapshot column schema");
      columnsSeen=true;
      continue;
    }
    const std::vector<std::string> values=Detail::SplitCsv(line);
    if (values.size()!=13)
      throw std::invalid_argument("SWMF snapshot row must contain exactly 13 columns");
    Cell cell;
    cell.blockId=Detail::ParseLong(values[0],"block_id");
    cell.i=static_cast<int>(Detail::ParseLong(values[1],"i"));
    cell.j=static_cast<int>(Detail::ParseLong(values[2],"j"));
    cell.k=static_cast<int>(Detail::ParseLong(values[3],"k"));
    for (int d=0;d<3;++d) cell.position_m[d]=Detail::ParseDouble(values[4+d],"position");
    for (int d=0;d<3;++d) cell.magneticField_T[d]=Detail::ParseDouble(values[7+d],"B");
    for (int d=0;d<3;++d) cell.plasmaVelocity_m_s[d]=Detail::ParseDouble(values[10+d],"u");
    snapshot.cells.push_back(cell);
  }
  if (!columnsSeen) throw std::invalid_argument("SWMF snapshot has no column header");

  snapshot.schema=Detail::RequireKey(metadata,"schema");
  snapshot.snapshotId=Detail::RequireKey(metadata,"snapshot_id");
  snapshot.contentFingerprint=Detail::RequireKey(metadata,"content_fingerprint");
  snapshot.meshRevision=Detail::RequireKey(metadata,"mesh_revision");
  snapshot.epochUTC=Detail::RequireKey(metadata,"epoch_utc");
  snapshot.simulationTime_s=Detail::ParseDouble(
      Detail::RequireKey(metadata,"simulation_time_s"),"simulation_time_s");
  snapshot.frame=Detail::RequireKey(metadata,"frame");
  snapshot.positionUnit=Detail::RequireKey(metadata,"position_unit");
  snapshot.magneticFieldUnit=Detail::RequireKey(metadata,"magnetic_field_unit");
  snapshot.plasmaVelocityUnit=Detail::RequireKey(metadata,"plasma_velocity_unit");
  snapshot.electricFieldConvention=Detail::RequireKey(metadata,"electric_field_convention");
  snapshot.electricFieldMode=Detail::RequireKey(metadata,"electric_field_mode");
  snapshot.usedLeafBlocks=Detail::ParseLong(
      Detail::RequireKey(metadata,"used_leaf_blocks"),"used_leaf_blocks");
  snapshot.blockCellsX=static_cast<int>(Detail::ParseLong(
      Detail::RequireKey(metadata,"block_cells_x"),"block_cells_x"));
  snapshot.blockCellsY=static_cast<int>(Detail::ParseLong(
      Detail::RequireKey(metadata,"block_cells_y"),"block_cells_y"));
  snapshot.blockCellsZ=static_cast<int>(Detail::ParseLong(
      Detail::RequireKey(metadata,"block_cells_z"),"block_cells_z"));
  snapshot.domain.enabled=true;
  snapshot.domain.minimum_m[0]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_min_x_m"),"domain_min_x_m");
  snapshot.domain.minimum_m[1]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_min_y_m"),"domain_min_y_m");
  snapshot.domain.minimum_m[2]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_min_z_m"),"domain_min_z_m");
  snapshot.domain.maximum_m[0]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_max_x_m"),"domain_max_x_m");
  snapshot.domain.maximum_m[1]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_max_y_m"),"domain_max_y_m");
  snapshot.domain.maximum_m[2]=Detail::ParseDouble(
      Detail::RequireKey(metadata,"domain_max_z_m"),"domain_max_z_m");

  // Unknown metadata are rejected.  This prevents a newer physical convention from
  // being silently read by an older executable that cannot implement it.
  if (metadata.size()!=22)
    throw std::invalid_argument("SWMF snapshot contains unknown or incomplete metadata");
  Validate(snapshot,true);
  return snapshot;
}

// PublicationQueue is the dependency-free model of the live field lifecycle.  The
// production AMR adapter uses the same state transitions around its compact arrays:
// one READY generation may be frozen for a product batch; a later complete receive is
// QUEUED and cannot replace the active generation until Release() verifies the batch
// identity. Invalid/incomplete publishes latch FAILED (or a pending failure behind the
// frozen batch), so a later product cannot fall back to an older/partial generation.
enum class PublicationState { Unavailable, Ready, Frozen, Queued, Failed };

class PublicationQueue {
public:
  PublicationState State() const {
    if (failed_ || pendingFailure_) return PublicationState::Failed;
    if (!hasActive_) return PublicationState::Unavailable;
    if (frozen_ && hasPending_) return PublicationState::Queued;
    if (frozen_) return PublicationState::Frozen;
    return PublicationState::Ready;
  }

  bool HasPending() const { return hasPending_; }

  void Publish(const Snapshot& snapshot) {
    try {
      Validate(snapshot,true);
    }
    catch (...) {
      // Never silently keep using an older READY generation after a corrupt receive.
      // A currently frozen batch may finish from its immutable active arrays, but its
      // release transitions to FAILED instead of falling back to that old generation.
      if (frozen_) pendingFailure_=true;
      else {
        hasActive_=false;
        failed_=true;
      }
      throw;
    }
    if (frozen_) {
      // Only the newest complete coupling receive is useful after a frozen batch.  A
      // replacement here is safe because neither pending snapshot has been exposed to
      // a trajectory; the active generation remains unchanged until Release().
      pending_=snapshot;
      hasPending_=true;
      pendingFailure_=false;
      return;
    }
    active_=snapshot;
    hasActive_=true;
    failed_=false;
  }

  const Snapshot& Freeze() {
    if (!hasActive_)
      throw std::runtime_error("Cannot freeze an unavailable SWMF snapshot");
    if (frozen_)
      throw std::runtime_error("SWMF snapshot is already frozen for a product batch");
    frozen_=true;
    return active_;
  }

  const Snapshot& Active() const {
    if (!hasActive_ || (failed_ && !frozen_))
      throw std::runtime_error("SWMF snapshot is unavailable or failed validation");
    return active_;
  }

  void Release(const std::string& expectedSnapshotId) {
    if (!frozen_)
      throw std::runtime_error("Cannot release an SWMF snapshot that is not frozen");
    if (active_.snapshotId!=expectedSnapshotId)
      throw std::runtime_error("Stale SWMF batch attempted to release another snapshot");
    frozen_=false;
    if (pendingFailure_) {
      hasActive_=false;
      hasPending_=false;
      pendingFailure_=false;
      failed_=true;
    }
    else if (hasPending_) {
      active_=pending_;
      hasPending_=false;
      failed_=false;
    }
  }

private:
  Snapshot active_;
  Snapshot pending_;
  bool hasActive_{false};
  bool hasPending_{false};
  bool frozen_{false};
  bool failed_{false};
  bool pendingFailure_{false};
};

struct Comparison {
  bool passed{false};
  std::size_t comparedCells{0};
  std::size_t failedComponents{0};
  double maxAbsPosition_m{0.0};
  double maxRelativeB{0.0};
  double maxRelativeVelocity{0.0};
  double maxRelativeElectric{0.0};
  std::string message;
};

// Absolute tolerances carry different physical units and therefore must not share one
// scalar in a release comparison.  For example, a useful velocity tolerance measured
// in m/s would be catastrophically large if reused for B in tesla.  Relative tolerance
// is dimensionless and is applied component by component in addition to the matching
// quantity-specific absolute floor.
struct ComparisonTolerances {
  double positionAbsolute_m{0.0};
  double magneticFieldAbsolute_T{0.0};
  double plasmaVelocityAbsolute_m_s{0.0};
  double electricFieldAbsolute_V_m{0.0};
  double relative{0.0};
};

inline double RelativeDifference(double a,double b,double absoluteFloor) {
  return std::fabs(a-b)/std::max(absoluteFloor,std::max(std::fabs(a),std::fabs(b)));
}

inline Comparison Compare(const Snapshot& left,const Snapshot& right,
                          const ComparisonTolerances& tolerance) {
  Validate(left,true);
  Validate(right,true);
  if (!(tolerance.positionAbsolute_m>=0.0) ||
      !(tolerance.magneticFieldAbsolute_T>=0.0) ||
      !(tolerance.plasmaVelocityAbsolute_m_s>=0.0) ||
      !(tolerance.electricFieldAbsolute_V_m>=0.0) ||
      !(tolerance.relative>=0.0) ||
      !std::isfinite(tolerance.positionAbsolute_m) ||
      !std::isfinite(tolerance.magneticFieldAbsolute_T) ||
      !std::isfinite(tolerance.plasmaVelocityAbsolute_m_s) ||
      !std::isfinite(tolerance.electricFieldAbsolute_V_m) ||
      !std::isfinite(tolerance.relative))
    throw std::invalid_argument("SWMF comparison tolerances must be finite and nonnegative");

  Comparison result;
  if (left.epochUTC!=right.epochUTC ||
      left.simulationTime_s!=right.simulationTime_s ||
      left.frame!=right.frame ||
      left.electricFieldMode!=right.electricFieldMode ||
      left.usedLeafBlocks!=right.usedLeafBlocks ||
      left.blockCellsX!=right.blockCellsX ||
      left.blockCellsY!=right.blockCellsY ||
      left.blockCellsZ!=right.blockCellsZ ||
      left.cells.size()!=right.cells.size()) {
    result.message="SWMF snapshots have different epoch/frame/topology dimensions";
    return result;
  }
  for (int d=0;d<3;++d) {
    if (left.domain.minimum_m[d]!=right.domain.minimum_m[d] ||
        left.domain.maximum_m[d]!=right.domain.maximum_m[d]) {
      result.message="SWMF snapshots have different validity domains";
      return result;
    }
  }

  const std::vector<Cell> a=CanonicalCells(left);
  const std::vector<Cell> b=CanonicalCells(right);
  result.comparedCells=a.size();
  const double bFloor=std::max(tolerance.magneticFieldAbsolute_T,
                               std::numeric_limits<double>::min());
  const double uFloor=std::max(tolerance.plasmaVelocityAbsolute_m_s,
                               std::numeric_limits<double>::min());
  const double eFloor=std::max(tolerance.electricFieldAbsolute_V_m,
                               std::numeric_limits<double>::min());
  for (std::size_t n=0;n<a.size();++n) {
    if (std::tie(a[n].blockId,a[n].i,a[n].j,a[n].k)!=
        std::tie(b[n].blockId,b[n].i,b[n].j,b[n].k)) {
      result.message="SWMF snapshots have different canonical cell keys";
      return result;
    }
    double ea[3],eb[3];
    DeriveElectricField(a[n].plasmaVelocity_m_s,a[n].magneticField_T,ea);
    DeriveElectricField(b[n].plasmaVelocity_m_s,b[n].magneticField_T,eb);
    for (int d=0;d<3;++d) {
      const double positionError=std::fabs(a[n].position_m[d]-b[n].position_m[d]);
      result.maxAbsPosition_m=std::max(result.maxAbsPosition_m,positionError);
      const double bRel=RelativeDifference(
          a[n].magneticField_T[d],b[n].magneticField_T[d],bFloor);
      const double uRel=RelativeDifference(
          a[n].plasmaVelocity_m_s[d],b[n].plasmaVelocity_m_s[d],uFloor);
      const double eRel=RelativeDifference(ea[d],eb[d],eFloor);
      result.maxRelativeB=std::max(result.maxRelativeB,bRel);
      result.maxRelativeVelocity=std::max(result.maxRelativeVelocity,uRel);
      result.maxRelativeElectric=std::max(result.maxRelativeElectric,eRel);
      if (positionError>tolerance.positionAbsolute_m) ++result.failedComponents;
      if (std::fabs(a[n].magneticField_T[d]-b[n].magneticField_T[d])>
          tolerance.magneticFieldAbsolute_T+
          tolerance.relative*std::max(std::fabs(a[n].magneticField_T[d]),
                                      std::fabs(b[n].magneticField_T[d])))
        ++result.failedComponents;
      if (std::fabs(a[n].plasmaVelocity_m_s[d]-b[n].plasmaVelocity_m_s[d])>
          tolerance.plasmaVelocityAbsolute_m_s+
          tolerance.relative*std::max(std::fabs(a[n].plasmaVelocity_m_s[d]),
                                      std::fabs(b[n].plasmaVelocity_m_s[d])))
        ++result.failedComponents;
      if (std::fabs(ea[d]-eb[d])>
          tolerance.electricFieldAbsolute_V_m+
          tolerance.relative*std::max(std::fabs(ea[d]),std::fabs(eb[d])))
        ++result.failedComponents;
    }
  }
  result.passed=(result.failedComponents==0);
  result.message=result.passed ? "snapshots agree within tolerance" :
                                 "snapshot component tolerance exceeded";
  return result;
}

// Compatibility overload for small exact/reference tests.  Production validation
// should use ComparisonTolerances so every absolute threshold has an explicit unit.
inline Comparison Compare(const Snapshot& left,const Snapshot& right,
                          double absoluteTolerance,double relativeTolerance) {
  ComparisonTolerances tolerance;
  tolerance.positionAbsolute_m=absoluteTolerance;
  tolerance.magneticFieldAbsolute_T=absoluteTolerance;
  tolerance.plasmaVelocityAbsolute_m_s=absoluteTolerance;
  tolerance.electricFieldAbsolute_V_m=absoluteTolerance;
  tolerance.relative=relativeTolerance;
  return Compare(left,right,tolerance);
}

// Fail-closed status JSON.  The live driver writes FAILED before doing any expensive
// work and overwrites it with PASS only after all requested products and identity gates
// complete.  A fatal exit therefore leaves an explicit non-passing artifact.
inline std::string BuildProductStatusJson(const std::string& status,
                                          const std::string& snapshotId,
                                          const std::string& epochUTC,
                                          const std::string& outputSuffix,
                                          const std::string& exportedSnapshot,
                                          bool cutoffRequested,
                                          bool fluxSpectrumRequested,
                                          const std::string& message) {
  if (status!="PASS" && status!="FAILED")
    throw std::invalid_argument("SWMF product status must be PASS or FAILED");
  std::ostringstream out;
  out << "{\n"
      << "  \"schema\": \"sep-in-geospace/swmf-product-status/v1\",\n"
      << "  \"status\": \"" << status << "\",\n"
      << "  \"snapshot_id\": \"" << Detail::JsonEscape(snapshotId) << "\",\n"
      << "  \"epoch_utc\": \"" << Detail::JsonEscape(epochUTC) << "\",\n"
      << "  \"output_suffix\": \"" << Detail::JsonEscape(outputSuffix) << "\",\n"
      << "  \"exported_snapshot\": \"" << Detail::JsonEscape(exportedSnapshot) << "\",\n"
      << "  \"cutoff_requested\": " << (cutoffRequested ? "true" : "false") << ",\n"
      << "  \"flux_spectrum_requested\": "
      << (fluxSpectrumRequested ? "true" : "false") << ",\n"
      << "  \"message\": \"" << Detail::JsonEscape(message) << "\"\n"
      << "}\n";
  return out.str();
}

inline void WriteProductStatus(const std::string& fileName,
                               const std::string& json) {
  std::ofstream output(fileName.c_str(),std::ios::out|std::ios::trunc);
  if (!output) throw std::runtime_error("Cannot write SWMF product status: "+fileName);
  output << json;
  output.flush();
  if (!output) throw std::runtime_error("Failed while writing SWMF product status: "+fileName);
}

} // namespace SWMFSnapshot
} // namespace Earth

#endif // _SRC_EARTH_UTIL_SWMF_SNAPSHOT_CONTRACT_H_
