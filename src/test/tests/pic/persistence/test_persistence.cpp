#include "test_harness.h"
#include "pic_callback_output_kernels.h"
#include "pic_deposition_kernels.h"
#include "pic_persistent_format.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace P=PIC::PersistentFormat;
static const std::uint64_t layout=P::LayoutIdentity({9,8,7,6});
static std::string path(const char* name) { return std::string(name); }
static void append_u64(std::vector<unsigned char>& out,std::uint64_t value) {
  const std::size_t offset=out.size(); out.resize(offset+8); P::Put64(out.data()+offset,value);
}
static void append_double(std::vector<unsigned char>& out,double value) {
  std::uint64_t bits; std::memcpy(&bits,&value,8); append_u64(out,bits);
}
static std::uint64_t take_u64(const std::vector<unsigned char>& in,std::size_t& offset) {
  if (offset+8>in.size()) throw AMPS::Test::Failure("truncated canonical payload");
  const auto value=P::Get64(in.data()+offset); offset+=8; return value;
}
static double take_double(const std::vector<unsigned char>& in,std::size_t& offset) {
  const std::uint64_t bits=take_u64(in,offset); double value; std::memcpy(&value,&bits,8); return value;
}
static P::Status read(const std::string& file,P::Kind kind,std::uint64_t expected_layout,
                      std::vector<unsigned char>& payload,bool legacy=false) {
  FILE* stream=nullptr; std::uint64_t size=0;
  const P::Status status=P::OpenValidated(file,kind,expected_layout,legacy,stream,size);
  if (status==P::Status::Success || status==P::Status::LegacyAccepted) {
    payload.resize(size);
    if (size && std::fread(payload.data(),1,size,stream)!=size) {
      std::fclose(stream); return P::Status::IoError;
    }
    std::fclose(stream);
  }
  return status;
}

struct State {
  std::uint64_t mesh_id,owner,sample_count,rng;
  double x,v,weight,dt,internal_energy,chemistry,boundary;
  bool operator==(const State& other) const {
    return mesh_id==other.mesh_id && owner==other.owner && sample_count==other.sample_count && rng==other.rng &&
      x==other.x && v==other.v && weight==other.weight && dt==other.dt && internal_energy==other.internal_energy &&
      chemistry==other.chemistry && boundary==other.boundary;
  }
};
static std::vector<unsigned char> encode(const State& s) {
  std::vector<unsigned char> out; append_u64(out,s.mesh_id); append_u64(out,s.owner); append_u64(out,s.sample_count); append_u64(out,s.rng);
  append_double(out,s.x); append_double(out,s.v); append_double(out,s.weight); append_double(out,s.dt);
  append_double(out,s.internal_energy); append_double(out,s.chemistry); append_double(out,s.boundary); return out;
}
static State decode(const std::vector<unsigned char>& in) {
  std::size_t o=0; State s; s.mesh_id=take_u64(in,o); s.owner=take_u64(in,o); s.sample_count=take_u64(in,o); s.rng=take_u64(in,o);
  s.x=take_double(in,o); s.v=take_double(in,o); s.weight=take_double(in,o); s.dt=take_double(in,o);
  s.internal_energy=take_double(in,o); s.chemistry=take_double(in,o); s.boundary=take_double(in,o);
  if (o!=in.size()) throw AMPS::Test::Failure("trailing canonical state bytes");
  return s;
}
static void step(State& s) {
  s.x+=s.v*s.dt; s.internal_energy+=0.125*s.dt; s.chemistry*=std::exp(-0.03*s.dt);
  s.boundary+=s.x>=1.0?1.0:0.0; ++s.sample_count; s.rng=6364136223846793005ull*s.rng+1442695040888963407ull;
}

AMPS_TEST_CASE("PIC-PERSISTENCE-INVENTORY-001","persistent-format inventory records metadata and compatibility policy") {
  std::ifstream input(AMPS_PERSISTENT_INVENTORY); AMPS_EXPECT_TRUE(input.good());
  const std::string text((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
  for (const char* token:{"sampling_restart","particle_restart","amr_mesh","background_binary","trajectory_dataset","tecplot_ascii","compatibility_policy"})
    AMPS_EXPECT_TRUE(text.find(token)!=std::string::npos);
}

AMPS_TEST_CASE("PIC-RESTART-ROUNDTRIP-001","canonical restart state round-trips topology, owner, particles, sampling, RNG, chemistry and boundary state") {
  const State expected{42,3,19,0x123456789abcdef0ull,0.25,-7.0,8.5,0.02,11.0,0.75,-2.0};
  const std::string file=path("restart-roundtrip.bin");
  AMPS_EXPECT_TRUE(P::WriteAtomic(file,P::Kind::TestFixture,P::CanonicalLittleEndian,layout,encode(expected))==P::Status::Success);
  std::vector<unsigned char> payload; AMPS_EXPECT_TRUE(read(file,P::Kind::TestFixture,layout,payload)==P::Status::Success);
  AMPS_EXPECT_TRUE(decode(payload)==expected); std::remove(file.c_str());
}

AMPS_TEST_CASE("PIC-RESTART-CONTINUATION-001","one step after restart is bitwise equivalent to uninterrupted deterministic continuation") {
  State uninterrupted{8,0,0,7,0.1,2.0,1.0,0.05,4.0,0.9,0.0}; step(uninterrupted);
  const std::string file=path("restart-continuation.bin");
  AMPS_EXPECT_TRUE(P::WriteAtomic(file,P::Kind::TestFixture,P::CanonicalLittleEndian,layout,encode(uninterrupted))==P::Status::Success);
  std::vector<unsigned char> payload; AMPS_EXPECT_TRUE(read(file,P::Kind::TestFixture,layout,payload)==P::Status::Success);
  State restored=decode(payload); step(uninterrupted); step(restored); AMPS_EXPECT_TRUE(restored==uninterrupted); std::remove(file.c_str());
}

AMPS_TEST_CASE("PIC-RESTART-EMPTY-MAX-001","empty and maximum bounded records round-trip without sentinel ambiguity") {
  for (const std::size_t size:{std::size_t(0),std::size_t(1<<20)}) {
    std::vector<unsigned char> source(size); for (std::size_t i=0;i<size;++i) source[i]=static_cast<unsigned char>((i*37u+11u)&255u);
    const std::string file=size?path("restart-max.bin"):path("restart-empty.bin");
    AMPS_EXPECT_TRUE(P::WriteAtomic(file,P::Kind::TestFixture,P::CanonicalLittleEndian,layout,source)==P::Status::Success);
    std::vector<unsigned char> copy; AMPS_EXPECT_TRUE(read(file,P::Kind::TestFixture,layout,copy)==P::Status::Success); AMPS_EXPECT_TRUE(copy==source); std::remove(file.c_str());
  }
}

static std::vector<unsigned char> bytes(const std::string& file) {
  std::ifstream in(file,std::ios::binary); return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}
static void put_bytes(const std::string& file,const std::vector<unsigned char>& data) {
  std::ofstream out(file,std::ios::binary|std::ios::trunc); out.write(reinterpret_cast<const char*>(data.data()),data.size());
}

AMPS_TEST_CASE("PIC-RESTART-CORRUPTION-001","truncation, payload corruption, unknown version, wrong kind and layout are rejected") {
  const std::string base=path("restart-valid.bin"); const std::vector<unsigned char> payload={1,2,3,4,5,6,7,8};
  AMPS_EXPECT_TRUE(P::WriteAtomic(base,P::Kind::TestFixture,P::NativeLayout,layout,payload)==P::Status::Success);
  auto data=bytes(base); std::vector<unsigned char> result;
  auto truncated=data; truncated.pop_back(); put_bytes("restart-truncated.bin",truncated);
  AMPS_EXPECT_TRUE(read("restart-truncated.bin",P::Kind::TestFixture,layout,result)==P::Status::Truncated);
  auto corrupt=data; corrupt.back()^=0x80; put_bytes("restart-corrupt.bin",corrupt);
  AMPS_EXPECT_TRUE(read("restart-corrupt.bin",P::Kind::TestFixture,layout,result)==P::Status::CorruptPayload);
  auto version=data; P::Put16(version.data()+8,99); P::Put32(version.data()+52,P::Crc32(version.data(),52)); put_bytes("restart-version.bin",version);
  AMPS_EXPECT_TRUE(read("restart-version.bin",P::Kind::TestFixture,layout,result)==P::Status::UnknownVersion);
  AMPS_EXPECT_TRUE(read(base,P::Kind::ParticleRestart,layout,result)==P::Status::WrongKind);
  AMPS_EXPECT_TRUE(read(base,P::Kind::TestFixture,layout+1,result)==P::Status::IncompatibleLayout);
  for (const char* f:{"restart-valid.bin","restart-truncated.bin","restart-corrupt.bin","restart-version.bin"}) std::remove(f);
}

AMPS_TEST_CASE("PIC-RESTART-ENDIAN-WIDTH-001","native payload rejects endian/width mismatch while canonical payload is width-independent") {
  std::vector<unsigned char> result; const std::vector<unsigned char> payload={9,8,7};
  P::WriteAtomic("native-compat.bin",P::Kind::TestFixture,P::NativeLayout,layout,payload); auto native=bytes("native-compat.bin");
  native[28]=native[28]==1?2:1; P::Put32(native.data()+52,P::Crc32(native.data(),52)); put_bytes("native-endian.bin",native);
  AMPS_EXPECT_TRUE(read("native-endian.bin",P::Kind::TestFixture,layout,result)==P::Status::WrongEndian);
  native=bytes("native-compat.bin"); native[22]=native[22]+1; P::Put32(native.data()+52,P::Crc32(native.data(),52)); put_bytes("native-width.bin",native);
  AMPS_EXPECT_TRUE(read("native-width.bin",P::Kind::TestFixture,layout,result)==P::Status::IncompatibleWidths);
  native=bytes("native-compat.bin"); native[20]=native[20]==4?8:4; P::Put32(native.data()+52,P::Crc32(native.data(),52)); put_bytes("native-pointer-width.bin",native);
  AMPS_EXPECT_TRUE(read("native-pointer-width.bin",P::Kind::TestFixture,layout,result)==P::Status::IncompatibleWidths);
  P::WriteAtomic("canonical.bin",P::Kind::TestFixture,P::CanonicalLittleEndian,layout,payload); auto canonical=bytes("canonical.bin");
  canonical[20]=4; canonical[22]=16; canonical[23]=16; P::Put32(canonical.data()+52,P::Crc32(canonical.data(),52)); put_bytes("canonical-other-build.bin",canonical);
  AMPS_EXPECT_TRUE(read("canonical-other-build.bin",P::Kind::TestFixture,layout,result)==P::Status::Success); AMPS_EXPECT_TRUE(result==payload);
  for (const char* f:{"native-compat.bin","native-endian.bin","native-width.bin","native-pointer-width.bin","canonical.bin","canonical-other-build.bin"}) std::remove(f);
}

AMPS_TEST_CASE("PIC-RESTART-ATOMIC-001","interrupted temporary publication leaves the prior complete restart readable") {
  const std::string file="atomic-restart.bin"; const std::vector<unsigned char> first={1,1,1},second={2,2,2,2};
  AMPS_EXPECT_TRUE(P::WriteAtomic(file,P::Kind::TestFixture,P::CanonicalLittleEndian,layout,first)==P::Status::Success);
  put_bytes(file+".amps-publish-interrupted.tmp",{0,1,2}); std::vector<unsigned char> value;
  AMPS_EXPECT_TRUE(read(file,P::Kind::TestFixture,layout,value)==P::Status::Success); AMPS_EXPECT_TRUE(value==first);
  AMPS_EXPECT_TRUE(P::WriteAtomic(file,P::Kind::TestFixture,P::CanonicalLittleEndian,layout,second)==P::Status::Success);
  AMPS_EXPECT_TRUE(read(file,P::Kind::TestFixture,layout,value)==P::Status::Success); AMPS_EXPECT_TRUE(value==second);
  std::remove(file.c_str()); std::remove((file+".amps-publish-interrupted.tmp").c_str());
}

AMPS_TEST_CASE("PIC-RESTART-LEGACY-POLICY-001","unversioned files require an explicit same-build legacy opt-in") {
  put_bytes("legacy.bin",{4,3,2,1}); std::vector<unsigned char> out;
  AMPS_EXPECT_TRUE(read("legacy.bin",P::Kind::TestFixture,layout,out,false)==P::Status::Truncated);
  AMPS_EXPECT_TRUE(read("legacy.bin",P::Kind::TestFixture,layout,out,true)==P::Status::LegacyAccepted); AMPS_EXPECT_TRUE(out==std::vector<unsigned char>({4,3,2,1})); std::remove("legacy.bin");
}

AMPS_TEST_CASE("PIC-DATAFILE-BACKGROUND-001","versioned background fixture round-trips and affine interpolation is exact") {
  std::vector<unsigned char> source; for (int c=0;c<8;++c) append_double(source,2.0+3.0*(c&1)-4.0*((c>>1)&1)+5.0*((c>>2)&1));
  P::WriteAtomic("background.bin",P::Kind::BackgroundData,P::CanonicalLittleEndian,layout,source); std::vector<unsigned char> data;
  AMPS_EXPECT_TRUE(read("background.bin",P::Kind::BackgroundData,layout,data)==P::Status::Success);
  const double xi[3]={0.2,0.65,0.4}; double w[8],value=0; PIC::DepositionKernels::TrilinearWeights(xi,w); std::size_t offset=0;
  for (int c=0;c<8;++c) value+=w[c]*take_double(data,offset);
  AMPS_EXPECT_NEAR(value,2+3*xi[0]-4*xi[1]+5*xi[2],2e-15); std::remove("background.bin");
}

AMPS_TEST_CASE("PIC-DATAFILE-DECOMPOSITION-001","concatenated rank chunks can be consumed under a different read decomposition") {
  put_bytes("rank0.part",{0,1,2}); put_bytes("rank1.part",{3,4,5,6,7});
  AMPS_EXPECT_TRUE(P::WriteAtomicConcat("rank-data.bin",P::Kind::BackgroundData,P::CanonicalLittleEndian,layout,{"rank0.part","rank1.part"})==P::Status::Success);
  std::vector<unsigned char> data; AMPS_EXPECT_TRUE(read("rank-data.bin",P::Kind::BackgroundData,layout,data)==P::Status::Success);
  std::vector<unsigned char> new_rank0(data.begin(),data.begin()+4),new_rank1(data.begin()+4,data.end());
  AMPS_EXPECT_TRUE(new_rank0==std::vector<unsigned char>({0,1,2,3})); AMPS_EXPECT_TRUE(new_rank1==std::vector<unsigned char>({4,5,6,7}));
  for (const char* f:{"rank0.part","rank1.part","rank-data.bin"}) std::remove(f);
}

AMPS_TEST_CASE("PIC-DATAFILE-COMPATIBILITY-001","data-file kind and layout identity prevent accidental cross-format reads") {
  P::WriteAtomic("typed-data.bin",P::Kind::BackgroundData,P::CanonicalLittleEndian,layout,{1,2}); std::vector<unsigned char> out;
  AMPS_EXPECT_TRUE(read("typed-data.bin",P::Kind::AssociatedNodeData,layout,out)==P::Status::WrongKind);
  AMPS_EXPECT_TRUE(read("typed-data.bin",P::Kind::BackgroundData,0,out)==P::Status::IncompatibleLayout); std::remove("typed-data.bin");
}

AMPS_TEST_CASE("PIC-OUTPUT-TECPLOT-LAYOUT-001","trajectory Tecplot variables retain order, units, dimensions, and finite records") {
  AMPS_EXPECT_TRUE(std::string(PIC::CallbackOutputKernels::TrajectoryBaseVariables())==
    "VARIABLES=\"x\", \"y\", \"z\", \"spec\", \"Speed\", \"vx\", \"vy\", \"vz\", \"Kinetic Energy [eV]\"");
  FILE* f=std::tmpfile(); const double x[3]={1,2,3},v[3]={4,5,6};
  AMPS_EXPECT_TRUE(PIC::CallbackOutputKernels::WriteTrajectoryBaseRecord(f,x,0,std::sqrt(77.0),v,12.0)); std::fflush(f); std::rewind(f);
  double values[8]; int spec; AMPS_EXPECT_EQ(std::fscanf(f,"%lf%lf%lf%d%lf%lf%lf%lf%lf",&values[0],&values[1],&values[2],&spec,&values[3],&values[4],&values[5],&values[6],&values[7]),9);
  for (double value:values) AMPS_EXPECT_TRUE(std::isfinite(value)); AMPS_EXPECT_EQ(spec,0); std::fclose(f);
}

AMPS_TEST_CASE("PIC-OUTPUT-NONFINITE-001","nonfinite output values are rejected before publication") {
  FILE* f=std::tmpfile(); const double x[3]={0,0,std::numeric_limits<double>::infinity()},v[3]={0,0,0};
  AMPS_EXPECT_TRUE(!PIC::CallbackOutputKernels::WriteTrajectoryBaseRecord(f,x,0,0,v,0)); std::fclose(f);
}
