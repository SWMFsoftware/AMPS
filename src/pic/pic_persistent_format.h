#ifndef _AMPS_PIC_PERSISTENT_FORMAT_H_
#define _AMPS_PIC_PERSISTENT_FORMAT_H_

/*
 * AMPS persistent-file envelope, version 1.
 *
 * The 64-byte header is always little-endian and contains magic, version,
 * format kind, compatibility flags, writer integer/floating widths, payload
 * byte order, generated-layout identity, payload length, CRC-32, and a header
 * CRC-32. Native-layout payloads are accepted only when byte order and widths
 * match the reader and layout_identity matches the caller's generated layout.
 * Canonical payloads use explicit fixed-width little-endian fields and may be
 * read across builds; the payload codec, not this envelope, owns conversion.
 *
 * Publication writes a sibling temporary file, flushes and fsyncs it, then
 * renames it over the destination. A crash before rename leaves the previous
 * complete destination intact. The parent directory fsync is intentionally
 * not attempted because several supported HPC filesystems reject it; this is
 * atomic-name publication, not a guarantee against total storage loss.
 */

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace PIC { namespace PersistentFormat {

constexpr std::size_t HeaderSize=64;
constexpr unsigned char Magic[8]={'A','M','P','S','P','F','1','\0'};
constexpr std::uint16_t Version=1;
enum class Kind : std::uint16_t {
  SamplingRestart=1,ParticleRestart=2,BackgroundData=3,
  AssociatedNodeData=4,TrajectoryData=5,TestFixture=65535
};
enum Flag : std::uint16_t { NativeLayout=1u,CanonicalLittleEndian=2u };
enum class Status {
  Success,LegacyAccepted,OpenFailed,IoError,BadMagic,Truncated,
  UnknownVersion,WrongKind,WrongEndian,IncompatibleWidths,
  IncompatibleLayout,CorruptPayload,TrailingData,InvalidArgument,
  PublishFailed
};

inline std::uint32_t Crc32Update(std::uint32_t crc,const unsigned char* data,
                                 std::size_t size) {
  crc=~crc;
  for (std::size_t i=0;i<size;++i) {
    crc^=data[i];
    for (int bit=0;bit<8;++bit) crc=(crc>>1)^
      (0xedb88320u & static_cast<std::uint32_t>(-static_cast<int>(crc&1u)));
  }
  return ~crc;
}
inline std::uint32_t Crc32(const unsigned char* data,std::size_t size) {
  return Crc32Update(0u,data,size);
}
inline bool HostLittleEndian() {
  const std::uint16_t value=1;
  return *reinterpret_cast<const unsigned char*>(&value)==1;
}
inline void Put16(unsigned char* p,std::uint16_t v) { p[0]=v&255u; p[1]=(v>>8)&255u; }
inline void Put32(unsigned char* p,std::uint32_t v) { for (int i=0;i<4;++i) p[i]=(v>>(8*i))&255u; }
inline void Put64(unsigned char* p,std::uint64_t v) { for (int i=0;i<8;++i) p[i]=(v>>(8*i))&255u; }
inline std::uint16_t Get16(const unsigned char* p) { return p[0]|(std::uint16_t(p[1])<<8); }
inline std::uint32_t Get32(const unsigned char* p) { std::uint32_t v=0; for(int i=0;i<4;++i)v|=std::uint32_t(p[i])<<(8*i); return v; }
inline std::uint64_t Get64(const unsigned char* p) { std::uint64_t v=0; for(int i=0;i<8;++i)v|=std::uint64_t(p[i])<<(8*i); return v; }

inline std::uint64_t LayoutIdentity(std::initializer_list<std::uint64_t> values) {
  std::uint64_t hash=1469598103934665603ull;
  for (std::uint64_t value:values) for (int byte=0;byte<8;++byte) {
    hash^=(value>>(8*byte))&255u; hash*=1099511628211ull;
  }
  return hash;
}

inline void BuildHeader(unsigned char header[HeaderSize],Kind kind,
                        std::uint16_t flags,std::uint64_t layout,
                        std::uint64_t payload_size,std::uint32_t payload_crc) {
  std::memset(header,0,HeaderSize); std::memcpy(header,Magic,8);
  Put16(header+8,Version); Put16(header+10,HeaderSize);
  Put16(header+12,static_cast<std::uint16_t>(kind)); Put16(header+14,flags);
  Put32(header+16,0x01020304u);
  // Native payloads can contain offsets or records whose generated layout is
  // pointer-width dependent even when no raw pointer is intentionally
  // serialized. Record the ABI pointer width explicitly instead of assuming
  // every supported build is LP64; OpenValidated rejects a 32/64-bit mismatch
  // before any native bytes are interpreted.
  header[20]=sizeof(void*); header[21]=sizeof(short); header[22]=sizeof(int);
  header[23]=sizeof(long); header[24]=sizeof(long long); header[25]=sizeof(std::size_t);
  header[26]=sizeof(float); header[27]=sizeof(double);
  header[28]=HostLittleEndian()?1:2;
  Put64(header+32,layout); Put64(header+40,payload_size); Put32(header+48,payload_crc);
  Put32(header+52,Crc32(header,52));
}

inline Status WriteAtomicConcat(const std::string& destination,Kind kind,
                                std::uint16_t flags,std::uint64_t layout,
                                const std::vector<std::string>& parts) {
  if (destination.empty() || (flags!=NativeLayout && flags!=CanonicalLittleEndian))
    return Status::InvalidArgument;
  const std::string temporary=destination+".amps-publish-"+
      std::to_string(static_cast<long long>(::getpid()))+".tmp";
  FILE* output=std::fopen(temporary.c_str(),"wb+");
  if (!output) return Status::OpenFailed;
  unsigned char header[HeaderSize]={0};
  if (std::fwrite(header,1,HeaderSize,output)!=HeaderSize) {
    std::fclose(output); std::remove(temporary.c_str()); return Status::IoError;
  }
  std::uint64_t total=0; std::uint32_t crc=0; unsigned char buffer[1<<16];
  for (const std::string& part:parts) {
    FILE* input=std::fopen(part.c_str(),"rb");
    if (!input) { std::fclose(output); std::remove(temporary.c_str()); return Status::OpenFailed; }
    while (true) {
      const std::size_t count=std::fread(buffer,1,sizeof(buffer),input);
      if (count) {
        if (std::fwrite(buffer,1,count,output)!=count) { std::fclose(input); std::fclose(output); std::remove(temporary.c_str()); return Status::IoError; }
        crc=Crc32Update(crc,buffer,count); total+=count;
      }
      if (count<sizeof(buffer)) {
        if (std::ferror(input)) { std::fclose(input); std::fclose(output); std::remove(temporary.c_str()); return Status::IoError; }
        break;
      }
    }
    std::fclose(input);
  }
  BuildHeader(header,kind,flags,layout,total,crc);
  if (std::fseek(output,0,SEEK_SET)!=0 || std::fwrite(header,1,HeaderSize,output)!=HeaderSize ||
      std::fflush(output)!=0 || ::fsync(::fileno(output))!=0 || std::fclose(output)!=0) {
    std::remove(temporary.c_str()); return Status::IoError;
  }
  if (std::rename(temporary.c_str(),destination.c_str())!=0) {
    std::remove(temporary.c_str()); return Status::PublishFailed;
  }
  return Status::Success;
}

inline Status WriteAtomic(const std::string& destination,Kind kind,
                          std::uint16_t flags,std::uint64_t layout,
                          const std::vector<unsigned char>& payload) {
  const std::string part=destination+".amps-payload-"+
      std::to_string(static_cast<long long>(::getpid()))+".tmp";
  FILE* file=std::fopen(part.c_str(),"wb");
  if (!file) return Status::OpenFailed;
  const bool ok=payload.empty() || std::fwrite(payload.data(),1,payload.size(),file)==payload.size();
  if (std::fclose(file)!=0 || !ok) { std::remove(part.c_str()); return Status::IoError; }
  const Status status=WriteAtomicConcat(destination,kind,flags,layout,{part});
  std::remove(part.c_str()); return status;
}

inline Status OpenValidated(const std::string& path,Kind expected_kind,
                            std::uint64_t expected_layout,bool allow_legacy,
                            FILE*& stream,std::uint64_t& payload_size) {
  stream=nullptr; payload_size=0;
  FILE* file=std::fopen(path.c_str(),"rb");
  if (!file) return Status::OpenFailed;
  unsigned char header[HeaderSize];
  const std::size_t got=std::fread(header,1,HeaderSize,file);
  if (got<8 || std::memcmp(header,Magic,8)!=0) {
    if (!allow_legacy) { std::fclose(file); return got<8?Status::Truncated:Status::BadMagic; }
    if (std::fseek(file,0,SEEK_END)!=0) { std::fclose(file); return Status::IoError; }
    const long end=std::ftell(file);
    if (end<0 || std::fseek(file,0,SEEK_SET)!=0) { std::fclose(file); return Status::IoError; }
    payload_size=static_cast<std::uint64_t>(end); stream=file; return Status::LegacyAccepted;
  }
  if (got!=HeaderSize) { std::fclose(file); return Status::Truncated; }
  if (Get16(header+8)!=Version) { std::fclose(file); return Status::UnknownVersion; }
  if (Get16(header+10)!=HeaderSize || Get32(header+16)!=0x01020304u ||
      Get32(header+52)!=Crc32(header,52)) { std::fclose(file); return Status::CorruptPayload; }
  if (Get16(header+12)!=static_cast<std::uint16_t>(expected_kind)) { std::fclose(file); return Status::WrongKind; }
  const std::uint16_t flags=Get16(header+14);
  if (flags==NativeLayout) {
    if (header[28]!=(HostLittleEndian()?1:2)) { std::fclose(file); return Status::WrongEndian; }
    if (header[20]!=sizeof(void*) || header[21]!=sizeof(short) || header[22]!=sizeof(int) ||
        header[23]!=sizeof(long) || header[24]!=sizeof(long long) ||
        header[25]!=sizeof(std::size_t) || header[26]!=sizeof(float) ||
        header[27]!=sizeof(double)) { std::fclose(file); return Status::IncompatibleWidths; }
  }
  else if (flags!=CanonicalLittleEndian) { std::fclose(file); return Status::InvalidArgument; }
  if (Get64(header+32)!=expected_layout) { std::fclose(file); return Status::IncompatibleLayout; }
  payload_size=Get64(header+40); const std::uint32_t expected_crc=Get32(header+48);
  std::uint64_t remaining=payload_size; std::uint32_t crc=0; unsigned char buffer[1<<16];
  while (remaining) {
    const std::size_t request=remaining<sizeof(buffer)?static_cast<std::size_t>(remaining):sizeof(buffer);
    const std::size_t count=std::fread(buffer,1,request,file);
    if (count!=request) { std::fclose(file); return Status::Truncated; }
    crc=Crc32Update(crc,buffer,count); remaining-=count;
  }
  if (crc!=expected_crc) { std::fclose(file); return Status::CorruptPayload; }
  if (std::fgetc(file)!=EOF) { std::fclose(file); return Status::TrailingData; }
  if (std::fseek(file,HeaderSize,SEEK_SET)!=0) { std::fclose(file); return Status::IoError; }
  stream=file; return Status::Success;
}

}} // namespace PIC::PersistentFormat
#endif
