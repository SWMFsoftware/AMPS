#include "DataIO.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace AMPS { namespace ChemicalDatabase {
namespace {
bool Fail(std::string* error,const std::string& message) {if(error)*error=message;return false;}
// SHA-256 (FIPS 180-4), kept private and independent of PIC/OpenSSL. The
// incremental implementation bounds memory for publication-file identity
// checks. Local tests include published empty/"abc" digests and Python's
// independent hashlib over the actual native sources, not a self-oracle.
class Sha256 {
 public:
  void Add(const unsigned char* data,std::size_t length) {
    bytes_+=length;
    while (length) {
      const std::size_t n=std::min(length,64-used_);
      std::copy(data,data+n,buffer_.begin()+used_);
      used_+=n;data+=n;length-=n;
      if (used_==64) {Block(buffer_.data());used_=0;}
    }
  }
  std::string Finish() {
    const std::uint64_t bitLength=bytes_*8;
    const unsigned char marker=0x80,zero=0;
    Add(&marker,1);
    while (used_!=56) Add(&zero,1);
    unsigned char tail[8];
    for (unsigned i=0;i<8;++i) tail[7-i]=static_cast<unsigned char>(bitLength>>(8*i));
    Add(tail,8);
    std::ostringstream output;output<<std::hex<<std::setfill('0');
    for (auto value:state_) output<<std::setw(8)<<value;
    return output.str();
  }
 private:
  static std::uint32_t Rotate(std::uint32_t v,unsigned n) {return (v>>n)|(v<<(32-n));}
  void Block(const unsigned char* b) {
    static const std::uint32_t k[64]={
      0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
      0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
      0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
      0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
      0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
      0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
      0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
      0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    std::uint32_t w[64];
    for (unsigned i=0;i<16;++i)
      w[i]=(std::uint32_t(b[4*i])<<24)|(std::uint32_t(b[4*i+1])<<16)|
          (std::uint32_t(b[4*i+2])<<8)|b[4*i+3];
    for (unsigned i=16;i<64;++i) {
      const auto s0=Rotate(w[i-15],7)^Rotate(w[i-15],18)^(w[i-15]>>3);
      const auto s1=Rotate(w[i-2],17)^Rotate(w[i-2],19)^(w[i-2]>>10);
      w[i]=w[i-16]+s0+w[i-7]+s1;
    }
    auto a=state_[0],c=state_[2],b0=state_[1],d=state_[3];
    auto e=state_[4],f=state_[5],g=state_[6],h=state_[7];
    for (unsigned i=0;i<64;++i) {
      const auto s1=Rotate(e,6)^Rotate(e,11)^Rotate(e,25);
      const auto t1=h+s1+((e&f)^((~e)&g))+k[i]+w[i];
      const auto s0=Rotate(a,2)^Rotate(a,13)^Rotate(a,22);
      const auto t2=s0+((a&b0)^(a&c)^(b0&c));
      h=g;g=f;f=e;e=d+t1;d=c;c=b0;b0=a;a=t1+t2;
    }
    state_[0]+=a;state_[1]+=b0;state_[2]+=c;state_[3]+=d;
    state_[4]+=e;state_[5]+=f;state_[6]+=g;state_[7]+=h;
  }
  std::array<std::uint32_t,8> state_{{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19}};
  std::array<unsigned char,64> buffer_{};
  std::size_t used_=0;
  std::uint64_t bytes_=0;
};

} // namespace
namespace ReactionData {
bool DigestSyntax(const std::string& digest) {
  return digest.size()==64 && std::all_of(digest.begin(),digest.end(),[](char c) {
    return (c>='0' && c<='9') || (c>='a' && c<='f');});
}

bool ParseCsv(const std::string& bytes,const Row& header,Csv* rows,std::string* error) {
  // Parse quoted commas/escaped quotes/newlines without locale or shell
  // interpretation. Garbage following a closing quote is rejected, rather
  // than silently absorbed into a source reference or numerical value.
  Csv candidate;Row row;std::string field;
  bool quoted=false,closed=false,started=false;
  const auto finishField=[&]() {row.push_back(field);field.clear();closed=false;started=false;};
  const auto finishRow=[&]() {finishField();candidate.push_back(row);row.clear();};
  for (std::size_t i=0;i<bytes.size();++i) {
    const char c=bytes[i];
    if (c=='\0') return Fail(error,"D04 CSV contains a NUL byte");
    if (quoted) {
      if (c=='"') {
        if (i+1<bytes.size() && bytes[i+1]=='"') {field+='"';++i;}
        else {quoted=false;closed=true;}
      }
      else field+=c;
      continue;
    }
    if (c==',') {finishField();continue;}
    if (c=='\n' || c=='\r') {
      if (c=='\r') {
        if (i+1==bytes.size() || bytes[i+1]!='\n') return Fail(error,"bare CR in D04 CSV");
        ++i;
      }
      finishRow();continue;
    }
    if (closed) return Fail(error,"characters after closing CSV quote");
    if (c=='"') {
      if (started || !field.empty()) return Fail(error,"quote inside unquoted CSV field");
      quoted=true;started=true;
    }
    else {field+=c;started=true;}
  }
  if (quoted) return Fail(error,"unterminated D04 CSV quote");
  if (!field.empty() || started || closed || !row.empty()) finishRow();
  if (candidate.empty() || candidate.front()!=header) return Fail(error,"D04 CSV exact schema/header mismatch");
  for (std::size_t i=1;i<candidate.size();++i)
    if (candidate[i].size()!=header.size()) return Fail(error,"D04 CSV row has wrong field count");
  candidate.erase(candidate.begin());*rows=std::move(candidate);return true;
}

bool ReadPinnedCsv(const std::string& path,const std::string& expected,std::size_t size,
    const Row& header,Csv* rows,std::string* error) {
  // v1 has fixed published identities/sizes. Bound allocation before reading,
  // then hash the SAME bytes that are parsed, avoiding a hash/reopen race.
  std::ifstream file(path,std::ios::binary);
  if (!file) return Fail(error,"cannot open frozen D04 source '"+path+"'");
  file.seekg(0,std::ios::end);
  if (file.tellg()!=static_cast<std::streamoff>(size)) return Fail(error,"frozen D04 source byte count changed: "+path);
  file.seekg(0,std::ios::beg);std::string bytes(size,'\0');
  file.read(&bytes[0],static_cast<std::streamsize>(size));
  if (!file || file.peek()!=std::char_traits<char>::eof()) return Fail(error,"D04 source changed or failed during read: "+path);
  Sha256 hash;hash.Add(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
  if (hash.Finish()!=expected) return Fail(error,"frozen D04 SHA-256 mismatch: "+path);
  return ParseCsv(bytes,header,rows,error);
}

bool FileSha256(const std::string& path,std::string* digest,std::string* error) {
  if (!digest || !error) return false;
  error->clear();std::ifstream file(path,std::ios::binary);
  if (!file) return Fail(error,"cannot open source for SHA-256: "+path);
  Sha256 hash;std::array<unsigned char,65536> bytes;
  while (file) {
    file.read(reinterpret_cast<char*>(bytes.data()),bytes.size());
    const auto n=file.gcount();
    if (n>0) hash.Add(bytes.data(),static_cast<std::size_t>(n));
  }
  if (!file.eof() || file.bad()) return Fail(error,"I/O failure hashing source: "+path);
  *digest=hash.Finish();return true;
}
bool VerifyFileSha256(const std::string& path,const std::string& expected,std::string* error) {
  if (!error) return false;
  if (!DigestSyntax(expected)) return Fail(error,"source SHA-256 must be exactly64 lowercase hexadecimal characters");
  std::string actual;
  if (!FileSha256(path,&actual,error)) return false;
  return actual==expected || Fail(error,"source SHA-256 mismatch: "+path);
}
bool BytesSha256(const std::string& bytes,std::string* digest,std::string* error) {
  if(!digest || !error)return false;
  Sha256 hash;hash.Add(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
  *digest=hash.Finish();error->clear();return true;
}
bool ReadFileBytes(const std::string& path,std::size_t ceiling,std::string* bytes,std::string* error) {
  if(!bytes || !error)return false;
  std::ifstream file(path,std::ios::binary|std::ios::ate);
  if(!file)return Fail(error,"cannot open runtime source: "+path);
  const auto length=file.tellg();
  if(length<0 || static_cast<unsigned long long>(length)>std::numeric_limits<std::size_t>::max() ||
      length>std::numeric_limits<std::streamsize>::max())return Fail(error,"unrepresentable runtime source byte count: "+path);
  const auto size=static_cast<std::size_t>(length);
  if(ceiling && size>ceiling)return Fail(error,"resource limit: runtime source exceeds explicitly requested byte ceiling: "+path);
  std::string candidate(size,'\0');file.seekg(0,std::ios::beg);
  if(size)file.read(&candidate[0],static_cast<std::streamsize>(size));
  if(!file || file.peek()!=std::char_traits<char>::eof())return Fail(error,"runtime source changed or failed during read: "+path);
  *bytes=std::move(candidate);error->clear();return true;
}
bool ValidUtf8(const std::string& bytes) {
  // Reject overlong encodings, surrogate scalars and values above Unicode's
  // maximum. Bytes are not normalized: hashes and source spelling stay exact.
  for(std::size_t i=0;i<bytes.size();) {
    const unsigned char first=static_cast<unsigned char>(bytes[i++]);
    if(first<0x80)continue;
    unsigned count=0;std::uint32_t value=0,minimum=0;
    if(first>=0xc2 && first<=0xdf){count=1;value=first&0x1f;minimum=0x80;}
    else if(first>=0xe0 && first<=0xef){count=2;value=first&0x0f;minimum=0x800;}
    else if(first>=0xf0 && first<=0xf4){count=3;value=first&7;minimum=0x10000;}
    else return false;
    if(bytes.size()-i<count)return false;
    for(unsigned j=0;j<count;++j) {
      const unsigned char c=static_cast<unsigned char>(bytes[i++]);
      if((c&0xc0)!=0x80)return false;
      value=(value<<6)|(c&0x3f);
    }
    if(value<minimum || value>0x10ffff || (value>=0xd800 && value<=0xdfff))return false;
  }
  return true;
}
} // namespace ReactionData

} } // namespace AMPS::ChemicalDatabase
