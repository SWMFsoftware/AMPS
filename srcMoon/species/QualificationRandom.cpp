#include "QualificationRandom.h"
#include "../../src/models/chemical_database/DataIO.h"
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
namespace Moon { namespace Species { namespace QualificationRandom {
namespace {
namespace D=AMPS::ChemicalDatabase::ReactionData;
bool Reject(std::string* error,const std::string& message) {if(error)*error=message;return false;}
void U64(std::uint64_t value,std::string* bytes) {
  for(int shift=56;shift>=0;shift-=8)bytes->push_back(static_cast<char>((value>>shift)&255U));
}
bool Decode(const std::string& hex,std::string* bytes,std::string* error) {
  if(!D::DigestSyntax(hex))return Reject(error,"maintained SHA256 returned a malformed digest");
  std::string value;value.reserve(32);
  for(std::size_t i=0;i<64;i+=2) {
    auto digit=[](char c){return c<='9'?c-'0':c-'a'+10;};
    value.push_back(static_cast<char>((digit(hex[i])<<4)|digit(hex[i+1])));
  }
  *bytes=std::move(value);return true;
}
}
bool DigestHex(const Key& key,Purpose purpose,std::uint64_t index,std::string* output,std::string* error) {
  if(!output)return Reject(error,"qualification RNG digest output is null");
  const auto code=static_cast<std::uint64_t>(purpose);
  if(key.cohort.empty()||!D::ValidUtf8(key.cohort)||key.cohort.find('\0')!=std::string::npos||code<1||code>4)
    return Reject(error,"qualification RNG requires nonempty non-NUL UTF8 cohort and declared purpose");
  if(key.cohort.size()>std::numeric_limits<std::uint64_t>::max()||
      key.genealogy.size()>std::numeric_limits<std::uint64_t>::max())
    return Reject(error,"qualification RNG semantic key length is not representable");
  try {
    // Standard HMAC inner/outer construction (RFC2104); the seed is a public
    // reproducibility key, NOT a secret authentication credential. The stream
    // only claims a deterministic finite-grid PRF/statistical construction.
    std::string seed;U64(key.campaignSeed,&seed);seed.resize(64,'\0');
    std::string inner=seed,outer=seed;
    for(std::size_t i=0;i<64;++i) {inner[i]^=0x36;outer[i]^=0x5c;}
    std::string message="AMPS-MOON-M3-HMAC-SHA256-COUNTER-v1";
    U64(key.cohort.size(),&message);message+=key.cohort;
    U64(key.originalParent,&message);U64(key.genealogy.size(),&message);
    for(const auto ordinal:key.genealogy)U64(ordinal,&message);
    U64(code,&message);U64(index,&message);
    std::string digest,raw;
    if(!D::BytesSha256(inner+message,&digest,error)||!Decode(digest,&raw,error)||
        !D::BytesSha256(outer+raw,&digest,error))return false;
    *output=std::move(digest);return true;
  } catch(const std::bad_alloc&) {
    return Reject(error,"qualification RNG key allocation failed; reduce genealogy/output concurrency or provide memory; no draw published");
  } catch(const std::length_error&) {
    return Reject(error,"qualification RNG key exceeds string capacity; reduce genealogy/key size; no draw published");
  }
}
bool OpenUniform(const Key& key,Purpose purpose,std::uint64_t index,double* output,std::string* error) {
  if(!output)return Reject(error,"qualification RNG uniform output is null");
  std::string digest;if(!DigestHex(key,purpose,index,&digest,error))return false;
  std::uint64_t bits=0;
  // Thirteen hex digits contain exactly52 bits. Both k+1/2 and multiplication
  // by2^-52 are exact binary64 operations, including the two extreme bins.
  for(std::size_t i=0;i<13;++i)bits=(bits<<4)|static_cast<std::uint64_t>(digest[i]<='9'?digest[i]-'0':digest[i]-'a'+10);
  const double value=std::ldexp(static_cast<double>(bits)+0.5,-52);
  if(!std::isfinite(value)||!(value>0&&value<1))return Reject(error,"qualification RNG finite-grid conversion violated strict open interval");
  *output=value;return true;
}
} } }
