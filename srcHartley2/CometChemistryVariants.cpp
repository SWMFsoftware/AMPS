// Rate tables of alternative reaction definitions; see CometChemistryVariants.h
// and srcHartley2/CHEMISTRY_VARIANTS.md.

#include "CometChemistryVariants.h"

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {

// ---------------------------------------------------------------- SHA-256
// Straightforward FIPS 180-4 implementation; only used for provenance of
// small input files, so clarity is preferred over speed.
const std::uint32_t K256[64]={
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};

inline std::uint32_t RotateRight(std::uint32_t x,int n) {
  return (x>>n)|(x<<(32-n));
}

bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

std::string Trim(const std::string &text) {
  const char *space=" \t\r\n";
  const std::size_t first=text.find_first_not_of(space);
  if (first==std::string::npos) return std::string();
  const std::size_t last=text.find_last_not_of(space);
  return text.substr(first,last-first+1);
}

// Strict decimal parse: the whole field must be one finite number.
bool ParseNumber(const std::string &field,double *value) {
  const std::string text=Trim(field);
  if (text.empty()) return false;
  char *end=NULL;
  errno=0;
  const double parsed=std::strtod(text.c_str(),&end);
  if (errno!=0 || end==text.c_str() || *end!='\0' || !std::isfinite(parsed))
    return false;
  *value=parsed;
  return true;
}

// Split "a,b" into exactly two fields.
bool SplitPair(const std::string &line,std::string *first,std::string *second) {
  const std::size_t comma=line.find(',');
  if (comma==std::string::npos || line.find(',',comma+1)!=std::string::npos)
    return false;
  *first=Trim(line.substr(0,comma));
  *second=Trim(line.substr(comma+1));
  return true;
}

}  // namespace

std::string Comet::ChemistryVariants::Sha256Hex(const std::string &bytes) {
  std::uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
                      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  std::string message=bytes;
  const std::uint64_t bitLength=static_cast<std::uint64_t>(bytes.size())*8u;
  message.push_back(static_cast<char>(0x80));
  while (message.size()%64!=56) message.push_back('\0');
  for (int shift=56;shift>=0;shift-=8)
    message.push_back(static_cast<char>((bitLength>>shift)&0xffu));

  for (std::size_t chunk=0;chunk<message.size();chunk+=64) {
    std::uint32_t w[64];
    for (int i=0;i<16;i++) {
      const unsigned char *p=
          reinterpret_cast<const unsigned char*>(message.data()+chunk+4*i);
      w[i]=(static_cast<std::uint32_t>(p[0])<<24)|
           (static_cast<std::uint32_t>(p[1])<<16)|
           (static_cast<std::uint32_t>(p[2])<<8)|
            static_cast<std::uint32_t>(p[3]);
    }
    for (int i=16;i<64;i++) {
      const std::uint32_t s0=RotateRight(w[i-15],7)^RotateRight(w[i-15],18)^(w[i-15]>>3);
      const std::uint32_t s1=RotateRight(w[i-2],17)^RotateRight(w[i-2],19)^(w[i-2]>>10);
      w[i]=w[i-16]+s0+w[i-7]+s1;
    }
    std::uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],k=h[7];
    for (int i=0;i<64;i++) {
      const std::uint32_t S1=RotateRight(e,6)^RotateRight(e,11)^RotateRight(e,25);
      const std::uint32_t ch=(e&f)^((~e)&g);
      const std::uint32_t t1=k+S1+ch+K256[i]+w[i];
      const std::uint32_t S0=RotateRight(a,2)^RotateRight(a,13)^RotateRight(a,22);
      const std::uint32_t maj=(a&b)^(a&c)^(b&c);
      const std::uint32_t t2=S0+maj;
      k=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=k;
  }

  char hex[65];
  for (int i=0;i<8;i++) std::snprintf(hex+8*i,9,"%08x",h[i]);
  return std::string(hex,64);
}

bool Comet::ChemistryVariants::ParseLifetimeTable(const std::string &text,
    const std::string &name,cLifetimeTable *table,std::string *error) {
  if (table==NULL) return Fail("internal error: no lifetime table output",error);
  cLifetimeTable result;
  std::istringstream stream(text);
  std::string raw;
  int lineNumber=0;
  int record=0;  // 0: expect tau_fuv_s, 1: expect header, 2: rows
  while (std::getline(stream,raw)) {
    lineNumber++;
    const std::string line=Trim(raw);
    if (line.empty() || line[0]=='#') continue;
    const std::string where=name+":"+std::to_string(lineNumber)+": ";
    std::string first,second;
    if (!SplitPair(line,&first,&second))
      return Fail(where+"expected two comma-separated fields, found '"+line+"'",error);

    if (record==0) {
      if (first!="tau_fuv_s")
        return Fail(where+"the first record must be 'tau_fuv_s,<seconds>'",error);
      if (!ParseNumber(second,&result.TauFuvS) || !(result.TauFuvS>0.0))
        return Fail(where+"tau_fuv_s must be a finite positive number of seconds",error);
      record=1;
    }
    else if (record==1) {
      if (first!="v_h_km_s" || second!="ohtaus_1e5_s")
        return Fail(where+"expected the column header 'v_h_km_s,ohtaus_1e5_s'",error);
      record=2;
    }
    else {
      double velocity=0.0,value=0.0;
      if (!ParseNumber(first,&velocity))
        return Fail(where+"v_h_km_s is not a finite number",error);
      if (!ParseNumber(second,&value) || !(value>0.0))
        return Fail(where+"ohtaus_1e5_s must be a finite positive number",error);
      if (!result.VelocityKmS.empty() && !(velocity>result.VelocityKmS.back()))
        return Fail(where+"v_h_km_s must increase strictly from row to row",error);
      result.VelocityKmS.push_back(velocity);
      result.Ohtaus.push_back(value);
    }
  }
  if (record<2)
    return Fail(name+": missing the tau_fuv_s record or the column header",error);
  if (result.VelocityKmS.size()<2)
    return Fail(name+": the table needs at least two velocity rows",error);

  result.Path=table->Path;
  result.Sha256=table->Sha256;
  *table=result;
  return true;
}

bool Comet::ChemistryVariants::ReadLifetimeTable(const std::string &path,
    cLifetimeTable *table,std::string *error) {
  if (table==NULL) return Fail("internal error: no lifetime table output",error);
  std::ifstream file(path,std::ios::in|std::ios::binary);
  if (!file) return Fail("cannot open the rate table file '"+path+"'",error);
  std::ostringstream content;
  content<<file.rdbuf();
  if (file.bad()) return Fail("cannot read the rate table file '"+path+"'",error);
  const std::string bytes=content.str();

  cLifetimeTable parsed;
  parsed.Path=path;
  parsed.Sha256=Sha256Hex(bytes);
  if (!ParseLifetimeTable(bytes,path,&parsed,error)) return false;
  *table=parsed;
  return true;
}

Comet::ChemistryVariants::cInterpolation
Comet::ChemistryVariants::InterpolateOhtaus(const cLifetimeTable &table,
    double velocityKmS) {
  cInterpolation result;
  const std::vector<double> &v=table.VelocityKmS;
  const std::vector<double> &y=table.Ohtaus;
  if (velocityKmS<=v.front()) {
    result.Value=y.front();
    result.Clamped=velocityKmS<v.front();
    return result;
  }
  if (velocityKmS>=v.back()) {
    result.Value=y.back();
    result.Clamped=velocityKmS>v.back();
    return result;
  }
  std::size_t i=0;
  while (i+1<v.size() && !(velocityKmS<v[i+1])) i++;
  const double fraction=(velocityKmS-v[i])/(v[i+1]-v[i]);
  result.Value=y[i]+fraction*(y[i+1]-y[i]);
  return result;
}

bool Comet::ChemistryVariants::EvaluateSA1982(const cLifetimeTable &table,
    double velocityKmS,double scaleS,double ionFraction,cSA1982Rates *rates,
    std::string *error) {
  if (rates==NULL) return Fail("internal error: no SA1982 rate output",error);
  if (table.VelocityKmS.size()<2 || !(table.TauFuvS>0.0))
    return Fail("the OH lifetime table is empty or incomplete",error);
  if (!std::isfinite(velocityKmS))
    return Fail("the heliocentric radial velocity is not finite",error);
  if (!(scaleS>0.0) || !std::isfinite(scaleS))
    return Fail("the OH lifetime scale S must be finite and positive",error);
  if (!(ionFraction>=0.0) || !std::isfinite(ionFraction))
    return Fail("the OH ionization fraction must be finite and non-negative",error);

  const cInterpolation ohtaus=InterpolateOhtaus(table,velocityKmS);
  cSA1982Rates result;
  result.VelocityKmS=velocityKmS;
  result.Ohtaus=ohtaus.Value;
  result.Clamped=ohtaus.Clamped;
  result.Dissociation1AU=1.0/table.TauFuvS+1.0/(scaleS*ohtaus.Value);
  result.Ionization1AU=ionFraction*result.Dissociation1AU;
  *rates=result;
  return true;
}
