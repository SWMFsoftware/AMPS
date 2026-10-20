#include "Provenance.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <ctime>

#include <omp.h>

namespace {
namespace fs=std::filesystem;

//Small self-contained SHA-256 implementation.  Keeping hashing local avoids a
//shell command (and its quoting/injection problems) and permits streaming very
//large meshes and kernel sets with fixed memory.
class SHA256 {
public:
  SHA256() {Reset();}

  void Update(const void* source,size_t count) {
    const unsigned char* data=(const unsigned char*)source;
    Total+=count;
    while (count>0) {
      size_t take=std::min(count,sizeof(Buffer)-Buffered);
      memcpy(Buffer+Buffered,data,take);
      Buffered+=take,data+=take,count-=take;
      if (Buffered==sizeof(Buffer)) Transform(Buffer),Buffered=0;
    }
  }

  std::string Final() {
    unsigned long long bits=Total*8;
    unsigned char one=0x80,zero=0;
    UpdatePadding(&one,1);
    while (Buffered!=56) UpdatePadding(&zero,1);
    unsigned char length[8];
    for (int i=0;i<8;i++) length[7-i]=(unsigned char)(bits>>(8*i));
    UpdatePadding(length,8);
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned int value : State) out << std::setw(8) << value;
    return out.str();
  }

private:
  unsigned int State[8];
  unsigned char Buffer[64];
  size_t Buffered;
  unsigned long long Total;

  static unsigned int R(unsigned int x,int n) {return (x>>n)|(x<<(32-n));}
  static unsigned int Ch(unsigned int x,unsigned int y,unsigned int z) {return (x&y)^(~x&z);}
  static unsigned int Maj(unsigned int x,unsigned int y,unsigned int z) {return (x&y)^(x&z)^(y&z);}

  void Reset() {
    const unsigned int initial[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
        0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    memcpy(State,initial,sizeof(State));
    Buffered=0,Total=0;
  }

  void UpdatePadding(const void* source,size_t count) {
    const unsigned char* data=(const unsigned char*)source;
    while (count>0) {
      size_t take=std::min(count,sizeof(Buffer)-Buffered);
      memcpy(Buffer+Buffered,data,take);
      Buffered+=take,data+=take,count-=take;
      if (Buffered==sizeof(Buffer)) Transform(Buffer),Buffered=0;
    }
  }

  void Transform(const unsigned char* block) {
    static const unsigned int K[64]={
      0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
      0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
      0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
      0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
      0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
      0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
      0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
      0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
    unsigned int w[64];
    for (int i=0;i<16;i++) w[i]=((unsigned int)block[4*i]<<24)|((unsigned int)block[4*i+1]<<16)|
        ((unsigned int)block[4*i+2]<<8)|block[4*i+3];
    for (int i=16;i<64;i++) {
      unsigned int s0=R(w[i-15],7)^R(w[i-15],18)^(w[i-15]>>3);
      unsigned int s1=R(w[i-2],17)^R(w[i-2],19)^(w[i-2]>>10);
      w[i]=w[i-16]+s0+w[i-7]+s1;
    }
    unsigned int a=State[0],b=State[1],c=State[2],d=State[3],e=State[4],f=State[5],g=State[6],h=State[7];
    for (int i=0;i<64;i++) {
      unsigned int s1=R(e,6)^R(e,11)^R(e,25),t1=h+s1+Ch(e,f,g)+K[i]+w[i];
      unsigned int s0=R(a,2)^R(a,13)^R(a,22),t2=s0+Maj(a,b,c);
      h=g,g=f,f=e,e=d+t1,d=c,c=b,b=a,a=t1+t2;
    }
    State[0]+=a,State[1]+=b,State[2]+=c,State[3]+=d;
    State[4]+=e,State[5]+=f,State[6]+=g,State[7]+=h;
  }
};

std::string HashFile(const fs::path& path) {
  std::ifstream input(path,std::ios::binary);
  if (!input) throw std::runtime_error("cannot hash missing/unreadable file: "+path.string());
  SHA256 hash;
  //The buffer is deliberately heap-backed: an 8 MiB automatic array can
  //exhaust the default process stack precisely when a run is being finalized.
  std::vector<char> buffer(8*1024*1024);
  while (input) {
    input.read(buffer.data(),buffer.size());
    std::streamsize n=input.gcount();
    if (n>0) hash.Update(buffer.data(),(size_t)n);
  }
  if (!input.eof()) throw std::runtime_error("read error while hashing: "+path.string());
  return hash.Final();
}

std::string Escape(const std::string& value) {
  std::ostringstream out;
  for (unsigned char c : value) {
    if (c=='"') out << "\\\"";
    else if (c=='\\') out << "\\\\";
    else if (c=='\n') out << "\\n";
    else if (c=='\r') out << "\\r";
    else if (c=='\t') out << "\\t";
    else if (c<0x20) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)c << std::dec;
    else out << c;
  }
  return out.str();
}

void FileJSON(std::ostream& out,const fs::path& path) {
  fs::path absolute=fs::absolute(path).lexically_normal();
  out << "{\"path\":\"" << Escape(absolute.string()) << "\",\"size_bytes\":"
      << fs::file_size(absolute) << ",\"sha256\":\"" << HashFile(absolute) << "\"}";
}

template<class T> bool Read(std::ifstream& in,T& value) {
  return (bool)in.read((char*)&value,sizeof(value));
}

bool ReadTablePoints(const fs::path& path,std::vector<int>& points,std::string& error) {
  std::ifstream input(path);
  std::string line;
  std::set<int> unique;
  while (std::getline(input,line)) if (!line.empty()&&line[0]!='#') {
    std::istringstream row(line);
    int point;
    if (!(row>>point)) {error="invalid data row in "+path.string(); return false;}
    if (!unique.insert(point).second) {error="duplicate point in "+path.string(); return false;}
    points.push_back(point);
  }
  if (!input.eof()||points.empty()) {error="empty/unreadable table "+path.string(); return false;}
  return true;
}
}

std::string ProvenanceUTCNow() {
  auto now=std::chrono::system_clock::now();
  std::time_t value=std::chrono::system_clock::to_time_t(now);
  std::tm utc;
  gmtime_r(&value,&utc);
  char text[32];
  strftime(text,sizeof(text),"%Y-%m-%dT%H:%M:%SZ",&utc);
  return text;
}

std::string ProvenanceFileSHA256(const std::string& filename) {
  return HashFile(filename);
}

cRunOutputAudit AuditRunOutputs(const std::string& outputDirectory,
    const std::vector<int>& expectedPoints,bool writeKernels,int kernelVersion) {
  cRunOutputAudit result;
  result.KernelVersion=writeKernels ? kernelVersion : 0;
  try {
    fs::path root(outputDirectory),table=root/"COPS_model.dat";
    std::vector<int> tablePoints;
    if (!ReadTablePoints(table,tablePoints,result.Error)) return result;
    if (tablePoints!=expectedPoints) {result.Error="COPS_model.dat point IDs do not match the requested sequence"; return result;}
    result.nPoints=(int)tablePoints.size();
    if (!writeKernels) {result.CacheVerified=true; return result;}

    fs::path elements=root/"elements.bin";
    std::ifstream elementInput(elements,std::ios::binary);
    int version,nfaces,nbasis;
    if (!Read(elementInput,version)||!Read(elementInput,nfaces)||!Read(elementInput,nbasis)||
        version!=1||nfaces<=0||nbasis!=25||fs::file_size(elements)!=(uintmax_t)(12LL+272LL*nfaces)) {
      result.Error="invalid elements.bin after writing";
      return result;
    }
    result.nFaces=nfaces;
    std::set<std::string> expectedNames;
    SHA256 setHash;
    for (int point : expectedPoints) {
      char name[64];
      snprintf(name,sizeof(name),"point=%04d.bin",point);
      expectedNames.insert(name);
      fs::path path=root/"kernels"/name;
      std::ifstream kernel(path,std::ios::binary);
      int savedVersion,savedPoint,nRecords,nSpecies=2;
      if (!Read(kernel,savedVersion)||!Read(kernel,savedPoint)||savedVersion!=kernelVersion||savedPoint!=point) {
        result.Error="invalid kernel prefix: "+path.string(); return result;
      }
      if (kernelVersion==1) {
        double header[10];
        if (!kernel.read((char*)header,sizeof(header))||!Read(kernel,nRecords)) {
          result.Error="truncated kernel header: "+path.string(); return result;
        }
      }
      else if (kernelVersion==2) {
        double header[13];
        if (!kernel.read((char*)header,sizeof(header))||!Read(kernel,nSpecies)||!Read(kernel,nRecords)||nSpecies!=2) {
          result.Error="invalid exact-kernel header: "+path.string(); return result;
        }
      }
      else {result.Error="unsupported kernel version in post-write audit"; return result;}
      if (nRecords<0||nRecords>nfaces) {result.Error="invalid kernel record count: "+path.string(); return result;}
      const long long headerBytes=(kernelVersion==1) ? 92 : 120;
      const long long recordBytes=(kernelVersion==1) ? 40 : 56;
      if (fs::file_size(path)!=(uintmax_t)(headerBytes+recordBytes*nRecords)) {
        result.Error="kernel length mismatch: "+path.string(); return result;
      }
      std::set<int> faces;
      for (int record=0;record<nRecords;record++) {
        int face,shadow;
        if (!Read(kernel,face)||face<0||face>=nfaces||!faces.insert(face).second) {
          result.Error="invalid/duplicate kernel face: "+path.string(); return result;
        }
        double values[6];
        int nDouble=(kernelVersion==1) ? 4 : 6;
        if (!kernel.read((char*)values,nDouble*sizeof(double))||!Read(kernel,shadow)||(shadow!=0&&shadow!=1)) {
          result.Error="invalid kernel record: "+path.string(); return result;
        }
      }
      std::string digest=HashFile(path);
      setHash.Update(name,strlen(name));
      setHash.Update(digest.data(),digest.size());
      result.KernelRecords+=nRecords;
      result.KernelBytes+=(long long)fs::file_size(path);
    }
    for (const auto& entry : fs::directory_iterator(root/"kernels"))
      if (entry.is_regular_file()&&entry.path().filename().string().rfind("point=",0)==0&&
          !expectedNames.count(entry.path().filename().string())) {
        result.Error="unexpected stale kernel: "+entry.path().string(); return result;
      }
    result.KernelSetSHA256=setHash.Final();
    result.CacheVerified=true;
  }
  catch (const std::exception& exception) {result.Error=exception.what();}
  return result;
}

void WriteRunManifest(const std::string& filename,const std::string& createdUTC,
    const std::string& updatedUTC,const std::string& processState,
    const std::string& runKind,const std::string& inputFile,
    const std::string& executable,const std::string& cspiceVersion,
    const cInput& input,
    const std::vector<cProvenanceFile>& dependencies,
    const std::string& legacyActivityIdentity,
    const ActivityEnhancement::Configuration& enhancement,
    const std::string& outputDirectory,const cRunOutputAudit& audit,
    const std::string& error) {
  fs::path target(filename),temporary=target.string()+".tmp";
  std::ofstream out(temporary);
  if (!out) throw std::runtime_error("cannot write run manifest: "+temporary.string());
  out << "{\n  \"schema\": \"liouville_cops.run_manifest\",\n  \"schema_version\": 1,\n"
      << "  \"created_utc\": \"" << Escape(createdUTC) << "\",\n"
      << "  \"updated_utc\": \"" << Escape(updatedUTC) << "\",\n"
      << "  \"run_kind\": \"" << Escape(runKind) << "\",\n"
      << "  \"states\": {\"process\": \"" << Escape(processState)
      << "\", \"cache\": \"" << (audit.KernelVersion==0 ? "not_applicable" : (audit.CacheVerified ? "verified" : "not_verified"))
      << "\", \"optimizer\": \"not_applicable\", \"cpp_prediction_agreement\": \"not_applicable\", "
         "\"scientific_review\": \"not_reviewed\"},\n";
  out << "  \"command\": [\"" << Escape(executable) << "\", \""
      << Escape(fs::absolute(inputFile).string()) << "\"],\n"
      << "  \"output_directory\": \"" << Escape(fs::absolute(outputDirectory).string()) << "\",\n";
  out << "  \"input_file\": "; FileJSON(out,inputFile); out << ",\n";
  out << "  \"software\": {\"executable\": "; FileJSON(out,executable);
  const unsigned int endianTest=1;
  out << ", \"compiler\": \"" << Escape(__VERSION__)
      << "\", \"cspice_toolkit\": \"" << Escape(cspiceVersion)
      << "\", \"openmp_max_threads\": " << omp_get_max_threads()
      << ", \"cache_platform\": {\"endianness\": \""
      << ((*(const unsigned char*)&endianTest)==1 ? "little" : "big")
      << "\", \"sizeof_int\": " << sizeof(int) << ", \"sizeof_double\": " << sizeof(double) << "}},\n";
  out << "  \"dependencies\": [\n";
  for (size_t i=0;i<dependencies.size();i++) {
    out << "    {\"role\": \"" << Escape(dependencies[i].Role) << "\", \"file\": ";
    FileJSON(out,dependencies[i].Path);
    out << "}" << ((i+1<dependencies.size()) ? "," : "") << "\n";
  }
  out << "  ],\n  \"effective_settings\": {\n";
  size_t keyIndex=0;
  for (const auto& item : input.EffectiveValues()) {
    out << "    \"" << Escape(item.first) << "\": [";
    for (size_t i=0;i<item.second.size();i++)
      out << (i ? ", " : "") << "\"" << Escape(item.second[i]) << "\"";
    out << "]" << ((++keyIndex<input.EffectiveValues().size()) ? "," : "") << "\n";
  }
  out << "  },\n  \"activity_enhancement\": {\n"
      << "    \"enabled\": " << (enhancement.Enabled() ? "true" : "false") << ",\n"
      << "    \"model\": \"" << Escape(enhancement.Model) << "\",\n"
      << "    \"legacy_activity_identity\": \"" << Escape(legacyActivityIdentity) << "\"";
  if (enhancement.Enabled()) {
    out << ",\n    \"file_format\": \"" << Escape(enhancement.Format) << "\",\n"
        << "    \"file_format_version\": " << enhancement.FormatVersion << ",\n"
        << "    \"basis_id\": \"" << Escape(enhancement.BasisId) << "\",\n"
        << "    \"species\": \"" << Escape(enhancement.Species) << "\",\n"
        << "    \"body_frame\": \"" << Escape(enhancement.BodyFrame) << "\",\n"
        << "    \"angle_convention\": \"" << Escape(enhancement.AngleDefinition) << "\",\n"
        << "    \"ordering\": \"" << Escape(enhancement.CoefficientOrdering) << "\",\n"
        << "    \"clipping\": \"" << Escape(enhancement.Clipping) << "\",\n"
        << "    \"minimum_degree\": " << enhancement.MinimumDegree << ",\n"
        << "    \"maximum_degree\": " << enhancement.MaximumDegree << ",\n"
        << "    \"coefficient_count\": " << enhancement.Coefficients.size() << ",\n"
        << "    \"facet_residual_count\": " << enhancement.FacetResiduals.size() << ",\n"
        << "    \"scale\": " << std::setprecision(17) << enhancement.Scale << ",\n"
        << "    \"scale_interpretation\": \"" << Escape(enhancement.ScaleMeaning) << "\",\n"
        << "    \"regularization\": \"" << Escape(enhancement.Regularization) << "\",\n"
        << "    \"lambda\": " << std::setprecision(17) << enhancement.Lambda << ",\n"
        << "    \"coefficient_file\": ";
    FileJSON(out,enhancement.FilePath);
  }
  out << "\n  },\n  \"outputs\": {\n";
  fs::path root(outputDirectory),model=root/"COPS_model.dat",elements=root/"elements.bin",snapshot=root/"input.snapshot.in",
      enhancementFacets=root/"ActivityEnhancementFacets.dat";
  if (fs::is_regular_file(model)) {out << "    \"model_table\": "; FileJSON(out,model);}
  else out << "    \"model_table\": null";
  if (fs::is_regular_file(elements)) {out << ",\n    \"elements\": "; FileJSON(out,elements);}
  if (fs::is_regular_file(snapshot)) {out << ",\n    \"input_snapshot\": "; FileJSON(out,snapshot);}
  if (fs::is_regular_file(enhancementFacets)) {out << ",\n    \"activity_enhancement_facets\": "; FileJSON(out,enhancementFacets);}
  out << ",\n    \"kernels\": {\"version\": " << audit.KernelVersion << ", \"count\": " << audit.nPoints
      << ", \"records\": " << audit.KernelRecords << ", \"size_bytes\": " << audit.KernelBytes
      << ", \"set_sha256\": \"" << Escape(audit.KernelSetSHA256) << "\"}\n  }";
  if (!audit.Error.empty()) out << ",\n  \"cache_error\": \"" << Escape(audit.Error) << "\"";
  if (!error.empty()) out << ",\n  \"error\": \"" << Escape(error) << "\"";
  out << "\n}\n";
  out.close();
  if (!out) throw std::runtime_error("failed while writing run manifest: "+temporary.string());
  fs::rename(temporary,target);
}
