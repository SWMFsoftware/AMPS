#include "sep_transactional_output.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace SEP {
namespace Output {
namespace {

Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,text);
}

Transport::Status EnsureDirectories(const std::string& path) {
  if (path.empty()) return Error("output root is empty");
  std::string current;
  std::size_t start=0;
  if (path[0]=='/') { current="/"; start=1; }
  while (start<=path.size()) {
    const std::size_t slash=path.find('/',start);
    const std::string component=path.substr(start,slash-start);
    if (!component.empty()) {
      if (current.size()>1 && current[current.size()-1]!='/') current+='/' ;
      current+=component;
      if (::mkdir(current.c_str(),0775)!=0 && errno!=EEXIST)
        return Error("cannot create output directory '"+current+"': "+
                     std::strerror(errno));
      struct stat info;
      if (::stat(current.c_str(),&info)!=0 || !S_ISDIR(info.st_mode))
        return Error("output path component is not a directory: "+current);
    }
    if (slash==std::string::npos) break;
    start=slash+1;
  }
  return Transport::Status::Ok();
}

Transport::Status WriteSyncedFile(const std::string& path,
                                  const std::string& bytes) {
  const int descriptor=::open(path.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0664);
  if (descriptor<0) return Error("cannot open temporary output '"+path+"': "+
                                 std::strerror(errno));
  std::size_t offset=0;
  while (offset<bytes.size()) {
    const ssize_t count=::write(descriptor,bytes.data()+offset,bytes.size()-offset);
    if (count<=0) {
      const std::string message="cannot write temporary output '"+path+"': "+
          std::strerror(errno);
      ::close(descriptor); ::unlink(path.c_str()); return Error(message);
    }
    offset+=static_cast<std::size_t>(count);
  }
  if (::fsync(descriptor)!=0 || ::close(descriptor)!=0) {
    const std::string message="cannot flush temporary output '"+path+"'";
    ::unlink(path.c_str()); return Error(message);
  }
  return Transport::Status::Ok();
}

std::string Parent(const std::string& path) {
  const std::size_t slash=path.find_last_of('/');
  return slash==std::string::npos ? std::string(".") : path.substr(0,slash);
}

}  // namespace

Transport::Status ValidateRelativePath(const std::string& path) {
  if (path.empty() || path[0]=='/' || path.find('\0')!=std::string::npos)
    return Error("output artifact path must be non-empty and relative");
  std::size_t start=0;
  while (start<=path.size()) {
    const std::size_t slash=path.find('/',start);
    const std::string component=path.substr(start,slash-start);
    if (component.empty() || component=="." || component=="..")
      return Error("output artifact path contains an unsafe component");
    if (slash==std::string::npos) break;
    start=slash+1;
  }
  return Transport::Status::Ok();
}

Transport::Status EnsureOutputDirectory(const std::string& root,
                                        const std::string& relative) {
  const Transport::Status valid=ValidateRelativePath(relative);
  if (!valid.ok() || root.empty()) return valid.ok() ? Error("output root is empty") : valid;
  return EnsureDirectories(root+(root[root.size()-1]=='/' ? "" : "/")+relative);
}

std::string Checksum64(const std::string& payload) {
  // FNV-1a is used as an accidental-corruption checksum, not as a cryptographic
  // signature.  Its specified byte order makes sidecars portable across ranks.
  std::uint64_t hash=UINT64_C(14695981039346656037);
  for (std::size_t i=0;i<payload.size();++i) {
    hash^=static_cast<unsigned char>(payload[i]);
    hash*=UINT64_C(1099511628211);
  }
  std::ostringstream out;
  out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;
  return out.str();
}

WriteResult WriteTransactional(const std::string& root,
                               const std::string& relative,
                               const std::string& payload,
                               const ArtifactMetadata& metadata) {
  WriteResult result;
  result.status=ValidateRelativePath(relative);
  if (!result.status.ok() || root.empty() || metadata.schemaVersion.empty() ||
      metadata.configurationFingerprint.empty()) {
    if (result.status.ok()) result.status=Error("output metadata/root is incomplete");
    return result;
  }
  result.finalPath=root+(root[root.size()-1]=='/' ? "" : "/")+relative;
  result.status=EnsureDirectories(Parent(result.finalPath));
  if (!result.status.ok()) return result;
  struct stat existing;
  if (::lstat(result.finalPath.c_str(),&existing)==0 ||
      ::lstat((result.finalPath+".manifest").c_str(),&existing)==0) {
    result.status=Error("refusing to overwrite an existing output artifact");
    return result;
  }
  if (errno!=ENOENT) {
    result.status=Error("cannot inspect output destination: "+
                        std::string(std::strerror(errno)));
    return result;
  }
  const std::string temporary=result.finalPath+".tmp";
  result.status=WriteSyncedFile(temporary,payload);
  if (!result.status.ok()) return result;
  if (::rename(temporary.c_str(),result.finalPath.c_str())!=0) {
    ::unlink(temporary.c_str());
    result.status=Error("cannot commit output artifact: "+
                        std::string(std::strerror(errno)));
    return result;
  }
  result.checksum=Checksum64(payload);
  std::ostringstream manifest;
  manifest<<"SEP_ARTIFACT 1\n"
      <<"schema "<<metadata.schemaVersion<<'\n'
      <<"records "<<metadata.recordCount<<'\n'
      <<"bytes "<<payload.size()<<'\n'
      <<"checksum "<<result.checksum<<'\n'
      <<"configuration "<<metadata.configurationFingerprint<<'\n'
      <<"complete true\n";
  const std::string manifestFinal=result.finalPath+".manifest";
  const std::string manifestTemporary=manifestFinal+".tmp";
  result.status=WriteSyncedFile(manifestTemporary,manifest.str());
  if (!result.status.ok()) return result;
  if (::rename(manifestTemporary.c_str(),manifestFinal.c_str())!=0) {
    ::unlink(manifestTemporary.c_str());
    result.status=Error("cannot commit output manifest: "+
                        std::string(std::strerror(errno)));
    return result;
  }
  result.status=Transport::Status::Ok();
  return result;
}

Transport::Status ValidateCompletedArtifact(const std::string& path,
                                            std::string* payload) {
  if (!payload) return Error("artifact payload destination is null");
  std::ifstream data(path.c_str(),std::ios::binary);
  std::ifstream manifest((path+".manifest").c_str());
  if (!data || !manifest) return Error("artifact or completion manifest is missing");
  std::ostringstream bytes; bytes<<data.rdbuf(); *payload=bytes.str();
  std::string magic,key,schema,checksum,configuration,complete;
  int version=0; std::uint64_t records=0,size=0;
  if (!(manifest>>magic>>version) || magic!="SEP_ARTIFACT" || version!=1 ||
      !(manifest>>key>>schema) || key!="schema" ||
      !(manifest>>key>>records) || key!="records" ||
      !(manifest>>key>>size) || key!="bytes" ||
      !(manifest>>key>>checksum) || key!="checksum" ||
      !(manifest>>key>>configuration) || key!="configuration" ||
      !(manifest>>key>>complete) || key!="complete" || complete!="true" ||
      size!=payload->size() || checksum!=Checksum64(*payload))
    return Error("artifact completion manifest is corrupt or mismatched");
  return Transport::Status::Ok();
}

}  // namespace Output
}  // namespace SEP
