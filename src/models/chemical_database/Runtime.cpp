#include "Runtime.h"
#include "../../general/input/amps_input.h"
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace AMPS { namespace ChemicalDatabase {
namespace fs=std::filesystem;
namespace detail {
struct CatalogueStorage {
  BundleIdentity identity;
  std::vector<SourceRecord> sources;
  std::vector<TableRecord> tables;
  std::map<std::string,std::size_t> sourceById,tableById;
  std::map<std::string,std::string> sourcePaths,tablePaths;
  std::map<std::string,std::shared_ptr<const ReactionData::Csv>> tableRows;
};
}
namespace {
Status Error(Code code,const std::string& text){return {code,text};}
bool Has(const std::vector<std::string>& values,const std::string& key) {
  return std::find(values.begin(),values.end(),key)!=values.end();
}
bool Identifier(const std::string& value) {
  return !value.empty() && std::all_of(value.begin(),value.end(),[](unsigned char c) {
    return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='-'||c=='.'||c==':';
  });
}
bool List(const std::string& text,std::vector<std::string>* output) {
  if(text.empty()||text=="none"){output->clear();return true;}
  std::vector<std::string> value;std::set<std::string> seen;std::size_t start=0;
  while(start<=text.size()) {
    const auto end=text.find(';',start);const auto item=text.substr(start,end==std::string::npos?end:end-start);
    if(!Identifier(item)||!seen.insert(item).second)return false;
    value.push_back(item);if(end==std::string::npos)break;start=end+1;
  }
  *output=std::move(value);return true;
}
bool Count(const std::string& text,std::size_t* output) {
  if(text.empty())return false;
  std::size_t value=0;
  for(unsigned char c:text) {
    if(c<'0'||c>'9'||value>(std::numeric_limits<std::size_t>::max()-(c-'0'))/10)return false;
    value=value*10+(c-'0');
  }
  *output=value;return true;
}
bool AvailabilityValue(const std::string& text,Availability* value) {
  if(text=="available")*value=Availability::Available;
  else if(text=="unavailable")*value=Availability::Unavailable;
  else if(text=="metadata-only")*value=Availability::MetadataOnly;
  else if(text=="quarantined")*value=Availability::Quarantined;
  else return false;
  return true;
}
bool PermissionValue(const std::string& text,Permission* value) {
  if(text=="allowed")*value=Permission::Allowed;
  else if(text=="pending")*value=Permission::Pending;
  else if(text=="prohibited")*value=Permission::Prohibited;
  else if(text=="unknown")*value=Permission::Unknown;
  else if(text=="not-applicable")*value=Permission::NotApplicable;
  else return false;
  return true;
}
Status Available(Availability value,const std::string& id) {
  return value==Availability::Available?Status{}:Error(Code::Unavailable,id+": "+AvailabilityName(value)+"; no physical zero or repaired payload is substituted");
}
Status SourcePermission(const SourceRecord& value) {
  if(value.storagePermission!=Permission::Allowed)
    return Error(Code::PermissionRequired,value.id+": numerical managed storage permission is "+PermissionName(value.storagePermission));
  if(value.evaluationPermission!=Permission::Allowed)
    return Error(Code::PermissionRequired,value.id+": evaluation permission is "+PermissionName(value.evaluationPermission));
  return {};
}
Status SourceGate(const SourceRecord& source) {
  // Permission precedes payload availability: metadata-only restricted records
  // explain the actual rights hold instead of requesting a forbidden download.
  auto status=SourcePermission(source);if(!status.ok())return status;
  return Available(source.availability,source.id);
}
Status TableGate(const detail::CatalogueStorage& storage,const TableRecord& table) {
  if(table.evaluationPermission!=Permission::Allowed)
    return Error(Code::PermissionRequired,table.id+": evaluation permission is "+PermissionName(table.evaluationPermission));
  auto status=Available(table.availability,table.id);if(!status.ok())return status;
  for(const auto& id:table.sourceIds) {
    const auto found=storage.sourceById.find(id);
    if(found==storage.sourceById.end())return Error(Code::CorruptSource,table.id+": dangling source ID "+id);
    status=SourceGate(storage.sources[found->second]);if(!status.ok())return status;
  }
  return {};
}
Status RelativePath(const fs::path& root,const std::string& text,fs::path* output) {
  // Reject syntactic traversal even when lexical normalization would put it
  // back inside root. Canonical containment ALSO rejects symlink escapes.
  if(text.empty()||text.find('\\')!=std::string::npos||text.find('\0')!=std::string::npos)
    return Error(Code::InvalidInput,"runtime object needs a nonempty forward-slash relative path");
  const fs::path relative(text);
  if(relative.is_absolute()||relative.has_root_name()||relative.has_root_directory())
    return Error(Code::InvalidInput,"absolute runtime object path refused: "+text);
  for(const auto& part:relative)if(part=="..")return Error(Code::InvalidInput,"runtime path traversal refused: "+text);
  std::error_code ec;const auto candidate=fs::weakly_canonical(root/relative,ec);
  if(ec)return Error(Code::Unavailable,"cannot resolve runtime object: "+text+": "+ec.message());
  auto r=root.begin(),c=candidate.begin();
  for(;r!=root.end();++r,++c)if(c==candidate.end()||*r!=*c)
    return Error(Code::InvalidInput,"runtime object escapes archive root: "+text);
  if(candidate==root)return Error(Code::InvalidInput,"runtime object path names archive root itself");
  *output=candidate;return {};
}
struct Reader {
  fs::path root;ResourceLimits limits;std::size_t total=0;
  std::map<std::string,std::pair<std::string,std::uintmax_t>> seen;
  Status CheckSize(const fs::path& path,const std::string& digest,std::uintmax_t* size) {
    if(!ReactionData::DigestSyntax(digest))return Error(Code::CorruptSource,"runtime object has invalid SHA-256: "+path.string());
    std::error_code ec;
    if(!fs::is_regular_file(path,ec)||ec)return Error(Code::Unavailable,"required runtime object unavailable: "+path.string());
    const auto bytes=fs::file_size(path,ec);
    if(ec)return Error(Code::Unavailable,"cannot stat runtime object: "+path.string());
    if(bytes>std::numeric_limits<std::size_t>::max() || (limits.maximumFileBytes&&bytes>limits.maximumFileBytes))
      return Error(Code::ResourceLimit,"runtime object exceeds representable or explicitly requested file-byte ceiling: "+path.string()+"; reduce selected data or increase the explicit ceiling/resources");
    const auto old=seen.find(path.string());
    if(old!=seen.end()) {
      if(old->second.first!=digest||old->second.second!=bytes)return Error(Code::CorruptSource,"conflicting/changed runtime object identity: "+path.string());
    } else {
      if(bytes>std::numeric_limits<std::size_t>::max()-total ||
          (limits.maximumTotalBytes && bytes>limits.maximumTotalBytes-total))
        return Error(Code::ResourceLimit,"runtime bundle exceeds explicit total-byte ceiling; reduce selected data or increase the explicit ceiling/resources");
      total+=static_cast<std::size_t>(bytes);seen.emplace(path.string(),std::make_pair(digest,bytes));
    }
    *size=bytes;return {};
  }
  Status Read(const std::string& relative,const std::string& digest,std::string* bytes,fs::path* resolved=nullptr) {
    fs::path path;auto status=RelativePath(root,relative,&path);if(!status.ok())return status;
    std::uintmax_t size=0;status=CheckSize(path,digest,&size);if(!status.ok())return status;
    std::string candidate,error,actual;
    if(!ReactionData::ReadFileBytes(path.string(),limits.maximumFileBytes,&candidate,&error))return Error(Code::CorruptSource,error);
    if(candidate.size()!=size||!ReactionData::BytesSha256(candidate,&actual,&error)||actual!=digest)
      return Error(Code::CorruptSource,"runtime source SHA-256/length mismatch: "+relative);
    if(resolved)*resolved=path;
    *bytes=std::move(candidate);return {};
  }
  Status Verify(const std::string& relative,const std::string& digest,fs::path* resolved=nullptr) {
    fs::path path;auto status=RelativePath(root,relative,&path);if(!status.ok())return status;
    std::uintmax_t size=0;status=CheckSize(path,digest,&size);if(!status.ok())return status;
    std::string error;
    if(!ReactionData::VerifyFileSha256(path.string(),digest,&error))return Error(Code::CorruptSource,error);
    if(resolved)*resolved=path;
    return {};
  }
};
bool AnyScope(const std::vector<std::string>& record,const std::vector<std::string>& profile) {
  for(const auto& cap:record)if(Has(profile,cap))return true;
  return false;
}
bool PayloadPair(const std::string& path,const std::string& hash,Availability availability) {
  if(path.empty()!=hash.empty())return false;
  return path.empty()?availability!=Availability::Available:ReactionData::DigestSyntax(hash);
}
Status Csv(const std::string& bytes,const ReactionData::Row& header,ReactionData::Csv* output) {
  std::string error;
  if(!ReactionData::ValidUtf8(bytes))return Error(Code::CorruptSource,"runtime CSV is not valid UTF-8");
  if(!ReactionData::ParseCsv(bytes,header,output,&error))return Error(Code::CorruptSource,error);
  return {};
}
}

const char* CodeName(Code c) {
  switch(c) {
    case Code::Evaluated:return "Evaluated";case Code::Unavailable:return "Unavailable";
    case Code::InvalidInput:return "InvalidInput";case Code::OutsideDomain:return "OutsideDomain";
    case Code::CorruptSource:return "CorruptSource";case Code::UnsupportedSchema:return "UnsupportedSchema";
    case Code::PermissionRequired:return "PermissionRequired";case Code::NumericalFailure:return "NumericalFailure";
    case Code::ResourceLimit:return "ResourceLimit";
  }return "InvalidCode";
}
const char* AvailabilityName(Availability x) {
  switch(x){case Availability::Available:return "available";case Availability::Unavailable:return "unavailable";
    case Availability::MetadataOnly:return "metadata-only";case Availability::Quarantined:return "quarantined";}return "invalid";
}
const char* PermissionName(Permission x) {
  switch(x){case Permission::Allowed:return "allowed";case Permission::Pending:return "pending";
    case Permission::Prohibited:return "prohibited";case Permission::Unknown:return "unknown";
    case Permission::NotApplicable:return "not-applicable";}return "invalid";
}
const BundleIdentity& Catalogue::identity()const {static const BundleIdentity empty;return storage_?storage_->identity:empty;}
const std::vector<SourceRecord>& Catalogue::sources()const {static const std::vector<SourceRecord> empty;return storage_?storage_->sources:empty;}
const std::vector<TableRecord>& Catalogue::tables()const {static const std::vector<TableRecord> empty;return storage_?storage_->tables:empty;}
Status Catalogue::RequireCapability(const std::string& capability)const {
  if(!storage_)return Error(Code::Unavailable,"chemical runtime bundle is not loaded");
  if(!Identifier(capability))return Error(Code::InvalidInput,"capability identifier must be explicit");
  return Has(storage_->identity.profileCapabilities,capability)?Status{}:Error(Code::Unavailable,"profile does not declare capability "+capability);
}
Status Catalogue::GetSource(const std::string& id,SourceRecord* output)const {
  if(!output)return Error(Code::InvalidInput,"source metadata output is null");
  if(!storage_)return Error(Code::Unavailable,"chemical runtime bundle is not loaded");
  const auto found=storage_->sourceById.find(id);if(found==storage_->sourceById.end())return Error(Code::Unavailable,"unknown chemical source "+id);
  try {auto value=storage_->sources[found->second];*output=std::move(value);return {};}
  catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"source metadata allocation failed; reduce selected metadata or provide more memory");}
}
Status Catalogue::GetTable(const std::string& id,TableRecord* output)const {
  if(!output)return Error(Code::InvalidInput,"table metadata output is null");
  if(!storage_)return Error(Code::Unavailable,"chemical runtime bundle is not loaded");
  const auto found=storage_->tableById.find(id);if(found==storage_->tableById.end())return Error(Code::Unavailable,"unknown chemical table "+id);
  try {auto value=storage_->tables[found->second];*output=std::move(value);return {};}
  catch(const std::bad_alloc&){return Error(Code::ResourceLimit,"table metadata allocation failed; reduce selected metadata or provide more memory");}
}
Status Catalogue::SourcePath(const std::string& id,const std::string& capability,std::string* output)const {
  if(!output)return Error(Code::InvalidInput,"source path output is null");
  auto status=RequireCapability(capability);if(!status.ok())return status;
  const auto found=storage_->sourceById.find(id);if(found==storage_->sourceById.end())return Error(Code::Unavailable,"unknown chemical source "+id);
  const auto& record=storage_->sources[found->second];
  if(!Has(record.capabilities,capability))return Error(Code::Unavailable,id+": source does not declare requested capability "+capability);
  status=SourceGate(record);if(!status.ok())return status;
  const auto path=storage_->sourcePaths.find(id);if(path==storage_->sourcePaths.end())return Error(Code::Unavailable,id+": permitted source payload was not prepared");
  *output=path->second;return {};
}
Status Catalogue::TablePath(const std::string& id,const std::string& capability,std::string* output)const {
  if(!output)return Error(Code::InvalidInput,"table path output is null");
  auto status=RequireCapability(capability);if(!status.ok())return status;
  const auto found=storage_->tableById.find(id);if(found==storage_->tableById.end())return Error(Code::Unavailable,"unknown chemical table "+id);
  const auto& record=storage_->tables[found->second];
  if(!Has(record.capabilities,capability))return Error(Code::Unavailable,id+": table does not declare requested capability "+capability);
  status=TableGate(*storage_,record);if(!status.ok())return status;
  const auto path=storage_->tablePaths.find(id);if(path==storage_->tablePaths.end())return Error(Code::Unavailable,id+": permitted table payload was not prepared");
  *output=path->second;return {};
}
Status Catalogue::TableRows(const std::string& id,const std::string& capability,std::shared_ptr<const ReactionData::Csv>* output)const {
  if(!output)return Error(Code::InvalidInput,"table rows output is null");
  std::string path;auto status=TablePath(id,capability,&path);if(!status.ok())return status;
  const auto rows=storage_->tableRows.find(id);if(rows==storage_->tableRows.end())return Error(Code::UnsupportedSchema,id+": only CSV payloads expose parsed rows");
  *output=rows->second;return {};
}

Status OpenRuntimeBundle(const BundleRequest& request,Catalogue* output) {
  if(!output)return Error(Code::InvalidInput,"chemical catalogue output is null");
  if(request.root.empty()||!Identifier(request.expectedReleaseId)||!Identifier(request.expectedProfileId)||
     !ReactionData::DigestSyntax(request.expectedExportSha256)||
     !ReactionData::DigestSyntax(request.expectedReleaseSha256)||
     !ReactionData::DigestSyntax(request.expectedProfileSha256))
    return Error(Code::InvalidInput,"opening chemistry requires explicit root, release/profile IDs and three pinned lowercase SHA-256 digests");
  if(request.expectedCodec!="csv-block-v1")return Error(Code::UnsupportedSchema,"unsupported requested chemical runtime codec");
  try {
    std::error_code ec;const auto root=fs::canonical(fs::path(request.root),ec);
    if(ec||!fs::is_directory(root))return Error(Code::Unavailable,"chemical archive root unavailable: "+request.root);
    Reader reader;reader.root=root;reader.limits=request.limits;
    std::string manifest;auto status=reader.Read(request.exportPath,request.expectedExportSha256,&manifest);
    if(!status.ok())return status;
    if(!ReactionData::ValidUtf8(manifest)||manifest.find('\0')!=std::string::npos)
      return Error(Code::CorruptSource,"chemical runtime manifest is not non-NUL UTF-8");
    // An export is self-contained and already hash-pinned. Reject includes
    // BEFORE the common parser could read an unpinned external file. This is
    // an owner policy, not another implementation of AMPS grammar.
    std::istringstream lines(manifest);std::string line;
    while(std::getline(lines,line)) {
      const auto first=line.find_first_not_of(" \t\r");
      if(first!=std::string::npos&&line[first]=='#'&&
         line.compare(first,6,"#begin")!=0&&line.compare(first,4,"#end")!=0)
        return Error(Code::CorruptSource,"detached chemical runtime manifest forbids includes and other directives");
    }
    Input::Document document;
    const auto parsed=Input::ParseText(manifest,request.exportPath,root.string(),&document);
    if(!parsed.ok())return Error(Code::CorruptSource,parsed.message);
    if(!document.root.assignments.empty()||document.root.blocks.size()!=1||document.files.size()!=1)
      return Error(Code::CorruptSource,"export must contain exactly one self-contained chemical_database: runtime block");
    const auto& block=document.root.blocks.front();
    if(block.kind!="chemical_database"||block.name!="runtime"||!block.blocks.empty())
      return Error(Code::UnsupportedSchema,"runtime export owner/schema block is not chemical_database: runtime");
    const std::set<std::string> keys={"schema","codec","release_id","release_manifest_path","release_manifest_sha256",
      "profile_id","profile_path","profile_sha256","source_index_path","source_index_sha256",
      "table_index_path","table_index_sha256","profile_capabilities","topology_fingerprint"};
    std::map<std::string,std::string> values;
    for(const auto& item:block.assignments) {
      if(!keys.count(item.key))return Error(Code::UnsupportedSchema,"unknown runtime manifest key: "+item.key);
      std::string value=item.value;
      if(value.size()>=2&&((value.front()=='"'&&value.back()=='"')||(value.front()=='\''&&value.back()=='\'')))value=value.substr(1,value.size()-2);
      values.emplace(item.key,std::move(value));
    }
    if(values.size()!=keys.size())return Error(Code::CorruptSource,"runtime manifest is missing a required v1 key");
    if(values["schema"]!="amps-chemical-runtime-v1"||values["codec"]!=request.expectedCodec)
      return Error(Code::UnsupportedSchema,"unsupported chemical export schema/codec");
    if(values["release_id"]!=request.expectedReleaseId||values["release_manifest_sha256"]!=request.expectedReleaseSha256||
       values["profile_id"]!=request.expectedProfileId||values["profile_sha256"]!=request.expectedProfileSha256)
      return Error(Code::CorruptSource,"runtime export does not bind the requested release/profile identities");
    auto storage=std::make_shared<detail::CatalogueStorage>();auto& identity=storage->identity;
    identity.root=root.string();identity.exportPath=request.exportPath;identity.exportSha256=request.expectedExportSha256;
    identity.releaseId=values["release_id"];identity.releaseManifestPath=values["release_manifest_path"];
    identity.releaseSha256=values["release_manifest_sha256"];identity.profileId=values["profile_id"];
    identity.profilePath=values["profile_path"];identity.profileSha256=values["profile_sha256"];
    identity.codec=values["codec"];identity.topologyFingerprint=values["topology_fingerprint"];
    if(identity.topologyFingerprint!="none"&&!ReactionData::DigestSyntax(identity.topologyFingerprint))
      return Error(Code::CorruptSource,"topology fingerprint must be none or lowercase SHA-256");
    if(!List(values["profile_capabilities"],&identity.profileCapabilities))return Error(Code::CorruptSource,"invalid/duplicate profile capability list");
    for(const auto& cap:request.requiredCapabilities) {
      if(!Identifier(cap))return Error(Code::InvalidInput,"invalid required capability identifier");
      if(!Has(identity.profileCapabilities,cap))return Error(Code::Unavailable,"profile does not declare required capability "+cap);
    }
    // Detached export cites FINAL JSON hashes. Those JSON files do not cite
    // this export; their acyclic release/import graph is verified by tools.
    // C++ authenticates their exact bytes without interpreting arbitrary JSON.
    status=reader.Verify(identity.releaseManifestPath,identity.releaseSha256);if(!status.ok())return status;
    status=reader.Verify(identity.profilePath,identity.profileSha256);if(!status.ok())return status;
    std::string sourcesText,tablesText;
    status=reader.Read(values["source_index_path"],values["source_index_sha256"],&sourcesText);if(!status.ok())return status;
    status=reader.Read(values["table_index_path"],values["table_index_sha256"],&tablesText);if(!status.ok())return status;
    ReactionData::Csv sources,tables;
    status=Csv(sourcesText,{"source_id","path","sha256","media_type","availability","storage_permission","evaluation_permission","capabilities","source_uri","notes"},&sources);if(!status.ok())return status;
    status=Csv(tablesText,{"table_id","path","sha256","format","schema_id","source_ids","availability","evaluation_permission","capabilities","columns","row_count","notes"},&tables);if(!status.ok())return status;
    for(const auto& row:sources) {
      SourceRecord source;source.id=row[0];source.path=row[1];source.sha256=row[2];source.mediaType=row[3];source.sourceUri=row[8];source.notes=row[9];
      if(!Identifier(source.id)||source.mediaType.empty()||
          !AvailabilityValue(row[4],&source.availability)||!PermissionValue(row[5],&source.storagePermission)||
          !PermissionValue(row[6],&source.evaluationPermission)||!List(row[7],&source.capabilities)||
          !PayloadPair(source.path,source.sha256,source.availability)||storage->sourceById.count(source.id))
        return Error(Code::CorruptSource,"invalid/duplicate source index record: "+source.id);
      if(!source.path.empty()) {fs::path resolved;status=RelativePath(root,source.path,&resolved);if(!status.ok())return status;}
      storage->sourceById.emplace(source.id,storage->sources.size());storage->sources.push_back(std::move(source));
    }
    for(const auto& row:tables) {
      TableRecord table;table.id=row[0];table.path=row[1];table.sha256=row[2];table.format=row[3];table.schemaId=row[4];table.notes=row[11];
      if(!Identifier(table.id)||!Identifier(table.schemaId)||!List(row[5],&table.sourceIds)||table.sourceIds.empty()||
          !AvailabilityValue(row[6],&table.availability)||!PermissionValue(row[7],&table.evaluationPermission)||
          !List(row[8],&table.capabilities)||!PayloadPair(table.path,table.sha256,table.availability)||storage->tableById.count(table.id))
        return Error(Code::CorruptSource,"invalid/duplicate table index record: "+table.id);
      if(table.format!="csv"&&table.format!="native-ascii"&&table.format!="opaque")return Error(Code::UnsupportedSchema,"unsupported table format: "+table.format);
      if(table.format=="csv") {
        // Column labels are exact strings, not chemical IDs. Preserve spaces,
        // units and punctuation; semicolon/newline columns need a new codec.
        std::set<std::string> names;std::size_t start=0;
        do {
          const auto end=row[9].find(';',start);const auto name=row[9].substr(start,end==std::string::npos?end:end-start);
          if(name.empty()||name.find_first_of("\r\n")!=std::string::npos||!names.insert(name).second)return Error(Code::CorruptSource,"invalid/duplicate exact table column: "+table.id);
          table.columns.push_back(name);if(end==std::string::npos)break;start=end+1;
        }while(start<=row[9].size());
        if(!Count(row[10],&table.rowCount))return Error(Code::CorruptSource,"CSV table needs an explicit nonnegative row count: "+table.id);
        table.rowCountKnown=true;
      } else if(!row[9].empty()||row[10]!="unknown")return Error(Code::CorruptSource,"native-ascii/opaque rows must use empty columns and unknown row_count: "+table.id);
      for(const auto& id:table.sourceIds)if(!storage->sourceById.count(id))return Error(Code::CorruptSource,"table references absent source: "+id);
      if(!table.path.empty()) {fs::path resolved;status=RelativePath(root,table.path,&resolved);if(!status.ok())return status;}
      storage->tableById.emplace(table.id,storage->tables.size());storage->tables.push_back(std::move(table));
    }
    // Only permitted and selected payloads are opened/parsed. Restricted and
    // unresolved metadata stay inspectable; merely listing them never reads
    // their numerical bytes or grants a capability to a dependent table.
    for(const auto& source:storage->sources)if(SourceGate(source).ok()&&AnyScope(source.capabilities,identity.profileCapabilities)) {
      fs::path resolved;status=reader.Verify(source.path,source.sha256,&resolved);if(!status.ok())return status;
      storage->sourcePaths.emplace(source.id,resolved.string());
    }
    for(const auto& table:storage->tables)if(TableGate(*storage,table).ok()&&AnyScope(table.capabilities,identity.profileCapabilities)) {
      fs::path resolved;
      // A source may not expose its own payload capability in this profile;
      // dependency identity is nevertheless verified before its table is used.
      for(const auto& id:table.sourceIds) {
        const auto& source=storage->sources[storage->sourceById.at(id)];
        status=reader.Verify(source.path,source.sha256);if(!status.ok())return status;
      }
      if(table.format=="csv") {
        std::string bytes;status=reader.Read(table.path,table.sha256,&bytes,&resolved);if(!status.ok())return status;
        auto rows=std::make_shared<ReactionData::Csv>();status=Csv(bytes,table.columns,rows.get());if(!status.ok())return status;
        if(rows->size()!=table.rowCount)return Error(Code::CorruptSource,"runtime table row count mismatch: "+table.id);
        storage->tableRows.emplace(table.id,std::move(rows));
      } else {status=reader.Verify(table.path,table.sha256,&resolved);if(!status.ok())return status;}
      storage->tablePaths.emplace(table.id,resolved.string());
    }
    identity.verifiedBytes=reader.total;
    Catalogue candidate;candidate.storage_=std::move(storage);*output=std::move(candidate);return {};
  } catch(const std::bad_alloc&) {
    return Error(Code::ResourceLimit,"chemical catalogue allocation failed; reduce selected data or concurrency, raise explicit byte ceiling only when resources permit, or provide more memory; no partial catalogue was published");
  } catch(const std::length_error&) {
    return Error(Code::ResourceLimit,"chemical catalogue allocation exceeds representable storage; reduce selected data; no partial catalogue was published");
  } catch(const fs::filesystem_error& e) {
    return Error(Code::Unavailable,std::string("chemical archive filesystem error: ")+e.what());
  } catch(const std::exception& e) {
    return Error(Code::CorruptSource,std::string("chemical runtime manifest/index error: ")+e.what());
  }
}
} } // namespace AMPS::ChemicalDatabase
