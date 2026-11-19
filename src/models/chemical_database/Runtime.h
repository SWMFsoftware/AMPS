#ifndef AMPS_CHEMICAL_DATABASE_RUNTIME_H
#define AMPS_CHEMICAL_DATABASE_RUNTIME_H
#include "DataIO.h"
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace AMPS { namespace ChemicalDatabase {
// These are DATA-ACCESS outcomes, not scientific PASS/FAIL or authority to
// execute a reaction. Application/test owners retain those separate meanings.
enum class Code { Evaluated, Unavailable, InvalidInput, OutsideDomain,
  CorruptSource, UnsupportedSchema, PermissionRequired, NumericalFailure, ResourceLimit };
struct Status {
  Code code=Code::Evaluated;
  std::string message;
  bool ok()const{return code==Code::Evaluated;}
};
const char* CodeName(Code);
enum class Availability { Available, Unavailable, MetadataOnly, Quarantined };
enum class Permission { Allowed, Pending, Prohibited, Unknown, NotApplicable };
const char* AvailabilityName(Availability);
const char* PermissionName(Permission);

struct SourceRecord {
  std::string id,path,sha256,mediaType,sourceUri,notes;
  Availability availability=Availability::Unavailable;
  Permission storagePermission=Permission::Unknown,evaluationPermission=Permission::Unknown;
  std::vector<std::string> capabilities;
};
struct TableRecord {
  std::string id,path,sha256,format,schemaId,notes;
  std::vector<std::string> sourceIds,capabilities,columns;
  Availability availability=Availability::Unavailable;
  Permission evaluationPermission=Permission::Unknown;
  bool rowCountKnown=false;
  std::size_t rowCount=0;
};
struct ResourceLimits {
  // Explicit byte ceilings, not physics/data-selection defaults. Zero means
  // no artificial ceiling; size_t/stream/allocation limits still fail closed.
  // MaximumTotalBytes counts each actually read file once per open operation.
  std::size_t maximumFileBytes=0,maximumTotalBytes=0;
};
struct BundleRequest {
  std::string root; // Existing persistent chemistry archive, explicitly supplied.
  std::string exportPath; // ROOT-relative detached export.in (never a cwd path).
  std::string expectedExportSha256,expectedReleaseId,expectedReleaseSha256;
  std::string expectedProfileId,expectedProfileSha256;
  std::string expectedCodec="csv-block-v1";
  std::vector<std::string> requiredCapabilities;
  ResourceLimits limits;
};
struct BundleIdentity {
  std::string root,exportPath,exportSha256,releaseId,releaseManifestPath,releaseSha256;
  std::string profileId,profilePath,profileSha256,codec,topologyFingerprint;
  std::vector<std::string> profileCapabilities;
  std::size_t verifiedBytes=0;
};
namespace detail { struct CatalogueStorage; }
class Catalogue {
 public:
  bool loaded()const{return static_cast<bool>(storage_);}
  // These references are valid while this Catalogue (or a copy) stays alive.
  // The underlying storage is immutable and shared: old views cannot dangle
  // when another Catalogue is opened. No hot reload mutates existing storage.
  const BundleIdentity& identity()const;
  const std::vector<SourceRecord>& sources()const;
  const std::vector<TableRecord>& tables()const;
  Status GetSource(const std::string& id,SourceRecord* output)const;
  Status GetTable(const std::string& id,TableRecord* output)const;
  // Checks the PROFILE'S explicit allowed selection scope, not scientific
  // completeness. Every actual source/table access also checks its own flags.
  Status RequireCapability(const std::string& capability)const;
  Status SourcePath(const std::string& id,const std::string& capability,std::string* output)const;
  Status TablePath(const std::string& id,const std::string& capability,std::string* output)const;
  // CSV rows are the exact parsed strings (no science conversion/repair).
  // Owning shared_ptr views remain valid after the Catalogue is destroyed.
  Status TableRows(const std::string& id,const std::string& capability,
      std::shared_ptr<const ReactionData::Csv>* output)const;
 private:
  friend Status OpenRuntimeBundle(const BundleRequest&,Catalogue*);
  std::shared_ptr<const detail::CatalogueStorage> storage_;
};
// One transactional read-only open. The manifest's SHA is checked before its
// generic-grammar parse. JSON release/profile bytes are authenticated, never
// interpreted by a second home-grown JSON parser. Includes are forbidden in
// this already self-contained detached export; every child object has an
// explicit digest. No network, runtime writes, RNG, PIC/MPI or application API.
Status OpenRuntimeBundle(const BundleRequest&,Catalogue* output);
} }
#endif
