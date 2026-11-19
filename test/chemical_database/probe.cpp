// Thin read-only adapter over the maintained native Catalogue.  Metadata-only
// fixtures and real releases use the SAME loader, hash checks and permission
// gates.  This executable implements no rates, repairs, parser or observations.
#include "../../src/models/chemical_database/Runtime.h"
#include "../../srcMoon/species/Photochemistry.h"
#include "../../srcMoon/species/ElectronImpact.h"
#include "../../srcMoon/species/ChargeExchange.h"
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <atomic>
#include <thread>

namespace DB=AMPS::ChemicalDatabase;
// Compile-time proof that old application names are actual shared service
// types. Merely linking a second unused implementation would not establish
// that Moon calls the maintained service exercised by database clients.
static_assert(std::is_same<Moon::Species::Photochemistry::Registry,DB::Photochemistry::Registry>::value,"photon compatibility must alias shared provider");
static_assert(std::is_same<Moon::Species::ElectronImpact::MaxwellianVoronovTable,DB::ElectronImpact::MaxwellianVoronovTable>::value,"electron compatibility must alias shared provider");
static_assert(std::is_same<Moon::Species::ChargeExchange::Registry,DB::ChargeExchange::Registry>::value,"capture compatibility must alias shared provider");
static std::string Json(const std::string& value) {
  static const char hex[]="0123456789abcdef";
  std::string out="\"";
  for (unsigned char c:value) {
    if(c=='"'||c=='\\') {out+='\\';out+=static_cast<char>(c);}
    else if(c<32) {out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}
    else out+=static_cast<char>(c);
  }
  return out+'"';
}
static void Strings(const std::vector<std::string>& values) {
  std::cout<<'[';
  for(std::size_t i=0;i<values.size();++i) {if(i)std::cout<<',';std::cout<<Json(values[i]);}
  std::cout<<']';
}
int main(int argc,char** argv) {
  if(argc<10) {std::cerr<<"ROOT EXPORT HASH RELEASE RELEASEHASH PROFILE PROFILEHASH CAP OP [ID] [FILELIMIT]\n";return 2;}
  DB::BundleRequest request;
  request.root=argv[1];request.exportPath=argv[2];request.expectedExportSha256=argv[3];
  request.expectedReleaseId=argv[4];request.expectedReleaseSha256=argv[5];
  request.expectedProfileId=argv[6];request.expectedProfileSha256=argv[7];
  const std::string capability=std::string(argv[8])=="none"?"":argv[8];
  if(!capability.empty())request.requiredCapabilities.push_back(capability);
  if(argc>11)request.limits.maximumFileBytes=std::stoull(argv[11]);
  if(argc>12)request.limits.maximumTotalBytes=std::stoull(argv[12]);
  DB::Catalogue catalogue;
  const DB::Status opened=DB::OpenRuntimeBundle(request,&catalogue);
  std::cout<<"{\"open_code\":"<<Json(DB::CodeName(opened.code))
           <<",\"open_ok\":"<<(opened.ok()?"true":"false")
           <<",\"open_message\":"<<Json(opened.message)
           <<",\"loaded\":"<<(catalogue.loaded()?"true":"false");
  if(opened.ok()) {
    std::cout<<",\"release_id\":"<<Json(catalogue.identity().releaseId)<<",\"source_ids\":[";
    for(std::size_t i=0;i<catalogue.sources().size();++i) {if(i)std::cout<<',';std::cout<<Json(catalogue.sources()[i].id);}
    std::cout<<"],\"table_ids\":[";
    for(std::size_t i=0;i<catalogue.tables().size();++i) {if(i)std::cout<<',';std::cout<<Json(catalogue.tables()[i].id);}
    std::cout<<']';
    const std::string operation=argv[9],id=argc>10?argv[10]:"";
    if(operation=="table"||operation=="lifetime") {
      auto sentinel=std::make_shared<const DB::ReactionData::Csv>(DB::ReactionData::Csv{{"UNCHANGED"}});
      std::shared_ptr<const DB::ReactionData::Csv> rows=sentinel;
      const auto status=catalogue.TableRows(id,capability,&rows);
      std::cout<<",\"query_code\":"<<Json(DB::CodeName(status.code))
               <<",\"query_ok\":"<<(status.ok()?"true":"false")
               <<",\"query_message\":"<<Json(status.message)
               <<",\"output_preserved\":"<<(rows==sentinel?"true":"false");
      if(status.ok()) {
        // Owning immutable row storage MUST outlive all catalogue handles.
        if(operation=="lifetime")catalogue=DB::Catalogue();
        std::cout<<",\"rows\":[";
        for(std::size_t i=0;i<rows->size();++i) {if(i)std::cout<<',';Strings((*rows)[i]);}
        std::cout<<']';
      }
    } else if(operation=="path"||operation=="source-path") {
      std::string path="UNCHANGED";
      const auto status=operation=="path"?catalogue.TablePath(id,capability,&path):catalogue.SourcePath(id,capability,&path);
      std::cout<<",\"query_code\":"<<Json(DB::CodeName(status.code))
               <<",\"query_ok\":"<<(status.ok()?"true":"false")
               <<",\"output_preserved\":"<<(path=="UNCHANGED"?"true":"false")
               <<",\"path\":"<<Json(path);
    } else if(operation=="denial-sweep") {
      // Inspect EVERY held metadata entry through the production permission
      // gate in one authenticated Catalogue. No science values are sampled,
      // no payload permission is changed, and output sentinels prove refused
      // queries are transactional rather than returning a fake empty table.
      bool all=true,first=true;
      std::cout<<",\"denials\":[";
      for(const auto& source:catalogue.sources()) {
        if(source.storagePermission==DB::Permission::Allowed &&
           source.evaluationPermission==DB::Permission::Allowed)continue;
        std::string path="UNCHANGED";
        const auto status=catalogue.SourcePath(source.id,"source-access",&path);
        const bool preserved=path=="UNCHANGED";
        all=all && status.code==DB::Code::PermissionRequired && preserved;
        if(!first)std::cout<<',';
        first=false;
        std::cout<<"{\"kind\":\"source-path\",\"id\":"<<Json(source.id)
                 <<",\"query_code\":"<<Json(DB::CodeName(status.code))
                 <<",\"query_ok\":"<<(status.ok()?"true":"false")
                 <<",\"output_preserved\":"<<(preserved?"true":"false")<<'}';
      }
      for(const auto& table:catalogue.tables()) {
        if(table.evaluationPermission==DB::Permission::Allowed)continue;
        const auto sentinel=std::make_shared<const DB::ReactionData::Csv>(
          DB::ReactionData::Csv{{"UNCHANGED"}});
        auto rows=sentinel;
        const auto status=catalogue.TableRows(table.id,"table-access",&rows);
        const bool preserved=rows==sentinel;
        all=all && status.code==DB::Code::PermissionRequired && preserved;
        if(!first)std::cout<<',';
        first=false;
        std::cout<<"{\"kind\":\"table\",\"id\":"<<Json(table.id)
                 <<",\"query_code\":"<<Json(DB::CodeName(status.code))
                 <<",\"query_ok\":"<<(status.ok()?"true":"false")
                 <<",\"output_preserved\":"<<(preserved?"true":"false")<<'}';
      }
      std::cout<<"],\"denial_sweep_ok\":"<<(all?"true":"false");
    } else if(operation=="concurrent") {
      std::shared_ptr<const DB::ReactionData::Csv> expected;
      const auto first=catalogue.TableRows(id,capability,&expected);
      std::atomic<bool> valid(first.ok());std::vector<std::thread> readers;
      // Repeated concurrent reads copy an owning immutable catalogue handle;
      // no thread mutates/reloads a shared Catalogue object during these reads.
      for(int worker=0;worker<8;++worker)readers.emplace_back([catalogue,id,capability,expected,&valid]() {
        for(int iteration=0;iteration<100;++iteration) {
          std::shared_ptr<const DB::ReactionData::Csv> rows;
          if(!catalogue.TableRows(id,capability,&rows).ok()||rows!=expected)valid=false;
        }
      });
      for(auto& reader:readers)reader.join();
      std::cout<<",\"concurrent_reads_ok\":"<<(valid?"true":"false");
    } else if(operation=="sha-vectors") {
      std::string empty,abc,error;
      const bool valid=DB::ReactionData::BytesSha256("",&empty,&error)&&
        DB::ReactionData::BytesSha256("abc",&abc,&error)&&
        empty=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"&&
        abc=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
      std::cout<<",\"known_sha_vectors_ok\":"<<(valid?"true":"false");
    } else if(operation=="guards") {
      DB::SourceRecord source;source.id="UNCHANGED";
      DB::TableRecord table;table.id="UNCHANGED";
      const auto sourceMissing=catalogue.GetSource("not-declared",&source);
      const auto tableMissing=catalogue.GetTable("not-declared",&table);
      DB::Catalogue unloaded;std::string path="UNCHANGED";
      const bool guards=!DB::OpenRuntimeBundle(request,nullptr).ok() &&
        !catalogue.GetSource(id,nullptr).ok()&&!catalogue.GetTable(id,nullptr).ok()&&
        !catalogue.SourcePath(id,capability,nullptr).ok()&&!catalogue.TablePath(id,capability,nullptr).ok()&&
        !catalogue.TableRows(id,capability,nullptr).ok()&&!sourceMissing.ok()&&!tableMissing.ok()&&
        source.id=="UNCHANGED"&&table.id=="UNCHANGED"&&
        !unloaded.SourcePath(id,capability,&path).ok()&&path=="UNCHANGED";
      std::cout<<",\"guards_ok\":"<<(guards?"true":"false");
    } else if(operation=="rollback") {
      const DB::Catalogue copy=catalogue;
      const auto* identity=&catalogue.identity();
      request.expectedExportSha256=std::string(64,'0');
      const auto rejected=DB::OpenRuntimeBundle(request,&catalogue);
      std::cout<<",\"reopen_code\":"<<Json(DB::CodeName(rejected.code))
               <<",\"reopen_ok\":"<<(rejected.ok()?"true":"false")
               <<",\"old_identity_preserved\":"<<(&catalogue.identity()==identity?"true":"false")
               <<",\"copy_identity_shared\":"<<(&copy.identity()==&catalogue.identity()?"true":"false");
    }
  }
  std::cout<<"}\n";
  return 0; // Data-access outcome is explicit JSON, not a scientific exit code.
}
