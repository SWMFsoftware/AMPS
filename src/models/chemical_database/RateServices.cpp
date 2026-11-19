#include "RateServices.h"
#include <new>
#include <utility>
namespace AMPS { namespace ChemicalDatabase {
namespace {
Status Role(const Catalogue& catalogue,const std::string& id,const char* schema,const char* format,std::string* path) {
  TableRecord table;auto status=catalogue.GetTable(id,&table);if(!status.ok())return status;
  if(table.schemaId!=schema||table.format!=format)
    return {Code::UnsupportedSchema,id+": selected table does not have the required prepared-service schema/format"};
  return catalogue.TablePath(id,"scalar-rates",path);
}
Status Authenticated(const Catalogue& catalogue,const std::string& id,const std::string& path) {
  TableRecord record;auto status=catalogue.GetTable(id,&record);if(!status.ok())return status;
  std::string error;
  return ReactionData::VerifyFileSha256(path,record.sha256,&error)?Status{}:Status{Code::CorruptSource,error};
}
}
Status PreparePhotochemistry(const Catalogue& catalogue,const PhotochemistryTables& ids,Photochemistry::Registry* output) {
  if(!output)return {Code::InvalidInput,"photochemistry output registry is null"};
  try {
    auto source=Photochemistry::FrozenV1Sources("");
    const std::string tableIds[]={ids.channels,ids.rates,ids.thresholds,ids.crossSections};
    const char* schemas[]={"phidrates-channels-v1","phidrates-rates-v1","phidrates-thresholds-v1","phidrates-cross-sections-v1"};
    std::string* paths[]={&source.channelsPath,&source.ratesPath,&source.thresholdsPath,&source.crossSectionsPath};
    for(int i=0;i<4;++i) {auto status=Role(catalogue,tableIds[i],schemas[i],"csv",paths[i]);if(!status.ok())return status;}
    Photochemistry::Registry candidate;std::string error;
    // The original v1 kernel authenticates fixed-size SAME bytes before its
    // strict CSV parse. None of its source hashes/counts/numerics is replaced
    // by the catalogue's selected filenames or arbitrary caller hashes.
    if(!candidate.LoadFrozenD04(source,&error))return {Code::CorruptSource,error};
    *output=std::move(candidate);return {};
  }catch(const std::bad_alloc&){return {Code::ResourceLimit,"photochemistry preparation allocation failed; reduce concurrent prepared registries or provide more memory"};}
}
Status PrepareElectronImpact(const Catalogue& catalogue,const std::string& id,ElectronImpact::MaxwellianVoronovTable* output) {
  if(!output)return {Code::InvalidInput,"electron-impact output table is null"};
  try {
    std::string path;auto status=Role(catalogue,id,"voronov-cfit-native-v1","native-ascii",&path);if(!status.ok())return status;
    status=Authenticated(catalogue,id,path);if(!status.ok())return status;
    ElectronImpact::MaxwellianVoronovTable candidate;std::string error;
    if(!candidate.LoadNativeTable(path,&error))return {Code::CorruptSource,error};
    // Native files are immutable release objects. Check again around legacy
    // structural parsing; a changed retained file is corruption, never a new
    // fit selected by matching only a chemical species name.
    status=Authenticated(catalogue,id,path);if(!status.ok())return status;
    *output=std::move(candidate);return {};
  }catch(const std::bad_alloc&){return {Code::ResourceLimit,"electron-impact preparation allocation failed; reduce selected concurrent state or provide more memory"};}
}
Status PrepareTotalCapture(const Catalogue& catalogue,const TotalCaptureTables& ids,ChargeExchange::Registry* output) {
  if(!output)return {Code::InvalidInput,"capture output registry is null"};
  try {
    auto sources=ChargeExchange::FrozenV1References("");
    auto status=Role(catalogue,ids.heliumPublication,"total-capture-publication-v1","opaque",&sources.heliumPublicationPath);if(!status.ok())return status;
    status=Role(catalogue,ids.sodiumPublication,"total-capture-publication-v1","opaque",&sources.sodiumPublicationPath);if(!status.ok())return status;
    ChargeExchange::Registry candidate;std::string error;
    if(!candidate.LoadPublishedReferences(sources,&error))return {Code::CorruptSource,error};
    *output=std::move(candidate);return {};
  }catch(const std::bad_alloc&){return {Code::ResourceLimit,"capture preparation allocation failed; reduce concurrent prepared registries or provide more memory"};}
}
} }
