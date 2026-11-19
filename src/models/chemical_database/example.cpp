// Standalone HOST-only catalogue inspection example. No PIC/MPI/application
// dependency, numerical model activation, network acquisition or runtime write.
#include "Runtime.h"
#include <iostream>
int main(int argc,char** argv) {
  if(argc!=9) {
    std::cerr<<"usage: chemdb-inspect ROOT EXPORT_PATH EXPORT_SHA RELEASE_ID RELEASE_SHA PROFILE_ID PROFILE_SHA CAPABILITY\n";
    return 2;
  }
  AMPS::ChemicalDatabase::BundleRequest request;
  request.root=argv[1];request.exportPath=argv[2];request.expectedExportSha256=argv[3];
  request.expectedReleaseId=argv[4];request.expectedReleaseSha256=argv[5];
  request.expectedProfileId=argv[6];request.expectedProfileSha256=argv[7];
  request.requiredCapabilities={argv[8]};
  AMPS::ChemicalDatabase::Catalogue catalogue;
  const auto status=AMPS::ChemicalDatabase::OpenRuntimeBundle(request,&catalogue);
  std::cout<<AMPS::ChemicalDatabase::CodeName(status.code)<<": "<<status.message<<'\n';
  if(!status.ok())return 1;
  std::cout<<"release="<<catalogue.identity().releaseId<<" profile="<<catalogue.identity().profileId
      <<" sources="<<catalogue.sources().size()<<" tables="<<catalogue.tables().size()<<'\n';
  for(const auto& table:catalogue.tables())std::cout<<table.id<<' '<<table.schemaId<<' '
      <<AMPS::ChemicalDatabase::AvailabilityName(table.availability)<<'\n';
  return 0;
}
