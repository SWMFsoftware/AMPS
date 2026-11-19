#ifndef AMPS_CHEMICAL_DATABASE_RATE_SERVICES_H
#define AMPS_CHEMICAL_DATABASE_RATE_SERVICES_H
#include "Runtime.h"
#include "Photochemistry.h"
#include "ElectronImpact.h"
#include "ChargeExchange.h"
namespace AMPS { namespace ChemicalDatabase {
// Roles are EXPLICIT caller/profile table IDs. No filename, chemical formula
// or default latest-source search selects a numerical model. All roles need
// profile/table scalar-rates scope AND allowed dependent source permissions.
struct PhotochemistryTables {
  std::string channels,rates,thresholds,crossSections;
};
struct TotalCaptureTables {
  std::string heliumPublication,sodiumPublication;
};
Status PreparePhotochemistry(const Catalogue&,const PhotochemistryTables&,Photochemistry::Registry*);
Status PrepareElectronImpact(const Catalogue&,const std::string& nativeCfitTable,ElectronImpact::MaxwellianVoronovTable*);
Status PrepareTotalCapture(const Catalogue&,const TotalCaptureTables&,ChargeExchange::Registry*);
// Required table schema IDs are phidrates-channels-v1, phidrates-rates-v1,
// phidrates-thresholds-v1, phidrates-cross-sections-v1, voronov-cfit-native-v1,
// and total-capture-publication-v1. These are prepared service roles, not the
// canonical catalogue's exhaustive scientific schemas. Unknown roles refuse
// evaluation; metadata may still describe/retain their source objects.
} }
#endif
