// Legacy translation-unit path retained for application and standalone-build
// compatibility. Implementations are compiled exactly once from the shared
// src/models/chemical_database sources (ChemDb*.o in the normal core archive).
// Standalone adapters use chemical_database_build.py for the same source list;
// including another implementation here would create duplicate definitions.
#include "Photochemistry.h"
