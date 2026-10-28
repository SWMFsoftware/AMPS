#ifndef _GENERIC_COMET_PHYSICS_H_
#define _GENERIC_COMET_PHYSICS_H_

#include <string>

// Keep the reusable database in src/models.  The second include path supports
// direct source-tree unit builds; normal AMPS generated builds use the first.
#if __has_include("models/cometary_physics/CollisionDatabase.h")
#include "models/cometary_physics/CollisionDatabase.h"
#elif __has_include("../models/cometary_physics/CollisionDatabase.h")
#include "../models/cometary_physics/CollisionDatabase.h"
#else
#error "Cannot locate models/cometary_physics/CollisionDatabase.h"
#endif

namespace PIC {
namespace Mesh {
class cDataBlockAMR;
class cDataCenterNode;
}
}

namespace Comet {
namespace Physics {

// The generated table contains only data.  These functions provide the small
// application adapter needed by the AMPS neutral-collision callback.
const CometaryPhysics::cCollisionDatabase& GetCollisionDatabase();
const char *GetDatabaseName();
bool Validate(std::string *error=NULL);

// Return the published cross section for an unordered species pair.  The
// temperature is read only for a TemperaturePowerLaw entry.  A missing pair
// is a configuration error, rather than a physical zero, because Rubin et al.
// publish parent-parent values only.
double GetNeutralCrossSection(double *v0,double *v1,int species0,int species1,
    PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell);

}
}

#endif
