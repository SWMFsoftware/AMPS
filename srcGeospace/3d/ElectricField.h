#ifndef _SRC_GEOSPACE_3D_ELECTRICFIELD_H_
#define _SRC_GEOSPACE_3D_ELECTRICFIELD_H_

#include "../util/amps_param_parser.h"

namespace Geospace {
namespace Mode3D {

void EvaluateBackgroundMagneticFieldSI(double B[3],const double xGSM_SI[3],const GeospaceUtil::AmpsParam& prm);
void EvaluateElectricFieldSI(double E[3],const double xGSM_SI[3],const GeospaceUtil::AmpsParam& prm);
double EvaluateElectricPotential_V(const double xGSM_SI[3],const GeospaceUtil::AmpsParam& prm);

} // namespace Mode3D
} // namespace Geospace

#endif
