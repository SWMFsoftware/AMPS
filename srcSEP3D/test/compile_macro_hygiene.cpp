// Compile-only regression for the AMPS constants.h namespace collision.
// AMPS defines Pi as a global macro before including an application header.
// This test deliberately reproduces that include order without requiring PIC.
#define Pi 3.141592653589793238462643383279502884
#include "../core/sep3d_types.h"

static_assert(SEP3D::Core::Const::kPi > 3.14,
              "srcSEP3D constants must survive the AMPS Pi macro");
