// Compile-only/runtime-smoke probe for the relocated SWCME public include root.
// The test intentionally uses no relative path and no srcSEP header. Therefore
// it succeeds only when the runner supplies AMPS/src/models/swcme explicitly,
// which is the same dependency contract used by both production applications.

#include "swcme1d.hpp"

#include <cmath>

int main() {
  swcme1d::Params parameters;
  parameters.V_sw_kms = 400.0;
  parameters.n1AU_cm3 = 5.0;
  parameters.B1AU_nT = 5.0;

  const swcme1d::Model model(parameters);
  const swcme1d::StepState state = model.prepare_step(0.0);

  // This is a path/API smoke test, not a scientific tolerance test. Requiring
  // finite SI state nevertheless prevents an empty compatibility shim from
  // satisfying the build check.
  return std::isfinite(state.r_sh_m) && std::isfinite(state.V_sh_ms) &&
                 state.r_sh_m > 0.0
             ? 0
             : 1;
}
