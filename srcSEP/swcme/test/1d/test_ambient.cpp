#include "test_framework.hpp"

#include <swcme1d.hpp>

void test_1d_ambient_at_one_au(swcme_test::Context& context) {
  swcme1d::Params parameters;
  parameters.V_sw_kms = 475.0;
  parameters.n1AU_cm3 = 8.25;

  const swcme1d::Model model(parameters);
  const swcme1d::StepState state = model.prepare_step(0.0);
  const double radius_m[] = {swcme1d::AU};
  double density_m3[] = {0.0};
  double speed_ms[] = {0.0};

  model.evaluate_radii_fast(state, radius_m, density_m3, speed_ms, 1);

  context.expect_true(state.r_sh_m < radius_m[0],
                      "the sample must be upstream of the launch-time shock");
  context.expect_near(density_m3[0], parameters.n1AU_cm3 * 1.0e6, 1.0e-6,
                      "production density must preserve the configured 1-AU value");
  context.expect_near(speed_ms[0], parameters.V_sw_kms * 1.0e3, 1.0e-9,
                      "production speed must preserve the configured ambient wind");
}
