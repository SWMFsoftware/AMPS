// Coefficient-kernel certification matrix (roadmap stage PS05).
// One row per released (model, kernel) pair:
//   {ModelId, Kernel, maximum nu_max*dt, epsilon min, epsilon max, evidence}
// A row is added only after its acceptance tests in
// test_pitch_angle_montecarlo.cpp pass at the stated limit; the run record
// (command, measured errors, uncertainties) is in IMPLEMENTATION_STATUS.md.
// Rows below are the certified pairs of the recorded campaign.
//
// Campaign of 2026-10-13: acceptance binary SHA-256 8e79b97e79e947a6...
// (build/acceptance/binary.sha256), full record build/acceptance/results.txt
// (52 cases incl. the f = 0.001 refinement; first campaign of binary
// 98238c7b... kept as first-campaign-results.txt, identical where shared).  Milstein passes every case at
// nu_max*dt = 0.001 and fails at 0.003 (AP-FOC-01 xi = 5: 101% of the 1%
// variance budget) and 0.01 (AP-ISO-01 endpoint beam).  PA-EPREM has the
// shape of PA-ISO (h = 1; only its amplitude-lambda relation differs), so the
// kernel sees the same D(mu)/nu_max and the PA-ISO evidence applies to it.
// Not certified: Euler-Maruyama (AP-ISO-01 fails at every step), PA-EPS with
// Milstein (AP-ZERO-02 band occupation fails at eps = 0.01, 0.048 for every
// step; eps = 1 inconclusive), random rotation (no AP-FOC-01 case run; no
// srcSEP3D scheme uses it), PA-QFORM (unbounded Ito drift, refused in code).
{ModelId::Isotropic, Kernel::Milstein, 0.001, 0.0, 0.0, "AP-ISO-01.milstein.f0.001, AP-ISO-02.milstein.f0.001.t10/t50, AP-TIME-01.milstein, AP-FOC-01.milstein.xi{0,+-0.1,+-1.5,+-5}.f0.001"},
{ModelId::Eprem, Kernel::Milstein, 0.001, 0.0, 0.0, "PA-ISO shape: AP-ISO-01.milstein.f0.001, AP-ISO-02.milstein.f0.001.t10/t50, AP-TIME-01.milstein, AP-FOC-01.milstein.xi{0,+-0.1,+-1.5,+-5}.f0.001"},
