// ============================================================================
// srcSEP3D/test/individual-test/test_kernels.cpp
//
// Test group UTIL — shared-kernel frozen record (Step 3).
//
// WHY THIS GROUP EXISTS:
//   Step 3 moves the shared SEP kernels (sep_transport_common,
//   sep_coefficient_physics, sep_coefficient_registry, sep_test_registry, ...)
//   out of per-application compilation and into a single archive, sep_common.a,
//   built once with one flag set.  That is a build-system change, and the
//   danger of a build-system change is silent numerical drift:
//
//     A different optimisation level, a stray -ffast-math, or a different
//     NDEBUG state will not break the build and will not break a
//     tolerance-based test, but it can change the last few bits of an
//     adaptive quadrature (e.g. IntegrateSpatialDiffusion).  Months later
//     that shows up as an unexplained srcSEP <-> srcSEP3D disagreement.
//
//   UTIL02 catches this by freezing a BYTE-EXACT record of a fixed set of
//   kernel calls, written at full %.17g precision to
//   test/frozen/S03_kernels.txt.  The first run writes the record; every
//   later run compares against it byte for byte.  Byte identity — not a
//   tolerance — is used deliberately, because the question is "did anything
//   change", not "is it still approximately right".
//
// TESTS:
//   UTIL02 — Reproduce (or, on first run, write) the frozen kernel record.
//     The record captures, per the plan:
//       * MomentumFromSpeed / SpeedFromMomentum round-trips at six momenta
//       * SelectSubstep over a fixed set of composable limits
//       * EvaluateJokipiiSlab at twelve (input, mu) pairs
//       * IntegrateSpatialDiffusion at four speeds
//       * ConfigurationFingerprint for eight configurations
//
//   (UTIL01 — one definition per kernel symbol — is enforced by the
//   sep_common makefile's `make verify` target and, at link time, by the
//   linker; it is not a C++ test.)
//
// RE-FREEZING PROTOCOL:
//   If a change legitimately alters the kernel outputs, delete
//   test/frozen/S03_kernels.txt, run this test to regenerate it, and commit
//   the code change and the new frozen file together with a one-line
//   justification in the file header.  See test/README.md.
// ============================================================================

#include "sep3d_test_registry.h"        // SEP3D::Testing::*

// Shared-kernel headers (from the sibling sep_common library via -I).
#include "sep_transport_common.h"
#include "sep_coefficient_physics.h"
#include "sep_coefficient_registry.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using R = SEP3D::Testing::Result;
using S = SEP3D::Testing::Status;

R Pass(std::string msg) { R r; r.status = S::Pass; r.message = std::move(msg); return r; }
R Fail(std::string msg) { R r; r.status = S::Fail; r.message = std::move(msg); return r; }
R Err (std::string msg) { R r; r.status = S::Error; r.message = std::move(msg); return r; }

// Namespace aliases for the shared kernels.
namespace TC = SEP::Transport;
namespace CP = SEP::Transport::CoefficientPhysics;
namespace CR = SEP::Transport::Coefficient;

// Physical constants used by the record (SI).  These match the srcSEP3D core
// constants; GEO03 (Step 6) will assert bitwise agreement with SWCME.
constexpr double C_LIGHT = 2.99792458e+08;   // m/s
constexpr double M_PROTON = 1.67262192369e-27; // kg
constexpr double E_CHARGE = 1.602176634e-19;   // C

// ----------------------------------------------------------------------------
// Deterministic record builder.
//
// Appends fixed-format lines to `out`.  Every floating-point value is written
// at %.17g so the record is byte-exact and portable across platforms that
// use IEEE-754 double.  No timestamps, no addresses, no rank — nothing that
// varies between runs.
// ----------------------------------------------------------------------------
void AppendDouble(std::ostringstream& out, const char* label, double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  out << label << ' ' << buf << '\n';
}

std::string BuildRecord() {
  std::ostringstream out;
  out << "# srcSEP3D frozen kernel record (UTIL02, Step 3)\n";
  out << "# Byte-exact reference of shared-kernel outputs. Re-freeze only\n";
  out << "# alongside a deliberate change, with a justification in this header.\n";
  out << "# Format: <label> <value at %.17g>\n";

  // ---- (1) MomentumFromSpeed / SpeedFromMomentum round-trip ---------------
  // Six speeds spanning non-relativistic to mildly relativistic protons.
  out << "SECTION momentum_speed_roundtrip\n";
  const double speeds[6] = {
    1.0e5, 1.0e6, 5.0e6, 1.0e7, 5.0e7, 1.0e8
  };
  for (int i = 0; i < 6; ++i) {
    TC::ScalarResult p = TC::MomentumFromSpeed(speeds[i], M_PROTON, C_LIGHT);
    TC::ScalarResult v = TC::SpeedFromMomentum(p.value, M_PROTON, C_LIGHT);
    std::ostringstream lbl_p, lbl_v;
    lbl_p << "  p_from_v[" << i << "]";
    lbl_v << "  v_roundtrip[" << i << "]";
    AppendDouble(out, lbl_p.str().c_str(), p.value);
    AppendDouble(out, lbl_v.str().c_str(), v.value);
  }

  // ---- (2) SelectSubstep over a fixed limit set ---------------------------
  out << "SECTION select_substep\n";
  {
    std::vector<TC::StepLimit> limits;
    limits.emplace_back("cell-crossing", 3.0);
    limits.emplace_back("diffusion",     1.5);
    limits.emplace_back("cooling",       7.2);
    limits.emplace_back("focusing",      0.85);
    TC::StepDiagnostics diag;
    TC::ScalarResult step =
        TC::SelectSubstep(10.0, limits, 1.0e-12, &diag);
    AppendDouble(out, "  selected_dt", step.value);
    // The binding limit name is deterministic; record it too.
    if (!diag.limitingNames.empty())
      out << "  binding_name " << diag.limitingNames.back() << '\n';
    else
      out << "  binding_name (none)\n";
  }

  // ---- (3) EvaluateJokipiiSlab at twelve (input, mu) pairs -----------------
  out << "SECTION jokipii_slab\n";
  {
    CP::SpeciesProperties proton;
    proton.modelSpecies = 0;
    proton.name = "h+";
    proton.signedChargeC = E_CHARGE;
    proton.restMassKg = M_PROTON;
    proton.nucleonCount = 1.0;

    CP::SpectrumParameters spectrum;   // defaults from the header

    CP::LocalInputView input;
    input.source = "prescribed";
    input.representation = "slab";
    input.generation = 1;
    input.checksum = 0x5eed;                      // must be non-zero (validated)
    input.heliocentricRadiusM = 1.495978707e11;   // 1 AU
    input.magneticFieldT = 5.0e-9;                // ~5 nT at 1 AU
    input.deltaB2T2 = (0.3 * 5.0e-9) * (0.3 * 5.0e-9);
    input.deltaBPlus2T2 = 0.5 * input.deltaB2T2;
    input.deltaBMinus2T2 = 0.5 * input.deltaB2T2;
    input.alfvenSpeedMPerS = 5.0e4;

    const double speed = 3.0e7;  // ~0.1c proton
    const double mus[12] = {
      -0.95, -0.75, -0.5, -0.25, -0.1, -0.01,
       0.01,  0.1,   0.25, 0.5,   0.75, 0.95
    };
    for (int i = 0; i < 12; ++i) {
      CP::PitchAngleResult res =
          CP::EvaluateJokipiiSlab(input, spectrum, proton, speed, mus[i]);
      std::ostringstream ld, ldd;
      ld  << "  Dmumu[" << i << "]";
      ldd << "  dDmumu_dmu[" << i << "]";
      AppendDouble(out, ld.str().c_str(),  res.dMuMuPerS);
      AppendDouble(out, ldd.str().c_str(), res.dDmuMuDmuPerS);
    }
  }

  // ---- (4) IntegrateSpatialDiffusion at four speeds -----------------------
  out << "SECTION spatial_diffusion\n";
  {
    // A simple analytic Dmumu = D0 (1 - mu^2) so the integral is well defined
    // and the record isolates the quadrature, not the QLT kernel.
    const double D0 = 1.0e-4;   // 1/s
    CP::PitchAngleFunction dmumu = [D0](double mu) -> CP::PitchAngleResult {
      CP::PitchAngleResult r;
      r.status = TC::Status::Ok();
      r.valueState = CP::ValueState::Finite;
      r.dMuMuPerS = D0 * (1.0 - mu * mu);
      r.dDmuMuDmuPerS = -2.0 * D0 * mu;
      return r;
    };
    CP::SpatialQuadratureConfiguration cfg;   // defaults
    const double speeds4[4] = {1.0e6, 5.0e6, 1.0e7, 3.0e7};
    for (int i = 0; i < 4; ++i) {
      CP::SpatialDiffusionResult res =
          CP::IntegrateSpatialDiffusion(speeds4[i], dmumu, cfg);
      std::ostringstream lk;
      lk << "  kappa_parallel[" << i << "]";
      AppendDouble(out, lk.str().c_str(), res.kappaParallelM2PerS);
    }
  }

  // ---- (5) ConfigurationFingerprint for eight configurations --------------
  out << "SECTION configuration_fingerprint\n";
  {
    // Eight distinct configurations toggling the primary selectors.  The
    // exact fingerprint string is what a restart must match, so it belongs
    // in the frozen record.
    CR::Configuration base;    // all defaults

    std::vector<CR::Configuration> configs;
    configs.push_back(base);                                    // 0: defaults

    CR::Configuration c1 = base; c1.invalidPolicy = CR::InvalidPolicy::Ballistic;
    configs.push_back(c1);                                      // 1

    CR::Configuration c2 = base; c2.prescribedDeltaBOverB = 0.5;
    configs.push_back(c2);                                      // 2

    CR::Configuration c3 = base; c3.constantDmumuPerS = 1.0e-4;
    configs.push_back(c3);                                      // 3

    CR::Configuration c4 = base; c4.spectrum.spectralIndex = 2.0;
    configs.push_back(c4);                                      // 4

    CR::Configuration c5 = base; c5.correlationLengthAt1AuM = 0.02 * 1.495978707e11;
    configs.push_back(c5);                                      // 5

    CR::Configuration c6 = base;
    c6.resonanceGapPolicy = CP::ResonanceGapPolicy::Ballistic;
    configs.push_back(c6);                                      // 6

    CR::Configuration c7 = base;
    c7.amplitudePolicy = CR::TurbulenceAmplitudePolicy::LimitToMeanField;
    configs.push_back(c7);                                      // 7

    for (std::size_t i = 0; i < configs.size(); ++i) {
      std::string fp = CR::ConfigurationFingerprint(configs[i]);
      out << "  fingerprint[" << i << "] " << fp << '\n';
    }
  }

  return out.str();
}

// ----------------------------------------------------------------------------
// Read an existing file into a string.  Returns false if it does not exist.
// ----------------------------------------------------------------------------
bool ReadFile(const std::string& path, std::string* content) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  std::ostringstream ss;
  ss << f.rdbuf();
  *content = ss.str();
  return true;
}

bool WriteFile(const std::string& path, const std::string& content) {
  std::ofstream f(path, std::ios::binary | std::ios::trunc);
  if (!f) return false;
  f << content;
  return f.good();
}

// ----------------------------------------------------------------------------
// UTIL02 — frozen kernel record
// ----------------------------------------------------------------------------
R run_UTIL02() {
  // Path is relative to the srcSEP3D root, which is the CWD when the runner
  // invokes the binary.
  const std::string path = "test/frozen/S03_kernels.txt";

  const std::string produced = BuildRecord();

  std::string frozen;
  if (!ReadFile(path, &frozen)) {
    // First run on a clean tree: write the record and pass, noting that it
    // was created (so the developer knows to commit it).
    if (!WriteFile(path, produced))
      return Err("could not write frozen record to " + path
                 + " (is test/frozen/ writable, and is the CWD the "
                   "srcSEP3D root?)");
    return Pass("frozen record did not exist; wrote " + path
                + " (commit it).  Re-run to verify byte-for-byte.");
  }

  // Subsequent runs: byte-exact comparison.
  if (produced == frozen)
    return Pass("kernel outputs match the frozen record byte-for-byte ("
                + path + ").");

  // Mismatch: report the first differing line to make diagnosis quick.
  std::istringstream a(produced), b(frozen);
  std::string la, lb;
  int line = 0;
  while (std::getline(a, la)) {
    ++line;
    if (!std::getline(b, lb)) {
      return Fail("frozen record is shorter than the produced output; first "
                  "extra produced line " + std::to_string(line) + ": " + la
                  + ".  If this change is intended, delete " + path
                  + " and re-run to re-freeze (commit with justification).");
    }
    if (la != lb) {
      return Fail("kernel output changed at line " + std::to_string(line)
                  + ":\n    frozen:   " + lb + "\n    produced: " + la
                  + "\n  A shared kernel's numerical output changed.  If this "
                    "is intended (e.g. a deliberate algorithm change), delete "
                  + path + " and re-run to re-freeze, committing the code "
                    "change and the new record together with a justification. "
                    "If it is NOT intended, a build-flag change (optimisation "
                    "level, -ffast-math, NDEBUG) has likely perturbed an "
                    "adaptive quadrature — check the sep_common.a build flags.");
    }
  }
  if (std::getline(b, lb)) {
    return Fail("produced output is shorter than the frozen record (frozen "
                "has extra line " + std::to_string(line + 1) + ": " + lb
                + ").  Re-freeze if intended.");
  }
  return Fail("records differ in a way not localised to a single line.");
}

} // anonymous namespace


// ============================================================================
// RegisterKernelTests — called from test/stage1.cpp
// ============================================================================
std::vector<SEP3D::Testing::Descriptor> RegisterKernelTests() {
  using D  = SEP3D::Testing::Descriptor;
  using IL = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;

  D d;
  d.id                  = "UTIL02";
  d.name                = "Shared-kernel frozen record (byte-exact)";
  d.group               = "UTIL";
  d.description         =
      "Reproduces a fixed set of sep_common.a kernel calls (momentum/speed, "
      "SelectSubstep, Jokipii Dmumu, spatial-diffusion quadrature, "
      "ConfigurationFingerprint) and compares byte-for-byte with "
      "test/frozen/S03_kernels.txt.";
  d.initialization      = IL::None;
  d.supportedBuildModes = "all";
  d.runtime             = RC::Routine;
  d.seedPolicy          = "deterministic-no-rng";
  d.stateIsolation      = "reads/writes only test/frozen/S03_kernels.txt";
  d.callback            = run_UTIL02;

  return { d };
}
