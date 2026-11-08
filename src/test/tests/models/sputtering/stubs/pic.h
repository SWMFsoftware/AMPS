#ifndef AMPS_TEST_MODELS_SPUTTERING_STUB_PIC_H
#define AMPS_TEST_MODELS_SPUTTERING_STUB_PIC_H

// Dependency boundary for the focused sputtering-law tests.  Production
// Sputtering.h includes the full PIC header even though its ice laws need only
// species identifiers, SI masses, one uniform random variate, and AMPS's
// fatal-error hook.  Supplying those exact contracts here lets the test compile
// the production .cpp files without MPI, mesh, particle, or application state.
// Values are explicit test fixtures, not a second physical implementation.

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

#define _PIC_MODE_OFF_ 0
#define _PIC_MODE_ON_ 1

// Particle species identifiers need only be distinct in this fixture.  Their
// numeric values have no physical meaning and are deliberately local to the
// test executable.
#define _H_PLUS_SPEC_ 0
#define _O_PLUS_THERMAL_SPEC_ 1
#define _O_PLUS_HIGH_SPEC_ 2
#define _O_PLUS_SPEC_ 3
#define _O2_PLUS_SPEC_ 4
#define _O2_SPEC_ 5
#define _H2O_SPEC_ 6

// AMPS mass-table keys used by the production implementation.
#define _H_ 100
#define _O_ 101
#define _O2_ 102
#define _H2O_ 103

// Exact 2019-SI elementary charge, also the exact conversion of one eV to J.
// Molecular masses use the same conventional integer atomic-mass composition
// assumed by the legacy sputtering model.  The tests form speeds and oracle
// energies from separate literal constants, so table expectations do not call
// this helper.
static const double ElectronCharge = 1.602176634e-19;
static const double eV2J = 1.602176634e-19;
static const double kAtomicMassKilogram = 1.66053906660e-27;

inline double amps_sputtering_test_mass(int key) {
  switch (key) {
    case _H_: return 1.0 * kAtomicMassKilogram;
    case _O_: return 16.0 * kAtomicMassKilogram;
    case _O2_: return 32.0 * kAtomicMassKilogram;
    case _H2O_: return 18.0 * kAtomicMassKilogram;
    default: throw std::invalid_argument("unknown focused-test mass key");
  }
}

#define _MASS_(KEY) amps_sputtering_test_mass(KEY)

// Each test sets the variate explicitly.  This turns the stochastic sampler
// into a reproducible evaluation at a chosen quantile without replacing the
// production transformation or weight calculation.
extern double amps_sputtering_test_random;
inline double rnd() { return amps_sputtering_test_random; }

// AMPS uses a three-argument fatal hook for unsupported compiled species and
// configurations.  Convert it to a typed exception only inside this isolated
// fixture so a failure remains observable instead of terminating the runner.
// This does not claim that the production application itself is exception
// based; unsupported angle/source-normalization capabilities are reported by
// the suite runner rather than routed through this hook.
class AmpsSputteringTestFatal : public std::runtime_error {
 public:
  explicit AmpsSputteringTestFatal(const std::string& message)
      : std::runtime_error(message) {}
};

[[noreturn]] inline void exit(int line, const char* file, const char* message) {
  std::ostringstream text;
  text << file << ':' << line << ": " << message;
  throw AmpsSputteringTestFatal(text.str());
}

#endif
