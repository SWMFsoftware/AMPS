// Copyright (C) 2002 Regents of the University of Michigan.
//
// Linked verification entry points for the srcMoon M0 baseline.  These
// routines are compiled into the normal Moon executable.  They are deliberately
// limited to selecting deterministic fixtures and recording diagnostics: the
// particle mover, force callback, surface callback, source dispatcher, mesh,
// and species tables remain the production implementations initialized by
// amps_init().

#ifndef MOON_INTEGRATION_TESTS_H_
#define MOON_INTEGRATION_TESTS_H_

#include <string>

namespace Moon {
namespace IntegrationTests {

// Run one linked diagnostic after amps_init() has completed.  The output is a
// UTF-8 key/value receipt consumed by the Python I-series runner.  Returning 0
// means the diagnostic executed; scientific PASS/FAIL is decided independently
// by the runner from the frozen acceptance contract.
int Run(const std::string& testId,const std::string& outputFile);

}  // namespace IntegrationTests
}  // namespace Moon

#endif  // MOON_INTEGRATION_TESTS_H_
