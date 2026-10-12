// ============================================================================
// Suspension of SWCME-input tests while the SWCME model is not ready.
//
// User decision 2026-10-12: the canonical SWCME model is not ready to be used
// in srcSEP3D, and every '#begin physics: swcme' block stops the parser with
// "the SWCME model is not ready yet to be used".  The SWCME example decks are
// kept as rejected fixtures, and the tests that ran them now assert that
// typed stop.
//
// SuspendedWhileSwcmeNotReady(deck, body) wraps such a test:
//   * the deck stops with the not-ready diagnosis  -> SKIP, naming the deck
//     and the asserted diagnosis (the SWCME physics checks of the body are
//     suspended, so the result is not reported as a PASS);
//   * the deck fails for any other reason           -> FAIL;
//   * the deck parses (SWCME released again)        -> the original body runs
//     unchanged.
// CFG3D08 additionally checks the exact file:line of each rejection.
// ============================================================================

#ifndef SEP3D_TEST_SWCME_SUSPENSION_H
#define SEP3D_TEST_SWCME_SUSPENSION_H

#include "sep3d_test_registry.h"
#include "configuration_io.h"

#include <string>
#include <utility>

namespace SEP3D {
namespace Testing {

inline TestCallback SuspendedWhileSwcmeNotReady(std::string deck,
                                                TestCallback body) {
  return [deck = std::move(deck), body = std::move(body)]() {
    RuntimeModel::RunConfiguration3DOptions options;
    const Core::Status status = RuntimeModel::LoadConfigurationFile(deck, &options);
    if (status.ok()) return body();
    Result result;
    if (status.message.find("the SWCME model is not ready yet to be used") !=
        std::string::npos) {
      result.status = SEP3D::Testing::Status::Skip;
      result.message = "SWCME model not ready (user decision 2026-10-12): " +
          deck + " stops at its physics: swcme block with the typed "
          "not-ready diagnosis (asserted); this SWCME test is suspended";
    } else {
      result.status = SEP3D::Testing::Status::Fail;
      result.message = deck + " did not stop with the SWCME not-ready "
          "diagnosis: " + status.message;
    }
    return result;
  };
}

}  // namespace Testing
}  // namespace SEP3D

#endif  // SEP3D_TEST_SWCME_SUSPENSION_H
