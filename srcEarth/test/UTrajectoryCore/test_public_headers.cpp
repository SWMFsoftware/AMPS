#include "../../gridless/CutoffRigidityGridless.h"
#include "../../3d/CutoffRigidityMode3D.h"

#include <type_traits>

// Compile-time API smoke test.  It deliberately includes both public backend headers
// in one translation unit so alias drift, duplicate result definitions, and overload
// ambiguity fail before a full AMPS link is attempted.
static_assert(std::is_same<Earth::GridlessMode::TrajectoryRequest,
                           Earth::Trajectory::Request>::value,
              "gridless request must alias the common Step-4 contract");
static_assert(std::is_same<Earth::GridlessMode::TrajectoryResult,
                           Earth::Trajectory::Result>::value,
              "Mode3D and gridless must return the same result type");

int main() {
  Earth::GridlessMode::TrajectoryRequest request;
  Earth::GridlessMode::TrajectoryResult result;
  return request.captureExitState || result.exitState.valid ? 1 : 0;
}
