#ifndef SRCSEP_VALIDATION_CASES_CV01_MODEL_H
#define SRCSEP_VALIDATION_CASES_CV01_MODEL_H

#include <string>
#include <vector>

namespace SEP {
namespace Validation {
namespace CV01 {

// Execute one controlled ballistic realization inside the linked srcSEP/AMPS
// process. Arguments use the same ``--name value`` representation as the
// former dependency-light driver, but no second executable or parser is
// involved in validation runs. The caller supplies the output path explicitly
// so MPI-aware registry code can enforce single-writer ownership.
//
// All physical quantities are SI: distances [m], speeds [m/s], mass [kg], and
// times [s]. A false return denotes configuration, physics-kernel, boundary,
// or transactional-output failure; ``error`` contains the operator-facing
// diagnostic and must be propagated as a registry ERROR.
bool RunModel(const std::vector<std::string>& arguments,
              const std::string& outputPath,
              std::string* error);

}  // namespace CV01
}  // namespace Validation
}  // namespace SEP

#endif  // SRCSEP_VALIDATION_CASES_CV01_MODEL_H
