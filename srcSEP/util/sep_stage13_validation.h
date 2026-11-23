#ifndef SRCSEP_STAGE13_VALIDATION_H
#define SRCSEP_STAGE13_VALIDATION_H

#include <iosfwd>
#include <string>

namespace SEP::Stage13 {

// Standalone campaign entry. Return 77 when the requested physical comparison
// has unmet prerequisites. Parsing and listing never initialize PIC or MPI.
int RunValidationDeck(const std::string& path, const std::string& caseName,
                      bool listOnly, bool validateOnly, std::ostream& out,
                      std::ostream& error);

}  // namespace SEP::Stage13

#endif
