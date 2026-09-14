#include "sep_configuration_matrix.h"

#include <iostream>

int main() {
  // Documentation uses the exact runtime classifier instead of duplicating
  // 90 hand-maintained rows.  Redirect this output when a versioned evidence
  // artifact is required; no physics or support decision is reimplemented here.
  std::cout << SEP::ConfigurationMatrix::RenderMarkdown();
  return 0;
}
