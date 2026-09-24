#ifndef _SRC_EARTH_TEST_UDIRECTIONAL_ACCESS_SPECFUNC_STUB_H_
#define _SRC_EARTH_TEST_UDIRECTIONAL_ACCESS_SPECFUNC_STUB_H_

// The production CLI reports validation failures through AMPS' three-argument exit()
// helper.  This dependency-free unit suite substitutes an exception so negative CLI
// cases can be asserted without linking the AMPS runtime or terminating the test
// process.  The production translation unit itself is compiled unchanged.

#include <stdexcept>
#include <string>

inline void exit(int line,const char* file,const char* message) {
  throw std::runtime_error(
      std::string(file ? file : "")+":"+std::to_string(line)+": "+
      (message ? message : ""));
}

#endif
