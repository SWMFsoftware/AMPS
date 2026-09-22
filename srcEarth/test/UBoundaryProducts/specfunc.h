#ifndef TEST_SPECFUNC_STUB_H
#define TEST_SPECFUNC_STUB_H

// The production AMPS tree provides this three-argument fatal helper through
// specfunc.h.  The dependency-free unit executable translates it into an exception so
// invalid-input contracts can be tested without linking the AMPS runtime.
#include <stdexcept>
#include <string>

inline void exit(int line,const char* file,const char* message) {
  throw std::runtime_error(std::string(file ? file : "")+":"+
                           std::to_string(line)+": "+(message ? message : ""));
}

#endif
