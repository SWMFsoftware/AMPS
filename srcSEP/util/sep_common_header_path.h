#ifndef SRCSEP_UTIL_SEP_COMMON_HEADER_PATH_H
#define SRCSEP_UTIL_SEP_COMMON_HEADER_PATH_H

// Resolve canonical SEP-common headers in both AMPS compilation contexts.
//
// Application objects are compiled by srcSEP/makefile, which supplies
//   -I<AMPS_ROOT>/src/models/sep_common
// and therefore uses the short public names.  AMPS also includes build/main/
// sep.h transitively while compiling sibling libraries such as meshAMR and
// PIC.  Those sibling submakes inherit -I<AMPS_ROOT>, but they do not inherit
// variables appended inside the already-completed build/main submake.  In that
// context a bare "sep_transport_common.h" cannot be found even though the
// canonical source is installed correctly.
//
// The macro below selects between those two existing include roots. The
// AMPS-root form is preferred whenever available so a stale short-name header
// elsewhere in an old build directory cannot win by include order. It does not
// create an application-local forwarding header and it never changes the
// binary owner: every selected file is the one under src/models/sep_common.
// Header users write, for example,
//
//   #include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
//
// A missing header fails during preprocessing with the explicit diagnostic
// below instead of being accidentally satisfied by a stale private copy.
#if defined(__has_include)
#  if __has_include(<src/models/sep_common/sep_transport_common.h>)
#    define SRCSEP_SEP_COMMON_HEADER(name) <src/models/sep_common/name>
#  elif __has_include(<sep_transport_common.h>)
#    define SRCSEP_SEP_COMMON_HEADER(name) <name>
#  else
#    error "Cannot locate canonical src/models/sep_common headers"
#  endif
#else
// Supported AMPS production builds use C++17 compilers with __has_include.
// The fallback retains compatibility with older dependency-light compilers
// when their caller provides the documented SEP_COMMON_DIR include path.
#  define SRCSEP_SEP_COMMON_HEADER(name) <name>
#endif

#endif  // SRCSEP_UTIL_SEP_COMMON_HEADER_PATH_H
