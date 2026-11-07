#ifndef AMPS_INTERFACE_TEST_PIC_STUB_H
#define AMPS_INTERFACE_TEST_PIC_STUB_H

/* Keep production GeopackInterface.cpp on its non-debug dependency-light path. */
#define DIM 3
#define _PIC_DEBUGGER_MODE_ 0
#define _PIC_DEBUGGER_MODE_ON_ 1
#define _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ 0
#define _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_ 1

/* The test driver defines this overload as an exception-producing diagnostic. */
[[noreturn]] void exit(int line, const char* file, const char* message);

#endif
