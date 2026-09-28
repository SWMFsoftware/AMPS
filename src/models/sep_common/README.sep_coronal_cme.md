# Neutral status support for `sep_coronal_cme`

`sep_status.h` is the dependency-free status/result contract required by the
shared coronal-CME model. It is owned by `sep_common`, not by `srcSEP3D`, so a
3-D provider and a field-line-bundle consumer can report the same typed state
without creating an application-to-application dependency.

The header contains no AMPS, PIC, MPI, output, `sep_coronal_cme`, or SWCME
include. `sep_status.cpp` is intentionally an archive-registration translation
unit; all small status helpers remain inline.
