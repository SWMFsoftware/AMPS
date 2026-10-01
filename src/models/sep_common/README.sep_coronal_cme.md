# Neutral status support for `sep_coronal_cme`

`sep_status.h` is the dependency-free status/result contract required by the
shared coronal-CME model. It is owned by `sep_common`, not by `srcSEP3D`, so a
3-D provider and a field-line-bundle consumer can report the same typed state
without creating an application-to-application dependency.

The header contains no AMPS, PIC, MPI, output, `sep_coronal_cme`, or SWCME
include. `sep_status.cpp` is intentionally an archive-registration translation
unit; all small status helpers remain inline.

Stage 14 adds `sep_coherent_transport.h/.cpp` here because coherent transport
is neutral and must be shared by both applications. The public header includes
only neutral status/STL types. Parker antisymmetric-operator and focused
full-characteristic ownership are distinct. The implemented stationary
inertial Hamiltonian has explicit invalid-region and magnetization guards;
moving-frame/plasma-frame terms remain unsupported. The shared model makefile
compiles its neutral object, and the srcSEP3D application archive enumerates
that object directly. A full srcSEP build must register the same canonical
source in its owning neutral build; its makefile was absent from this upload.
See the [research guide](../sep_coronal_cme/docs/STAGE14_RESEARCH_EXTENSIONS.md).
