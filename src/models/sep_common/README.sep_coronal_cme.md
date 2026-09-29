# Neutral status support for `sep_coronal_cme`

`sep_status.h` is the dependency-free status/result contract required by the
shared coronal-CME model. It is owned by `sep_common`, not by `srcSEP3D`, so a
3-D provider and a field-line-bundle consumer can report the same typed state
without creating an application-to-application dependency.

The header contains no AMPS, PIC, MPI, output, `sep_coronal_cme`, or SWCME
include. `sep_status.cpp` is intentionally an archive-registration translation
unit; all small status helpers remain inline.

Stage 10 adds two further neutral interfaces:

- `sep_field_line_exchange.{h,cpp}` owns SI-valued `NodeState`, front,
  observer-mapping, connection-summary, line, and bundle records. It validates
  stable identities, strictly outward arc length, open/smooth topology,
  positive physical primitives, flux-derived measures, and one-sided bounded
  interpolation. Categorical topology is copied, never reconstructed from
  interpolated floating-point columns.
- `sep_field_line_bundle_io.{h,cpp}` writes a canonical JSON manifest and one
  tabular member per line. Each member has a SHA-256 checksum; bundle identity
  covers all canonical metadata and member bytes. A complete sibling temporary
  directory is atomically renamed only after every member and manifest write
  succeeds. Import verifies schema, safe IDs, hashes, member identities, and
  the recomputed bundle identity before returning a value.

These sources likewise have no upward dependency on `sep_coronal_cme` or an
application. Consequently `srcSEP` can read a bundle and publish it through
its existing background snapshot store without linking the 3-D model.
