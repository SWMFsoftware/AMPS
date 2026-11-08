# AMPS persistent-format policy

The machine-readable inventory is
[`src/test/persistent_formats.json`](../test/persistent_formats.json). Every
core binary family found in `src/pic`, `src/meshAMR`, and the reusable stack
serializers is listed there with its metadata and compatibility boundary.

Version 1 envelope
------------------

`pic_persistent_format.h` defines a 64-byte, canonical little-endian envelope:

- `AMPSPF1` magic, envelope version, and typed format kind;
- compatibility flags (`NativeLayout` or `CanonicalLittleEndian`);
- writer byte order, pointer width, and C/C++ integer/floating widths;
- a caller-owned 64-bit generated-layout identity;
- payload byte length and CRC-32, plus a header CRC-32.

Native payloads are accepted only when widths, byte order, kind, and layout
match. Canonical payloads use explicit fixed-width little-endian fields and can
cross builds after the same validation. Publication uses a sibling temporary
file, `fflush`, `fsync`, and atomic rename. An interrupted temporary file never
replaces the last complete destination.

The sampling and particle restart writers now publish this envelope. Their
payloads remain native AMPS layouts, so cross-build reads are deliberately
rejected. Existing unversioned files can still be read through the same-build
legacy path, with a warning. `PIC::Restart::AllowLegacyRestartFiles` defaults
to true to preserve established applications during migration and can be set
false before recovery when a run requires the corruption-detecting envelope.
Legacy files do not count as portable evidence.

The particle native-layout identity is at payload revision 2. In the supported
species-global one-rank MPI profile it serializes the global statistical-
weight table, global SI time-step table and initialization flag, scalar RNG
seed, AMR block/cell particle lists, and each complete generated particle
record. Complete records include configured internal energy, chemistry/
tracker/boundary fields and application extension bytes. Sampling restart
separately carries the sampling counters and complete center/corner associated
data. The native test perturbs the in-memory tables and deletes its source
particle before reading, so recovery cannot pass by observing pre-save state.

Hybrid per-thread RNG arrays, local time-step/weight modes not already present
in associated block state, changed-rank ownership, and user callback schema
identity do not yet have a portable payload contract. Applications using those
capabilities must not infer deterministic continuation from the one-rank MPI
receipt.

Legacy boundary
---------------

AMR mesh, raw background/associated-node data, trajectory scratch/dataset,
stack, and cut-cell scratch formats remain inventoried legacy formats. They are
not silently described as portable. Mesh section markers or record-size magic
provide useful corruption hints but are not substitutes for a version, byte
order, width table, and complete layout identity. T9 tests the new envelope and
pins those legacy policies; migrating each legacy producer is follow-up work
that requires its application-level continuation fixtures.

Run the T9 suites from the repository root:

```sh
python3 src/test/run_tests.py --test CORE-PIC-RESTART-001
python3 src/test/run_tests.py --test CORE-PIC-DATAFILE-001
python3 src/test/run_tests.py --test CORE-PIC-OUTPUT-001
python3 src/test/run_tests.py --test CORE-PIC-RESTART-NATIVE-001
```
