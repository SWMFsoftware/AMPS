# AMPS manual development status

Edition: 2026-11-10.

## Performance-test registry refresh

Refreshed the modern suite appendix from the expanded registry after adding
18 independent fast-wave optimization cases: 56 suites / 1,192 cases.
The manual remains 87 pages and its current compilation passes. This table
records registration rather than scientific qualification. Detailed algorithm
contracts, measurements, core revalidation and known limitations remain in
test/srcFastWave/STATUS.md and the adjacent core READMEs. The legacy
application catalog is unchanged and its source check still passes.

## Fast-wave performance controls documented

Added the measured fast-wave-only GNU optimization fallback and opt-out to
the individual regression chapter. Refreshed the legacy catalog source
hashes and phase locations after the maintained Table/generated recipe
change. Root make PDF and both catalog checks PASS; the PDF remains 87 pages
and compiles without warnings. Native performance validation (four passing
runs, unchanged bounds and byte-identical sampled output) is documented in
test/srcFastWave/README.md and src/test/STATUS.md, with the machine-readable
receipt at test_output/fast-wave-performance-20261110/summary.json.

## Architecture, application regressions and source-guide update

Current PDF: 87 pages, 20 chapters and three appendices. Added the layered
architecture illustration and its local vector PDF/PNG plus commented
generator. The architecture chapter explains interfaces, lifecycle and
ownership without prescribing a universal operator order.

Added the individual application-regression chapter, further end-to-end
examples, and a catalog derived from MakefileTest/Table and Makefile.test:
153 declarations / 150 distinct targets, including all 100 currently generated
host targets and the 50 additional declared names. Each generated entry
records its purpose, source locations, comparison method, explicit tolerance
overrides and detected prerequisites. The companion JSON and optional editorial
generator preserve the distinction between registration and successful runs.
Removed the repository instruction-file reference from the manual's source
guide; verified its filename is absent from all reachable TeX sources and
the rendered PDF. The repository instruction file itself remains in place.

Validation PASS:

- Exact root command make PDF; zero compiler warnings, unresolved references
  or overflowing/underfull boxes.
- Local structure check: 26 reachable TeX files, all local assets available.
- Both catalog --check commands; modern snapshot remains 55 suites / 1174 cases.
- Independent reconciliation against the generated compile-target names:
  all 100 owners cataloged, all 100 target names present in the PDF.
- Both legacy authority files retain the cataloged SHA-256 hashes.
- Python syntax and documentation whitespace checks.
- Visual inspection of architecture, application chapter and catalog pages
  (PDF pages 16, 50 and 74).

Receipt: build-manual/legacy-documentation-validation.json; build logs and
PDF hash are also retained in build-manual/. The generated legacy test_help
command has a pre-existing shell quoting failure (exit 2); this is documented,
not repaired or reported as passing. No application simulation, native build,
configuration or Git-writing operation was executed for this update.

The sections below retain the earlier 65-page compilation/rename checkpoints;
their page counts and archive inventories describe those earlier snapshots.

## Manual source filename cleanup

Removed the AMPS_ prefix from all 22 chapter/appendix/bibliography TeX source
filenames. Updated the main input list, suite-appendix generator and README.
AMPS.tex remains the main document; the supplied logo retains its requested
filename. No chapter content was changed by the rename.

Validation PASS: local include/asset/reference checks, exact suite snapshot,
and root make PDF. The rebuilt manual remains 65 pages, with zero TeX warnings
or overflowing boxes and byte-identical extracted text compared with the
pre-rename PDF. Every renamed file has the same content hash as before.
The complete mapping and verification receipt are in
build-manual/source-renames.json. No Git state was modified.

## Implementation

- Rebuilt AMPS.tex as the main user/developer manual with the supplied
  AMPS_clean_sweep_white.jpg on the title page.
- Added chapters spanning core/applications, configuration/runtime input,
  numerics/weights/sampling, mesh/movers, collisions/internal energy,
  photochemistry, dust, surfaces/arbitrary shapes, planetary/comet and
  heliospheric physics, SWMF, output/restart, testing, development and diagnosis.
- Added end-to-end test recipes with actual registered owners, independent
  acceptance criteria, outputs and explicit qualification boundaries.
- Added a generated complete suite-owner appendix; the individual global
  catalog remains the exact-case description authority.
- Replaced external shared make/style dependencies with a local Makefile,
  standard-library Python build/check tools, ordinary TeX packages and embedded
  bibliography. PDF and MANUAL aliases are available in this directory;
  the existing root make PDF delegates here.
- Preserved the four former operator chapters under legacy/ and retained
  the longer specialist research notes unchanged.

## Validation

Completed successfully with the exact command **make PDF**, executed from
/home/vtenishe/Mars2/AMPS. Output: doc/AMPS.pdf, 65 pages, containing 19 chapters,
two appendices, embedded bibliography, supplied title-page logo and an
analytical reference figure. The figure is labelled as analytical teaching
material and is not AMPS simulation output.

- Tectonic 0.17.0 compilation PASS: no TeX warnings, undefined citations or
  references, or overflowing boxes.
- Local structure check PASS: 24 reachable TeX files, local graphics,
  unique labels and resolved citations/references.
- Suite snapshot check PASS: 55 suites / 1174 exact registered cases.
- Standalone-copy build PASS: copied 36 local manual/support/graphic files
  and supplied a private offline standard TeX cache; no AMPS, SWMF, application,
  input or test source tree was copied. Normal compilation does not need the
  surrounding repository.
- Rendered title, analytical-reference figure, dust equations, end-to-end
  recipes and suite-table pages were visually inspected (PDF pages
  1, 21, 31, 48 and 59).
- Source-guide paths and Python syntax checks PASS. Original specialist
  notes/bibliography are byte-identical to the tracked versions; the four
  former operator chapters are preserved byte-for-byte under legacy/.
- Supplied logo SHA-256:
  45e0abc6da6c5746ca7e1784ee7a0a64b00edb55386b38828049d0ef54bd824a.

Receipts/logs: build-manual/build-receipt.json, validation-receipt.json,
AMPS.log, compiler.log, standalone-source-inventory.json and the standalone/
build's own build-manual/build-receipt.json. The optional compiler download
has its official release URL/archive/binary hashes in
.build-tools/download-receipt.json. Package/compiler caches are generated
tools, not dependencies on the AMPS source tree.

An initial private compiler preflight required correcting its relative HEADER
path, and the optional figure generator was adjusted to the installed older
matplotlib API. Final cached and root make builds pass. No native AMPS rebuild
or physics test was needed for this documentation-only work.

## Evidence boundaries

The core 2026-11-10 result quoted in the manual is retained historical
evidence, not a fresh run. Full-CME G1, live SWMF fixtures, some application
science campaigns and strict complete coverage remain open. Source availability
and registration are never relabelled as qualification. No Git state is
modified by this manual task.
