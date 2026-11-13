# AMPS manual development status

Edition: 2026-11-10.

## Unconditional standard-url formatting — 2026-11-13

Supersedes the optional-xurl implementation recorded below. HEADER.tex now
always loads standard `url` and applies the same ASCII letter/digit wrapping
rules. There is no xurl load, availability check, or alternate URL-formatting
branch. Existing punctuation breaks, displayed text and link targets remain
intact. README.md describes the resulting dependency requirements.

Validation PASS:

- Root `make PDF`: 87 pages; xurl does not appear in the input recorder even
  though it is installed on this host.
- A relocated source-only manual compiled using real pdfLaTeX with xurl hidden
  from TEXINPUTS and standard url available; no system package was modified.
- `pdftotext -layout` output is byte-identical across the previous xurl-based
  PDF, the new ordinary build and the build with xurl unavailable. All retain
  the same one existing overfull box (4.52621 pt), with no additional overflow.
- `make -C doc check` and documentation whitespace checks pass.

The preceding temporary four-way comparison found 25 overflow warnings with
plain url, four with its hyphens option, and one with either xurl or the retained
project wrapping rules. This motivated keeping those rules unconditionally.
Evidence: `build-manual/url-only-validation.json`; dependency-absence transcript:
`/tmp/amps-url-only-611unmr6/without-xurl-build.log`. The comparison remains under
`/tmp/amps-url-layout-plcbeddw/`. These are local TeX tests, not remote macOS runs.
Source revision is 183c8b2514a1cc5c46af05740afbcf9bf5149d82 plus the uncommitted
manual changes. Application sources/configuration and both Makefiles remain
unchanged. No commit or push was performed. Earlier sections retain their
historical implementation and validation context.

## Optional xurl dependency — 2026-11-13

A macOS TeX installation reported `File xurl.sty not found` in HEADER.tex.
Made the extended URL-wrapping package optional with `\IfFileExists`.
If absent, the header uses the standard `url` package already needed by
hyperref and extends its break opportunities to ASCII letters/digits. Long
source paths and registry hashes still wrap without changing displayed
characters or hyperlink targets. Core TeX dependencies remain required;
no package installation or system-wide changes are performed by this fallback.
The chosen fallback is reported in the compiler log.

Validation PASS:

- Root `make PDF` with installed xurl: 87-page manual compiled successfully.
- A relocated source-only manual compiled with real pdfLaTeX under a restricted
  TEXINPUTS tree containing its actual package dependencies except xurl.
  `kpsewhich xurl.sty` finds nothing; `kpsewhich url.sty` succeeds. The log
  confirms `using url fallback`, and the input recorder shows no xurl loaded.
- Extracted text/layout (`pdftotext -layout`) is byte-identical between the
  installed-xurl and fallback PDFs. Both have 87 pages and retain the existing
  single overfull box; the fallback adds no layout warnings.

No system TeX file was renamed/removed. Tests are local dependency-absence
checks, not a claim of execution on the reported macOS host. Transfer the
updated `doc/HEADER.tex` to that checkout and rerun `make PDF`. Required versus
optional dependencies are clarified in README.md. Evidence is preserved in
`build-manual/xurl-compatibility-validation.json` and
`/tmp/amps-manual-xurl-6wa5wvuc/without-xurl-build.log`.
Source revision: 183c8b2514a1cc5c46af05740afbcf9bf5149d82 plus these uncommitted
documentation changes. Application code, configuration, both Makefiles and
scientific references are untouched; no commit or push was performed.

## Compiler portability — 2026-11-13

Compiler selection is evaluated on the receiving machine. Normal PATH tools
are invoked as `pdflatex` or `tectonic`, rather than using a stored compiler
installation prefix. `LATEX_ENGINE` and `--engine` accept names; explicit paths
remain optional. The receipt separately records `engine_command` and the
actual absolute executable/hash for provenance. It is never read back as
compiler configuration. Explicit/relative tool paths retain their original
format-selecting symlink name and are made absolute only to survive the
compiler's working-directory change.

Validation PASS:

- Root `make PDF`: invokes `pdflatex` by name and produces the 87-page manual.
- A source-only copy under `/tmp/amps-manual-portability-c86323ux/relocated manual`
  builds with `make PDF LATEX_ENGINE=pdflatex`. No local compiler, TeX cache,
  AMPS sources or generated build tree was copied into that manual directory.
- The same copied source builds from another working directory with
  `LATEX_ENGINE=pdflatex` and a different compiler prefix on PATH, containing
  spaces. The receipt identifies the changed prefix; no source was edited.
- `python3 -B doc/build_manual.py --engine tectonic --only-cached`: local
  fallback and offline build remain supported.
- `make -C doc check`: local include/asset/reference checks pass.

These are local relocation/PATH tests, not execution on a second physical
machine. The current pdfLaTeX build retains its pre-existing one overfull-box
warning; no manual TeX content was changed by the portability work. The root
and documentation Makefiles are unchanged. Evidence and exact compiler
identities are in `build-manual/portability-validation.json` and the compiler
logs; source revision is be7ede4cce60c1cd98950768a7f9a43714ce9bed plus these
uncommitted documentation changes. No application configuration/build/run,
regression reference change, commit or push was performed.

## PDF build repair — 2026-11-13

The root `make PDF` target failed because `doc/build_manual.py` was absent.
The working driver in `~/Mars2/AMPS/doc` was ignored there by the broad
`*build*` Git rule and consequently missing from this checkout. Restored the
source driver and added a narrow `.gitignore` exception. Generated TeX engines,
package/font caches, build products and partial PDFs are explicitly ignored.
The driver also resolves an explicit `tectonic` name to the same local compiler
as auto mode, so the documented cached build works without changing PATH.

This host's compiler and TeX cache were copied locally from the working Mars
checkout. The compiler matches its existing download receipt SHA-256
`a98aa59ad5c1df39a6c9e56cbfc5088f2b11d6c179c0130b97998e4bd46a46da`.
The original source/Makefile/PDF backups are under
`/tmp/amps-pdf-fix-0xzsq6zr`; local-copy provenance is in
`.build-tools/local-copy-receipt.json`. The maintained driver is visible to
ordinary Git source selection, while generated dependencies remain ignored.

Validation PASS:

- `make PDF` from the AMPS root: rebuilt `doc/AMPS.pdf`, 87 pages, zero TeX
  warnings, unresolved references or overfull/underfull boxes.
- `make -C doc check`: 26 reachable TeX sources with all local assets,
  unique labels and resolved citations/references.
- `python3 -B doc/build_manual.py --engine tectonic --only-cached`: successful
  offline compilation using the locally installed engine/package cache.
- Source inclusion/ignore checks and `git diff --check`.

The final PDF hash, engine hash and compiler transcript are retained in
`build-manual/build-receipt.json` and `build-manual/compiler.log`.
No manual science content, application source, AMPS configuration, root
Makefile or regression reference was changed. No application was compiled or
run, and no Git commit/push was performed. The manual edition below is unchanged.

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
