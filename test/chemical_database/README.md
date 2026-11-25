# Shared chemical database software checks

This is the single global runner for the shared C++/Python data-access service.
It performs no AMPS configuration, particle simulation or observational
validation. Read [CONTRACT.md](CONTRACT.md) for the preimplementation inventory,
scientific boundaries and exact software acceptance rules.

From the repository root, choose a new output directory:

```bash
/usr/bin/python3 -B test/chemical_database/run_all_tests.py \
  --output-dir test_output/chemical-database/my-new-run \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/database-lock.json
```

The optional `--database-lock` selects a real prepared release for additional
Python/native exact-source agreement and actual installed CLI/C++ SDK checks.
The suite copies that release (never hardlinks), relocates the copy, runs its
own `bin/chemdb` and compiles all six shared CPPs from its own `sdk/` tree
without AMPS checkout includes/libraries. `installed-deployment/result.json`
records those deployed source/executable/lock identities. The original release
and canonical sources remain read-only. Omit the lock to run the native and
Python administrative metadata/integrity controls; dependent real-release
checks report SKIPPED. Use the actual
installation lock filename written by the builder. No latest-release search,
implicit download or alternative source is used. Python3 and a C++17 compiler
(`g++`, or the compiler argv named by `CXX`) are required. No MPI/PIC/SPICE is used.
The CORE adapter explicitly passes `--protocol-only` for the site-independent
software selection. It preserves its configured compiler argv through
`--cxx-json`; there is no site-data prerequisite hidden in that CORE case.

Nightly execution is owned by `CORE-CHEMICAL-DATABASE-001` in
[`src/test/manifest.json`](../../src/test/manifest.json), which belongs to the
`nightly` campaign. The `core` entry in
[`utility/nightly/runners.json`](../../utility/nightly/runners.json) selects
that campaign. The core report has one CORE suite/case; its nested
`artifacts/CORE-CHEMICAL-DATABASE-001/shared-protocol/result.json` and suite log
retain every child method and its counts. These checks are separate from the
26 lunar U-series rows displayed by the nightly `moon` runner.

On the nightly web page (since 2026-11-25) the suite and its checks form the
group "Core: chemical database":

- The suite row keeps the CORE verdict.
- Each executed child method is a row with its own status, reason and log.
- Each `outside_selected_scope` ID is a SKIPPED row saying that the nightly's
  `--protocol-only` selection leaves it out.
- A missing or inconsistent receipt shows as a FAILED row
  `CORE-CHEMICAL-DATABASE-001 checks`.

The nightly reader is `utility/nightly/amps_nightly/readers.py`
(`NESTED_CHECKS`). These rows only display the receipt: they add no CORE case
and no test-count contract.

To repeat the actual nightly adapter without running other CORE suites:

```bash
/usr/bin/python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --campaign nightly --test CORE-CHEMICAL-DATABASE-001 \
  --output-dir test_output/chemical_database/NEW_NIGHTLY_ADAPTER_RUN
```

`--protocol-only` deliberately excludes classes marked
`requires_authoritative_source=True` and the four `PythonProduction` methods.
It therefore does **not** run the explicit real-source, pinned-release,
installed/relocated SDK, or full acquired-release audits described below.
Excluded source-method IDs, and since 2026-11-25 the four `PythonProduction`
method IDs, appear in `outside_selected_scope`, rather than being reported as
executed or passed. Supply the documented original-file
selectors and lock to the full runner for that separate scope. The scheduled
nightly tests a committed `master` checkout; uncommitted tests in another
working tree are not included until published. The November 23 audit and
actual historical/current counts are recorded in `srcMoon/STATUS.md` section
5.42; those counts are dated evidence, not a fixed test-count contract.

The native controls compile the actual maintained `Runtime.cpp` and `DataIO.cpp`
and call them through a thin read-only adapter. Administrative fixtures contain
only labelled parser strings, never observations or numerical chemistry data.
They verify exact SHA-256 binding, published empty/`abc` SHA vectors,
schema/codec ownership, ordered identifiers, quoted/newline CSV strings and
exact padded column labels, capability and permission refusal, unavailable
metadata, resource bounds, path and symlink confinement, immutable shared
lifetime and concurrent readers, null/missing-output guards, failed-reload
rollback, row counts and relocation. Compile-time type identities prove that
the historical Moon namespaces alias the shared services. Real installed and
canonical source bytes are read-only. Altered fixtures remain inside this run.

`result.json` contains every evaluated test or fixture outcome and the executed,
PASS, FAIL, SKIPPED and ERROR totals. Every failed/error entry names the test,
runner and exact log path. `unittest.log`, `build.json`, compiler stdout/stderr
and `queries/query-*/{stdout.log,stderr.log,result.json}` preserve the actual
commands and responses. Prior outputs are never overwritten. Exit codes are
0=PASS,1=FAIL,2=ERROR,77=SKIPPED. A prerequisite fixture outcome is counted as
such; it does not pretend its dependent methods were executed.

When adding a new software test, add its method to `test_protocol.py` or a
maintained `test_*.py` module in this directory; the runner discovers its
`unittest.TestCase` methods automatically and binds their source hashes. Keep native
adapters thin and update the contract before implementation. Existing Moon
U15/U16/U17 numerical/statistical runners retain their original independent
oracles and acceptance criteria; the shared-service move changes their source
link lists and hashes, not their physical scope. No permanent U/I ID is added.

The November 19 source-acquisition checks cover strict packet admission,
NASA CEA fixed-width records and domains, IAEA native printed fit coefficients,
ExoMol source projections and TU Wien native CSV members. Administrative parser
fixtures are labelled as such. Real-source tests read the explicit selectors
below; absent inputs are SKIPPED and never substituted with generated science:

```bash
AMPS_CEA_SOURCE=/path/to/acquisition/raw/nasa-cea-data-thermo.inp \
AMPS_EXOMOL_PACKET=/path/to/exomol-native \
AMPS_SURFACE_SOURCE_ZIP=/path/to/As_published_in_Comm_Earth_Env.zip \
/usr/bin/python3 -B test/chemical_database/run_all_tests.py \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2/database-lock.json \
  --output-dir test_output/chemical_database/NEW_REAL_RELEASE_RUN
```

[`verify_acquired_release.py`](ACQUIRED_RELEASE_QA.md) separately checks the actually published new tables
against their acquired originals, preserving a result/manifest and per-check
logs. Invoke its `--help` for the required packet and lock arguments. This is
source reproduction and deployed API verification, not linked I21–I24 or
observational validation. Both old and new locks can be tested independently;
no test changes a selected source, threshold or default installation.

The r3 follow-up adds `test_native_chianti.py`, `test_native_apollo_tpd.py`,
`test_native_workbooks.py`, `test_model_definitions.py` and
`test_reduced_models.py` to the same automatic discovery. The first three keep
their actual-source classes marked `requires_authoritative_source=True`:
missing originals are SKIPPED in the full suite and explicitly outside the
CORE protocol selection. Generic parser controls are labelled administrative,
and analytic/conservation tests exercise the actual shared C++ helpers. Model
metadata controls prove old-lock compatibility and that definitions grant no
physical capability. No lunar U/I ID or acceptance threshold is changed.

Select CHIANTI originals with `AMPS_CHIANTI_LEVELS=/path/to/o_1.elvlc` and
`AMPS_CHIANTI_RADIATIVE=/path/to/o_1.wgfa`, and Apollo originals with
`AMPS_APOLLO_TPD_SOURCE_DIR=/path/to/jones2020-native`; the existing CEA/ExoMol
selectors above remain applicable. Workbook real-source controls require the
exact archived ZIPs documented in `NATIVE_WORKBOOK_SOURCES.md`. The global
runner records all actual selections and distinguishes source access from
physical agreement. Choose r3's own lock to exercise its contents.

[`verify_physical_model_release.py`](PHYSICAL_MODEL_RELEASE_QA.md) checks an
explicit r3-style release and selected packet family through Python and C++,
including exact rows, recipe members, held-source refusal and model-registry
agreement. Use `--help`; supply the compiled shared `catalogue-probe` with
`--native-probe`, the actual lock and repeated `--source-packet` arguments.
Its output directory must be new. This is shared source/API verification, not
formal linked I21–I24 qualification or observational validation.

[REDUCED_MODEL_NUMERICS.md](REDUCED_MODEL_NUMERICS.md) records the independent
Decimal/Fraction oracles, retained development failures and the proof behind
the binary floating-point zero-COM certificate. Mathematical roundoff criteria
are distinct from physical approximation or observational agreement limits.

The r4 physical-input follow-up adds automatically discovered
`test_new_molecular_exomol.py`, `test_new_surface_apollo15221.py`,
`test_new_surface_workbooks.py`, `test_new_surface_relaxation.py` and
`test_amend_exomol_metadata.py`. Administrative controls test exact tokens,
missing/duplicate metadata, dependency hashes, source gaps and refusal paths.
They are not observations. Molecular real-source controls require
`AMPS_NEW_MOLECULAR_SOURCE_ROOT` to name the acquired ExoMol directory. The
metadata-amendment real-source controls require
`AMPS_EXOMOL_AMENDMENT_OLD_PACKET` and
`AMPS_EXOMOL_AMENDMENT_EVIDENCE_PACKET`; these tests authenticate their frozen
packet identities before copying any science values. Missing actual inputs
remain SKIPPED in the full suite and outside CORE's protocol-only scope.

Example from the retained work package, using a new output directory:

```bash
AMPS_NEW_MOLECULAR_SOURCE_ROOT="$PWD/test_output/chemical_database/m3-physical-inputs-20261119T231417Z/molecular/database-originals-network/exomol" \
AMPS_EXOMOL_AMENDMENT_OLD_PACKET="$PWD/test_output/chemical_database/physical-input-integration-20261119T062501Z/exomol-native/source_packet_v2.json" \
AMPS_EXOMOL_AMENDMENT_EVIDENCE_PACKET="$PWD/test_output/chemical_database/m3-physical-inputs-20261119T231417Z/molecular/source_packet.json" \
/usr/bin/python3 -B test/chemical_database/run_all_tests.py \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r4/database-lock.json \
  --output-dir test_output/chemical_database/NEW_R4_RUN
```

Set the earlier CEA/CHIANTI/ExoMol selectors too when their real-source controls
are desired. `verify_physical_model_release.py` additionally verifies all
selected r4 packets and all actual available CSV cells through Python/C++,
including complete packet/recipe/dependency archives and held-source refusals.
Source-parser QA independently compares new author files with lexical views;
release QA compares published views with those originals. These are separate
layers, and neither substitutes for linked physical verification/validation.

The targeted2026-11-20 follow-up adds source-specific admission controls to
`test_source_packets.py`. The Figshare controls require the exact reviewed
deposit/version/licence/listed file and refuse unrelated or missing grants.
The NIST controls confine reuse to JRes121.005 and its native technical-series
notice. Their small administrative records exercise permission/identity
software only; actual native-source QA is recorded separately. The normal
automatic discovery includes these tests in full and CORE protocol scopes.

`test_native_nist_hxoy.py` also joins automatic discovery. Four administrative
parser controls test exact cell/spans/footnotes, malformed/truncated markup,
hash/article identity and preservation of existing outputs. The separate
authoritative-source class compares the actual article/CSV cells and declares
`requires_authoritative_source=True`. Select it explicitly with
`AMPS_NIST_HXOY_SOURCE=/path/to/PMC7339710.html` and
`AMPS_NIST_HXOY_SHA256=EXPECTED_SHA256`; otherwise it is SKIPPED in the full
suite and outside CORE's protocol-only scope. It preserves source units and
uncertainty strings and performs no thermochemical state/energy selection.

The November 24 helium archive adds a nightly packet-admission control for
unreviewed deposit bytes. Three opt-in tests verify the actual snapshot,
including all 507,648 source rows against original ZIP members, the exact
native C++/Python energy index, the live IAEA byte identity and the ten held
source permissions. Select the real packet and snapshot explicitly.

```bash
AMPS_M3_ARCHIVE_PACKET=/path/to/packet/source_packet.json \
python3 -B test/chemical_database/run_all_tests.py \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r6/database-lock.json \
  --output-dir test_output/chemical_database/NEW_ARCHIVE_TEST_RUN
```

The [source contract](../../utility/chemical_database/NATIVE_M3_ARCHIVE_SOURCES.md)
describes the producer references and the difference between exact table
reproduction and physical model qualification.

The nightly protocol also checks the reviewed article and permission
binding for the r7 publication. For the three independent sources in
the published r7 tables, set `AMPS_CHECK_M3_LITERATURE_RELEASE=1` and use
the r7 lock with the same shared runner. The tests compare published
tokens, preserve the missing helium measurement, exercise the actual C++
reader, and check that metadata-only records grant no payload access.

The r8 source checks run against the explicit r8 lock with
`AMPS_CHECK_M3_PARTIAL_RELEASE=1`. They compare every projected cross
section and blank against the original XLSX cells, verify the three published
He tables, exercise the C++ catalogue, and check held-payload behavior.
The nightly protocol verifies the new version and rights pinning policies.

For r9, set `AMPS_CHECK_CEX_SCATTERING_RELEASE=1` with the r9 lock to
verify the unchanged D119990 CSV, source-review link, PDF permission
hold and exact C++/Python definition records. The source-independent
nightly protocol checks the two new definitions and their activation
flags through the model API.

Set `AMPS_CHECK_M3_LAB_RELEASE=1` with the r10 lock to compare the
new OH and He rows with the printed tables and check the native
C++ view and held-source behavior. The nightly protocol verifies
the new publication and source-specific rights pins.

The r11 installed-source checks select
`AMPS_CHECK_M3_PHOTON_BENCHMARKS=1`, `AMPS_CHECK_OH_ICE_RELEASE=1`,
`AMPS_CHECK_M3_HELIUM_EMISSION=1` and `AMPS_CHECK_SODIUM_SOURCE_RELEASE=1`
with the explicit r11 lock. They check primary numerical cells, figure
lineage, source permissions, relative versus absolute quantities, and
shared-language table access. The nightly core database case exercises the
source-independent protocol. Retain the real-source suite result separately.

Select `AMPS_CHECK_WATER_YIELDS_RELEASE=1` to check the three
21.2 eV water fragment signals, denominator, ion-collection bias and
held proceedings metadata. [Water benchmark](../../utility/chemical_database/NATIVE_M3_WATER_BRANCH_SOURCES.md)

For the r12 real-source suite, use its explicit lock and set
`AMPS_CHECK_M3_PICKUP_PRECEDENT=1`, `AMPS_CHECK_WATER_PRECEDENT_RELEASE=1`,
and `AMPS_CHECK_LUNAR_SURFACE_PRECEDENT_RELEASE=1` along with the earlier
release selectors. Test discovery includes all three modules. They check
the deposited 1 au frequencies, conditional water branches, identified
source algorithms, and separate OH sticking, residence and launch
assumptions. The core nightly protocol runs without external acquisitions.
[Inputs, scope and query examples](../../utility/chemical_database/DATABASE_README.md#18-reproduce-published-exosphere-model-choices-in-r12)

On the r13 lock, add `AMPS_CHECK_M3_LUNAR_NOBLE_GAS=1`,
`AMPS_CHECK_LUNAR_WATER_COVERAGE_RELEASE=1` and
`AMPS_CHECK_EXOSPHID_ATOMS_RELEASE=1` to the previous release selectors.
The real-source checks verify the Table 1 rates, conversion to per-second
frequencies, coverage and light conditions of the water model, and the
thresholds and product labels in the atomic code. The nightly protocol
also rejects altered inputs presented under the Peschel publication's
permission policy. [r13 model sources](../../utility/chemical_database/DATABASE_README.md#19-lunar-neutral-losses-water-coverage-and-noble-gas-product-assumptions)

### Sodium pickup and ground-state capture scenarios (r14)

The publication test exercises the two named database models through the Python
reader and through the C++ CSV reader. The independent oracles are Poppe Table 1
(for the three Na median frequencies, rounded to paper precision), the deposited
Na row (for full precision), and Shingal Table 4 (4.0 keV/u and 0.22 Å²).
The tests require zero lunar-frame birth velocity for the bulk source, and the
pair-COM isotropic law and two-body recoil for the ground-state case. A release
omitting a case's native sources must omit that case. The nightly protocol checks
source selection without external files. For the four source-backed checks supply the
r14 lock and set `AMPS_CHECK_NA_REFERENCE_MODELS=1` in the shared suite. These
checks verify data access and the physical meaning of the inputs. Linked ion
trajectories, recoil spectra and return fractions require separate simulations.

From the AMPS root, the following runs the shared suite and requires both sodium
cases in the selected release. Other source-specific checks depend on the
acquisitions present and the corresponding release selections.

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
work=$(mktemp -d)
AMPS_CHECK_NA_REFERENCE_MODELS=1 python3 -B test/chemical_database/run_all_tests.py \
  --database-lock "$db/database-lock.json" --output-dir "$work/suite"
```
