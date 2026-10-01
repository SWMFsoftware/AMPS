# Complete coupled SEP+corona software test entry point

The current shared registry contains **222** canonical cases through Stage 14.
The linked AMPS `sep-corona` suite still contains **seven** generic native
initialization cases. Their complete aggregate contains **229** records:

```sh
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0
```

Run from the configured AMPS root after overlaying these sources and rebuilding
AMPS. The driver discovers both registries, builds shared tests/adapters, runs
all shared cases, then launches the native MPI suite. Current/future cases do
not need explicit names. Do not launch the Python driver itself under mpiexec.
For the portable subset use `--model-only`; for scheduler launchers use e.g.
`--launcher 'srun -n {ranks}'`. Discovery is available through `--list`.

The baseline Stage-13 gate has 209 shared cases. Stage 14 adds thirteen
independently selected research/protocol verification IDs. The new source
adds release profiles/checksummed packaging, explicit cross-application report
ownership, guarded research kernels, immutable offline products, and detailed
comments/READMEs. `make -C src/models/sep_coronal_cme test-stage13` runs the
baseline release machinery; `test-stage14` and `test` run all 222 shared cases.

Every evidence scope remains explicit. The active-tube deck selects analytic
Parker/SWCME (application schema 4). Its seven host passes do not establish a
PFSS/SCS coronal provider, production Stage-14 movers or observed-event
validation. Shared JSON and aggregate reports explicitly leave production and
observational qualification false. Stage-13 release qualification requires
actual clean owning-app/core builds, MPI/convergence/cross-model evidence,
D1--D10, validity/limitations and any preregistered event/campaign-C evidence.
The partial source upload has srcSEP adapters but not its full application
build/test tree; missing host/campaign evidence is not replaced by mock passes.
The research guide lists implemented domains and outstanding scientific gates.

The outside-Sun active corridor and solar absorbing boundary remain as before.
New research JSON uses schema 6; momentum-dependent families use bundle major
4. They do not add unsupported selectors to the current schema-5 input or
silently change baseline fixed-length first passage, shock calibration or
passive waves. No production last-pass record is fabricated.

Live progress includes phase start/end, completed counts, percentages, elapsed
time and quiet-phase heartbeats. The final failure summary gives each failing
ID, its diagnosis and diagnostic/execution log paths. Fresh per-run directories
retain all child logs and original reports. The latest JSON/JUnit and
`failures.txt` are refreshed; older runs remain available. Missing, duplicate,
extra, stale or exit-inconsistent reports become errors. Use
`--require-no-skips` when every selected native case must execute.

The ARCHSCCM01 source audit retains the earlier UTF-8 object fix: only genuine
C/C++ sources/headers enter text decoding. Native archive symbols and a public
C++17 consumer are still checked. New `.d` dependency files ensure public
header changes rebuild affected model/neutral/test/adapter objects. When
replacing an older pre-dependency build, clean the shared model once:

```sh
make -C src/models/sep_coronal_cme clean
```

Verification for this update: clean portable shared build and 222 selected cases: 219 PASS, 3 campaign SKIP, 0 FAIL;
application Python regressions 30/30; thirteen new research IDs include
manufactured independent numerical references and negative provenance/domain
fixtures. Complete linked AMPS/MPI execution and real observation campaigns
were unavailable in this partial upload and remain required on the configured
installation. Manufactured EVT/XMD/SLM tests verify protocols, not observations.

The source archive is an overlay, not a full AMPS checkout. It contains all
previous sources and incremental patches, the new `stages13_14.patch`, and an
embedded `RELEASE_SOURCE_MANIFEST.json` with required source/build/test/example
paths and SHA-256 hashes. Compiled objects, archives, executables, generated run
reports and caches are excluded. Use either the complete archive overlay or
apply `stages13_14.patch` to the immediately preceding ARCHSCCM01-fix source
package. Older updates retain their historical patch order: stages11_12,
aggregate suite, progress, failure summary, then architecture source-scan fix.

Details:

- `src/models/sep_coronal_cme/docs/STAGE13_RELEASE_QUALIFICATION.md`
- `src/models/sep_coronal_cme/docs/STAGE14_RESEARCH_EXTENSIONS.md`
- `src/models/sep_coronal_cme/examples/stage14/README.md`
- `srcSEP3D/validation/COUPLED_SUITE_CLI.md`

Campaign dispositions: EVT3D01, XMD3D01 and SLM3D01 run their synthetic
contract checks successfully, then report SKIP for absent actual campaign
evidence. Their rows retain `verification_passed=true` and a reason. This is
not a failed physics check. `--require-no-skips` makes an incomplete campaign
a nonzero exit. If all seven native host checks pass, the current aggregate
contains 229 records with 226 PASS and these three SKIP.
