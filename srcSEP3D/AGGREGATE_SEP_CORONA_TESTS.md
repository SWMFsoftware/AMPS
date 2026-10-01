# Complete shared-model plus native host test entry point

The direct `amps --test-suite sep-corona` command selects seven generic native
initialization checks. It does not select the separate 206 Stage-0--12
shared-model/architecture/documentation/offline preprocessing gates. This
source update adds one command that discovers and executes both registries:

```sh
python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0
```

Run once from the AMPS root, after overlaying the updated files and rebuilding
`amps`. The runner builds the shared test binary and adapters, runs all shared
gates, then launches MPI for the live native suite. Current discovery is 206
shared + 7 native = 213 canonical records. Future registry members are included
without naming IDs. This is an aggregate of two evidence scopes; shared tests
are not relabeled as live AMPS mover tests.

The native report now captures input schema, configured background/shock names
and actual completed steps. The supplied active-tube deck uses schema 4,
analytic Parker and SWCME. Its initialization passes do not establish an
installed PFSS/SCS or Stage-11 coronal runtime provider: the current host
selectors do not provide that integration. The existing corridor/solar-boundary
contracts are retained. A complete coronal host adapter is still required for
that production capability.

`test_output/coupled-sep-corona/summary.json` and `junit.xml` combine the
results. Fresh per-run directories retain original reports/logs/artifacts.
Missing/duplicate/extra/stale results, schema/provider/rank failures or exit
disagreement become ERROR. Native FAIL/SKIP remain explicit. Add
`--require-no-skips` to require every selected case to execute successfully.
Use `--list` for discovery, `--model-only` for an explicitly shared-only run,
`--launcher 'srun -n {ranks}'` for a suitable allocation, and `--jobs 16` for
shared compilation parallelism. Do not wrap this Python driver in `mpiexec`.

Verification: real shared aggregate 206/206; portable native C++ boundary/CLI
regression 45 checks plus provider/active-region JSON checks; application Python
runner regressions 26/26, including seven new process-level failure/discovery
tests. Synthetic executable/launcher fixtures test orchestration mechanics
only. Full linked AMPS/MPI execution was unavailable in this partial upload and
must be run on the configured installation.

The archive includes all previous Stage-11/12 sources. The new
`srcSEP3D/validation/aggregate_coupled_suite.patch` applies to the preceding
Stage-11/12 package. Use either the archive overlay or this incremental patch.
If starting from the earlier Stage-0--10 package, apply `stages11_12.patch`
first, then the aggregate patch. Generated binaries, run reports and caches
are excluded from the source archive.

Detailed runner/evidence documentation is in
`srcSEP3D/validation/COUPLED_SUITE_CLI.md`.
