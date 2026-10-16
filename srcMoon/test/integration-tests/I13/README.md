# I13 — controlled spherical-Moon baseline

I13 derives the exact finest receipt (60 km for the approved memory campaign,
5 km for legacy FULL) with the same transformation as I11 and requires its
hash to match the passing `I11/<finest>km/application.in`; it
never freezes the coarser human-authored template with finest-grid science.
The runner accepts the approved memory-profile settings or legacy FULL;
arbitrary CUSTOM/QUICK meshes are rejected. Missing/failed I11 evidence
stops I13 before any linked mesh allocation.
It first evaluates the registered production sphere against the independent
analytical area `4*pi*R^2`, with relative tolerance `1e-10`.  It also records
the realized AMR leaf/block counts, compiled species properties, particle
weight and time step, and the enabled Na source totals.

The complete roadmap gate additionally needs equilibrated surface-return flux,
profiles, column density, and synthetic brightness produced by a passing I11
campaign. Those values must be reviewed once and frozen as
`reference/I13_archived_baseline_metrics.csv`. The exact input hash, realized
mesh counts, compiled species properties, particle weights/time steps, and Na
source totals are frozen separately as
`reference/I13_archived_configuration.json`, since categorical species symbols
and hashes are not SI-valued science metrics. This runner never creates or
overwrites either reference. After I11 passes, it writes
`I13_candidate_baseline_metrics.csv` and `I13_candidate_configuration.json`
for independent review. Until both candidates are reviewed and installed at
their reference paths, a successful sphere-area check is `SKIPPED`, not a
baseline PASS. Once frozen, deterministic configuration must match exactly;
each sampled metric is compared using its I11 numerical tolerance plus current
and reference 95% relative CI half-widths.

Run from the repository root with a clean SPICE Moon executable:

```sh
python3 srcMoon/test/integration-tests/I13/test.py \
  --amps ./amps \
  --input srcMoon/examples/m1_memory200.in \
  --output-dir test_output/srcMoon/I13
```

The output directory contains the raw receipt, logs, manifest, `result.json`,
and `process_resource_metrics.csv` with one-second process-tree RSS and
host-memory observations. Monitoring imposes no time or memory cap. Sphere mode
and a complete `.in` receipt are mandatory.
