# I19 — Surface inventory and global mass-conservation accounting

## Supporting development before M1 closure

The user-approved November 12 development-order exception permits a separate
supporting run while formal I19 remains predecessor-gated:

```bash
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

This builds the dedicated Na/Ar executable, enabling individual weights for
partial finite-reservoir packets and compiling chemistry OFF. It runs local
U08-U13, read-only package QA, supporting I14/I18 adapters and the supporting
I19 budget fixture, then restores canonical Moon/SPICE. The normal formal
I14-I20 runners are also invoked to preserve their genuine I13-blocked SKIPPED
results. Supporting PASSes have their own counts and never become formal
I19 or M2 PASS. An existing M2 build may be explicitly supplied with
`--m2-amps /absolute/path/to/amps-moon-m2-na-ar` and matching generated tree.

`MoonM2Budget` observes the same production source, surface and mover paths
as ordinary runs. It independently sums owned gas particles, replicated
surface stocks counted once, and local pending adsorption-minus-release
fluxes. Surface deposit is distinct from whole-particle gas adsorption;
uncredited adsorption is a physical surface sink. External injection equals
accepted emission minus surface release. Signed prescribed replenishment is
separate from gas injection. Gas, surface, combined atom and mass identities
are retained by species and independently recomputed with Decimal arithmetic.
The floating-point error bound includes the actual double surface-array
updates as well as long-double sums; the 0.5% source criterion remains intact.

The fixture synchronizes both species to the smaller initialized global step,
performs one initialization step, adds 64 controlled incoming particles and
one outward high-energy particle per species, and seeds 32 model-weight
equivalents into one logical surface element before the initial ledger row.
These are explicit verification initial conditions, never observations or
production equilibrium claims. Twenty ordinary steps exercise actual sources,
adsorption, release and escape. The predeclared I09 classifier uses
`v_inertial=v_LSO+Omega_LSO cross r_LSO` and specific two-body lunar energy:
positive-energy exits are lunar escape; bound finite-box exits and known
mover removals are numerical losses. This is not a claim about ultimate
three-body escape. Both classifier branches are checked explicitly; numerical
loss during the ordinary fixture prevents a supporting PASS.

The `supporting/I19/` directory contains `m2_budget.csv`, `budget.png/.eps`,
input, receipt, source/build/executable/SPICE manifest, stdout/stderr and
resource history. `analyze_budget.py` preserves all rows and requires actual
Na/Ar emission, adsorption, deposit and escape plus Ar release, preventing a
vacuous zero-event PASS. The signed source/sink equations are documented in
`test/M2_IMPLEMENTATION_PLAN.md`; complete formal I19 remains unresolved until
its predecessor and complete scope pass.

Purpose: Close synchronized Na and Ar particle and mass budgets across adsorption, desorption, inventory change, escape, physical sinks, and numerical loss.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I19 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: No external observation data. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The approved M2 budget case uses neutral Na and Ar with photolytic chemistry
compiled off. Consequently `converted_in` and `converted_out` are required to
be exactly zero, not estimated from a noisy ensemble. The dedicated Na+Ar
build helper verifies the generated `_PIC_PHOTOLYTIC_REACTIONS_MODE_` macro and
records its SHA-256 before it preserves the executable.

The budget table must separately report physical atoms and kilograms for
external injection, adsorption, desorption, initial/final exosphere,
initial/final surface inventory, physical escape, and numerical domain loss.
Adsorption and desorption are internal transfers and cancel from the global
closure equation; they remain explicit diagnostic columns. Deterministic
bookkeeping must close to machine precision, and the full weighted Monte Carlo
window must have absolute imbalance below 0.5% of its dominant source.

Formal I19 is currently `SKIPPED / NOT VALIDATED` because same-campaign I13
is unavailable. The supporting eight-rank Na/Ar fixture now closes its actual
gas/surface identities and exercises weighted finite release, incoming
adsorption and outward escape. This scoped result does not qualify every
mover, source, chemistry conversion, hybrid layout or time/weight resolution.
Those required paths and the complete formal window must be established
before formal I19 PASS; existing surface flux arrays alone remain insufficient.

An executed formal comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
