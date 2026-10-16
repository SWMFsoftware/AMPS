# I19 — Surface inventory and global mass-conservation accounting

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

I19 is currently `SKIPPED / NOT VALIDATED`, because production-wide weighted
event hooks for every required boundary/source/domain-loss path are not yet
complete. Existing surface flux arrays do not distinguish physical escape
from numerical deletion and therefore cannot be relabelled as a closed global
budget. This is a missing production capability, not a reason to weaken or
approximate the acceptance criterion.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
