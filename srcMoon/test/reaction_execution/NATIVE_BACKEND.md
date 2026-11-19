# Real PIC reaction adapter controls

Read [the pre-implementation contract](../M3_PARTICLE_BACKEND_CONTRACT.md).
This runner clean-builds an eight-species Moon configuration and links a thin
fixture to the **actual generated production library**, including
`MoonReactionAdapter.cpp` and `ReactionExecution.cpp`. It initializes real
PIC mesh topology, MPI block ownership, cell lists and a64-record particle
buffer. It does not initialize lunar ephemerides/surfaces, score observations,
enable physical particle chemistry or replace a kinetic distribution.

Run from a disposable/isolated AMPS checkout with copied machine dependencies:

```bash
/usr/bin/python3 -B srcMoon/test/run_m3_particle_backend.py \
  --isolated-checkout \
  --output-dir /absolute/path/to/new/m3-particle-backend-evidence \
  --wall-seconds 3300 --build-jobs 8 \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels
```

The script removes that checkout's `build/`, runs normal `Config.pl` and
`make -j8`, then executes the same retained probe at **1 and8 MPI ranks**,
OMP=1. `--isolated-checkout` explicitly permits changes to the selected
checkout's generated configuration; it does not create a clone itself. Do
not run this on an active build tree. The one monotonic budget includes
configuration/build/probe/MPI attempts; the caller must also budget clone
setup/cleanup so the complete development package remains below3600s.

Individual weights and species-local weight/time-step modes are enabled only
in unique derived input files. Normal input/reference bytes are retained.
Parent dt.5s and product dt.125/.0625s deliberately differ. For represented
parent number1.5, each declared product must still represent1.5: no daughter
birth-dt multiplier is permitted. Products retain the original physical
start/duration/elapsed interval and their fresh pending clocks.
The canonical velocity-distribution sampler is disabled in this fixture:
its Moon sampling coordinates are outside the analytical unit box. This
changes only that unrelated diagnostic and retains normal PIC
`Init_AfterParser`, actual mesh allocation and production particle services.

The topology is explicitly an **analytical stationary zero-Q nominal-heavy-
mass storage control**: Na + bathHe+ → Na+ + He. It uses actual compiled
mass/charge values and refuses unexpected neutral/ion nominal-mass
inequality. It supplies no cross section or lunar CEX/recoil authority.
Incoming bath nuclei/charge/mass are explicit and conservation is checked
by the production planner. This is an adapter test, not ionization or CEX
kinetic qualification.

A separate one-dimensional elastic nominal-mass control gives the Na parent
velocity+1m/s and incoming He+ bath velocity−1m/s. Analytical two-body
conservation gives distinct daughter velocities, which are read back from
actual PIC records. The bath moments come from its declared incoming state;
they are not fitted residuals. This remains a numerical storage/field test,
with no physical CEX angular distribution, rate or threshold claim.

The real allocator is filled until one free slot remains for a two-product
packet. Deferral must preserve parent bytes, public list, count and accepted
ledger; after real filler release the same immutable plan must commit.
Additional cases cover private staging/rollback, parent-byte and catalog
mutation refusal, actual compiled child indices/state/clocks, list links,
repeated commit refusal and detached products with mover-owned parent
deletion. Detached ownership does not by itself qualify MPI transport or a
trajectory-event schedule.
Two staged transactions also share a real accepted-event ledger; a one-shot
C++ allocation fault checks that renewed capacity reservation fails before
physical publication. Rollback and same-plan retry must preserve the earlier
accepted event. The PIC particle allocator itself is never mocked.

Outputs are `result.json`, `run_manifest.json`, configure/build/probe logs,
exact compile-input and generated-header snapshots and source/library/binary
identities. `np1/` and `np8/` each retain native stdout/stderr, result and every
rank's checks/products/species CSV. PASS means completed declared native
controls. Evaluated invariant misses are FAIL; build/launcher/runtime/partial
evidence failures ERROR; absent explicit isolated-build authorization is
SKIPPED. Exit codes:0 PASS,1 FAIL,2 ERROR,77 SKIPPED. The global runner must
invoke this supporting runner and retain its structured status; it does not
assign a new permanent test ID or close formal chemistry gates.

After a successful normal application build followed by a thin-probe or
launcher infrastructure failure, `--retained-build /path/to/run_manifest.json`
may relink the test and rerun **all**1/8MPI controls without recompiling
unchanged production code. It requires the actual prior make phase to have
exited0 and verifies source, library, generated ABI and compile-input snapshot
hashes. The earlier native ERROR stays retained; it is never relabelled PASS.
Thin mesh probes explicitly use the production `-DMPI_ON` preprocessing flag:
linking an MPI library alone cannot enable MPI branches in inline templates.
