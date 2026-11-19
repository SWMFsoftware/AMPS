# Conditional normal Na cohort

This owner uses normal `amps_init`, actual furnished SPICE, the configured sphere
or LOLA mesh, the maintained physical-state registry, accepted-trajectory exposure,
conservative controller/transaction, source-mean photon composer and semantic RNG.
It is **EXPERIMENTAL**. A data/preflight PASS is not a physical run or formal I21
qualification. Current authoritative profile lacks `23Na+:unresolved` and the
source/threshold crosswalk, so physical operation remains `SKIPPED / NOT QUALIFIED`.

## Source availability without invented physics values

Generate the maintained source profile, then project it into the common AMPS
input grammar. Use a NEW output directory for every command:

```sh
/usr/bin/python3 -B utility/chemical_database/state_profiles/run.py \
  --output-dir /tmp/moon-state-profile-review

/usr/bin/python3 -B srcMoon/test/sodium_cohort/prepare_source_preflight.py \
  --profile-json /tmp/moon-state-profile-review/profile/profile.json \
  --application-input srcMoon/examples/m3_common_time_second_order.in \
  --parent-state-id '23Na:ground' \
  --ion-state-id '23Na+:unresolved' \
  --electron-state-id 'electron:free' \
  --output-dir /tmp/moon-Na-source-admission

# Select the SAME disposable source checkout recorded by this producer.
MOON_NATIVE_CHECKOUT=/ABSOLUTE/PATH/TO/DISPOSABLE-AMPS
/usr/bin/python3 -B "$MOON_NATIVE_CHECKOUT/srcMoon/test/run_m3_sodium_cohort.py" \
  --producer-manifest /ABSOLUTE/PATH/TO/FRESH/PRODUCER/run_manifest.json \
  --input /tmp/moon-Na-source-admission/source-preflight.in \
  --mode source-preflight --mpi-ranks 8 --wall-seconds 3600 \
  --output-dir /tmp/moon-Na-native-source-check
```

The producer is the retained ordinary native producer from the explicit M3
campaign, with the CURRENT production source and generated ABI. It must be built
in an isolated checkout through normal Config.pl and make. This wrapper verifies
its library, generated snapshots and production identity before linking the thin
adapter. It does not reconfigure the active checkout or infer a different binary.
Invoke the native wrapper and qualification executor FROM THAT SOURCE CHECKOUT
(or by their absolute script paths there). A producer compiled in a private
checkout is deliberately rejected by the active checkout's wrapper even when
the two trees initially have similar sources. Changing only the shell working
directory is insufficient: the script's own location establishes its source
root. Keep absolute source-preflight/profile/input paths when crossing checkouts.
On this machine `/usr/bin/python3` supplies the maintained scientific packages.

Expected current result: native exit77, `SKIPPED`, with the missing ion/state/rest
mass reason, before SPICE/mesh/particle initialization. Lightweight MPI setup is
used for one output writer and coherent status; it is not a physical startup.
`source-preflight` checks the real registry and all three state/rest-mass entries;
it deliberately selects no model-error bounds just to report missing data.

## Profile and input binding

`profile_projection.py PROFILE.json NEW.in` produces a deterministic
`#begin physical_state_profile: native` block. It pins original JSON bytes,
isotope/constants/states CSV paths+hashes and every ordered source ID/path/hash/
role. Production `Prepared::InspectProfile` compares the complete projection to
the supplied request, then calls the SAME source-bound registry. JSON metadata
never grants activation. Files must remain unmodified; corruption is ERROR,
missing authoritative data is SKIPPED, and true unavailable state rows stay gaps.
Paths in each native .in assignment resolve against their owning included file.
Unknown/duplicate keys and blocks are errors.

`native.cpp` owns one named `#begin sodium_qualification: COHORT` block, schema
`srcMoon-native-sodium-cohort-v1`. `--mode preflight` additionally prepares the
actual source/state/composer; `--mode run` performs the explicit finite physical
experiment. Neither mode assigns absent values. The full phase requires:

- application input; original profile and native projection paths+SHA; all pinned
  CSV/source files; explicit parent/ion/electron state IDs and compiled mapping;
- a reviewed complete channel/state/common-energy-reference crosswalk CSV and
  pinned reduced transport policy evidence;
- separate quiet or active source endpoint, model/channel/common-energy IDs;
- each state's inertial-mass, acceleration and relativistic model-error limits;
- numerical mass/momentum/energy/threshold/vector/rotation and planner bounds;
  separate represented-number and source-rate bounds for independent checking;
- positive explicit physical horizon, equal representable step and step count,
  original parent count and represented number, source position/velocity in LSO;
- seed/cohort/full genealogy capacity, exposure/root controls, SUN-origin inertial
  frame, Earth-occultation choice, source reference distance and source identity.

No positive physical full-input example is supplied by substituting toy masses,
ion formation energies or guessed tolerances into current missing source rows.
The canonical source IDs above document the actual gap, rather than manufacturing
a runnable channel.

## Reference problem versus transport

`homogeneous-fixed-position` declares zero parent velocity and requires the actual
normal generated acceleration macro to select existing core
`PIC::Mover::TotalParticleAcceleration_default`, plus explicit
`_MOON_COHORT_ZERO_FORCE_REFERENCE_=1`. A marker alone is refused; a normal Moon
force binary cannot serve that problem. The ordinary compile recipe must retain
effective macro proof and source identities. No duplicate test-only force law is
introduced. This reference problem is appropriate for independent exponential
survival and parent-keyed first-removal times; it is not full lunar dynamics.

`trajectory-integrated-hazard` calls the SAME extracted normal coefficient stage
before each particle interval, without unrelated source injection/sampling.
Its physical source/state/model-error and independent per-parent integrated-H
requirements remain mandatory. Existing actual force and variable-exposure
supporting tests retain their separate scopes. Never fit a constant lifetime to
a spatially changing case.
Unexpected surface/escape loss of an initial neutral is censored data: this
owner refuses to manufacture a complete homogeneous binomial cohort.

The native explicit-ET field factory queries Sun/Earth/body transformations at
one actual furnished SPICE epoch. It does not borrow uninitialized first-step
caches or advance global ET. Frozen environmental basis and advancing common
transport time are recorded separately. Source reference distance is explicit;
legacy core `_AU_` is not silently relabelled exact SI.

At an accepted event the controller supplies the preserved positive LEFT-limit
rate and exposure identity. A second static shadow query must not reject that
same event on a discontinuity. The composer uses actual position/ray/range and
full affine inertial-to-LSO velocity transform, including dR/dt and Moon velocity;
it records accepted affine rate, exact source law and numerical depth bounds
separately. This numerical endpoint convention is not a universal physical
visibility-discontinuity law.

## Output and independent checks

The wrapper retains compile/native/resource logs, executable/source/generated
identities, input snapshot, `run_manifest.json` and `result.json`. Native output:

- `input_receipt.json`: exact parsed source/profile/controls plus actual geometry
  and actual acceleration macro after normal startup;
- `parent_records.csv`: every original parent, survivor or first reaction;
- `events.rankR.csv`: event identity, time and source;
- `event_metrics.rankR.csv` and `.hex.csv`: actual native initial PIC product
  snapshots, parent, signed bath, source and weighted ledger moments. Hex
  significands are exact and never narrowed through Python binary64;
- result77 for unavailable source/capability, result2 for corruption/runtime/
  infrastructure failure. Formal status remains separate.

The optional opaque initializer is checked against live PIC species, position,
velocity, weight, state/clock/interval and epoch BOTH after staging and immediately
before atomic publication. A bad offset cannot corrupt a product while leaving
an apparently conserved precomputed ledger. Actual read-back snapshots accompany
accepted ledgers; independent verification recomputes source state/nuclei/charge,
weighted moments and signed photon/electron changes rather than trusting flags.

`verify_ledger.py` and its maintained controls describe the independent exact
Fraction/Decimal checks. Qualification additionally requires the frozen same-
campaign predecessors, numerical-bias evidence, original-parent statistics and
approved simultaneous branch criteria. Missing Na physics or an incomplete gate
never becomes PASS because the parser, source audit or analytical fixture passes.

The normal force/environment preparation seam is now factored into
`Moon::PrepareProductionParticleEnvironment`, called by BOTH the ordinary step
and finite native owner after common Begin and before particle execution.
See [../coefficient_stage/README.md](../coefficient_stage/README.md) for the new
fresh-build proof. A physical trajectory remains experimental and requires its
own independent integrated-hazard reference; the constant-lifetime scorer never
admits a spatially varying case.

Before physical startup every rank participates in the same admission reduction:
ERROR outranks missing prerequisites, which outrank availability. Real per-rank
reasons are retained and rank0 reports a collective refusal if its own local
copy was available. No rank may finalize while another enters `amps_init`.
The actual common transport origin is recorded separately from the frozen
geometry epoch; source-only preflight invents neither value.

## Opt-in physical Na profile and explicit case producer — November 19

The default partial profile retains its unavailable Na ion record. The additive
`utility/chemical_database/state_profiles/prepare_atomic.py` producer now derives
a separate sodium-23 ground-ion state from the retained isotope, CODATA, NIST
ground-state and ionization evidence. It preserves the original unavailable
entries and the source-reported uncertainty semantics. This is an authenticated
research profile, not an automatic release or production activation.

For the approved experimental source-mean/isotropic model, the opt-in
`photon_energy_policy = source-effective-photon` binds the eight-column CSV in
[the energy contract](../M3_SOURCE_EFFECTIVE_PHOTON_CONTRACT.md). The frozen D04
threshold and mean still define the representative incoming photon; the resolved
physical formation gap enters the complete electron/photon energy ledger.
The physical-minus-source threshold is derived from pinned sources, never a
user-tuned correction. Default `source-equals-state` behavior and its independent
strict representation guard remain unchanged. The independent ledger verifier
authenticates the selected crosswalk before reconstructing either convention.

[PREPARATION.md](PREPARATION.md) describes the new input producer and exact
commands. It requires all 64 scientific/numerical controls explicitly, preserves
the actual native parser grammar and refuses missing controls without publishing
a runnable case. It chooses no P02 error budget. `reference_weight_mode` can be
`production` or explicitly `explicit`; the latter requires both positive species
reference weights. Native startup assigns those real block weights before
creating parents, and the existing represented-number and committed-particle
readback checks still apply. Weight correction must be enabled through a fresh
normal input/build; the isolated source-build helper proves the generated ABI.

The new preparation and source-build protocol controls are registered with U15
and preserve result.json, unittest logs and source hashes. A real eight-rank
source-preflight PASS checks available states/masses and compiled bindings only.
It is not a physical cohort, numerical-convergence gate or I21/M3 qualification.

## Approved P02 commissioning workflow

The user approved the four mechanical limits in
[M3_P02_NUMERICAL_BUDGET.json](../M3_P02_NUMERICAL_BUDGET.json) on November 19.
They are explicit per-state inputs, not changes to conservation, numerical-bias
or observational thresholds. [The controlled protocol](../M3_P02_CONTROLLED_NA_PROTOCOL.json)
and `srcMoon/examples/m3_sodium_p02_commissioning.{in,json}` freeze the first
128-parent quiet-source case. It is a zero-force analytical chemistry experiment
with actual source states, normal startup/dispatch/mover/transaction and frozen
SPICE geometry. It does not model a complete force-driven lunar exosphere.

The case uses an exterior sunward LSO point at `(10000000,0,0)` m, zero initial
LSO velocity, eight equal 8192 s steps and 2.625 represented atoms per independent
parent. Both reference weights are 7, giving exactly representable correction
0.375. Its 65536 s horizon is an explicit commissioning choice, not an assertion
that lambda*T equals log(2). Actual source range/visibility and endpoint rate
come from the retained startup receipt and authenticated D04 data. The explicit
149597870700 m normalization is IAU 2012 Resolution B2, not a relabelled legacy
AMPS `_AU_`. Arithmetic tolerances remain separate from approved P02 model limits.

In global timestep modes the public core setter keeps a minimum. The isolated
native owner now resets that reduction accumulator before applying this declared
reference interval and checks actual getter values before creating its source.
It does not change the core setter, normal application tables or common-clock
admission. A fresh source-matched producer is mandatory after native code edits.

[COMMISSIONING_AUDIT_README.md](COMMISSIONING_AUDIT_README.md) describes the
independent original-owner, complete parent/event, ledger, source-rate and
single-level first-hit checks, with CSV/PNG/EPS artifacts.
[P02_README.md](P02_README.md) describes the separate physical/reference mass,
acceleration, actual inertial-speed and relativistic-kinetic audit. Neither
claims statistical qualification or three-level convergence. Run commands are
also in [srcMoon/RUN.md](../../RUN.md). Both verifiers consume pinned original
owner manifests/results rather than accepting an arbitrary plausible event CSV.
