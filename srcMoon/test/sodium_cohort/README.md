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
