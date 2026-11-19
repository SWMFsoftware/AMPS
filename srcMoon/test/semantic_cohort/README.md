# Normal startup and semantic PIC composition

This is supporting U15/U16 software verification. It does not qualify a real
sodium channel, manufacture an ion state, or close formal I21–I24. See the
[contract](CONTRACT.md) for the exact inventory, independent references,
configuration, units, artifacts and status rules.

The maintained `run_m3_development.py --isolated-checkout` worker invokes the
SI mechanics owner after all consumers of its legacy producer finish, then
restores ordinary Moon once. Its normal registration owner reuses the same
authenticated restored archive after the common-clock startup checks. Every
child receives the original package's remaining budget; the cleanup/report
reserve is retained. Without explicit isolation, the worker never selects the
derived SI build. The seven administrative receipt checks are loaded by the
registered U16 orchestration owner, so they remain reachable from the global
test campaign.

`m3_normal_semantic_registration_probe.cpp` calls the normal `amps_init()`.
Its optional callback registers real opaque semantic storage at the actual
pre-buffer boundary. It checks normal initialization invokes it once before
allocation, refuses replacement/reentry, and leaves a usable real PIC buffer.
The callback constructs the production photon field at the normal SPICE ET;
independent direct CSPICE **state** queries verify the actual Sun/Earth
positions and epoch. The production contract uses `spkezr`'s six-component
state, not bitwise equivalence with the separate `spkpos` API: dynamic-frame
evaluation can differ by a few binary64 roundings between those algorithms.
The oracle checks state mapping/units exactly and introduces no geometry
tolerance. A diagnostic comparing both APIs is retained for the original
overstrict-oracle error; it is not a physical validation residual.
Five named checks run at one/eight MPI ranks and both tracing orders.

`m3_native_semantic_cohort_probe.cpp` reuses the distributed controller fixture
and production Moon movers. Its 21 named checks per layout exercise complete
staged rollback, immutable-plan retry, source/key reuse, full daughter genealogy,
fresh child-clock counter one, actual ownership clones, MPI continuation and
capacity refusal. It also checks actual compiled SI charges and refusal when a
transaction or collector supplies the other charge unit. An actual sphere
impact records the same source-SI charge and preserves mover deletion ownership.
Represented-number reconstruction
has a documented two-operation binary64 rounding bound; it is not an adjustable
science acceptance tolerance. A four-generation identity storage chain reaches
the declared capacity without asserting any repeated physical sodium reaction.

The SI mechanics build is explicitly analytical. It keeps repository nominal
masses, uses the exact SI elementary charge in the ordinary species parser's
numeric token convention, and independently compiles the emitted expression.
It changes no shared constant or canonical species file. The normal application
build path, archive and generated headers are preserved in the run manifest.
The existing legacy backend belongs to a separate build and must not be run
against SI charges and then relabelled.

Run the inexpensive administrative receipt controls from the repository root:

```sh
python3 -B -m unittest discover -s srcMoon/test -p test_m3_normal_cohort.py -v
```

Run startup only after an authentic successful ordinary-SPICE build has produced
its `normal_build_evidence.json`. This command does not configure the application:

```sh
python3 -B srcMoon/test/run_m3_normal_cohort.py \
  --normal-startup-only \
  --build-evidence /absolute/path/to/normal_build_evidence.json \
  --output-dir test_output/srcMoon/semantic-startup-NEW \
  --wall-seconds 600
```

The analytical SI build requires a disposable checkout. From that checkout:

```sh
python3 -B srcMoon/test/run_m3_normal_cohort.py \
  --isolated-checkout \
  --output-dir /absolute/path/to/semantic-mechanics-NEW \
  --wall-seconds 1800 --build-jobs 8
```

The runner removes that checkout's `build` before normal `Config.pl` and
`make -j8`. Its enclosing campaign must subsequently restore ordinary Moon with
another clean normal configuration/build, count setup/restoration in the same
original wall budget, and preserve the active user's checkout. Do not invoke
this build mode in an active working checkout.

The Python API is `run(output, isolated_checkout=True, wall_seconds=...,
build_jobs=...)` and `run_startup(output, wall_seconds=...,
build_evidence=...)`. Result/manifest and native logs are retained in the chosen
new output directory. Mechanics also produces `key_comparison.csv/.png/.eps`
against an independent Python HMAC oracle. These figures are analytical storage
diagnostics, not observations. Missing normal build evidence is SKIPPED;
malformed proofs, incomplete check sets, signals or arbitrary exits are ERROR;
an evaluated invariant miss is FAIL. Native exit zero alone never means PASS.
