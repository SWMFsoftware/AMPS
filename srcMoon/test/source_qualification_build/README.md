# Build a source-bound Na/Na+ producer in a disposable checkout

This is a build/preparation runner, not a particle run or an I21 PASS. Read
[CONTRACT.md](CONTRACT.md). It reuses the maintained source-profile compile
recipe, normal clean Config/make path and actual normal-producer proof function.
All resources, staging, compiler children and reporting share one allowance
of at most3600s and a200GB process-tree memory ceiling; make uses16jobs and
the existing -O3/-g setting. It never configures the active checkout.

First prepare the explicit source-bound profile with
`utility/chemical_database/state_profiles/prepare_atomic.py`. The old partial
profile keeps Na ion unavailable and correctly causes SKIPPED when selected.
Create the exact mapping (canonical compiled symbol is **NA_PLUS**, as in
`input/species.input`; the generated species macro is `_NA_PLUS_SPEC_`):

```sh
cat > /tmp/moon-source-mapping.json <<'JSON'
{"Na":"23Na:ground","NA_PLUS":"23Na+:ground"}
JSON
```

Create/prepare a separate disposable checkout outside the active code tree,
including the current working source bytes, workflow tools and private/rebased
site configuration. The runner checks exact source identity (including extra
or stale native files), workflow hashes and unre-based active-root paths before
any build. It does not create a clone, overwrite a stale checkout or infer its
dependency/site configuration. A source mismatch is ERROR before configuration.

Bootstrap a freshly cloned disposable checkout with `./Config.pl -install`
inside that checkout before invoking this builder. The normal installation
generates `Makefile.test`, which site `Makefile.def` includes even for an ordinary
application build. Missing generated prerequisites are refused before staging;
do not manufacture or patch a generated makefile. Shared-library/site paths must
point to private dependencies or explicitly read-only external inputs, never
back to writable files in the active checkout.

Invoke from the active checkout with absolute paths to the new profile and
private checkout; OUTPUT must be new:

```sh
python3 -B srcMoon/test/build_m3_source_qualification.py \
  --profile /absolute/path/to/atomic-profile/profile.json \
  --mapping /tmp/moon-source-mapping.json \
  --checkout /tmp/disposable-AMPS \
  --output-dir test_output/srcMoon/NEW-source-qualification-build \
  --reference-problem homogeneous-zero-acceleration \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels \
  --wall-seconds 3600
```

`normal-template` retains the supplied canonical force binding instead. The
explicit homogeneous option uses the existing zero-acceleration function and
marker through the normal AMPS compile input; it does not invent a force law.
`lunar-point-mass-gravity` builds the gravity-only reference for the I21
decision-1 force-driven path: only the production lunar point-mass term is
applied ([I21 decision-1 contract](../integration-tests/I21/DECISION1_CONTRACT.md)).
Neither option selects error bounds, a radiation spectrum, source kinetics,
event execution or observational acceptance.

The supervised private worker calls the production resolver-backed compile
recipe and copies its three normal application companions only into the private
`input/` directory. It removes that checkout's `build`, configures
`moon_source_profile` with TESTMODE off, runs make-j16, and retains the executable,
archive/generated-header identities and exact input/species/profile recipe.
The disposable checkout remains in this source-bound configuration for the
next explicitly authorized cohort runner. No active/canonical file or baseline
reference is restored/replaced because none was changed there.

`normal_build_evidence.json` has schema
`srcMoon-actual-normal-restore-producer-v1` and is produced by the actual existing
`normal_restore_evidence` function only after the two Config/make commands
actually returned0. It includes the source-before/after identity, library,
retained executable and current/preserved generated ABI hashes; additional
`source_qualification_recipe` metadata binds exact derived inputs/mapping and
source-SI charge policy. This is usable by the existing cohort wrapper **from
the SAME private checkout**, not an unrelated active generated tree.

Inspect `result.json`, `worker.log`, `resources/`, `build_logs/moon-spice.log`,
`source-recipe/run_manifest.json`, `generated/` and the producer receipt. PASS
means actual build/provenance criteria evaluated and met; absent source state
or dependency is SKIPPED77; wrong source/configuration, compiler/runtime/resource
failure is ERROR2. No physical model or formal predecessor is qualified here.

Protocol guards can be checked without a configure/build:

```sh
python3 -B srcMoon/test/source_qualification_build/test_protocol.py
```

These labelled administrative fixtures test refusal semantics only and never
manufacture successful production build receipts or physical profiles.
# Conservative product storage

The builder explicitly stages the normal recipe's
`packet_layout='conservative-products'` option. This enables individual weight
correction through the normal input parser and requires a fresh clean build;
default recipe callers still retain their template layout. The actual generated
ON macro is checked after compilation and recorded in the producer evidence.
No timestep policy, core constant, species mass, or source budget is changed by
this packet-layout choice.

The custom species include preserves the canonical complete symbol table and
changes only the selected states' masses/charges. Inactive hydrogen/ENA symbols
are needed by existing core compilation even when only Na and NA_PLUS are active.
Keeping those declarations does not add injected or transported species. The
actual normal species parser and clean compile exercise this compatibility
requirement; no core compatibility macro is invented by the test.

### Build the I24 sodium capture case

Use the ordered [Na, Na⁺, H mapping](../integration-tests/I24/reference/compiled_species.json)
and a source-authenticated profile containing all three ground states. The
projectile is a prescribed proton bath with explicit mass and momentum
exchange. The outgoing H particle is tracked on the mesh. The builder checks
the profile, source checkout and generated individual-weight ABI. Prepare
the disposable checkout as described above. From the source checkout run

```bash
/usr/bin/python3 -B srcMoon/test/build_m3_source_qualification.py \
  --checkout "$private_checkout" --profile "$profile" \
  --mapping srcMoon/test/integration-tests/I24/reference/compiled_species.json \
  --reference-problem homogeneous-zero-acceleration \
  --spice-path /home/vtenishe/SPICE/cspice \
  --spice-kernels /home/vtenishe/SPICE/Kernels \
  --wall-seconds 1200 --output-dir "$build_receipt"
```

Supply `$build_receipt/normal_build_evidence.json` to the I24 aggregate
while the private checkout and its build artifacts are available. Putting
a build directory on temporary storage means a machine reboot may require
rebuilding it. [Run the linked tests](../integration-tests/I24/README.md#running-both-the-mpi-cohort-and-the-linked-pic-movers)
