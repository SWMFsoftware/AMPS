# Step 12 cross-path validation and release gate

This package closes the Phase-1 release without changing a trajectory, cutoff,
unresolved-support, convergence, or observation threshold. It verifies that the
standalone and SWMF-coupled products differ only in field acquisition—not in particle
physics, units, integration controls, boundary distribution, or output semantics.

The released scope is `INSTANTANEOUS_QUASI_STATIC_MAGNETIC`. Dynamic E/B
characteristics, long-duration trapping across snapshots, local acceleration, physical
loss processes, and forward trapped-flux prediction remain unsupported and must be
listed that way in the capability table.

## What must exist before release

1. Build manifests for `STANDALONE` and `SWMF_COUPLED`. Both must contain the tag
   `sep-in-geospace-phase1-static-characteristics-v1`, the same SHA-256 of the common
   physics source set, the same clean source revision, compiler/dependency records, and
   their own executable SHA-256. The executable hashes need not be equal.
2. PASS evidence for the fixed Phase-1 matrix: U-F01–U-F10, U-F12–U-F13;
   I-F01–I-F07 and I-F10; C1–C19; F1–F17; and O1–O4. Dynamic and forward-only U/I
   gates are not silently claimed by the static release.
3. A sampled-field campaign comparing the direct analytic/phenomenological field with
   the same field sampled on the AMR mesh.
4. An SWMF replay campaign comparing a live frozen SWMF snapshot with its exported
   standalone offline replay.
5. Both campaigns must compare all seven roles: field samples, individual trajectory
   results, directional access `A(E,Omega)`, cutoff, spectrum, density, and detector
   products. CSV key sets must match; terminal/identity columns match exactly; every
   remaining column must be explicitly numeric and finite.
6. A preregistered, frozen O4 result with `retuned=false`. A failed O4 is reported as a
   failed gate; it is never absorbed by a new scale, exclusion, or threshold.
7. Hash-pinned capability and measured resource-estimate tables, plus concrete CCMC,
   snapshot-export, replay, standalone, and coupled commands.

Every evidence record also carries revision/dirty state, compiler/dependencies,
executable and input hashes, field and snapshot identity, boundary and response hashes,
species, epoch, frame, mover, integrator tolerances, complete grid history, MPI/thread/
scheduler settings, random seeds, wall time, termination counts, and closed artifact
hashes. Observation-facing F8/F9/F10/F17 and O1–O4 require a nonzero comparison count
and the shared-event/no-platform-scale normalization policy. A dry run cannot be used
as PASS evidence.

## Thresholds are ceilings

The manifest may make a comparison stricter, never looser. `EXACT` is byte/numeric
exact; `KERNEL` is at most `rtol=1e-10`, `atol=1e-12`; `INTEGRATED` is at most 2%;
and `MESH`/`DETECTOR` are at most 5%. Frozen live-versus-replay values are additionally
limited to kernel precision. Sampled analytic-versus-mesh values may use the 5% mesh
representation ceiling, but discrete termination/access outcomes remain exact.

The comparator uses

```text
abs(candidate-reference) <= atol + rtol*max(abs(reference),abs(candidate))
```

and rejects NaN/Inf. It also rejects duplicate/missing keys, ignored columns, duplicate
roles, missing roles, and a reference file reused as its own candidate.

## Prepare a campaign

Create the full registry skeleton:

```bash
python3 srcEarth/release_validation/make_manifest_skeleton.py \
  --output release/step12_release.json
```

The skeleton is intentionally not runnable: every `REPLACE_...` value and all-zero
digest must be replaced with real evidence. Compute a resource digest with
`sha256sum FILE` and keep paths relative to the release manifest when possible. Use the
examples in `srcEarth/examples/step12_release/` for build, gate, capability, and
resource records.

The archived command map must contain both supported standalone workflows. A cutoff only
example is `mpirun -np 8 ./amps -mode gridless -i standalone_cutoff.in -nt 16`.
A cutoff + flux/spectrum example is
`mpirun -np 8 ./amps -mode 3d -i standalone_products.in -nt 16`. The analogous SWMF
cutoff-only and combined commands use the site's coupled launcher, and offline replay
uses `-mode 3d` with the exported `SWMF_SNAPSHOT` input. Record resolved inputs rather
than leaving template tokens in an actual release manifest.

Each parity CSV must have a header and at least one row. Declare every header exactly
once across `key_columns`, `exact_columns`, and `numeric_columns`. Recommended exact
columns include snapshot/control fingerprints, access state, termination reason,
retry state, channel name, and units; numerical columns hold coordinates, fields,
phase-space states, access, cutoffs, intensities, density, flux, and rates.

## Run and review

Run only the cross-path comparison:

```bash
python3 srcEarth/release_validation/compare_cross_path.py \
  --manifest release/step12_release.json \
  --output release/cross_path_comparison.json
```

Run the complete release decision:

```bash
python3 srcEarth/release_validation/run_release.py \
  --manifest release/step12_release.json \
  --output release/step12_release_summary.json
```

Both commands print `RESULT: PASS` only for a complete successful evaluation and
return zero. A malformed, missing, mismatched, relaxed, retuned, or failed input prints
`RESULT: FAIL`, writes a structured error, and returns 2. Configuration preflight is
available as `--validate-only`; it prints `RESULT: VALID_NOT_RELEASE` and deliberately
does not claim physics PASS.

## Restart behavior

Add `--restart` to the complete command only after an earlier PASS. The evaluator
rehashes the manifest and every referenced build, evidence, artifact, capability,
resource, and parity file. It reuses the prior summary only if the resulting input
fingerprint is identical, printing `RESTART: SKIPPED_VERIFIED_PASS`. Any modification,
missing file, failed hash, previous non-PASS, or changed release ID forces validation
or failure. The verified prior report is not rewritten during a skip.

## Portable tests versus release evidence

Run the evaluator tests with:

```bash
./srcEarth/test/UStep12Release/run_test.sh
```

They use manufactured bundles to prove that the contract accepts a complete exact
case and rejects all listed negative cases. They do not claim that a linked AMPS/SWMF
campaign or any public observation comparison has passed. Only archived production
U/I/C/F/O evidence may populate a real release manifest.
