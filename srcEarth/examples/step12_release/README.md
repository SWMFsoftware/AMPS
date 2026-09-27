# Step 12 release-record examples

These files document the records consumed by `release_validation/run_release.py`.
They are examples, not scientific evidence and not a pre-approved resource estimate.

Start by generating the complete manifest (all gates and all fourteen parity pairs):

```bash
python3 srcEarth/release_validation/make_manifest_skeleton.py \
  --output release/step12_release.json
```

Then copy `build_provenance.template.json` once for the standalone build and once for
the SWMF-coupled build. Replace every `REPLACE_...` field and hash the exact executable
and common physics source set. Copy `gate_evidence.template.json` for every fixed gate;
the observation fields are mandatory for F8, F9, F10, F17, and O1–O4, and the O4
freeze fields are mandatory for O4.

`capabilities.json` is the released Phase-1 capability boundary. Add site-specific
capabilities if useful, but do not relabel the listed unsupported physics as supported
without implementing and passing the later dynamic/trapping/loss validation program.
Copy `resource_estimates.template.json`, replace illustrative numbers with measured
values for the intended production grid/event, and describe the host and measurement
basis.

The manifest's command strings must be concrete, reproducible commands. Typical
standalone forms are:

```bash
# Cutoff only; CALC_TARGET=CUTOFF_RIGIDITY in standalone_cutoff.in
mpirun -np 8 ./amps -mode gridless -i standalone_cutoff.in -nt 16

# Cutoff plus flux/spectrum; CALC_TARGET=CUTOFF_RIGIDITY+DENSITY_SPECTRUM
mpirun -np 8 ./amps -mode 3d -i standalone_products.in -nt 16
```

For SWMF, use the site's normal SWMF launch command with the Step-10 cutoff-only or
Step-11 combined input. Export the frozen snapshot, then replay the exact CSV offline:

```bash
mpirun -np 8 ./amps -mode 3d \
  -i srcEarth/examples/standalone_step9_swmf_replay.in.template -nt 16
```

The input templates must be copied and filled with real paths, epochs, snapshot IDs,
boundary spectra, and responses. Commands containing placeholders are rejected by the
release evaluator.
