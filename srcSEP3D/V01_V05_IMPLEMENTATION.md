# Improvements V01–V05

## V01 — controlled perpendicular diffusion and drift

The spatial tensor is `K = kappa_perp I + (kappa_parallel-kappa_perp) bb`.
`transport/perpendicular_transport.*` is the single owner of tensor assembly,
its Itô drift, a deterministic transverse basis, and relativistic first-order
gradient-B and curvature drifts. Supported closures are a constant SI
`kappa_perpendicular` and a constant ratio to the local parallel coefficient.
For the ratio closure, `d(kappa_perpendicular)/ds` is the ratio times the
available field-aligned parallel derivative; no transverse gradient is
invented.

The drift uses signed charge and `p*v`. Focused transport uses the particle's
instantaneous pitch; Parker transport uses the isotropic averages
`<mu^2>=1/3` and `<(1-mu^2)/2>=1/3`. Two new keyed purposes own the transverse
Wiener processes, so enabling V01 cannot shift released parallel/pitch random
histories. The timestep limiter uses the largest diffusion-tensor eigenvalue.
Current-sheet drift remains unsupported because no background contract defines
a sheet surface, thickness, or regularization.

Input fields under `[transport]` are `perpendicular_diffusion =
none|constant|constant-ratio`, `constant_kappa_perpendicular_m2_per_s`,
`kappa_perpendicular_to_parallel_ratio`, and `drifts =
none|gradient-b|curvature|gradient-curvature`. Inert legacy `false` remains
readable; ambiguous `true` is rejected. Drift freezes magnetic-gradient
storage before mesh allocation. `V1D01`–`V1D05` test the algebra, moments,
charge/polarity reversal, focused invariants, and timestep.

## V02 — genuine 1-D/3-D parity

`test/run_cross_model_parity.py` compiles different executables. The 1-D
producer links the real `srcSEP` Parker/focused cores; the 3-D producer links
the real `srcSEP3D` cores. Both write canonical JSON under matched SI inputs.
V2D01 requires distinct executable hashes, hashes both records, and compares
Parker and focused moments against predeclared tolerances. A shared helper
compared with itself cannot satisfy the gate.
The 1-D root is supplied explicitly with `--sep1d-source`; srcSEP3D never
searches for or assumes a sibling application tree.

## V03 — native MPI qualification

`validation/native_profiles.json` defines small, medium, and production
rank/thread matrices. `run_native_matrix.py` launches the real linked AMPS
binary, recording argv, thread count, executable SHA-256, wall time, native
reports, and output tail. It never emulates MPI and errors when a prerequisite
is absent.

## V04 — scientific validation ladder

`validation/v04_campaign.json` orders analytic, native, cross-model,
observational, and live-SWMF evidence. Each rung names required IDs and uses
immutable SHA-256/SI provenance. Later evidence cannot replace an earlier
rung. `XM3D03` is now an evidence-backed perpendicular-transport diagnostic.
The live-SWMF rung is explicitly blocked by deferred R8; analytic or replayed
snapshots cannot be reported as live-coupling success.

## V05 — release governance

`release/capabilities.json` separates coded, component-tested, integrated,
validated, and blocked states. `release/profiles.json` defines development,
native-integration, MPI-qualification, and scientific-release gates. The
generator reads normal test/validation reports and emits JSON plus Markdown.
Any required MISSING, SKIP, FAIL, ERROR, or structural block gives INCOMPLETE
and nonzero exit. The scientific-release profile remains blocked by R8.
