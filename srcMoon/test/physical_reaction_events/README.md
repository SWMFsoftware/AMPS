# Source-backed event-composer numerical controls

Run the maintained production composer and independent Decimal reference:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_physical_reaction_events.py \
  --output-dir test_output/srcMoon/physical-events-NEW
```

This is a PIC-independent source/kinematics/plan test; no AMPS configuration,
MPI launch or physical lunar state profile is enabled. C++17 and Python3 are
required; Matplotlib supplies optional PNG/EPS while numerical CSV is retained
without it. The complete preimplementation contract is [CONTRACT.md](CONTRACT.md).
The runner is registered as `physical-reaction-events` in the existing global
M3 supporting services and U15 owner. Full U15/I21–I24 remain unqualified when
their actual physical/source/predecessor scope is unavailable.

The input rates, thresholds and quiet/active means are the REAL frozen D04
PHIDRATES rows. The runner reads package README/provenance/QA, requires native
extraction QA PASS and checks every processed-output hash before compiling.
The production registry additionally checks its original frozen file identities
and complete Na family. Source extraction integrity is distinct from a physical
Na state/threshold crosswalk; the latter remains unavailable and is not invented.

All numerical states/masses in `ANALYTICAL_FIXTURE_NOT_PHYSICAL_PROFILE.json`
are explicitly fictional mechanics inputs: neutral reference mass4e-26kg,
electron mass1e-30kg, additive ion mass, and an ion formation gap set equal to
the source threshold ONLY for the analytical closure control. These are not
isotope/electron measurements or approved actual Na populations. The caller's
planner admission bit is labelled analytical storage/closure scope; a separate
negative control requires refusal without it. No Boolean can create a physical
profile or formal qualification from these tests.

Six predeclared cases cover quiet/active endpoints and rest, moving and exact
affine-frame parents. The independent70-digit Decimal evaluator reconstructs
the inertial parent from `v_event=R*v_inertial+offset`, computes the representative
photon Doppler transform and inverse-boost momentum, solves the two-body COM
momentum/energy equations directly and independently forms outgoing-electron
minus incoming-photon signed bath moments. It then applies the same declared
affine coordinate change, including offset momentum/work terms. The reference
does not call production kinematics or use a stored production output as truth.
All comparisons use the predeclared64-binary64-epsilon arithmetic bound and
explicit vector scales; no failed result changes that bound.

The incoming representative photon is transformed consistently, but D04's
frozen rate/mean source law is not a Doppler-shifted solar spectrum or a measured
differential product PDF. The composition is explicitly EXPERIMENTAL, classical
material kinematics with supplied error limits, a source-mean/isotropic closure
and an electron/radiation reservoir. It is not full relativistic transport,
measured daughter escape/brightness or physical validation. Complete molecular,
electron-impact and all-state CEX products are not silently substituted.

Native controls cover exact nuclei/charge/represented number, original unused
time, fresh clocks, source/state identities, caller admission, independent error
limits, improper/missing frame derivatives, physical epoch consistency,
nonopen draws, shadow/distance/domain rejection, unsupported families and
transactional outputs. The same admitted ion definition is checked at the
optional ideal collector; this remains an artificial boundary inventory transfer,
not a lunar neutralization, heating, implantation or release law.

The new output directory retains compiler/command/stdout/stderr, executable and
source hashes, source-data identities, explicit analytical fixture, complete
`native.checks.csv`, raw `native.metrics.csv`, independent `comparison.csv`,
`comparison.png`, `comparison.eps`, `run_manifest.json` and `result.json`.
PASS means the declared software/numerical control passed; FAIL means an actual
criterion missed; ERROR is build/runtime/parser/corruption; SKIPPED is a genuine
unavailable source/compiler. Exit codes are0/1/2/77 respectively. Prior outputs
are never overwritten. These files are never observations or a frozen physical
reference solution.
