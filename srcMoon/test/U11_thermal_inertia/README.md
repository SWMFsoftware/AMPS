# U11 — Thermal inertia and eclipse response

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
/usr/bin/python3 -B srcMoon/test/U11_thermal_inertia/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The November 12 radiative-equilibrium case initializes every depth at 250 K,
applies exactly balancing absorbed and emitted surface flux, and advances the
same production column solver at 1, 60 and 3600 s. Uniform temperature and
zero conductive gradient are an analytical steady solution; each result must
remain within 1e-8 K. This tests a thermal invariant, without observational
inputs, spin-up assumptions or a claim of transient time-step convergence.

## Integrated direct-flux entry point — 2026-11-16

`Moon::Thermal::Column::AdvanceWithAbsorbedFlux(dt_s, absorbed_W_m2, error)`
now accepts direct solar absorption already integrated over the apparent disc.
It enters the same conservative nonlinear conduction/enthalpy solve as
`Advance`; the legacy scalar-cosine/fraction API computes its original flux
and delegates to that solver. No second thermal equation or test-only boundary
law was introduced.

The finite-disc production path evaluates the existing directional Bond albedo
at each visible ray's local incidence, integrates the joint terrain/Earth mask
and radiometric Jacobian, then supplies nonnegative W m^-2. Replacing that
integral with a central cosine or multiplying separate mask fractions would
change the physical absorbed energy. PSD uses a separate unabsorbed photon
factor and does not reuse this thermal-albedo-weighted flux.

The U11 production probe registers an additional analytical 250 K radiative
equilibrium case through `AdvanceWithAbsorbedFlux`, at a declared 60 s step and
exact flux `epsilon*sigma*250^4`. Its surface-temperature limit remains
`1e-8 K`. A negative flux must be rejected and leave the temperature unchanged.
These checks supplement the existing whole-column energy identity and legacy
1/60/3600 s equilibrium checks; they do not relax any earlier tolerance.
Their current outcome must come from the retained native campaign, not from
the presence of this entry point or the probe source.

No external observation is used in these fixtures. Passing local analytical
checks does not establish an actual-epoch lunar initial state, synodic spin-up,
SPICE wiring, spatial/temporal convergence or I17 Diviner agreement. Formal
M1/M2 and every data/operator/hold-out predecessor remain independently gated.

## Applied initialization and export controls

The `thermal-initialization` production-probe mode exercises a 2x4 logical
surface and closed octahedral terrain with the actual compiled SPICE kernels.
It checks exact uniform 250 K export of all boundaries/layers, ordered indices
and centre depths, and refusal to overwrite existing evidence. It then selects
automatic `periodic-spinup` with an explicitly declared 1 ms UTC replay ending
at simulation start. It preserves the approved 0.01 K and 10000 J m^-2
successive-cycle criteria and checks both actual completed metrics against
those limits. The 1 ms interval remains a dispatch fixture; convergence on
that short replay cannot qualify a physical lunar periodic initial state.

Actual method, TDB epoch, completed cycles and source-BEGIN epoch must agree.
An immediate explicit spin-up call must reuse exact unchanged bytes/cycles.
Positive evolution ending at the same numerical epoch must invalidate reuse
and initial-state export; a new spin-up must recompute changed state. A
mismatched cycle-end/start initialization must fail with no available proof or
export and must preserve caller outputs. These checks guard against treating
a uniform fallback or evolved columns as a separately qualified initial state.

Fixture CSVs stay below `U11/initialization_fixture/`; actual furnished-kernel
hashes are retained in `U11/kernel_hashes.json`, including failed metric runs
when the inventory is available. See the
[native initialization contract](../NATIVE_THERMAL_INITIALIZATION_CONTRACT.md).
Passing these local controls establishes neither a physical periodic state nor
I19 numerical/statistical qualification. Newly added native checks require a
freshly rebuilt production library and retained campaign evidence.
