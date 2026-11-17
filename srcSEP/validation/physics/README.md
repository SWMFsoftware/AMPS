# Controlled wave-physics studies

These tools verify controlled wave-physics references: linear parallel
proton–electron Maxwellian dispersion and finite-time relativistic Lorentz
response in prescribed circular waves. They do not set production domain
limits, adopt an amplitude cutoff, infer missing plasma temperatures, or perform
observational validation. The nonlinear particle/wave,
background-work, reflection, cascade and moving-tube contracts remain separate.

[FINITE_WAVE_ENERGY_COLLISION_DESIGN.md](FINITE_WAVE_ENERGY_COLLISION_DESIGN.md)
records the separate CP05/Q04 finite-budget collision derivation and proposed
independent benchmarks. It is a design, not an implemented finite-reservoir
solver: fixed-W versus constrained drift, positivity, sharp spectral
interfaces, finite-weight noise, combined event timing and physical versus
canonical momentum must all be resolved before production qualification.

The source-owned finite-amplitude experiment is documented in
[FINITE_AMPLITUDE_RESPONSE_STUDY.md](FINITE_AMPLITUDE_RESPONSE_STUDY.md).
It reuses the same equilibrium/constants helper, adds explicit monoenergetic
particle/wave controls, independently compares Lorentz orbits with finite-time
unperturbed response, and separates orbit/phase/mode/representation errors.
Its default action prepares only; trajectory runs are opt-in. Neither a
single-wave variance nor a mixture of separate single-wave experiments is
reported as a Markov diffusion coefficient. Its bounded gate is
`srcSEP/test/run_finite_amplitude_study_tests.sh`.

`controlled_ep_maxwellian_baseline.json` records both temperatures and densities
explicitly. Its ambient n/B/T values are documented in the SWCME defaults; its
Maxwellian distributions are controlled-study inputs. The production
`ProtonOnly` thermodynamic closure does not establish an electron distribution
or downstream electron temperature. All numerical constants are explicit;
the file does not assert identity with every compiled AMPS constant set.

## Prepare without calculating roots

From the repository root:

```sh
python3 srcSEP/validation/physics/run_warm_parallel_study.py \
  --case srcSEP/validation/physics/controlled_ep_maxwellian_baseline.json \
  --output-directory /tmp/sep-warm-prepared --prepare-only
```

Preparation compiles `export_wave_grid.cpp` against the actual
`ResolveCenteredLogGrid` API. The native default is 128 logarithmic centers
from 1e-8 to 1e-2 m^-1, with geometric half-width outer edges. Native doubles
are retained as hexadecimal strings. Explicit fresh-study grids can use
`--grid N FIRST_CENTER LAST_CENTER`; this never resizes or remaps live arrays.

The output includes raw input bytes, requested decimal coordinates, native
grid/coordinate mapping, source archives and hashes, compiler identity, and
`protocol.json`. The default `native-grid` sampling uses actual native double
values. `--sampling nominal-grid` uses the high-precision geometric reference
grid of the recorded configured decimal bounds. Their differences remain
separate from frequency numerical errors. `--sampling explicit-kdi
--points-json FILE` accepts an explicit ordered decimal-string list without
claiming association with native bins. Grid sampling also accepts recorded
`--points-json` values, but verifies their exact decimal rounding to the
declared grid; it does not apply a physical comparison tolerance.

The capture uses the production value API under the recorded exporter compiler
and flags. It is not a live PIC-state dump. Binding a study to a particular
running application additionally requires its resolved configuration and actual
snapshot grid; this standalone tool does not claim that production binding.

Archived inputs and sources are immutable. Choose a fresh output directory for
a new invocation. Generated artifacts, temporary compiler tools
and calculation outputs belong outside the source tree.

## Explicit calculation or import

Calculation is intentionally opt-in and can be slow:

```sh
python3 srcSEP/validation/physics/run_warm_parallel_study.py \
  --case srcSEP/validation/physics/controlled_ep_maxwellian_baseline.json \
  --output-directory /tmp/sep-warm-calculated --run --plot
```

`--helicities Right` or `Left` can select one study branch. `--digits` and
`--verification-digits` control numerical convergence, with defaults 40 and 70.
The exact solver command, raw requested points, solver implementation and input
hashes are recorded. Numerical root isolation is a refined winding/integral
check, not an interval-arithmetic proof or full global mode-family certificate.
Unresolved numerical damping remains unavailable; it is not replaced with zero.

Existing studies can be verified and plotted without solving any roots:

```sh
python3 srcSEP/validation/physics/run_warm_parallel_study.py \
  --case /path/to/archived/raw-input.json \
  --output-directory /tmp/sep-warm-imported \
  --sampling nominal-grid --points-json /path/to/archived/requested-points.json \
  --import-studies /path/to/right.json /path/to/left.json --plot
```

Import requires the exact raw-input hash, current solver-source hash and actual
implementation archive, numerical settings, requested coordinate spellings,
reported-coordinate agreement within the declared numeric serialization bound,
complete requested row coverage, and nonoverlapping helicities. It preserves
raw input studies byte-for-byte. A changed solver requires reviewed regeneration
or execution of the archived protocol, rather than relabeling old evidence.

Exit 0 means prepared, or evaluated with all requested roots numerically
resolved. Exit 1 means an evaluated study retains explicit `UNRESOLVED` rows.
Exit 2 means protocol/input/provenance failure. No unresolved result selects a
physical cutoff, permits interpolation over an unknown interval, or enables
production. A numerical failure retains its reason.

## Plotting and interpretation

```sh
python3 srcSEP/validation/physics/plot_warm_parallel_study.py \
  --study /path/to/study.json --output-prefix /tmp/sep-warm-figure
```

The source-owned plotter produces PNG, PDF and a `.plot.json` provenance record.
`--plot` executes the pinned content-addressed plotter archive, then verifies
its actual source and study hashes. That plotter is standalone: copying the
bundle does not require importing mutable sibling scripts from the checkout.
`protocol.json` records both the original absolute solver command and portable
`calculation_command_relative` / `plot_command_relative` argument lists, run
with the bundle as working directory under the recorded dependency versions.
Replaying the calculation remains an explicit slow action; plotting existing
evidence calculates no new roots.
Labels use the study's actual input quantities. Explicit unresolved rows break
curves; missing rows are rejected. Damping below numerical resolution is omitted,
and exact zero error is omitted from logarithmic error panels without adding a
floor. The 1%, 5%, and 10% lines are reference percentages, not adopted budgets.
Positive values outside Matplotlib's floating-point range produce an explicit
error instead of underflowing to zero or disappearing. The reported frequency
convergence estimate is not a separately certified per-point error bar for the
group derivative or Hermitian response.
Sampling brackets never bridge an unresolved row and are not interpolated
cutoffs. For strongly damped roots, the Hermitian response is not an exact
conserved wave energy and the real group derivative need not be an energy speed.

Dependencies: Python, mpmath and a C++11 compiler for preparation/calculation;
Matplotlib only for plotting. No external dataset or network call is required.
Run bounded mathematical and protocol checks with:

```sh
srcSEP/test/run_warm_parallel_study_tests.sh
```

The susceptibility and weak-damping response use the primary
[NHDS description](https://verscharen.com/papers/NHDS.pdf),
[author susceptibility implementation](https://github.com/danielver02/NHDS/blob/master/src/suscep.f90),
and [author wave-energy implementation](https://github.com/danielver02/NHDS/blob/master/src/waveenergy.f90).
The Landau-contour and cold-limit tests are independent of the production MHD
and scattering kernels. Manufactured protocol records test software handling
only and are never scientific reference solutions.
