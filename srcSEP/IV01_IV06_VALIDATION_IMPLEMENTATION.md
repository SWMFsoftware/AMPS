# IV01–IV06 integrated manufactured validation

## Common execution contract

All six cases are compiled into the selected srcSEP/AMPS application and
advertised by its native test registry. `test/run_tests.py` first checks
`--list-tests`, then supplies the reviewed SI configuration through the strict
native argument manifest. The C++ callback publishes raw numerical evidence;
case-local Python programs independently generate references; scoring reloads
both files before creating authoritative JSON/JUnit metrics and PNG/EPS review
figures. Missing or stale executables are errors, never local fallbacks.

## Implemented cases

- **IV01:** analytic Parker geometry, arclength inversion, focused streaming,
  weak scattering, two wind speeds, three timesteps, and reversed vertex order.
- **IV02:** `D_mumu=D0(1-mu^2)` plus constant focusing for three ratios and two
  initial distributions versus the exact zero-flux exponential equilibrium.
- **IV03:** a positive manufactured distribution in `s`, `mu`, `ln p`, and
  time, with exposed exact/numerical residuals and coordinate refinement.
- **IV04:** production conservative overlap remapping under rigid translation,
  stretching, compression, and sinusoidal grid motion, including cell-integrated
  free-stream and ledger diagnostics.
- **IV05:** production shock-trajectory endpoints plus continuously interpolated
  crossings for slow/fast moving shocks, stationary/moving frame controls,
  exact-node cases, and the CV09 DSA slope.
- **IV06:** frozen, one-way, and two-way resonant turbulence controls linking
  streaming, production-ledger wave growth, wave-dependent scattering,
  saturation, and a closed particle/wave energy reservoir.

## Scope and interpretation

IV03 verifies the integrated manufactured residual assembly used by the
field-line operators; it is not a grid-based replacement for production Monte
Carlo transport. IV05 deliberately isolates shock crossing and frame logic;
validation against a resolved time-dependent heliospheric shock remains a
separate native/event campaign. IV06 uses a controlled prescribed early growth
law so causal feedback and synchronization are testable independently of a
particular physical instability closure.

## Commands

```sh
make test-iv01-iv06-unit
make test-iv01-iv06-unit SEP_EXECUTABLE=/absolute/path/to/amps

python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case IV01 --validation-case IV02 \
  --validation-case IV03 --validation-case IV04 \
  --validation-case IV05 --validation-case IV06 \
  --output-dir /absolute/path/to/evidence/IV01-IV06
```

The first command is only a warning-clean source/registry gate and reports the
linked part as `SKIP` without `SEP_EXECUTABLE`. The other commands validate the
actual linked application and require every case to produce its complete
model/reference/report/figure evidence set.
