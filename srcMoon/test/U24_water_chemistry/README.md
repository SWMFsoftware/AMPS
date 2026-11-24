# U24 — H2O/OH migration and chemistry

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U24_water_chemistry/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

### Conditional residence on a smooth surface

The U24 runner selects the Laferriere et al. Figshare v1 model from the
locked r14 catalogue. For a molecule **already on the surface**, the
model defines `tau = exp(E/(k_B T))/(f A)` with `f = 0.01`,
`A = 10^13 s^-1` and `E = 0.6 eV` for OH and H₂O. The probe passes
these parameters to the lunar `ResidenceTimeS` and `DesorptionProbability`
kernels. At fixture temperatures of 150, 250 and 400 K it compares both
outputs with a 60-digit decimal solution over 0.01, 1 and 5 lifetimes.
The predeclared relative tolerance is 1e-12 for these finite arithmetic
checks. There are no random draws, spatial coordinates or observation
epochs. AMPS uses the exact SI joule/eV conversion and Boltzmann constant.
The archived hopping code uses rounded values of 1.60218e-19 J/eV and
1.38064e-23 J/K. The comparison CSV reports the author's resulting lifetime
separately. That difference does not enter the arithmetic tolerance.
[Author code and model](https://doi.org/10.6084/m9.figshare.28548857.v1)

Set `AMPS_CHEMISTRY_LOCK` to the r14 lock on another machine. By default
the test selects the r14 snapshot in `/data/vtenishe/CHEMICAL_DATABASE`.
The `residence-*/` child under the U24 output contains the actual command,
compiler log, numerical comparison CSV, PNG and EPS comparison figures,
source hashes and child result.
A missing catalogue or model skips the child. Changed evidence produces
ERROR and evaluated arithmetic disagreement produces FAIL. No incident
sticking, launch-speed distribution, migration trajectory or D04 event
is evaluated by this child. U24 retains its U11/U12/U13/U15 predecessors.
