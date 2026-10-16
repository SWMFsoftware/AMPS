# AMPS pitch-angle implementation supplement — specification v1.2

This package accompanies the revised `PITCH_ANGLE_SCATTERING_MODEL_AMPS_v1_2.md`.
It supplies an implementation review, a detailed roadmap, a compile-checked API
contract and independent mathematical checks. It is not an implemented AMPS
library and does not replace the missing canonical v1.1 data archive.

| File | Purpose |
|---|---|
| `CHANGE_REVIEW.md` | Findings, improvements and evidence boundary |
| `detailed_roadmap.md` | Complete 20-stage AMPS work plan and 21 acceptance cases |
| `normative_additions.md` | Sections 16.12/17: physical execution and interface contract |
| `pitch_angle_contract.hpp` | C++17 declarations, scalar/tensor separation and typed diagnostics |
| `compile_contract.cpp` | Compile-only interface smoke check; no physics implementation |
| `check_physics.py` | Independent mathematical/reference checks; never overwrites supplied fixtures |
| `verification_results.json` | Actual supplementary run: 262 passed, zero failed |
| `physical_contract_checks.png` | Exact/reference model predictions; no AMPS simulation data |
| `make_figure.py` | Reproduces the reference figure with numpy/scipy/matplotlib |
| `review_manifest.json` | Input/output digests, evidence and scope |
| `SHA256SUMS` | Digests of package contents |

Run the supplementary checks with Python 3.9+ and numpy/scipy/mpmath:

```bash
python check_physics.py --output verification_results.local.json
g++ -std=c++17 -Wall -Wextra -Wpedantic -Werror -fsyntax-only compile_contract.cpp
```

To repeat the archive-version evidence check as well, add
`--bundle /absolute/path/to/extracted/PITCH_ANGLE_SCATTERING_MODEL_DATA`.
The archive inspected here is v1.0; the v1.1 fixture/registry claims were not
certified. Its earlier 567-value/253-check results are historical checks of
v1.0, not new evidence for v1.1.

The tests evaluate analytical/mathematical contracts and independent references.
They are not the implementation of the header's declared functions, a pusher
convergence campaign, observational validation or an AMPS integration test.
The roadmap specifies those outstanding gates and leaves optional blocked
physics explicit.

The first production scope should be a selected, certified prescribed-scattering
subset. Data/specification pinning, scalar/tensor measure, frames/clocks,
trajectory streaming, physical-time/RNG ownership and geometry lifecycle are
part of correctness, not optional performance details.
