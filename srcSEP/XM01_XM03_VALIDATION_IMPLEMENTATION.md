# XM01–XM03 cross-model and observational validation

## Evidence boundary

XM01 is a linked cross-solver benchmark: the selected AMPS executable advances
production focused-transport characteristics and an independent finite-volume
PDE solver supplies the reference. XM02 is a publication-informed controlled
first-passage comparison for three M-FLAMPA mean free paths. XM03 is now a
linked observational comparison: an event-informed one-field-line srcSEP model
is scored against Earth measurements in Liu et al. Figure 12.

All three numerical stages execute through the `--test XM0N` registry of the
binary passed with `--amps`. Standalone callback builds are compile/runtime
smoke tests only and cannot produce scientific PASS evidence.

## XM03 reference extraction

Liu et al., *Physics-based Simulation of the 2013 April 11 Solar Energetic
Particle Event*, ApJ 985:82 (2025), DOI `10.3847/1538-4357/adc4e3`, provide
Figure 12 as a vector PDF in arXiv source `2412.07581v2`.
`validation/reference/digitize_liu_figure12.py` converts that PDF to SVG and
reads the actual Matplotlib path coordinates. It extracts:

- 80 ACE/EPAM, GOES-13/EPEAD, and SOHO/ERNE Earth measurements from panels
  (a–c), including each channel's low/high energy bounds; and
- 97 points from the Earth-connected shock thermal-energy-density trace in
  panel (d).

Only the first product is a scored observational reference. The second is a
source-time input because the publication explicitly makes particle injection
proportional to shock thermal energy density. STEREO-B measurements belong to
a different field line, while plotted SOFIE curves would make the case a
model-to-model comparison, so both are excluded.

The digitizer checks the reviewed observation count and records figure/PDF
hashes, vector calibration, selections, units, and limitations in
`validation/cases/XM03/reference/provenance.json`.

## XM03 linked model

The paper does not publish its evolving Earth-connected AWSoM line, complete
SOFIE restart, or shock surface. XM03 therefore makes a reduced calculation
whose substitutions are explicit rather than fabricating missing global-model
arrays:

1. The paper's 363 km/s Earth wind and the Carrington 25.38-day sidereal rate
   define a constant-speed equatorial Parker spiral from 2.5 solar radii to
   Earth at 1 AU.
2. A shock moves radially from the inner radius at the paper's 675 km/s EEGGL
   input speed. Release begins after the reported 15-minute Earth connection.
3. Release times are sampled from the exact piecewise-linear Figure 12(d)
   Earth trace. Injection follows the reported 10 keV threshold,
   `f(p) proportional to p^-5`, and factor 1.2.
4. `SEP::Transport::AdvanceParker` advances each proton with
   `kappa_parallel=lambda_parallel*v/3` and the paper's
   `lambda_parallel=0.3 AU*(r/AU)*(pc/GeV)^(1/3)` law.
5. The background supplies Parker-tangent wind projection and spherical
   `div(U)=2U/r`, so the production core applies adiabatic momentum loss.
6. Sunward particles crossing 2.5 solar radii are absorbed; outward first
   passage at 1 AU supplies the modeled flux spectrum in two-hour windows
   around 4, 12, and 36 hours after the 07:24 UTC flux-rope launch. The Figure
   12(d) source table instead retains its 06:00 UTC civil-time origin; the
   native callback applies the 5040 s offset explicitly.

The IAU nominal solar radius and exact astronomical unit are used only as unit
conversions. The Parker geometry, constant-speed shock, and boundary policies
are reduced-model choices—not claims about unpublished SOFIE values. Their
rationale and the missing artifacts are machine-readable in
`validation/cases/XM03/publication_input.json`.

## XM03 score and plots

The global simulation's shock area and the one-dimensional connected flux-tube
area are not published, so absolute pfu normalization is not identifiable.
The scorer applies one log-least-squares multiplicative amplitude across every
scored time, energy, and instrument. It never fits per-panel, per-energy, or
per-instrument factors. Formal scoring begins at 1 MeV because the paper warns
about low-energy contamination; all 80 measurements remain in evidence and in
the plot.

Outputs include the raw relative model spectrum, copied observation CSV,
pointwise scaled comparison, resolved input, native manifest/JSON/JUnit/log,
provenance, and PNG/EPS overlays. The registered gates cover point coverage,
global log-RMSE, median absolute log error, and log-intensity correlation.
They are provisional reduced-model validation criteria, not a statement that
srcSEP reproduces the original global SOFIE run.

## Execution

```sh
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir /evidence/XM01-XM03
```

XM02 and XM03 reject `--case-input` so archived commands always refer to the
one source-reviewed reconstruction. Neither accepts an external production
CSV. A missing linked executable or missing registered source/reference is an
error; XM03 no longer skips for a missing `model/srcsep_output.csv`.
