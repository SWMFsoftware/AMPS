# Supporting M3 photon registry and total-capture rate checks

This runner exercises the **production** PIC-independent libraries
`src/models/chemical_database/Photochemistry.*` and `ChargeExchange.*`, through
the compatible `Moon::Species` namespaces. Its contract is in
[CONTRACT.md](CONTRACT.md). It configures no AMPS application, changes no generated
build tree and activates no particle chemistry. A local PASS does not qualify
molecular/CEX daughter transport or replace linked U/I predecessors.

From the AMPS repository root, choose a new output directory:

```bash
/usr/bin/python3 -B srcMoon/test/common/m3_rates/test.py \
  --output-dir test_output/srcMoon/m3-rates-my-new-run \
  --data-root /data/vtenishe/moon_validation_data \
  --cex-reference-dir test_output/srcMoon/m3-charge-exchange-research-20261117
```

The runner compiles a thin adapter, the compatibility translation units and
the shared `DataIO`, `Photochemistry` and `ChargeExchange` production CPPs
with C++17, `-Wall -Wextra -Werror -pedantic`, then records every command/response.
The source receipt includes the shared headers and dependency recipe; rate
formulas, source hashes, domains and declared tolerances are unchanged.
`--cxx` selects another C++ compiler executable. A native AMPS/SPICE/MPI build is
unnecessary. Actual linked/global registration remains with the existing parent
U15/U17 orchestration; this directory creates no new permanent test ID.

`result.json` distinguishes PASS/FAIL/ERROR/SKIPPED. `checks.log`, `build.json`,
`call-*.json`, `photon_inventory.json`, `primary_reference.json`, independent
comparison results and the He source-table diagnostic retain the evidence.
Exit codes are 0/1/2/77 respectively. Missing authoritative D04/CEX sources skip
the corresponding prerequisite group and do not become PASS. Group-level
`setUpClass` outcomes and the count of started methods are reported explicitly.
Canonical data are read-only; all source-corruption fixtures stay in the owned
output directory. Existing output directories are rejected.

## U15 photon-selection statistical extension

The independent photon-only runner adds the predeclared statistical branch
control from [PHOTON_SELECTION_CONTRACT.md](PHOTON_SELECTION_CONTRACT.md):

```bash
/usr/bin/python3 -B srcMoon/test/common/m3_rates/test_photon_selection.py \
  --output-dir test_output/srcMoon/u15-photon-selection-my-new-run \
  --data-root /data/vtenishe/moon_validation_data
```

It compiles only the actual production photon registry and a thin uniform-draw
adapter. No CEX publication is needed. SciPy, NumPy and Matplotlib are explicit
optional dependencies: missing modules produce SKIPPED, without installation.
Each of12 parent/solar families gets100000 draws, master MT seed731993 and a
separate documented stream. All30 positive channel claims use frozen exact
binomial bounds with familywise alpha0.01 allocated by Bonferroni. Native data
QA/provenance/raw/processed identities are checked before execution. The
complete bounds are written before sampled output is generated.

`branch_statistics.csv/.json`, `branch_statistics.png/.eps`, `frozen_bounds.json`,
`data_hashes.json`, native stdout/stderr and build/call manifests retain all
comparisons. These are local numerical selector comparisons, not observations.
This can satisfy U15's local table/sampler statistical requirement when its
companion gates pass; daughter transport and linked I22 remain separate gates.

## Supporting supervision and global status integrity

The maintained protocol controls in
[`test_m3_supporting_protocol.py`](../../test_m3_supporting_protocol.py) exercise
the actual supporting supervisor and global result loader. They verify full
child-tree cleanup, one shared wall interval, separate resource receipts,
status/exit agreement, exact counters and visible child runner/log diagnostics.
The [protocol contract](../../M3_SUPPORTING_PROTOCOL_CONTRACT.md) explains why
formal SKIPPED can retain local PASS, but cannot conceal local FAIL or ERROR.

```bash
/usr/bin/python3 -B srcMoon/test/test_m3_supporting_protocol.py \
  --output-dir test_output/srcMoon/m3-supporting-protocol-my-new-run
```

These checks use tiny temporary process/metadata controls, not AMPS builds,
plasma data or observations. U16 orchestration/global registration imports the
same unittest module; the protocol adds no permanent test ID or physics PASS.

## Frozen photon registry and its limits

The D04 v1 schema pins the four canonical CSV identities **and exact schemas**:
`reaction_channels.csv`, `photo_rates_1au.csv`, `thresholds.csv` and
`cross_sections.csv`. Their original package provenance/QA remains a separate
formal-run prerequisite. A numerically plausible edit cannot be self-approved
by updating a caller hash under the same v1 schema. Hashing and CSV parsing use
the same bytes. Failed loading is transactional; old loaded state and output
arguments survive a refusal.

The full 17-row catalogue retains four atomic, seven H2O and four experimental-OH
reactive channels, and two **nonreactive** total-photoabsorption rows. `OH_exp`
must be selected explicitly; it preserves the archive's experimental `OHe`
compilation and is not a different isotope or a chosen rovibrational population.
Products retain charge, nuclear multiplicity and source electronic-state labels.
The O+H+H branch has two H atoms. O(3P), O(1D) and O(1S) are distinct labels;
plain O and unspecified parent/internal states remain unspecified. Missing
molecular threshold metadata is not changed into a zero threshold. Source tails,
duplicate wavelengths, empty references and printed rounding remain intact.

Rates are the separate quiet/active endpoints at 1 AU and optical depth zero,
in s^-1. Mean excess energies in eV are **not** sampled daughter distributions.
The only rate scaling is inverse-square distance and a supplied normalized gas
photon-visibility fraction. The production geometry adapter must evaluate that
operator at the particle's actual position/epoch; facet incidence, albedo and
temperature cannot stand in for it. Tabulated wavelength/energy extrema describe
archived samples, not a continuous or extrapolatable cross-section model.

Branch fractions condition on the **complete** reactive family. No difficult
channel is dropped/renormalized. Local checks compare every source row and exact
stoichiometry, source domains, whole-family rates, scaling and inverse-CDF
selection. A 100000-point deterministic uniform grid verifies the piecewise-CDF
rounding identity; this does not replace U15's declared stochastic 99% gate or
the linked 95% branching gate.

## CEX sources, units and independent reference

The original research PDFs must be supplied with the exact reviewed identities:

| Original source | Filename and SHA-256 |
|---|---|
| [Barnett 1990 ORNL-6086/V1 A-32](https://www.osti.gov/servlets/purl/6570226), PDF page48; Appendix I | `ornl_redbook_v1_osti.pdf`, `3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030` |
| [Morgan et al. 1985](https://srd.nist.gov/jpcrdreprint/1.555752.pdf?download=1), printed p984/PDF15, Table2 Na column | `nist_morgan_1985_retry2.pdf`, `7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e` |

These are source-supported **total-capture rate laws**, not a newly READY data
package or measured plasma driver. The original numerical tables/coefficients
were checked against the original PDF pages; no plotted curve was digitized.
`primary_reference.py` commits that explicit source-number transcription and
generates a reproducible reference artifact after verifying the original PDFs:

```bash
/usr/bin/python3 -B srcMoon/test/common/m3_rates/primary_reference.py \
  --research-directory test_output/srcMoon/m3-charge-exchange-research-20261117 \
  --output-directory test_output/srcMoon/m3-primary-reference-my-new-run
```

For He, the native kernel uses the original nine-term A0/2 Chebyshev series over
99..11000000 eV/u; the independent oracle constructs T_n directly with 70-digit
Decimal arithmetic. Cross sections convert cm2 to m2 by 1e-4. The approximate
fit's published 10% RMS/18.7% maximum difference from its printed table is retained
as a diagnostic, separate from the <=1% implementation gate. The source quotes
40% accuracy below5keV/u and20% above; those are not Gaussian1sigma errors.
No physical-accuracy claim is upgraded by near-machine-precision code agreement.

For Na, the selected channel ID explicitly names log-log interpolation between
the eleven positive **all-H-state** Table2 knots from1..200keV/u. Exact nodes are
returned unchanged; blanks below1keV/u are missing, not zero. Source units are
10^-16cm2, hence a1e-20 conversion to m2, and energy converts keV/u to eV/u once.
The independent reference evaluates linear interpolation of logarithms with
Decimal arithmetic. No Table3 2s-only values or combined capture+ionization data
are substituted. Invalid/mixed source metadata and out-of-domain queries refuse
without clamping, zero replacement or extrapolation.

The pure API accepts target-rest kinetic energy per projectile u directly.
`ColdBeamRateS` is only the dimensional n*sigma(E)*g wrapper for an explicitly
declared, consistent source-energy/relative-speed pair. It chooses no projectile
mass or classical/relativistic conversion. The test uses the rounded E/g pair
printed in ORNL; its density is an analytic fixture, not solar-wind observations.
There is no Maxwellian integration/tail policy or species-resolved measured bath.

Both CEX metadata records explicitly report that particle execution and daughter
kinetics are unavailable. Total capture does not determine electronic branching,
signed energy defects, angular distributions, ENAs or pickup-ion velocities.
Ne/Ar/alpha channels remain unsupported. Formal full U17/I24 and molecular gates
remain SKIPPED until those production/data/predecessor contracts are satisfied.
