# I17 — Time-dependent thermal inertia and eclipse response

## Status and scope

I17 is the linked numerical/observational gate for the production
depth-resolved thermal state. It is not established by U11 solver verification,
successful RDR preprocessing, or I16 climatological-map agreement. A formal
`PASS` requires a same-campaign qualified I13 `PASS`, a READY D03 package, the normal
terrain/illumination/thermal path, time-step convergence, and the frozen
December hold-out comparison.

The test keeps two physical quantities distinct:

- `physical_surface_temperature_K` is the kinetic boundary temperature solved
  by the Hayne-type conduction column.
- `bolometric_brightness_temperature_K` is the blackbody temperature having
  the same spectrally integrated emitted flux as the observed/modelled surface.

For the production grey thermal state, `M = epsilon sigma T^4` and therefore
`T_bol = epsilon^(1/4) T`. Comparing kinetic `T` directly with Diviner `T_bol`
would silently set emissivity to unity and is prohibited.

## Frozen assumptions, physics, and implementation decisions

The complete machine-readable selection is
[`reference/selection_contract.json`](reference/selection_contract.json). Its
hash is written into every processed-data and linked-run receipt. These choices
are frozen before any December AMPS residual is examined.

1. **Evidence.** Use NASA PDS4 calibrated Diviner RDR products identified by
   the Diviner team as RDR processing version 3.1. Native XML labels control
   fields, units, UTC, coordinates, missing constants, activity, and quality
   flags. Directory presence is not readiness. The processed receipt records
   every consumed ZIP/XML SHA-256 and product LID, the processor SHA-256, the
   frozen-contract SHA-256, and the output SHA-256.
2. **Training/hold-out.** 2011-06-15 is training for operator/QA and threshold
   qualification only; it is not used to fit the Hayne thermal parameters.
   2011-12-10 is an untouched one-pass hold-out. Windows extend at least 30
   minutes before/after NASA GSFC P1/P4 contacts and are rounded outward to
   ten-minute RDR boundaries.
3. **Footprints.** Each channel uses its native FOV-centre latitude, east
   longitude, UTC, detector, and emission angle. The two focal planes are only
   approximately boresighted. Co-registration means membership in the same
   predeclared 0.5-degree site and 600-second UTC bin; channel-specific mean
   footprints/times remain in the product.
4. **Quality.** Retain only activity flag `110` and emission angle below 15
   degrees. Require `qca=0`. Retain geometry bits 2/3 because they mean
   definitive but not reprocessed pointing/ephemeris; reject predictive or
   missing geometry. Retain only the physical `qmi` eclipse bit; reject
   transient, abnormal thermal state, drift, noise, saturation, and moving
   records. Missing/invalid constants are never interpolated.
5. **Aggregation.** Average valid calibrated radiances before conversion to
   brightness temperature. Require at least three real records in every C3–C9
   channel. Preserve gaps; never smooth, extrapolate, gap-fill, or independently
   renormalize cases.
6. **Native calibration.** RDR supplies paired calibrated radiance and
   brightness temperature. The processor reconstructs and QA-checks that native
   monotone inverse. It does not substitute a band centre or use AMPS output.
7. **Bolometric operator.** Use Paige et al. (2010) C3–C9 partitions at 0,
   8.075, 8.40, 13, 25, 50, 100, and 1000 micrometres. Sensitivity floors for
   C3–C8 are 190, 190, 180, 95, 60, and 40 K. Replace a below-floor channel by
   the next longer-wavelength channel above its own floor; C9 is required.
8. **Uncertainty.** Use Diviner Radiometric Test Report equations 1–3 and
   Tables 4/10. Combine the documented calibration terms with real within-bin
   radiance dispersion, then propagate through native radiance conversion and
   nonlinear T_bol with an independent-channel first-order Jacobian. Numerical
   derivatives remain on the current published sensitivity/positivity branch:
   a substituted below-floor C3–C8 channel has zero local derivative, and an
   exact-threshold value uses the valid-branch right derivative. This prevents
   a cold C9 one-sigma interval from being misinterpreted as a non-physical
   negative Planck temperature. No assumed fixed kelvin error is allowed. The
   unavailable cross-channel covariance remains an explicit QA limitation and
   is not replaced by invented correlations.
9. **Diagnostics.** Preserve C6–C9 radiance and brightness-temperature
   residuals, footprints, times, detector coverage, and uncertainties for
   anisothermality/emissivity diagnosis. The model has one grey emissivity, so
   its documented rectangular-half-maximum channel approximation is diagnostic
   only and excluded from the primary score.
10. **Thermal state.** All thermophysical parameters remain explicit in the
    `.in` receipt and are not adjusted to improve either eclipse. Spin-up uses
    no Diviner temperature and must demonstrate periodic convergence before an
    observational score can be promoted from `SKIPPED`.
11. **Direct eclipse irradiance.** The first-comparison model computes the
    analytical overlap of circular apparent Sun/Earth discs at every
    body-fixed surface point and SPICE epoch. The Sun has uniform projected
    radiance; Earth is an opaque sphere. Input radii are 695,700,000 m (IAU
    nominal solar radius) and 6,371,000 m (AMPS mean solid-Earth radius).
    The unocculted area fraction scales the direct flux evaluated with the
    solar-centre incidence cosine; incidence is not integrated across the
    0.27-degree solar disc. Terrain shadowing remains a separate all-or-none
    BVH ray to the solar centre, so terrain penumbrae are not modeled. Solar
    limb darkening and refracted, scattered, or transmitted
    terrestrial-atmosphere light are absent, not fitted. NASA eclipse material
    establishes that the latter is physically nonzero, so plots produced under
    this approximation are explicitly `EXPERIMENTAL / NOT VALIDATED`.
12. **Periodic initialization.** The input interval
    2011-05-17T04:05:57.024Z--2011-06-15T16:50:00.000Z has an explicit
    2,551,442.976 s mean-synodic duration and ends at the first June training
    bin. Replaying one ephemeris interval assumes the forcing is periodic even
    though real lunations differ. Every one of the 60x100 logical-surface
    columns retains all 80 finite-volume layers. Replay stops only after the
    maximum surface/layer boundary change is <=0.01 K and the maximum
    whole-column sensible-enthalpy change is <=10,000 J m^-2. These controls
    are initialization tolerances, not observational acceptance limits. There
    is no maximum cycle count; allocation/job exhaustion is ERROR.
13. **First-plot isolation.** Until a June threshold is reviewed and frozen,
    the campaign stages a training-only history. The linked executable rejects
    any non-training row and records `holdout_rows_accessed=0`. It evaluates
    the production kinetic temperature at each declared 0.5-degree bin centre
    and mean UTC, applies `epsilon^(1/4) T`, and retains C6--C9 rectangular-band
    diagnostics. Channel-specific native coordinates/times/detectors remain in
    the immutable processed history; the primary model value uses the
    predeclared aggregate bin centre/time and does not pretend to reconstruct a
    wavelength-dependent sub-footprint temperature field.

## Complete test specification

### Approved numerical error budget (2026-10-14)

Before computing any June model-minus-observation residual, run the production
history at the input thermal step, half that step, and one quarter of that
step (currently 300, 150, and 75 s). Each run performs its own full-column
spin-up with unchanged physical parameters and spin-up tolerances. For every
identical site/time bin, require the fine pair's T_bol difference to be strictly
below 10% of its documented observational one-sigma uncertainty and smaller
than the coarse pair's difference. Two exactly zero differences count as zero
measured sensitivity. Require <50 K per production update at all three steps.
An evaluated failure is FAIL; missing or mismatched bins are ERROR. Neither
condition produces observational residuals or comparison figures.

The driver preserves `dt_div_2/` and `dt_div_4/` inputs, manifests, logs, model
tables, and resource records. `thermal_time_step_convergence.csv` records each
bin's differences and error budget. After this gate passes, the finest run's
comparison CSV/PNG/EPS and threshold candidate are written in `dt_div_4/`.
The overall I17 scientific status remains SKIPPED pending the separate
observational threshold freeze and hold-out validation. This time-step gate
does not demonstrate spatial convergence or remove the documented physical
approximations.

**ID / name:** I17, time-dependent thermal inertia and eclipse response.

**Purpose:** verify linked thermal evolution through changing illumination and
Earth eclipse, establish update-size convergence/no discontinuity, and validate
the frozen hold-out bolometric history without hold-out fitting.

**Production capability / entry points:** normal srcMoon initialization;
terrain/Earth-eclipse forcing in `MoonSurfaceTemperature.cpp`;
`Moon::Thermal::Column::Advance`; `Moon::DivinerObservation`; normal
synchronized epoch advancement.

**Configuration:** SPICE-enabled Moon build, eight MPI ranks, approved-campaign
I13-qualified input, `surface_temperature_mode=thermal-model`, explicit Hayne
parameters, and paired history/selection-contract paths. MPI/OpenMP, seed,
revision, executable, input, kernels, data, and contract hashes are retained.

**Inputs / data:** same-campaign I13 result; D03 PDS4 RDR v3.1 ZIP/XML products
for the frozen windows; RDR SIS; Paige supplementary algorithm; radiometric
calibration report. Raw files are immutable.

**Independent reference:** U11 analytical conservation/convergence is the
solver oracle. Native calibrated Diviner radiance plus an independent published
Paige implementation is observational evidence. Production output is never
treated as its own physical truth.

**Units / frames / time:** model SI; RDR radiance W m^-2 sr^-1; temperature K;
native planetocentric latitude and positive-east longitude; observation SCET
UTC, internal SPICE TDB, model MOON_ME_DE421. Converted epochs and footprint
centroids are saved.

**Deterministic controls:** thermal evolution and preprocessing are
deterministic. MPI layout, OpenMP layout, and AMPS seed are still recorded.

**Procedure:** verify I13 and D03; run the committed processor and verify hashes;
advance the unmodified production state after model-only periodic spin-up;
repeat at smaller declared thermal steps; apply `epsilon sigma T^4` at actual
bin footprints/epochs; qualify thresholds with June; freeze them; score
December once; preserve row-level CSV and render unchanged PNG/EPS plots.

**Metrics:** maximum per-update temperature change; successive-dt differences;
June/December weighted bias, RMSE, and uncertainty coverage for
`bolometric_brightness_temperature_K`; normalized residuals; C6–C9 diagnostics;
native-flag exclusions; bin/detector counts; resource history.

**Acceptance:** no unphysical discontinuity; less than 50 K per update;
time-step convergence; phase/amplitude within a Diviner tolerance derived and
frozen from June plus documented uncertainty before December is exposed to the
scorer. Until the threshold artifact is frozen, scoring is `SKIPPED`, not PASS.

**Expected failures / status:** missing predecessor, data, channel bins,
spin-up, or threshold freeze is `SKIPPED`; corrupted ZIP/label/hash or build,
parser, and runtime failures are `ERROR`; only an executed frozen numerical or
physical mismatch is `FAIL`.

**Predecessors:** applicable U09–U11 and same-campaign qualified I13.

**Artifacts:** logs, exact input, receipt, result/manifest, D03 QA, contract
snapshot, per-dt histories, comparison CSV/PNG/EPS, C6–C9 diagnostics, hashes,
and resource metrics. A skipped test creates no fake/empty comparison plot.

## Commands

### Launcher and plotting infrastructure contract (2026-10-14)

Before a long linked run, `test_infrastructure.py` checks that both the MPI
command and its provenance use the derived `application.in`, rather than the
original template whose observation paths are `none`. It also checks that the
plot reader accepts the explicit `training-experimental` selection and rejects
unknown selections. These are isolated infrastructure checks, not I17 physics
or observational evidence. Their temporary CSV records are parser fixtures;
they are never rendered or retained as comparison figures.

The overview uses unconnected markers for experimental training samples:
consecutive times can belong to different lunar sites and must not be joined
into an apparent cooling curve. Site/time identities remain in the linked
model CSV for subsequent site-specific figures.

```bash
python3 -B srcMoon/test/integration-tests/I17/test_infrastructure.py
```

```bash
python3 -B srcMoon/test/data_processing/prepare_d03_rdr_eclipse.py \
  --output-dir /data/vtenishe/moon_validation_data/diviner/processed/i17

python3 srcMoon/test/run_all_test.py \
  --campaign-mode qualification --mpi-ranks 8 \
  --only I11 --only I12 --only I13 --only I17
```

If a predecessor or D03 qualification is missing, no linked model or plot is
launched. Once those gates pass, the first execution may produce the June
experimental CSV/PNG/EPS and threshold candidate, but it still reports
`SKIPPED / NOT VALIDATED` until that candidate is reviewed and frozen.
The processor command is an explicit data-maintenance action; `run_all_test.py`
does not generate or alter canonical processed data.  Once qualified, the
runner verifies the QA hashes and stages byte-identical history, contract, QA,
and processor logs under the unique linked campaign directory. The processor
writes the science CSV under a temporary name, atomically publishes it only
after all bins finish, and publishes the QA receipt last; a killed run cannot
leave a new file/receipt pair that appears READY.

## Current processed-evidence receipt

The 2026-10-13 processor execution examined 77,075,901 native records from 87
immutable ZIP/XML archive pairs, retained 25,848,333 records under the frozen
cuts, and produced 6,353 complete C3--C9 bins: 2,039 training and 4,314
hold-out.  It rejected 287 incomplete site/time bins and rejected no bin for a
non-positive mean C9.  The resulting temperature range is 18.5316--306.4471 K
and the propagated one-sigma range is 0.4231--31.6302 K; every published value
is finite, every temperature is positive, and every uncertainty is
non-negative.  These ranges are integrity diagnostics, not acceptance limits.

The canonical receipt is `READY` for this deterministic transformation and
records these SHA-256 values:

- history CSV: `870577c82b430b098d816d210b528c714ec6cfa3a20ae9be3f6482a3639ae182`;
- frozen contract: `88fc5196fee5863691937a505873d420990b98056db8a1e8167255cd736189b0`;
- processor: `5a5d29d0f6dba88cd49c899e541a980eb84ce24717481b0cad471f8012459a65`.

`READY` here applies only to the processed I17 transformation. D03 as a
dataset package is still not READY because its package-level README,
provenance, and QA records are absent.  I17 also still lacks the linked
June execution evidence and a reviewed/frozen numerical threshold artifact.
The production spin-up and first-plot code now exist, but no plot is created
while I13 or D03 readiness blocks launch. When those prerequisites pass, June
plots are experimental training diagnostics and I17 remains
`SKIPPED / NOT VALIDATED`; the real processed observations are not a substitute
for predecessor, threshold-freeze, or untouched-holdout gates.
