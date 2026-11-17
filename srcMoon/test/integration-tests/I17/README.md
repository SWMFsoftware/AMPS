# I17 — Time-dependent thermal inertia and eclipse response

**Current implementation checkpoint — 2026-11-16.** The isolated June path is
preserved. An explicit, review-gated December path now exists, with independent
observation-free periodic spin-up, selected-only native input/output, three-dt
convergence, and an independent frozen-criterion scorer. The source/data/native
build identities, scientific model-scope evidence and investigator approval are
checked before a December CSV is parsed for staging. No actual approved
threshold file or December interval is created by the runner. Current formal
I13/threshold/operator prerequisites remain unavailable; these implementation
and protocol tests are not an I17 PASS or Diviner validation. Selected D03
package inputs became READY in the November 12 qualification; the October
package-metadata limitations retained below are historical.

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

## Historical processed-evidence receipt — 2026-10-13

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

`READY` in that historical receipt applied only to the processed I17
transformation. At that checkpoint D03 package README/provenance/QA were
absent; the November qualification subsequently made the selected package
inputs READY. Data readiness still does not supply I13, linked thermal
qualification, an investigator-reviewed acceptance artifact or an untouched
hold-out result. June plots remain experimental training diagnostics until
all corresponding scientific gates exist.

## Reviewed December execution — implementation completed 2026-11-16

The user approved **independent December spin-up**. Its interval must be
declared in the reviewed December `.in`; the runner neither copies June's
final state nor invents another lunation. Keep the material, surface geometry,
forcing and species settings identical to June. Supply these explicit keys:

```text
diviner_eclipse_selection = holdout
diviner_eclipse_freeze_file = /absolute/path/to/committed/reviewed_holdout.json
thermal_spinup_cycle_start_utc = <reviewed UTC start>
thermal_spinup_cycle_end_utc = <reviewed UTC end on 2011-12-10>
simulation_start_utc = <reviewed December initialization UTC>
```

The angle brackets are documentation placeholders, not accepted input values.
The ordinary `thermal_time_step_s`, full-column temperature/enthalpy criteria,
and all original production parameters remain explicit. CSPICE validates UTC
syntax and ordering; the native adapter refuses an observation preceding the
spin-up end. Allocation/runtime failure remains ERROR, without a shortened
spin-up or looser physical/numerical criterion.

[`i17_holdout.py`](../../i17_holdout.py) defines the release schema
`srcMoon-I17-reviewed-holdout-v1`. An investigator creates and separately
commits the artifact after reviewing actual June evidence. The script has no
approve/freeze command. Required content is:

| Field | Meaning and required evidence |
|---|---|
| `status`, `approved_by`, `approved_utc` | Exactly `APPROVED_FROZEN`, a named reviewer and explicit UTC timestamp. A candidate or missing record stays SKIPPED. |
| `frozen_before_holdout_residuals` | Must be true; membership/parameters/limits are frozen before December residual inspection. |
| `holdout_campaign_id` | Basename of the one planned aggregate output directory. Another directory cannot reuse that release automatically. |
| `initialization` | Exactly `independent-observation-free-periodic-spinup`, the user-approved method. |
| `processed_history_sha256`, `selection_contract_sha256`, `processor_sha256` | Exact canonical processed history, unchanged selection contract and producing script. |
| `executable_sha256`, `production_source_identity_sha256` | Exact retained executable and compact/sorted JSON hash of authoritative native build-input hashes. The June/December binary and native source must match. |
| `training_evidence` | Hash-bound `input`, `candidate`, `result`, `manifest`, `model_csv`, and `convergence_csv` records, each with `path` and `sha256`. June must actually pass the unchanged numerical gate. |
| `reviewed_holdout_input` | `path`/`sha256` of the exact independently reviewed December receipt, including its explicit spin-up interval. |
| `model_scope_qualification` | READY plus hash-bound evidence of the physical approximation and observation-operator scope. A scalar threshold cannot qualify absent atmospheric transmission or an unresolved terrain/footprint operator. |
| `acceptance` | Explicit reviewed limits listed below; never filled from December output or assigned by the script. |
| `phase_amplitude_groups` | Frozen named groups of exact `site_time_bin_ids`, at least three unique bins per group, established without December model residuals. |

The freeze source must be tracked and byte-identical to its own HEAD Git blob;
the enclosing repository, revision and blob identity are retained. The review
may remain in its original committed checkout during private native builds.
All Git operations are read-only. An untracked/edited candidate cannot become
approved merely because its JSON says READY. The complete reviewed hold-out
input is byte-bound. A second identity check permits only declared epoch,
dt-refinement and file-transport differences from the June receipt; changes to
material/geometry/forcing remain errors. JSON file paths resolve relative to
the review artifact, while `.in` paths keep the native input-file convention.

Supported acceptance fields are `weighted_rmse_k_max`,
`absolute_weighted_bias_k_max`, `normalized_residual_rms_max`,
`coverage_sigma`, `minimum_coverage_fraction`,
`sampled_phase_error_s_max` and `sampled_amplitude_error_k_max`. All must be
explicit, finite and nonnegative; coverage fraction lies in [0,1]. They have
**no defaults**. A reviewer must also explicitly select
`phase_amplitude_operator = sampled-extrema-on-frozen-bin-groups` before that
implemented operator can be used; its scientific suitability/group membership
has not been approved by implementing the backend.

For each group the scorer requires the same native model latitude/longitude
and distinct actual epochs. Sampled amplitude is max(T_bol)-min(T_bol) on
those exact bins; sampled phase is the epoch of the unique sampled minimum.
The model/reference differences are compared to the predeclared limits. No
interpolation, gap filling, curve fitting or post-hoc site selection occurs.
These metrics describe sampled extrema, not an assumed continuously observed
eclipse minimum. An ambiguous observed minimum is SKIPPED because phase is
unidentifiable; a model with no unique sampled minimum FAILs that declared
test. Missing/mismatched bins are ERROR. Absent reviewed operator, group,
scope or threshold evidence remains SKIPPED before launch.

The release is claimed once under `_m2_validation_data/D03-RDR-I17/` before
December science staging. The claim covers the predeclared dt, dt/2 and dt/4
executions, not three independently tuned observational experiments. A failed
or interrupted campaign preserves the claim and artifacts; another execution
requires investigator review and a newly authorized campaign/release. The
runner does not append seeds, modify limits or silently inspect December
again. Fine-step model CSV, comparison CSV/PNG/EPS and channel diagnostics are
preserved; no observational comparison is produced if numerical convergence
fails first.
Reviewing an infrastructure retry does not permit changing frozen membership,
physics or acceptance limits after hold-out inspection. An evaluated scientific
FAIL and the original observations remain evidence, rather than an invitation
to select a more agreeable case.

After review/hash qualification and selected-only staging, the runner writes
`holdout_native_authorization.txt` and supplies its generated path through
`diviner_eclipse_authorization_file`. This optional input transport key is
added to the derived receipt if absent. Do not hand-author the sidecar as an
approval substitute. Production `MoonI17Authorization::Validate` checks its
exact schema, selection/status, staged history/review/validation paths,
SHA-token shape and committed-review token before the native adapter opens
December observations. Full JSON approval and actual file SHA-256 checks
remain in Python and the run manifest: this is an accident guard, not a
cryptographic signature or separate qualification. Missing, duplicate, extra,
wrong-selection/path/token or unavailable-file records stop the native run.
The same production header is exercised by a compiled thin protocol adapter
inside `test_holdout.py`; its temporary files contain no science values.

Once all scientific predecessors, data and review evidence actually exist:

```bash
python3 -B srcMoon/test/integration-tests/I17/test.py \
  --amps /absolute/path/to/qualified/amps \
  --input /absolute/path/to/reviewed-december.in \
  --output-dir test_output/srcMoon/<approved-campaign-id> --mpi-ranks 8
```

The same-campaign I13 PASS must already be present in that qualification
directory. Use the normal aggregate runner to construct the predecessor chain;
this individual command does not bypass it. With unavailable prerequisites
the command emits a truthful SKIPPED result and does not launch the model.

### Protocol tests registered in the individual/global runner

`test.py` always runs both `test_infrastructure.py` and `test_holdout.py`
before delegating the scientific ID, so the global runner includes the new
coverage without assigning another U/I number. Run their inexpensive checks
directly with:

```bash
python3 -B srcMoon/test/integration-tests/I17/test_infrastructure.py
python3 -B srcMoon/test/integration-tests/I17/test_holdout.py
```

They cover pre-CSV/pre-plot rejection, source/data/build/input hashes, approval
and campaign binding, independent explicit initialization, unchanged physical
parameters, one-time release, selected-only staging, tamper detection, matched
three-dt bins, and frozen arithmetic phase/amplitude/coverage decisions.
All tiny CSV/JSON values live in temporary protocol fixtures. They are never
observations, comparison figures or thermal/Diviner validation evidence. Actual
December observations were not scored to implement or test this control path.
