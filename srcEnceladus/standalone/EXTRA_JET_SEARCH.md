# E17 extra-jet search

## Purpose and status

The E17 INMS mass-44 profile contains a narrow feature near closest approach
that the fitted Porco 98-jet shared spatial pattern does not reproduce. The
search asks a limited inverse question: *where would one additional
Porco-like gas source have to be placed to reduce that residual under the
current steady, collisionless plume model?* It does not assume beforehand
that the residual is a new vent, and it does not establish a geological jet
detection.

The selected definition is stored in `e17_extra_jet_candidate.json`. It is a
development alternative to `porco2014_98`, not part of any default source
case. The model contains 98 Porco sources unless the user explicitly passes

```bash
--extra-jet e17_extra_jet_candidate.json
```

The candidate is at planetocentric latitude -88.6656758 degrees and west
longitude 296.3659692 degrees on the modeled Baghdad Sulcus polyline. Its axis
is radial. It lies 0.545856 degrees from Porco source 42, so even a real local
enhancement would not necessarily imply a spatially distinct vent.

## Why the search was performed

The baseline profile used one relative shape fitted jointly to E3, E5, E14,
E17, and E18 CO2 data. A structured E17 closest-approach residual remained.
Because Cassini crossed a narrow spatial region rapidly, a compact source
whose ballistic signature intersects the E17 track is a plausible diagnostic
basis function. Other plausible explanations remain: time-variable strength
of a cataloged source, imperfect source direction, grains or heavy-molecule
fragmentation in mass 44, inlet response, missing collisions, trajectory-scale
model error, or statistical fluctuation.

## Screening inverse problem

For a trial candidate geometry `theta`, the modeled E17 density is

```text
n_i = A S_i + B J_i(theta) + C + D (r_i/R)^alpha,
```

where `S_i` is the fixed unit-rate signature of the established 98-source
shared shape, `J_i` is the unit-rate signature of the trial source, and the
four nonnegative coefficients are re-fitted at every trial location. The
objective in the declared closest-approach window is

```text
chi2(theta) = sum_i [(n_i - d_i)/sigma_i]^2.
```

Holding the 98-source *relative* pattern fixed during this screening step
prevents its 98 ill-conditioned amplitudes from rearranging separately for
every candidate position. The search therefore measures whether one new
signature complements the established pattern; it is not yet the final
99-source inversion.

The checked search performed five related tests:

1. A global radial-source search over the south-polar latitude range occupied
   by the Porco catalog and all west longitudes.
2. A radial search constrained to the modeled tiger-stripe polylines. This is
   the declared primary search because it incorporates the geological prior
   that plume vents occur on the active fissure system.
3. Repetition with 1.5, 2.5, and 4.0 s half-windows around closest approach.
4. Free-direction searches with tilt limited to the Porco catalog envelope,
   used to diagnose the location/direction degeneracy.
5. Fits at the selected location to E3, E5, E14, E17, and E18, both with an
   independently fitted candidate rate and with the E17 candidate-to-baseline
   ratio transferred unchanged.

The location was selected using E17. Its E17 improvement is consequently a
training result affected by the look-elsewhere search and is not independent
validation.

## Complete 98-versus-99 refit

After selecting the location, the candidate was appended as source 99 and the
same complete procedures used for the original Porco case were repeated:

- independent and shared-shape H2O fits for E3 and E5;
- independent and shared-shape UVIS fits for the 2007 and 2010 occultations;
- independent, E3+E5-shape, and all-CO2-shape fits for E3, E5, E14, E17, and
  E18; and
- a direct comparison of every 98-source and 99-source target fit.

The comparison contains 31 observation/strategy combinations. In independent
fits, E17 CO2 improves by Delta chi-square = 3.40154. E3, E5, E14, and E18 CO2
and E3/E5 H2O are unchanged to numerical precision. UVIS 2010 improves by
3.89155, whereas UVIS 2007 worsens by 0.09343. Shared-shape fits assign
essentially zero production fraction to the candidate. The complete result
does not require a persistent 99th source.

Negative improvements occasionally reported for nested shared-shape fits are
optimizer-sensitivity diagnostics. In exact arithmetic a 99-source
nonnegative training fit cannot be worse than the nested 98-source model,
because zero candidate amplitude recovers the original model. These cases are
not physical evidence against the candidate.

## Reproduce the search on another machine

Set the validation-data root once. The tree must contain
`data/processed/co2_fit_v1`, `data/derived/porco2014_jet_catalog.json`, and the
baseline study result named below.

```bash
cd /path/to/AMPS/srcEnceladus/standalone
export ENCELADUS_VALIDATION_DATA_ROOT=/path/to/cassini_validation_data

python3 extra_jet_e17_study.py \
  --config config/tenishev2014_stage0.json \
  --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" \
  --co2-baseline-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/co2_h2o_tidal_20261015/porco2014_98/study_results.json" \
  --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/my_e17_extra_jet"
```

`--jet-catalog` and `--co2-data-root` may be supplied explicitly to override
their conventional locations beneath the root. Output directories must not
already exist.

To run the complete 99-source H2O/UVIS and CO2 refits, pass the local candidate
file to both fitting programs:

```bash
python3 joint_fitting_study.py \
  --config config/tenishev2014_stage0.json \
  --source-case porco2014_98 \
  --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" \
  --extra-jet e17_extra_jet_candidate.json \
  --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/my_99_source_h2o_uvis"

python3 co2_fitting_study.py \
  --config config/tenishev2014_stage0.json \
  --source-case porco2014_98 \
  --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" \
  --extra-jet e17_extra_jet_candidate.json \
  --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/my_99_source_h2o_uvis/study_results.json" \
  --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/my_99_source_co2"
```

Omitting `--extra-jet` always restores the published 98-source basis.

## Interpretation limits

- The search assumes steady, collisionless, gravity-free propagation.
- The candidate temperature and speed, 180 K and 1.5 km/s, are controlled
  comparison values rather than an independent measurement at the candidate.
- One flyby track cannot uniquely determine latitude, longitude, tilt, and
  azimuth; the window and free-direction tests demonstrate this degeneracy.
- Counting uncertainties do not include the common calibration scale or model
  discrepancy, so Delta chi-square is diagnostic rather than a discovery
  significance.
- Mass 44 is treated as apparent CO2 where high-speed fragmentation can
  contribute. A localized mass-44 residual is not uniquely compositional.
- The candidate should remain optional unless independent observations or a
  better physical model provide reproducible evidence for it.
