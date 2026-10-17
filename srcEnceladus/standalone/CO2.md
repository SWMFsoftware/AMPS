# Cassini INMS CO₂ calibration and Enceladus plume-source study

## Scope and scientifically defensible result

This document describes the mass-44/CO₂ data, calibration, physics, fitting,
uncertainties, and reproducible results used by the standalone Enceladus plume
application.  Five Cassini INMS flybys are fitted: E3, E5, E14, E17, and E18.

Two materially different data regimes must remain separate.  E14, E17, and E18
use the Perry et al. (2015)/SwRI calibrated mass-44 density product and can be
treated as CO₂ under that product's interpretation.  E3 and E5 are high-speed
encounters for which impacts of heavy ice/organic grains can fragment into mass
44.  Their calibrated channel is therefore labeled **apparent CO₂**.  The
fitted apparent source rate is bracketed to an actual CO₂ rate by dividing by
the published enhancement range 2–5; this is a systematic range, not Gaussian
point noise.

The strongest conclusion is methodological.  A fixed relative source shape
fitted to all five observations permits a controlled temporal comparison, but
the absolute total and even some apparent temporal behavior depend strongly on
the assumed source map and gas kinematics.  Independent 41- and 98-jet
inversions can produce (10^{30}–10^{34}) molecules s⁻¹ while obtaining smaller
residuals.  Those totals are numerical manifestations of a severely
underdetermined inverse problem and are not credible plume production rates.

## Reproducible products

- Calibrated profiles and metadata:
  [`co2_fit_v1`](/data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1).
- Five-flyby three-geometry study:
  [`co2_h2o_tidal_20261015`](/data/vtenishe/cassini_validation_data/studies/co2_h2o_tidal_20261015).
- Release builder:
  [`prepare_co2_fit_data.py`](/data/vtenishe/cassini_validation_data/scripts/prepare_co2_fit_data.py).
- CO₂ fitter:
  [`co2_fitting_study.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/co2_fitting_study.py).
- Joint H₂O reference fit:
  [`joint_fitting_study.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/joint_fitting_study.py).
- CO₂/H₂O and orbital-proxy figures:
  [`plot_co2_h2o_tidal_study.py`](/data/vtenishe/cassini_validation_data/scripts/plot_co2_h2o_tidal_study.py).

All observational paths enter through `--validation-data-root`; changing
machines does not require editing the physics code.

## 1. Instrument observable and data selection

INMS does not count neutral molecules directly.  In closed-source neutral mode,
incoming gas is ionized by electron impact and the quadrupole records detector
counts in mass-to-charge channels.  The mass-44 channel is dominated by CO₂
under the conditions and calibration adopted here, but a channel assignment is
not a proof that every count was created by an intact ambient CO₂ molecule.
Fragmentation and impact chemistry must be treated explicitly.

| flyby | released interpretation | fitted rows | interval/provenance | common calibration scale |
|---|---|---:|---|---:|
| E3 | apparent CO₂ | 44 | reliable closed-source mass-44 samples, ±50 s | 15% plus 2–5 fragmentation bracket |
| E5 | apparent CO₂ | 44 | reliable closed-source mass-44 samples, ±50 s | 15% plus 2–5 fragmentation bracket |
| E14 | CO₂ density | 371 | `Flybys2.mat`, matched to PDS Level-1A rows | 15% |
| E17 | CO₂ density | 377 | `Flybys2.mat`, matched to PDS Level-1A rows | 15% |
| E18 | CO₂ density | 366 | `Flybys2.mat`, matched to PDS Level-1A rows | 15% |

Archive rows marked unreliable are excluded.  No undocumented automated spike
removal is applied.  Possible remaining grain events are therefore an explicit
limitation rather than a hidden preprocessing choice.

## 2. Count-to-density calibration

### 2.1 E3 and E5 apparent CO₂

For a reliable mass-44 sample containing (C_i) high-sensitivity counts, the
adopted Perry et al. conversion is

\[
n^{\rm app}_{{\rm CO_2},i}=C_i
(660 {m cm^{-3} count^{-1}})
=C_i(6.60\times10^8 {m m^{-3} count^{-1}}).
\]

Poisson counting uncertainty is

\[
\sigma_{n,i}=6.60\times10^8
\max(\sqrt{C_i},1) {m m^{-3}}.
\]

The one-count floor prevents a zero-count sample from receiving zero formal
uncertainty.  A separate 15% common multiplicative calibration uncertainty is
not inserted as independent noise at every point.

At the E3/E5 encounter speeds, impact fragmentation of molecules or grains
heavier than 100 u can increase mass-44 counts by approximately a factor
(f_{44}=2–5).  The pointwise channel cannot determine (f_{44}).  Accordingly,

\[
n_{\rm CO_2}={n^{\rm app}_{\rm CO_2}\over f_{44}},\qquad
Q_{\rm CO_2}={Q^{\rm app}_{\rm CO_2}\over f_{44}},\qquad 2\le f_{44}\le5.
\]

This interval is carried as a systematic bracket.  Choosing its midpoint and
calling it a 1σ uncertainty would have no statistical basis.

### 2.2 E14, E17, and E18 calibrated density

The later-flyby density is read from the published SwRI/Perry `Flybys2.mat`
mass-44 `Dens` field.  The associated count sequence is matched exactly to
unique official PDS Level-1A mass-44 rows, which supplies spacecraft position.
The fitted time alignment and maximum absolute match residual are:

| flyby | constant product-to-PDS lag | maximum residual | exact count matches |
|---|---:|---:|---:|
| E14 | +0.286 s | 0.008 s | 371 |
| E17 | −0.018 s | 0.010 s | 377 |
| E18 | +0.194 s | 0.008 s | 366 |

For each nonzero-count sample, the local density-per-count factor is

\[
k_i={n_i\over C_i},
\]

and the point uncertainty is propagated as

\[
\sigma_{n,i}=k_i\max(\sqrt{C_i},1).
\]

This retains the published calibration's sample-dependent response while using
the archived counting statistic.  The 15% absolute calibration uncertainty is
again common-mode.  It scales a fitted total but does not alter the normalized
profile shape.

## 3. Neutral-plume physics

The CO₂ fits call the same production-normalized collisionless plume model as
H₂O, but with the CO₂ molecular mass

\[
m_{\rm CO_2}=44.0095m_u=7.308\times10^{-26} {m kg}.
\]

For source (j),

\[
v_{{\rm th},j}=\sqrt{{2k_BT_j\over m_{\rm CO_2}}},\qquad
s_j={u_j\over v_{{\rm th},j}}.
\]

At point (mathbf x), define

\[
\mathbf R=\mathbf x-\mathbf x_j,quad R=|\mathbf R|,quad
\hat{\mathbf r}=\mathbf R/R,quad
\mu=\hat{\mathbf r}\cdot\hat{\mathbf n}_j,quad
c=\hat{\mathbf r}\cdot\hat{\mathbf a}_j .
\]

With

\[
K_p(c,s)=\int_0^\infty t^p\exp[-t^2+2sct-s^2],dt
\]

and outward-flux normalization

\[
\mathcal N_j=\int_{\mu\ge0}L(\mu)K_3(c,s_j),d\Omega,
\]

the density from a source producing (Q_j) CO₂ molecules s⁻¹ is

\[
n_j(\mathbf x)=
\begin{cases}
\displaystyle {Q_jL(\mu)K_2(c,s_j)\over
\mathcal N_jv_{{\rm th},j}R^2},&\mu\ge0,\\[5pt]
0,&\mu<0.
\end{cases}
\]

The default law is the Tenishev (2010) point-source law (L=1).  The selectable
surface-flux Maxwellian uses (L=\max(\mu,0)), but it was not used for the
reported study.  Tiger-stripe line sources are length-weighted quadratures of
normalized point kernels.

CO₂ is assigned the same source temperature and bulk speed specified for H₂O
in each geometry, while its own mass changes the thermal spreading.  This is a
controlled comparison assumption, not an INMS determination that both species
have identical vent thermodynamics.

| source case | geometry | adopted gas kinematics |
|---|---|---|
| `tenishev2014_hybrid_8` | 8 jets plus 7 tiger-stripe components | source-specific values from the legacy configuration |
| `porco2014_98` | 98 image-derived jet axes | 180 K, 1500 m s⁻¹ |
| `portyankina2022_solar_41` | 41 local-normal gas sources | 180 K, 1500 m s⁻¹ |

The calculation is steady, collisionless, gravity-free, and chemically
loss-free over a flyby.  Collisions and vent/nozzle expansion can affect the
near-field spatial kernel; the inferred (Q) is therefore conditional on this
ballistic model.

## 4. Fitting mechanics

### 4.1 Local-density inverse problem

For flyby sample (i) and source (j), let
(G_{ij}=n_j(\mathbf x_i)/Q_j).  The model is

\[
\hat n_i=\sum_jG_{ij}Q_j+c_0+c_2(r_i/R_E)^{-2},
\qquad Q_j,c_0,c_2\ge0,
\]

where (R_E=252.1) km.  The two broad nuisance backgrounds are not counted in
the vent source rate.  The nonnegative weighted least-squares objective is

\[
\chi^2=\sum_i\left({n_i-\hat n_i\over\sigma_{n,i}}\right)^2 .
\]

Columns are scaled by their weighted norms before solution.  This improves
numerical conditioning but cannot resolve physically indistinguishable jets.

### 4.2 Three fitting strategies

1. `independent_each_observation`: every flyby independently fits all source
   amplitudes and its backgrounds.
2. `shared_shape_e3_e5`: E3 and E5 jointly determine one relative source shape;
   all five flybys then refit one total amplitude for that shape.
3. `shared_shape_all_co2`: all five flybys jointly determine one relative source
   shape; every flyby then refits its own total amplitude.

For a combined shape solution (mathbf q_*),

\[
w_j={q_{*,j}\over\sum_kq_{*,k}},\qquad \sum_jw_j=1,
\]

and target flyby (d) fits

\[
\hat n_{di}=Q_d\sum_jw_jG_{dij}+c_{0,d}+c_{2,d}(r_{di}/R_E)^{-2}.
\]

Each dataset's point uncertainties are multiplied by (sqrt{N_d}) in the
combined shape solve.  The objective is thus a sum of dataset mean-square
normalized residuals rather than being dominated by the longer E14–E18 tables.
The all-five shape is the primary controlled temporal comparison because every
epoch contributes to the map.  It is not necessarily the best fit to each
individual flyby.

## 5. CO₂ results

### 5.1 Primary all-five shared-shape comparison

Rates are apparent CO₂ molecules s⁻¹ for E3/E5 and CO₂ molecules s⁻¹ for
E14/E17/E18.  Numbers in parentheses are χ² per sample.

| source case | E3 | E5 | E14 | E17 | E18 |
|---|---:|---:|---:|---:|---:|
| hybrid 8 + stripes | 5.356e25 (7.23) | 1.013e26 (84.17) | 8.227e25 (31.74) | 6.774e25 (31.11) | 6.385e25 (75.91) |
| Porco 98 | 1.966e27 (10.19) | 2.906e27 (98.93) | 2.400e27 (10.73) | 2.535e27 (28.56) | 2.558e27 (14.28) |
| Portyankina 41 | 9.961e26 (21.20) | 1.557e27 (153.47) | 1.223e27 (39.78) | 1.586e27 (61.26) | 1.212e27 (63.58) |

The high χ²/N values reject counting noise as a complete description of the
residuals.  They also show that one fixed source shape does not reproduce all
five profiles in the current collisionless model.

Within the shared-all shape, later-flyby variation is small for Porco
(E14–E18 span 6.6%) and larger for the hybrid and Portyankina cases (about 29%
and 31%).  Because the relative level and trend depend on the source map, this
cannot yet be claimed as a geometry-independent CO₂ temporal variation.

### 5.2 Independent fits

| source case | E3 | E5 | E14 | E17 | E18 |
|---|---:|---:|---:|---:|---:|
| hybrid 8 + stripes | 2.534e25 | 1.103e26 | 1.939e27 | 2.981e25 | 9.314e26 |
| Porco 98 | 4.009e26 | 1.252e29 | 1.455e32 | 4.971e34 | 7.934e28 |
| Portyankina 41 | 5.005e25 | 7.358e27 | 4.788e30 | 8.766e33 | 3.067e32 |

These totals illustrate why an independent multi-jet solution cannot be read as
a physical global source rate.  The solver can place enormous production into
sources having very weak sensitivity along the trajectory and use a small
fraction of their broad tails to match the data.  Nonnegativity does not prevent
this.  The associated profile residual can improve while the total becomes
physically meaningless.

### 5.3 E3/E5-trained shared shape

| source case | E3 | E5 | E14 | E17 | E18 |
|---|---:|---:|---:|---:|---:|
| hybrid 8 + stripes | 5.303e25 | 1.010e26 | 1.169e26 | 8.095e25 | 6.140e25 |
| Porco 98 | 1.941e31 | 2.156e31 | 1.731e27 | 4.428e28 | 1.592e27 |
| Portyankina 41 | 1.589e26 | 2.864e26 | 2.332e25 | 2.836e26 | 1.006e25 |

Poor transfer to E14–E18 demonstrates that E3/E5 alone do not determine a
portable detailed vent shape.  The E3/E5 high-speed fragmentation issue further
limits physical interpretation of that training strategy.

## 6. H₂O/CO₂ relation and the H₂O proxy

### 6.1 Paired E3/E5 source-rate ratio

Only E3 and E5 provide paired H₂O and mass-44 constraints in this study.  For
paired source rates (H_d) and apparent CO₂ rates (C_d^{\rm app}), the
through-origin proportional ratio is

\[
R_{\rm app}={\sum_d C_d^{\rm app}H_d\over
\sum_d(C_d^{\rm app})^2},\qquad d\in\{E3,E5\}.
\]

Using the H₂O INMS-shape rates and the all-five CO₂-shape rates gives:

| source case | E3 (H/C^{app}) | E5 (H/C^{app}) | (R_{app}) | intrinsic (H/C) bracket after 2–5 correction |
|---|---:|---:|---:|---:|
| hybrid 8 + stripes | 656.0 | 457.3 | 500.7 | 1001–2503 |
| Porco 98 | 1023.3 | 557.7 | 703.9 | 1408–3519 |
| Portyankina 41 | 105.3 | 110.0 | 108.6 | 217–543 |

Because (C^{\rm app}=f_{44}C), the intrinsic ratio is

\[
R={Q_{\rm H_2O}\over Q_{\rm CO_2}}=f_{44}R_{\rm app},
\qquad 2R_{\rm app}\le R\le5R_{\rm app}.
\]

The 23% H₂O and 15% CO₂ common calibration scales would give a 27.5%
relative scale uncertainty if independent,
(sqrt{0.23^2+0.15^2}=0.275), in addition to the much larger fragmentation and
source-geometry systematics.  Two epochs are insufficient for a meaningful
correlation coefficient; none is reported.

### 6.2 CO₂-derived H₂O proxy

For later flyby (d), the conditional water proxy is

\[
Q_{{\rm H_2O},d}^{\rm proxy}=R,Q_{{\rm CO_2},d},
\]

with the full (2R_{app}) to (5R_{app}) bracket retained.  The geometric
midpoint (sqrt{(2R_{app})(5R_{app})}) is used only to place a marker on a
logarithmic plot; it is not a maximum-likelihood estimate.

This proxy is justified only as a conditional transfer within one source case:
CO₂ and H₂O are both plume species, E3/E5 constrain their ratio, and the
same ballistic machinery is applied with species-specific molecular mass.  It
requires the additional assumptions that their source composition and relative
vent partition remain stable and that the E3/E5 impact-fragment enhancement
bracket is applicable.  Enceladus chemistry, phase partition, grain release,
and species-dependent transport can violate those assumptions.  Therefore the
proxy is displayed as a bracket and must not be called a direct H₂O
measurement.

## 7. Figures and machine-readable results

- Cross-case study summary:
  [`study_summary.md`](/data/vtenishe/cassini_validation_data/studies/co2_h2o_tidal_20261015/study_summary.md).
- Hybrid results:
  [`tenishev2014_hybrid_8`](/data/vtenishe/cassini_validation_data/studies/co2_h2o_tidal_20261015/tenishev2014_hybrid_8).
- Porco results:
  [`porco2014_98`](/data/vtenishe/cassini_validation_data/studies/co2_h2o_tidal_20261015/porco2014_98).
- Portyankina results:
  [`portyankina2022_solar_41`](/data/vtenishe/cassini_validation_data/studies/co2_h2o_tidal_20261015/portyankina2022_solar_41).

Each case directory contains observed/model profiles, `co2_source_rates.csv`,
H₂O/CO₂ ratio and proxy tables, diagnostics, and full-precision JSON.  The
study root also contains plots versus mean anomaly, true anomaly, the
second-focus angle, and the adopted crack-opening/tidal proxy.

## 8. Uncertainty and interpretation limits

- **Counting uncertainty:** included sample by sample through the propagated
  Poisson term.
- **Absolute CO₂ scale:** 15% common-mode; it should scale a fitted rate rather
  than be divided by the square root of the number of samples.
- **E3/E5 impact fragmentation:** factor 2–5, retained as a systematic range.
- **Residual grain events:** possible because no undocumented clipping is used.
- **Species identity:** E3/E5 remain an apparent mass-44 CO₂ product.  E14–E18
  use the published CO₂ interpretation but remain tied to that calibration.
- **Kinematics:** identical H₂O/CO₂ vent temperature and bulk speed is a model
  assumption.  Molecular mass is treated correctly in thermal spreading.
- **Source geometry:** the largest uncertainty in absolute totals.  The spread
  among the three discrete source hypotheses is not a Gaussian error bar.
- **Rank deficiency:** many jet kernels are nearly collinear on one flyby track;
  independent individual-jet amplitudes are not physically identifiable.
- **Model discrepancy:** large χ²/N for shared shapes shows missing spatial or
  temporal structure beyond formal counting errors.
- **H₂O transfer:** the later-flyby H₂O series is a composition proxy, not a
  water observation.  It must remain visually and numerically distinct from
  direct H₂O measurements.

## 9. Conclusions

1. E14, E17, and E18 provide calibrated, fit-ready CO₂ density profiles.  E3
   and E5 provide useful apparent-CO₂ profiles but require a factor 2–5
   impact-fragmentation bracket.
2. A shared all-five source shape is the appropriate current product for
   comparing epoch amplitudes under one controlled geometry.  Its large
   residuals show that it is not a complete physical description.
3. Independent 98- and 41-jet totals are non-identifiable and can become
   orders of magnitude too large.  They must not be used as plume production
   measurements merely because their residuals are smaller.
4. The magnitude and temporal variation of CO₂ source rate are source-map
   dependent.  The present data do not select one of the three maps or establish
   a geometry-independent orbital modulation.
5. CO₂ can extend the temporal sampling of H₂O only as a conditional proxy.
   The required E3/E5 composition transfer spans a broad fragmentation and
   geometry-dependent range and should be propagated, not hidden.

## 10. Commands

From `/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone`, after producing the
matching H₂O studies:

```bash
export ENCELADUS_VALIDATION_DATA_ROOT=/data/vtenishe/cassini_validation_data
export H2O_RUN=source_shape_tidal_20261014
export CO2_RUN=co2_h2o_tidal_20261015
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$H2O_RUN/tenishev2014_hybrid_8/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$CO2_RUN/tenishev2014_hybrid_8"
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$H2O_RUN/porco2014_98/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$CO2_RUN/porco2014_98"
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case portyankina2022_solar_41 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$H2O_RUN/portyankina2022_solar_41/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$CO2_RUN/portyankina2022_solar_41"
```

Use a new output directory for a new run.  To rebuild the calibrated release,
run `scripts/prepare_co2_fit_data.py` from the validation-data root as documented
in `README.md`; do not overwrite a released study while changing calibration or
source physics.

## References

1. Waite, J. H. et al. (2004), “The Cassini Ion and Neutral Mass Spectrometer
   (INMS) Investigation,” *Space Science Reviews* 114, 113–231,
   <https://doi.org/10.1007/s11214-004-1408-2>.
2. Perry, M. E. et al. (2015), “Cassini INMS measurements of Enceladus plume
   density,” *Icarus* 257, 139–162,
   <https://doi.org/10.1016/j.icarus.2015.04.037>.
3. Teolis, B. D. et al. (2010), “Detection and measurement of ice grains and gas
   distribution in the Enceladus plume by Cassini's Ion Neutral Mass
   Spectrometer,” *JGR Space Physics* 115, A09222,
   <https://doi.org/10.1029/2009JA015192>.
4. Teolis, B. D. et al. (2015), “A Revised Sensitivity Model for Cassini INMS:
   Results at Titan,” *Space Science Reviews* 190, 47–84,
   <https://doi.org/10.1007/s11214-014-0133-8>.
5. Tenishev, V. et al. (2010), “An approach to numerical simulation of the gas
   distribution in the atmosphere of Enceladus,” *JGR Space Physics* 115,
   A09302, <https://doi.org/10.1029/2009JA015223>.
6. Tenishev, V. et al. (2014), “Effect of the Tiger Stripes on the water vapor
   distribution in Enceladus' exosphere,” *JGR Planets* 119, 2658–2667,
   <https://doi.org/10.1002/2014JE004700>.
7. Porco, C. C., D. DiNino, and F. Nimmo (2014), “How the Geysers, Tidal
   Stresses, and Thermal Emission across the South Polar Terrain of Enceladus
   Are Related,” *Astronomical Journal* 148, 45,
   <https://doi.org/10.1088/0004-6256/148/3/45>.
8. Portyankina, G. et al. (2022), “Modeling the complete set of Cassini's UVIS
   occultation observations of Enceladus' plume,” *Icarus* 374, 114918,
   <https://doi.org/10.1016/j.icarus.2022.114918>.
9. NASA PDS Atmospheres Node, Cassini INMS archive and calibration guidance,
   <https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/inst-inms.html>.
