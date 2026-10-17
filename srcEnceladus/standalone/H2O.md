# Cassini H₂O calibration and Enceladus plume-source study

## Scope and scientifically defensible result

This document describes the H₂O observations, physical model, inverse method,
uncertainties, and reproducible results used by the standalone Enceladus plume
application.  It covers two independent measurement types:

1. local, sticking-corrected Cassini INMS H₂O number-density profiles from
   E3 and E5; and
2. line-of-sight H₂O constraints from the 24 October 2007 stellar and 18 May
   2010 solar UVIS occultations.

The calculations constrain effective source rates conditional on a specified
source map and collisionless gas kernel.  They do not uniquely recover the
activity of every vent.  E3 and E5 independent fits with the legacy eight-jet
plus tiger-stripe geometry give (2.06\times10^{28}) and
(4.52\times10^{28}) H₂O molecules s⁻¹.  The independently fitted UVIS
events give (7.93\times10^{27}) and (4.80\times10^{27}) s⁻¹.  Those
differences are model- and observation-conditioned: the current data do not
support interpreting them as an instrument-independent measurement of global
plume variability.

The 98-jet independent INMS inversions produce rates as high as
(10^{30}–10^{31}) s⁻¹ while fitting the profiles well.  This is an inverse-
problem degeneracy, not evidence for such large production.  A one-dimensional
track cannot identify tens of nearly collinear source kernels.  Shared-shape
fits, rank diagnostics, and comparisons among source maps must therefore be
examined together.

## Reproducible products

- Processed observations and metadata:
  [`plume_fit_v1`](/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1).
- Three-geometry fit outputs:
  [`source_shape_tidal_20261014`](/data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014).
- INMS profile construction:
  [`build_teolis_profiles.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/tools/build_teolis_profiles.py).
- Data release assembly:
  [`prepare_plume_fit_data.py`](/data/vtenishe/cassini_validation_data/scripts/prepare_plume_fit_data.py).
- Joint fitter:
  [`joint_fitting_study.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/joint_fitting_study.py).
- Source physics:
  [`plume_model.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/plume_model.py).

The processed release carries source hashes and metadata.  The paths above are
not hard-coded by the fit: the validation-data root is a command-line input.

## 1. Observations and calibration lineage

### 1.1 INMS E3 and E5 local density

The two fit tables contain 79 outbound points each, sampled every 5 s from
+10 to +400 s relative to closest approach.  Spacecraft positions are linearly
interpolated from official PDS INMS Level-1A trajectory columns in the
IAU_ENCELADUS frame.

The released density shape is a digitization of the red, sticking-corrected
trace in Figure 5 of Tenishev et al. (2010), which reproduces the treatment of
Teolis et al. (2010).  Page 7 of the source PDF was rendered at 200 dpi, red
curve pixels were located in the two panels, and the vertical pixel coordinate
was converted through the plotted logarithmic density axis.  Each curve was
then normalized to the reported peak density:

| flyby | rows | interval | published peak (n_{\rm H_2O}) | independent digitization σ/n | common scale uncertainty |
|---|---:|---:|---:|---:|---:|
| E3 | 79 | +10 to +400 s | (7.6\times10^{13} {m m^{-3}}) | 5.157% | 23% |
| E5 | 79 | +10 to +400 s | (3.4\times10^{14} {m m^{-3}}) | 5.174% | 23% |

The independent uncertainty is the density change corresponding to two pixels
of line thickness on the logarithmic axis,

\[
{\sigma_{n,i}\over n_i}=\left|2a\ln 10\right|,
\]

where (a=d\log_{10}n/dy_{\rm pixel}).  It is used in the pointwise fit.  The
23% revised-INMS sensitivity uncertainty from Teolis et al. (2015) is a common
multiplicative scale and is deliberately not added independently to every
sample.

These are calibrated, publication-derived profiles, not a new inversion of
raw Level-1A water counts.  That qualification matters: digitization cannot
recover covariance introduced by the original sticking correction, and the
common normalization does not supply 79 independent absolute measurements.

### 1.2 UVIS 2007 stellar and 2010 solar occultations

The UVIS input is broadband attenuation along occultation rays:

| event | date | samples | reference counts (C_0) | published peak H₂O column used for scale |
|---|---|---:|---:|---:|
| ζ Orionis stellar occultation | 2007-10-24 | 97 | 37,597 | ((1.4\pm0.14)\times10^{20} {m m^{-2}}) |
| E10 solar occultation | 2010-05-18 | 1,133 | 16,816 | ((9.0\pm2.3)\times10^{19} {m m^{-2}}) |

For observed counts (C_i), the transmission and raw optical-depth proxy are

\[
T_i={C_i\over C_0},\qquad \tau_i=-\ln T_i .
\]

Ignoring uncertainty in (C_0), Poisson propagation gives

\[
\sigma_{\tau,i}=\left|{\partial\tau\over\partial C}\right|\sqrt{C_i}
= {1\over\sqrt{C_i}}.
\]

The release stores both raw optical depth and a nonnegative normalized shape.
The fitter uses raw τ, including negative noise excursions, and fits one signed
reference-level offset.  Clipping the fitted data would bias a weak absorption
signal upward.

For a monochromatic, optically thin absorption measurement,

\[
I_\lambda=I_{0,\lambda}e^{-\tau_\lambda},\qquad
\tau_\lambda=\int n_{\rm H_2O}(\mathbf x)\sigma_{\rm H_2O}(\lambda)\,d\ell .
\]

A spectrally integrated detector measures an instrument-response-weighted
integral of (I_{0,\lambda}e^{-\tau_\lambda}), not a single τ/σ.  The archived
broadband series does not contain the wavelength-resolved response needed for a
new sample-by-sample absolute column inversion.  The implemented, auditable
procedure therefore fits the attenuation **shape** and applies the separately
published event column once to set its absolute scale.  No fluorescence
g-factor appears: a g-factor converts emission brightness to column abundance,
whereas these data are absorption occultations.

The tangent-point Cartesian coordinate follows the archived latitude,
longitude, altitude transformation with (R_E=252.1) km.  The 2007 stellar
ray is treated as an effectively infinite chord.  The 2010 solar event uses
SPICE-derived Cassini-to-Sun direction and finite observer range when the
tangent point lies ahead of the observer; outward rays are integrated from the
observer.  This avoids extending a ray through material that cannot lie
between Cassini and the Sun.

## 2. Neutral-plume physics

### 2.1 Production-normalized point source

Source (j) is at (mathbf x_j), has outward surface normal
(hat{\mathbf n}_j), source axis (hat{\mathbf a}_j), temperature (T_j),
bulk speed (u_j), and total production (Q_j).  With
(m_{\rm H_2O}=18.01528m_u=2.9915\times10^{-26}) kg, define

\[
v_{{\rm th},j}=\sqrt{{2k_BT_j\over m_{\rm H_2O}}},\qquad
s_j={u_j\over v_{{\rm th},j}}.
\]

At field point (mathbf x), let

\[
\mathbf R=\mathbf x-\mathbf x_j,\quad R=|\mathbf R|,\quad
\hat{\mathbf r}=\mathbf R/R,
\]

\[
\mu=\hat{\mathbf r}\cdot\hat{\mathbf n}_j,\qquad
c=\hat{\mathbf r}\cdot\hat{\mathbf a}_j .
\]

Velocity integration of the shifted Maxwellian introduces

\[
K_p(c,s)=\int_0^\infty t^p
\exp[-t^2+2sct-s^2],dt .
\]

The angular-flux normalization is

\[
\mathcal N_j=\int_{\mu\ge0}L(\mu)K_3(c,s_j),d\Omega .
\]

The density per source is then

\[
n_j(\mathbf x)=
\begin{cases}
\displaystyle {Q_jL(\mu)K_2(c,s_j)\over
\mathcal N_jv_{{\rm th},j}R^2},&\mu\ge0,\\[5pt]
0,&\mu<0.
\end{cases}
\]

The normalization makes (Q_j) the total outward molecular flux.  The default
Tenishev (2010) law is (L(\mu)=1).  The selectable surface-flux Maxwellian is
(L(\mu)=\max(\mu,0)).  All reported results in this document use the default.

### 2.2 Tiger-stripe line sources and source cases

A line source is evaluated as a length-weighted quadrature of normalized point
sources along its polyline.  It has one fitted production coefficient, so its
production is not multiplied by the number of quadrature nodes.

| source case | physical components | comparison kinematics |
|---|---|---|
| `tenishev2014_hybrid_8` | 8 discrete jets and 7 modeled tiger-stripe polylines | source-specific jet values; stripe (T=180) K and (u=350) or 950 m s⁻¹ |
| `porco2014_98` | 98 image-derived jet locations and axes | (T=180) K, (u=1500) m s⁻¹ |
| `portyankina2022_solar_41` | 41 sources selected for the 2010 solar-occultation model, using local-normal gas axes | (T=180) K, (u=1500) m s⁻¹ |

Porco et al. inferred particulate-jet geometry.  Treating those axes as gas
sources is a test hypothesis, not a statement that the gas and particles share
the same angular distribution.  The Portyankina 41-source selection is tied to
the 2010 UVIS geometry and is not an independent global gas catalog.

### 2.3 Line-of-sight column

For occultation ray origin (mathbf x_0), unit direction (hat{\mathbf d}),
and physical integration interval ([\ell_a,\ell_b]), the model column is

\[
N_{ij}=\int_{\ell_a}^{\ell_b}
n_j(\mathbf x_0+\ell\hat{\mathbf d}),d\ell .
\]

The code evaluates this linear operator with Gauss quadrature.  Convergence is
checked by repeating the fit at higher quadrature order.  An (r^{-2}) diffuse
density background can be integrated because its infinite-chord column is
finite.  A spatially constant local-density background is not integrated on an
infinite chord because it would diverge; UVIS instead has a signed additive
optical-depth/reference offset.

### 2.4 Physical approximations

The model is steady, collisionless, gravity-free, and loss-free over each
observation.  It neglects vent/nozzle collisions, plume-plume collisions,
surface re-impact, rotation during flight, and time evolution during a profile.
Photochemical loss is safely slow relative to the seconds-to-minutes flight
times represented here, but collisions and near-vent acceleration can alter
angular spreading.  Thus the result is an effective source rate for the stated
kernel, not a substitute for the later collisional AMPS calculation.

## 3. Fitting mechanics

### 3.1 INMS local-density fit

For sample (i), let (G_{ij}=n_j(\mathbf x_i)/Q_j).  The model is

\[
\hat n_i=\sum_jG_{ij}Q_j+c_0+c_2(r_i/R_E)^{-2},
\qquad Q_j,c_0,c_2\ge0 .
\]

The two background terms absorb broad unresolved plume/exosphere signal and
instrumental baseline.  They are nuisance parameters and are excluded from the
reported vent production.  Design columns are internally expressed per
(10^{28}) molecules s⁻¹ and scaled by their weighted norms.  The solver
minimizes

\[
\chi^2_{\rm INMS}=\sum_i\left({n_i-\hat n_i\over\sigma_{n,i}}\right)^2
\]

subject to nonnegative physical coefficients.

### 3.2 UVIS shape and absolute scale

Let (A_{ij}=N_{ij}/Q_j) be the source column kernel.  The shape fit minimizes

\[
\chi^2_{\rm UVIS}=\sum_i
\left[{\tau_i-b-\sum_jA_{ij}a_j\over\sigma_{\tau,i}}\right]^2,
\qquad a_j\ge0,quad b\in\mathbb R .
\]

The signed (b) accounts for uncertainty in the reference level; it is not a
negative molecular column.  If (S_i=\sum_jA_{ij}a_j), the physical scaling is

\[
f={N_{\rm pub,peak}\over\max_iS_i},\qquad Q_j=f,a_j,qquad
N_{{\rm model},i}=fS_i .
\]

Thus the profile shape and published peak column enter once each.  Optical
depth is never mislabeled as column density.

For the all-observation shape diagnostic only, each UVIS event is mapped to a
pseudo-column,

\[
N_i^*=\max(\tau_i-b,0)
{N_{\rm pub,peak}\over\max_k\max(\tau_k-b,0)},
\]

with its uncertainty scaled by the same factor.  This construction permits a
balanced joint linear fit; it is not an independently calibrated UVIS column at
every sample.

### 3.3 Shape strategies and dataset balancing

Four declared strategies are computed:

1. `independent_each_observation`: every event determines its own relative
   source shape and total rate;
2. `shared_shape_all_observations`: all four events determine one shape;
3. `shared_shape_inms_only`: E3 and E5 determine the shape;
4. `shared_shape_uvis_only`: the two UVIS occultations determine the shape.

After a shared shape is found, it is normalized to fractions
(w_j=Q_j/\sum_kQ_k).  Each target then fits one total source amplitude plus
target-specific nuisance backgrounds:

\[
\hat y_i=Q_{\rm target}\sum_jw_jG_{ij}+\hbox{backgrounds}.
\]

In a combined shape objective, the point uncertainty for dataset (d) is
multiplied by (sqrt{N_d}).  The minimized objective is consequently the sum
of per-dataset mean squared normalized residuals, preventing the 1,133-point
solar event from dominating merely by sample count.

## 4. Results

### 4.1 Total H₂O source rates

The table gives molecules s⁻¹.  `Ind.`, `All`, `INMS`, and `UVIS` denote the
four shape strategies above.  These are conditional estimates, not four
independent measurements of the same quantity.

| source case | strategy | E3 | E5 | UVIS 2007 | UVIS 2010 |
|---|---|---:|---:|---:|---:|
| hybrid 8 + stripes | Ind. | 2.061e28 | 4.520e28 | 7.930e27 | 4.801e27 |
| hybrid 8 + stripes | All | 1.272e28 | 2.607e28 | 7.415e27 | 5.635e27 |
| hybrid 8 + stripes | INMS | 3.513e28 | 4.631e28 | 1.630e27 | 1.069e27 |
| hybrid 8 + stripes | UVIS | 8.046e27 | 1.745e28 | 6.449e27 | 4.231e27 |
| Porco 98 | Ind. | 1.044e31 | 2.305e30 | 1.685e28 | 1.314e28 |
| Porco 98 | All | 1.931e28 | 2.294e28 | 1.643e28 | 1.169e28 |
| Porco 98 | INMS | 2.012e30 | 1.621e30 | 3.895e27 | 1.755e27 |
| Porco 98 | UVIS | 1.934e28 | 4.960e28 | 1.684e28 | 1.387e28 |
| Portyankina 41 | Ind. | 2.263e28 | 5.331e28 | 8.338e27 | 1.405e28 |
| Portyankina 41 | All | 1.692e28 | 2.369e28 | 8.417e27 | 1.034e28 |
| Portyankina 41 | INMS | 1.048e29 | 1.712e29 | 7.495e26 | 1.255e27 |
| Portyankina 41 | UVIS | 1.698e28 | 2.603e28 | 8.495e27 | 1.221e28 |

For mass loss in kg s⁻¹, multiply the molecular rate by
(2.9915\times10^{-26}) kg molecule⁻¹; (10^{28}) H₂O molecules s⁻¹ is
299.15 kg s⁻¹.

### 4.2 Residual quality for independent fits

| source case | E3 χ²/N | E5 χ²/N | UVIS 2007 χ²/N | UVIS 2010 χ²/N |
|---|---:|---:|---:|---:|
| hybrid 8 + stripes | 5.115 | 1.051 | 1.533 | 3.278 |
| Porco 98 | 6.045 | 0.954 | 0.362 | 0.664 |
| Portyankina 41 | 8.900 | 1.198 | 3.029 | 1.687 |

Values below unity for the most flexible Porco fits do not establish physical
superiority.  The 98-source model has at least as many source degrees of freedom
as the 97-point stellar profile, and many source columns are linearly dependent.
Conversely, χ²/N well above unity shows that point uncertainty alone does not
describe model error.

For the hybrid geometry, removing all stripe components increases χ² from
404.1 to 536.7 (E3), 83.0 to 117.7 (E5), 142.3 to 375.6 (2007), and 3857.0 to
7133.4 (2010) in the corresponding ablation calculation.  Removing discrete
jets has a much smaller effect.  These observations support a distributed
tiger-stripe contribution, but they do not identify unique individual stripe
or jet strengths.

### 4.3 Figures and machine-readable tables

- Hybrid comparison:
  [`tenishev2014_hybrid_8`](/data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014/tenishev2014_hybrid_8).
- Porco comparison:
  [`porco2014_98`](/data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014/porco2014_98).
- Portyankina comparison:
  [`portyankina2022_solar_41`](/data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014/portyankina2022_solar_41).
- Cross-case summary:
  [`source_case_comparison.md`](/data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014/source_case_comparison.md).

Each case directory contains observed/model profiles, source-rate tables,
cross-predictions, rank and condition diagnostics, quadrature checks, and
`study_results.json` with full precision.

## 5. Uncertainty and interpretation limits

- **Digitization:** about 5.2% per INMS point; included in the local-density
  objective.
- **INMS common scale:** 23%; it shifts an event coherently and should be
  propagated to its fitted total, not treated as 79 independent errors.
- **UVIS counting noise:** the pointwise Poisson approximation is included.
  Reference-count and published-column errors remain common/systematic terms.
- **UVIS absolute column:** 10% for 2007 and 25.6% for 2010 from the adopted
  published constraints; source rates scale linearly with them.
- **Occultation spectral response:** the archived broadband trace supports
  shape fitting, not a new wavelength-resolved absorption inversion.
- **Source geometry and kinematics:** dominant for many fits.  The spread among
  source cases is not a statistical confidence interval.
- **Identifiability:** nonnegativity reduces numerical ambiguity but does not
  create information absent from one-dimensional tracks or chords.
- **Steady-state assumption:** each observation is treated independently in
  time.  A shared shape is an explicit hypothesis tested through residuals,
  not an imposed physical fact.
- **Model inadequacy:** collisionless kernels and two smooth backgrounds cannot
  reproduce every fine feature; ordinary covariance based only on point noise
  would therefore understate uncertainty.

The INMS and UVIS fitted totals should not be averaged into one “best” number
without a hierarchical treatment of these common scales and model discrepancy.

## 6. Conclusions

1. The current release supports reproducible H₂O fitting of E3/E5 local
   density and two UVIS occultation profiles, with distinct and documented
   calibration paths.
2. UVIS optical depth is used as a profile-shape observable and is normalized
   by a published H₂O column.  It is not directly relabeled as column density,
   and no emission g-factor is involved.
3. A distributed tiger-stripe component improves all four hybrid fits.  The
   available data do not uniquely determine individual vent amplitudes.
4. The extremely large independent Porco-INMS totals are clear evidence of
   ill-conditioning; residual minimization alone is insufficient for physical
   model selection.
5. No single time-invariant relative source map reproduces all INMS and UVIS
   observations at their point-error level.  This can reflect real source
   redistribution, different event geometries, gas/particle catalog mismatch,
   omitted collisional physics, calibration systematics, or a combination.
6. These results are appropriate as Stage-0 constraints and sensitivity tests.
   Publication-level global rates should be reported conditional on source map,
   kinematics, and observation/calibration strategy.

## 7. Commands

From `/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone`:

```bash
export ENCELADUS_VALIDATION_DATA_ROOT=/data/vtenishe/cassini_validation_data
export RUN_TAG=source_shape_tidal_20261014
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/tenishev2014_hybrid_8" --quadrature-order 192
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco2014_98" --quadrature-order 192
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case portyankina2022_solar_41 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/portyankina2022_solar_41" --quadrature-order 192
```

Use a new output directory for a new run.  Existing result directories are not
silently overwritten.  Data preparation, when intentionally rebuilding the
release, is described in `README.md` and must be run with the same external
validation-data root.

## References

1. Waite, J. H. et al. (2004), “The Cassini Ion and Neutral Mass Spectrometer
   (INMS) Investigation,” *Space Science Reviews* 114, 113–231,
   <https://doi.org/10.1007/s11214-004-1408-2>.
2. Teolis, B. D. et al. (2010), “Detection and measurement of ice grains and gas
   distribution in the Enceladus plume by Cassini's Ion Neutral Mass
   Spectrometer,” *JGR Space Physics* 115, A09222,
   <https://doi.org/10.1029/2009JA015192>.
3. Teolis, B. D. et al. (2015), “A Revised Sensitivity Model for Cassini INMS:
   Results at Titan,” *Space Science Reviews* 190, 47–84,
   <https://doi.org/10.1007/s11214-014-0133-8>.
4. Hansen, C. J. et al. (2011), “The composition and structure of the Enceladus
   plume,” *Geophysical Research Letters* 38, L11202,
   <https://doi.org/10.1029/2011GL047415>.
5. Hansen, C. J. et al. (2020), “The composition and structure of Enceladus'
   plume from the complete set of Cassini UVIS occultation observations,”
   *Icarus* 344, 113461,
   <https://doi.org/10.1016/j.icarus.2019.113461>.
6. Tenishev, V. et al. (2010), “An approach to numerical simulation of the gas
   distribution in the atmosphere of Enceladus,” *JGR Space Physics* 115,
   A09302, <https://doi.org/10.1029/2009JA015223>.
7. Tenishev, V. et al. (2014), “Effect of the Tiger Stripes on the water vapor
   distribution in Enceladus' exosphere,” *JGR Planets* 119, 2658–2667,
   <https://doi.org/10.1002/2014JE004700>.
8. Porco, C. C., D. DiNino, and F. Nimmo (2014), “How the Geysers, Tidal
   Stresses, and Thermal Emission across the South Polar Terrain of Enceladus
   Are Related,” *Astronomical Journal* 148, 45,
   <https://doi.org/10.1088/0004-6256/148/3/45>.
9. Portyankina, G. et al. (2022), “Modeling the complete set of Cassini's UVIS
   occultation observations of Enceladus' plume,” *Icarus* 374, 114918,
   <https://doi.org/10.1016/j.icarus.2022.114918>.
10. NASA PDS Atmospheres Node, Cassini INMS archive and calibration guidance,
    <https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/inst-inms.html>.
