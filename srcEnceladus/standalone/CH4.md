# Cassini INMS CH₄ calibration and Enceladus plume-source study

## Scope and scientifically defensible result

This study reconstructs a time-resolved methane proxy from Cassini INMS Level-1A counts and fits it with the three source geometries implemented in the standalone Enceladus plume model. It is reproducible from the archived inputs in this data tree.

The principal result is deliberately narrower than “a unique CH₄ production rate.” E14, E17, and E18 contain statistically strong, plume-correlated mass-15 signals that are usable for studying CH₄ spatial structure. Their background-subtracted, time-integrated CH₄/CO₂ density ratios over ±50 s are 0.233, 0.250, and 0.291. Absolute CH₄ source rates, however, change by more than two orders of magnitude among the three assumed source maps and gas-kinematic prescriptions. The data therefore constrain a CH₄ profile and conditional source rates; they do **not** select a unique vent distribution or geometry-independent global CH₄ production rate.

E14, E17, and E18 are the primary set because the accumulated low-velocity spectra used by Waite et al. (2017) are the least altered compositional measurements and explicitly constrain CH₄ through masses 15 and 13. E7 is secondary. E21 has only 32 usable primary-filament mass-15 samples inside ±100 s and reaches a 3.42° ram angle, so it is also secondary. E3 and E5 are retained only as impact-chemistry diagnostics: ice-grain impacts and the higher encounter speed make a simple gas-composition interpretation unsafe. E13 is a non-detection/control (12 mass-15 counts in 87 samples).

No UVIS CH₄ fit is made. The UVIS products in this archive constrain H₂O absorption; there is no CH₄ occultation product with a validated CH₄ cross section and sufficient signal in the selected data. Consequently, fluorescence *g*-factors and line-of-sight brightness equations are not part of this analysis. Applying an H₂O UVIS calibration to CH₄ would be physically invalid.

## Reproducible products

- Official calibration: [`2004001_CAL_01.CSV`](/data/vtenishe/cassini_validation_data/data/raw/inms/calibration/2004001_CAL_01.CSV), its PDS label, and `CALINFO.TXT`.
- Calibrated profiles and release manifest: [`ch4_fit_v1`](/data/vtenishe/cassini_validation_data/data/processed/ch4_fit_v1).
- Three-geometry fits: [`ch4_20261016`](/data/vtenishe/cassini_validation_data/studies/ch4_20261016).
- Profiles including the one cross-channel-rejected E14 event, retained only for sensitivity testing: [`ch4_20261016_including_flagged`](/data/vtenishe/cassini_validation_data/studies/ch4_20261016_including_flagged).
- Processing: [`prepare_ch4_fit_data.py`](/data/vtenishe/cassini_validation_data/scripts/prepare_ch4_fit_data.py).
- Fitting: [`run_ch4_fitting_study.py`](/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/run_ch4_fitting_study.py).
- Cross-geometry summary: [`summarize_ch4_study.py`](/data/vtenishe/cassini_validation_data/scripts/summarize_ch4_study.py).
- Verification checks: [`test_ch4_release.py`](/data/vtenishe/cassini_validation_data/tests/test_ch4_release.py).

The release manifest records SHA-256 hashes of the calibration table and every processed profile. The source CSVs retain the original Level-1A file name for every sample.

## 1. Instrument observable and data selection

INMS does not archive a calibrated CH₄ density channel. The PDS Level-1A products contain counts and count rates. The official INMS user guide likewise states that the user must combine the count rate, species sensitivity, fragmentation pattern, and inlet response to obtain ambient density.

The selected records satisfy all of the following:

1. target `ENCELADUS`;
2. closed-source neutral mode (`source = csn`);
3. primary closed-source filament in its high-electron-energy state (`csp_fil_3_status = High-EV`), corresponding to the 70 eV calibration;
4. C1 high-sensitivity detector;
5. archived `data_reliability = 0`;
6. mass-to-charge 15 u/e for fitting, and 13 u/e for an independent fragmentation check;
7. |t − t_CA| ≤ 100 s for the fitted profile.

The off-plume detector background for each flyby and mass is the mean C1 counts per integration over 500 ≤ |t − t_CA| ≤ 1000 s. This is measured separately for every flyby rather than imposed from another encounter.

| Flyby | mass-15 samples | mass-15 counts | off-plume counts/integration | maximum ram angle | role |
|---|---:|---:|---:|---:|---|
| E3 | 86 | 4,727 | 1.201 | 0.167° | impact diagnostic |
| E5 | 87 | 14,460 | 2.264 | 0.192° | impact diagnostic |
| E7 | 85 | 1,226 | 0.903 | 0.314° | secondary |
| E13 | 87 | 12 | 0.155 | 0.306° | non-detection/control |
| E14 | 173 | 3,401 | 1.008 | 0.304° | primary; one point flagged |
| E17 | 173 | 3,332 | 0.931 | 0.303° | primary |
| E18 | 173 | 4,007 | 1.373 | 0.305° | primary |
| E21 | 32 | 161 | 1.219 | 3.416° | sparse secondary |

The median integration time derived from `c1counts/c1rate` is 0.031103 s. The archived `c1error` equals the Poisson count uncertainty propagated to count rate.

## 2. Why mass 15 is used

Electron-impact ionization fragments CH₄. The official 70 eV flight-model, primary-filament, closed-source calibration row gives

\[
s_{16}=7.47\times10^{-4}\ {\rm cm^3\,s^{-1}},\qquad
p_{15}=0.723,\qquad p_{13}=0.0205,
\]

where \(s_{16}\) is the major-peak sensitivity and \(p_m\) is the mass-\(m\) intensity relative to the mass-16 CH₄ peak. Thus the pure-CH₄ count ratio is

\[
\left({X_{13}\over X_{15}}\right)_{\rm CH_4}={p_{13}\over p_{15}}=0.02835.
\]

Mass 16 itself cannot isolate CH₄ in an Enceladus plume because H₂O fragmentation dominates the neighboring 16–18 u/e spectrum. Waite et al. (2017, supplementary material) report that CH₄ supplies about 90% of the accumulated mass-15 and mass-13 signals in the combined low-velocity spectra. The nominal time-resolved conversion therefore assigns a fraction \(f_{15}^{\rm CH_4}=0.90\) of mass 15 to methane. Because an accumulated spectral fraction need not be constant through every spatial sample, 0.80–1.00 is retained as a sensitivity range, not called a confidence interval.

Mass 13 is too sparse for a sample-by-sample density fit. After subtracting the independently measured detector backgrounds, the accumulated mass-13/mass-15 count ratios are 0.0441, 0.0472, and 0.0491 for E14, E17, and E18. Their excess over 0.02835 demonstrates that mass 13 also contains background and/or other fragments; it cannot simply be combined with mass 15 as a second pure-CH₄ channel. It is used as a diagnostic only.

### Cross-channel event rejection

One E14 sample at +52.948 s contains 536 mass-15 counts. The nearest mass-13 measurement, 0.102 s earlier, contains zero counts, whereas the CH₄ flight fragmentation pattern predicts about 13.6 mass-13 counts after background and the nominal 90% assignment. This is incompatible with CH₄ at far beyond ordinary Poisson scatter and is consistent with a grain/other-species event of the type documented by Teolis et al. (2010) and Perry et al. (2015).

The algorithm flags a point only if (i) the nearest mass-13 scan is within 0.2 s, (ii) the expected CH₄ mass-13 signal is at least 10 counts, and (iii) the observed mass-13 signal is at most one count. Only that E14 point is flagged. It remains in the released profile with `fit_usable = 0`; it is not deleted. Including it changes the fitted E14 source rate by −0.34% to +1.66%, depending on source case and fitting strategy ([table](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/e14_flagged_event_sensitivity.csv)), so it does not drive the main source-rate result.

## 3. Count-rate to density calibration

Let \(X_m\) be the archived C1 count rate at fragment mass \(m\), \(T_0=293\) K the calibration-instrument temperature, \(T_a=150\) K the nominal ambient neutral temperature used in the inlet response, \(u\) the archived spacecraft/target relative speed, \(\theta\) the archived ram angle, and \(m_{\rm CH_4}=16.04246\,m_u=2.664\times10^{-26}\) kg. Define

\[
S=u\sqrt{m_{\rm CH_4}\over 2k_B T_a}\cos\theta,
\]

\[
F(S)=e^{-S^2}+\sqrt{\pi}S\,[1+\operatorname{erf}(S)].
\]

For the archived preflight sensitivity and ideal on-ram response, the density associated with a fragment rate is proportional to

\[
n_{\rm old}={X_m\over s_{16}p_m}\sqrt{T_0\over T_a}\,{1\over F(S)}.
\]

Teolis et al. (2015) revised the closed-source leakage and detector-gain model. For the special near-on-ram case, they prescribe multiplying densities obtained with the original calibration by \(2.2\pm23\%\). The nominal mass-15 methane-equivalent density used here is therefore

\[
\boxed{
n_{15}^{\rm CH_4}=f_{15}^{\rm CH_4}\,2.2\,
{X_{15}\over (7.47\times10^{-10}\ {\rm m^3\,s^{-1}})(0.723)}
\sqrt{293\over150}\,{1\over F(S)} }
\]

with \(f_{15}^{\rm CH_4}=0.90\). The primary flybys have \(|\theta|<0.305^\circ\), which makes use of the published on-ram shortcut justified. E21's 3.42° maximum angle is one reason it is not in the primary set.

At 7.5 km s⁻¹, \(S=19.0\) and the nominal conversion is \(7.60\times10^7\ {\rm m^{-3}}\) per count s⁻¹. At 17.7 km s⁻¹ it is \(3.22\times10^7\ {\rm m^{-3}}\) per count s⁻¹. In this high-\(S\) limit, \(F\simeq2\sqrt\pi S\), so the factor \(\sqrt{T_0/T_a}/F\) becomes independent of \(T_a\). Direct evaluation from 85 to 210 K changes it by less than \(4\times10^{-16}\) relative at these speeds; ambient-temperature uncertainty is therefore negligible for this conversion. This cancellation does not mean temperature is irrelevant to plume propagation, where it controls angular spreading.

For a sample containing \(C_i\) counts in integration time \(\Delta t_i\), the gross fitted value, background-subtracted diagnostic value, and counting uncertainty are

\[
y_i=K_i{C_i\over\Delta t_i},\qquad
y_{i,{\rm net}}=K_i{C_i-\bar B\over\Delta t_i},
\]

\[
\sigma_i=K_i\max\left({\sqrt{C_i}\over\Delta t_i},{1\over\Delta t_i}\right),
\]

where \(K_i\) is the full density-per-count-rate factor above and \(\bar B\) is the mean off-plume counts per integration. The one-count floor prevents a zero-count sample from receiving zero uncertainty.

The fitter uses the nonnegative **gross** density \(y_i\), with free instrumental/physical nuisance-background columns. It does not clip negative values from \(y_{i,{\rm net}}\), a practice that would bias a low-count signal upward. The net product is used only for figures and integral diagnostics.

## 4. Neutral plume physics

The fit calls the existing standalone collisionless, steady, gravity-free model. For source \(j\), the thermal speed and bulk-to-thermal speed ratio are

\[
v_{{\rm th},j}=\sqrt{2k_BT_j\over m_{\rm CH_4}},\qquad s_j={u_j\over v_{{\rm th},j}}.
\]

For a field point \(\mathbf x\), define the source displacement \(\mathbf R=\mathbf x-\mathbf x_j\), \(R=|\mathbf R|\), ray \(\hat{\mathbf r}=\mathbf R/R\), surface-normal cosine \(\mu=\hat{\mathbf r}\cdot\hat{\mathbf n}_j\), and source-axis cosine \(c=\hat{\mathbf r}\cdot\hat{\mathbf a}_j\). The velocity moments are

\[
K_p(c,s)=\int_0^\infty t^p\exp[-t^2+2sct-s^2]\,dt.
\]

The total-flux normalization over the outward hemisphere is

\[
N_j=\int_{\mu\ge0}L(\mu)K_3(c,s_j)\,d\Omega,
\]

and the density from production rate \(Q_j\) is

\[
n_j(\mathbf x)=
\begin{cases}
\displaystyle {Q_j L(\mu)K_2(c,s_j)\over N_jv_{{\rm th},j}R^2},&\mu\ge0,\\[6pt]
0,&\mu<0.
\end{cases}
\]

The default Tenishev (2010) point-source law uses \(L(\mu)=1\). The alternative surface-flux Maxwellian, selectable in the standalone configuration but not used in this nominal study, uses \(L(\mu)=\max(\mu,0)\). Tiger-stripe line sources are length-weighted sums of these normalized point kernels along each published polyline.

CH₄ uses the same source temperatures and bulk speeds as the corresponding H₂O source configuration, but its own molecular mass in \(v_{\rm th}\). This is a controlled modeling assumption, not an INMS measurement of species-specific CH₄ kinematics. It permits like-for-like source-map comparison while retaining the physically required mass dependence of thermal spreading.

The three cases are:

| Case | Geometry | gas kinematics used here |
|---|---|---|
| `tenishev2014_hybrid_8` | 8 discrete jets plus the modeled tiger-stripe components | source-specific 170–210 K and 580–1338 m s⁻¹ jets; stripes at 180 K and 350 or 950 m s⁻¹ |
| `porco2014_98` | all 98 Porco et al. image-derived jet axes | 180 K, 1500 m s⁻¹ controlled-comparison setting |
| `portyankina2022_solar_41` | 41 jets selected by the 2010 solar-occultation solution, local-normal axes | 180 K, 1500 m s⁻¹ |

Gravity, neutral-neutral collisions, photochemistry, surface re-impact after launch, and time variability during an individual flyby are omitted. The first two can matter close to dense vents; therefore the inferred rates are effective rates conditional on a collisionless kernel, not a replacement for a collisional AMPS calculation. Chemical loss is negligible over the seconds-to-minutes sampled by a flyby, but this analysis does not extrapolate the kernels to the long-lived neutral torus.

## 5. Fitting mechanics

For sample \(i\), source \(j\), and model kernel \(G_{ij}=n_j(\mathbf x_i)/Q_j\), the fitted model is

\[
\hat y_i=\sum_j G_{ij}Q_j+c_0+c_2\left({r_i\over R_E}\right)^{-2},
\qquad Q_j,c_0,c_2\ge0,
\]

where \(R_E=252.1\) km. Source design-matrix columns are stored per \(10^{28}\) molecules s⁻¹. The constant and \(r^{-2}\) terms are nuisance backgrounds; they absorb detector background and broad unresolved plume/exosphere signal and are not counted as vent production.

The solver minimizes the weighted nonnegative least-squares objective

\[
\chi^2=\sum_i\left({\hat y_i-y_i\over\sigma_i}\right)^2,
\qquad \mathbf q\ge0,
\]

after scaling every design column by its weighted norm. Three strategies are reported:

1. `independent_each_observation`: every flyby fits all source amplitudes and its rate independently;
2. `shared_shape_primary_low_impact`: one relative source shape is fitted jointly to E14, E17, and E18, then each flyby refits one total source amplitude plus its backgrounds;
3. `shared_shape_all_low_speed`: the relative shape is fitted jointly to E7, E14, E17, E18, and E21, then each flyby refits its total amplitude.

For a combined shape fit, the uncertainty of each dataset is multiplied by \(\sqrt{N_d}\). The combined minimized objective is therefore the sum of the datasets' mean squared normalized residuals, preventing a dataset from dominating only because it contains more samples. If \(\mathbf q_*\) is the combined solution, its source portion is normalized,

\[
w_j={q_{*,j}\over\sum_k q_{*,k}},
\]

and a target flyby fits

\[
\hat y_i=Q_{\rm target}\sum_jw_jG_{ij}+\hbox{target-specific backgrounds}.
\]

Individual jet amplitudes are not interpreted. Many source columns are nearly dependent along a one-dimensional spacecraft trajectory, and the 98-jet solution can fit small residuals with very large, poorly observed rates. Nonnegativity stabilizes the numerical solution but does not create missing spatial information.

## 6. Results

### Primary conditional source rates

Rates below are molecules s⁻¹. The three rightmost values are \(\chi^2/N\). These formal residual statistics use counting noise only and are greater than unity in most cases, demonstrating unresolved structure and model/systematic error.

| source case | strategy | E14 Q | E17 Q | E18 Q | E14 χ²/N | E17 χ²/N | E18 χ²/N |
|---|---|---:|---:|---:|---:|---:|---:|
| hybrid 8 + stripes | independent | 3.020e25 | 9.582e24 | 1.492e25 | 5.90 | 5.10 | 7.10 |
| hybrid 8 + stripes | shared E14/E17/E18 shape | 1.053e25 | 1.149e25 | 1.557e25 | 6.19 | 5.45 | 7.80 |
| hybrid 8 + stripes | shared all-low-speed shape | 5.034e24 | 5.287e24 | 8.101e24 | 6.41 | 5.52 | 8.04 |
| Porco 98 | independent | 9.564e25 | 1.488e29 | 5.029e26 | 2.94 | 1.41 | 3.30 |
| Porco 98 | shared E14/E17/E18 shape | 1.801e27 | 1.851e27 | 1.652e27 | 4.12 | 3.88 | 3.72 |
| Porco 98 | shared all-low-speed shape | 8.970e24 | 1.023e25 | 1.105e25 | 6.02 | 5.28 | 7.58 |
| Portyankina 41 | independent | 3.132e26 | 9.065e27 | 2.865e28 | 6.88 | 4.93 | 8.57 |
| Portyankina 41 | shared E14/E17/E18 shape | 2.348e26 | 3.248e26 | 3.721e26 | 7.51 | 6.75 | 9.38 |
| Portyankina 41 | shared all-low-speed shape | 7.123e24 | 9.266e24 | 1.238e25 | 7.92 | 6.84 | 10.67 |

The independent Porco and Portyankina fits illustrate the inverse problem's instability: a low residual does not imply a physically determined total when the fitted shape places most sources where the trajectory has little sensitivity. Even the shared-primary totals span approximately \(10^{25}\) to \(2\times10^{27}\) s⁻¹ across source cases. Therefore none of these cross-geometry differences should be interpreted as methane variability.

Within one fixed shared-primary shape, E14-to-E18 variation is modest for Porco (about −9%) and larger for the hybrid and Portyankina cases (about +48% and +59%). Because the sign and magnitude depend on the source assumption and because only three primary epochs exist, the analysis does not establish a secular or orbital CH₄ trend.

### Geometry-independent CH₄/CO₂ profile check

For each primary flyby, positive background-subtracted CH₄ density is integrated from −50 to +50 s. The calibrated CO₂ profile is integrated over the same interval after subtracting the median CO₂ density at 40–50 s on both sides. This comparison does not assume a vent map:

| flyby | ∫n(CH₄)dt (m⁻³ s) | ∫n(CO₂)dt (m⁻³ s) | ratio |
|---|---:|---:|---:|
| E14 | 6.343e12 | 2.722e13 | 0.233 |
| E17 | 7.908e12 | 3.165e13 | 0.250 |
| E18 | 9.262e12 | 3.187e13 | 0.291 |

These ratios are compatible with the broad CH₄/CO₂ abundance envelope implied by the low-velocity spectral solutions in Waite et al. (2017). They are not identical to source-rate ratios: species-dependent thermal spreading and possibly species-dependent vent partitioning mean that a local trajectory integral need not equal the global production ratio.

### Figures

- [Cross-geometry primary source rates](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/ch4_primary_source_rates.png)
- [CH₄/CO₂ profile-integral ratios](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/ch4_co2_profile_integral_ratios.png)
- [Hybrid-profile comparisons](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/tenishev2014_hybrid_8/ch4_fit_profiles.png)
- [Porco-98 profile comparisons](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/porco2014_98/ch4_fit_profiles.png)
- [Portyankina-41 profile comparisons](/data/vtenishe/cassini_validation_data/studies/ch4_20261016/portyankina2022_solar_41/ch4_fit_profiles.png)

## 7. Uncertainty and limitations

The following must not be collapsed into one misleading error bar:

- **Counting uncertainty:** sample-specific Poisson C1 uncertainty, included in \(\chi^2\).
- **INMS absolute scale:** ±23% common-mode uncertainty on the revised near-on-ram neutral correction. It moves all points in a flyby coherently and is not independent point noise.
- **Mass-15 allocation:** the 0.80–1.00 sensitivity range around the nominal 0.90 changes rates by −11.1% to +11.1%. This is a bracket, not a Gaussian 1σ uncertainty.
- **Cross-channel event classification:** the one flagged E14 point changes fitted rates by no more than 1.66% in the tested fits.
- **Source geometry and kinematics:** the dominant uncertainty. Shared-primary total rates vary by over two orders of magnitude across the three cases.
- **Model inadequacy:** \(\chi^2/N>1\) shows that Poisson noise alone does not explain the residuals. Collisionless propagation, two smooth backgrounds, and steady sources are incomplete descriptions.
- **Identifiability:** a single flyby is a one-dimensional cut through a three-dimensional plume. Near-collinear jet kernels make individual source coefficients and ordinary covariance estimates unreliable.

No formal Gaussian confidence interval is quoted for the global rate because it would omit the dominant discrete geometry/kinematic uncertainty and would be physically misleading. The CSV tables expose the count-scale and mass-15-allocation sensitivities separately.

## 8. Conclusions

1. Calibrated mass 15 provides a useful **CH₄-dominated spatial proxy** for E14, E17, and E18 after revised INMS on-ram calibration, explicit background treatment, and cross-channel screening.
2. The mass-13 data validate the presence of additional spectral/background contributions but are too sparse and impure for a time-resolved two-mass inversion.
3. The most robust quantitative result is the trajectory-integrated CH₄/CO₂ density ratio, 0.233–0.291 for the three primary flybys.
4. A unique absolute CH₄ source rate is **not** determined by these data. Conditional totals depend strongly on the assumed source map and gas kinematics. The hybrid shared-primary result is of order \(10^{25}\) s⁻¹, the Portyankina result several \(10^{26}\) s⁻¹, and the Porco result about \(2\times10^{27}\) s⁻¹.
5. The present data do not establish temporal or tidal modulation of CH₄ independently of source geometry. Such a claim requires contemporaneous H₂O/CO₂ constraints, a collisional multi-species model, and regularized or ensemble source-map inference.

## 9. Commands

Prepare the calibrated profiles from `/data/vtenishe/cassini_validation_data`:

```bash
python3 scripts/prepare_ch4_fit_data.py --validation-data-root /data/vtenishe/cassini_validation_data
python3 -m pytest -q tests/test_ch4_release.py
```

Then run the plume fits from
`/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone`:

```bash
.venv/bin/python run_ch4_fitting_study.py --validation-data-root /data/vtenishe/cassini_validation_data --source-case tenishev2014_hybrid_8 --output-dir /data/vtenishe/cassini_validation_data/studies/ch4_20261016/tenishev2014_hybrid_8
.venv/bin/python run_ch4_fitting_study.py --validation-data-root /data/vtenishe/cassini_validation_data --source-case porco2014_98 --output-dir /data/vtenishe/cassini_validation_data/studies/ch4_20261016/porco2014_98
.venv/bin/python run_ch4_fitting_study.py --validation-data-root /data/vtenishe/cassini_validation_data --source-case portyankina2022_solar_41 --output-dir /data/vtenishe/cassini_validation_data/studies/ch4_20261016/portyankina2022_solar_41
```

Finally, generate the cross-geometry summary from the validation-data root:

```bash
python3 scripts/summarize_ch4_study.py
```

If `pytest` is unavailable, the checks contain plain assertions and can be run without installing anything:

```bash
python3 -c "import runpy; n=runpy.run_path('tests/test_ch4_release.py'); [n[k]() for k in sorted(n) if k.startswith('test_')]"
```

An output directory must not already exist; this prevents accidental replacement of a prior study. Use a new dated directory for a new run.

## References

1. Waite, J. H. et al. (2004), “The Cassini Ion and Neutral Mass Spectrometer (INMS) Investigation,” *Space Science Reviews* 114, 113–231, <https://doi.org/10.1007/s11214-004-1408-2>.
2. Teolis, B. D. et al. (2015), “A Revised Sensitivity Model for Cassini INMS: Results at Titan,” *Space Science Reviews* 190, 47–84, <https://doi.org/10.1007/s11214-014-0133-8>. The correction is explicitly stated to apply to prior Titan, Enceladus, and other Saturn-system neutral results.
3. NASA PDS Atmospheres Node, “Cassini INMS: Ion and Neutral Mass Spectrometer,” including the data-analysis and calibration guidance, <https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/Cassini/inst-inms.html>.
4. Teolis, B. D. et al. (2010), “Detection and measurement of ice grains and gas distribution in the Enceladus plume by Cassini's Ion Neutral Mass Spectrometer,” *JGR Space Physics* 115, A09222, <https://doi.org/10.1029/2009JA015192>.
5. Perry, M. E. et al. (2015), “Cassini INMS measurements of Enceladus plume density,” *Icarus* 257, 139–162, <https://doi.org/10.1016/j.icarus.2015.04.037>.
6. Waite, J. H. et al. (2009), “Liquid water on Enceladus from observations of ammonia and 40Ar in the plume,” *Nature* 460, 487–490, <https://doi.org/10.1038/nature08153>.
7. Waite, J. H. et al. (2017), “Cassini finds molecular hydrogen in the Enceladus plume: Evidence for hydrothermal processes,” *Science* 356, 155–159, <https://doi.org/10.1126/science.aai8703>, especially the supplementary mass-spectrum analysis.
8. Tenishev, V. et al. (2010), “An approach to numerical simulation of the gas distribution in the atmosphere of Enceladus,” *JGR Space Physics* 115, A09302, <https://doi.org/10.1029/2009JA015223>.
9. Tenishev, V. et al. (2014), “Effect of the Tiger Stripes on the water vapor distribution in Enceladus' exosphere,” *JGR Planets* 119, 2658–2667, <https://doi.org/10.1002/2014JE004700>.
10. Porco, C. C., D. DiNino, and F. Nimmo (2014), “How the Geysers, Tidal Stresses, and Thermal Emission across the South Polar Terrain of Enceladus Are Related,” *Astronomical Journal* 148, 45, <https://doi.org/10.1088/0004-6256/148/3/45>.
11. Portyankina, G. et al. (2022), “Modeling the complete set of Cassini's UVIS occultation observations of Enceladus' plume,” *Icarus* 374, 114918, <https://doi.org/10.1016/j.icarus.2022.114918>.
