ASK_MAVEN.txt -- 3I/ATLAS MAVEN/IUVS H Lyman-alpha and OH(0,0) profiles:
model-data differences (case 1 above all) and questions for the MAVEN/IUVS team
================================================================================

Date: 2026-11-09
Model: AMPS DSMC, 1D spherical comet application, I3 "campaign 2"
Paths in Appendix B are relative to the AMPS repository root.


1. SUMMARY
----------

- OH(0,0) setup. We scaled the water production rate of each of the six
  cases so that the modeled OH(0,0) brightness matches the IUVS profile.
  OH then agrees within 6% in every case.
- H Lyman-alpha results. With OH matched, H Lyman-alpha agrees in the
  perihelion cases (4, 5, 6). It is too faint in the pre-perihelion, Mars
  closest-approach cases:
  - case 1 (Sep 28-29): 1.7x too faint;
  - case 2 (Oct 3): 1.3x too faint.
- Shape of the case-1 deficit. It is a uniform multiplicative factor
  (k = 1.72) over 0.4-4.3e5 km, with no additive offset. The modeled H/OH
  brightness ratio is about half the observed ratio.
- Independent check. The red "Model" curves in the IUVS figures show the same
  H deficit in cases 1 and 2 while fitting OH. The problem is therefore not
  specific to our code.
- Trend. The deficit follows the CO+CO2-to-H2O production ratio. It also
  follows heliocentric distance, approach velocity, and the observing
  geometry near Mars. These quantities all change together across the six
  cases, so the data alone cannot separate them.
- Excluded model-side causes. Domain size, collisions, optical depth,
  convergence, a constant background offset, time-variable production, and
  an icy-grain water source cannot explain it.
- Remaining explanations:
  (a) hydrogen-bearing volatiles other than water;
  (b) fluorescence or solar-input errors at r_dot = -40 km/s;
  (c) observation-side calibration or background subtraction.
  For (c) we need information from the MAVEN team (Section 7).
- Public archive. The case 1-3 observations are in the NASA PDS, labelled
  with the generic target "MARS". A quick look confirms that the scans cross
  the comet and include off-comet frames about 6 deg away, with a large,
  time-variable Lyman-alpha background. Imaging for cases 4-6 was not found
  (Section 9).


2. THE MODEL (what the comparison assumes)
------------------------------------------

Geometry and transport:
- 1D spherical DSMC of H2O, CO, CO2 and their photoproducts H, H2, O, OH.
- Nucleus radius 2.5 km.
- Uniform surface source at 175 K.
- Neutral-neutral collisions.
- Outer boundary at 3e9 m (3e6 km). Particles that cross it are removed.

Photochemistry:
- Per-channel 1-au rates, scaled by the solar-activity index F10.7P at the
  time of each exposure (Moore & Mapaye 2025 law), and by (1 au / r_h)^2.
- Products receive the full channel excess energy as translation. For
  H2O -> H + OH (4.04 eV) this gives H at 27 km/s.
- Lifetimes at the case-1 distance, r_h = 1.75 au:
  - H2O: 2.0e5 s
  - OH: 2.1e5 s
  - H: 2.8e7 s

Brightness:
- Particle-resolved line-of-sight integration, with a per-particle Doppler
  g-factor (Swings and Greenstein effects included).
- H Lyman-alpha g-factor: IKL2020 solar-profile model, driven by the 1-au
  solar Lyman-alpha photon flux from the LASP/LISIRD Version-4 daily
  composite (Machol et al. 2019). This record is referenced to Earth. In
  cases 1-3 the comet saw a solar hemisphere far from the one facing Earth.
  F10.7P (shifted to the Earth-equivalent time) is used only for the
  photochemical rates.
- OH(0,0) g-factor: Schleicher & A'Hearn (1988) Lowell table, -60 to
  +60 km/s, unquenched.
- Emission is optically thin. Rays with H line-centre optical depth above
  0.1 are flagged but not attenuated.

Edge correction for H:
- The 3e9-m boundary removes the far H coma. The H scale length is
  4-8e11 m, so this matters.
- We correct H with a collisionless reference model built from the same
  inputs (factor N_full/N_trunc; srcHartley2/EDGE_CORRECTION.md).
- The correction is 1.01-1.09 at the case-1 points, rising to 1.30-1.73 at
  the outermost points of cases 4-6.
- OH needs no correction.

Inputs per case:

              exposure (UTC)        r_h    r_dot   res.    Q(H2O) s^-1               CO+CO2
  case        start -> end          (au)   (km/s)  (1e5km) table   red   AMPS(OH)  /H2O
  1 pre-Mars  09-28 17:45 -> 09-29 05:43  1.75 -39.7  0.39   8.0e27  1.07e28 3.16e28  1.28
  2 Mars CA   10-03 00:06 -> 10-03 12:04  1.66 -36.2  0.28   1.4e28  2.34e28 6.73e28  0.81
  3 post-Mars 10-09 09:38 -> 10-09 21:36  1.54 -29.6  0.53   5.3e28  6.81e28 1.57e29  0.25
  4 pre-peri  10-24 05:57 -> 10-27 21:40  1.36  -5.8  1.6    2.6e29  2.38e29 4.52e29  0.07
  5 peri      10-28 04:59 -> 10-31 06:05  1.36   0.4  1.8    2.8e29  2.48e29 3.50e29  0.06
  6 post-peri 11-01 04:03 -> 11-03 21:51  1.37   7.1  2.1    2.6e29  1.96e29 3.98e29  0.07

  - "table" is the Model_Table_6cases.pdf value.
  - "red" is the Q printed on the red model curve in the IUVS figures.
  - "AMPS(OH)" is our OH-matched value.
  - CO and CO2 were scaled by the same factor as H2O, so the CO+CO2/H2O
    column is the table ratio.

g-factors at r_h (photons/s/molecule) and the 1-au solar Lyman-alpha flux used:

  case   g_OH(0,0)   g_H(Lya)    F_Lya(1 au), photons cm^-2 s^-1
  1      1.80e-4     9.66e-4     5.10e11
  2      1.88e-4     1.07e-3     5.28e11
  3      2.00e-4     1.20e-3     5.44e11
  4      1.44e-4     1.31e-3     5.24e11
  5      1.26e-4     1.37e-3     5.45e11
  6      2.24e-4     1.36e-3     5.45e11

  These are the g-factors of an atom or molecule at rest relative to the
  nucleus. Scaled to 1 au, g_OH at -39.7 km/s (case 1) is 2.4 times
  g_OH at +0.4 km/s (case 5), a strong Swings dependence.


3. RESULTS: MODEL / OBSERVATION
-------------------------------

Medians over the IUVS points. Case 3 was still running when this was
written.

  case  state                 OH     H raw   H edge-corr   H/OH model/obs
  1     converged (70k it.)   0.95   0.55    0.60          0.63
  2     converged (70k it.)   0.99   0.72    0.75          0.76
  3     running (40k it.)     0.97   0.94    1.02          1.05
  4     converged (80k it.)   1.02   1.01    1.16          1.14
  5     converged (70k it.)   1.00   0.73    0.89          0.89
  6     converged (80k it.)   1.06   0.91    1.19          1.12

Case 1 in detail. The model is edge-corrected; the observations are
digitized from the IUVS figure.

  H Lyman-alpha
    b (1e5 km)  obs (kR)  model (kR)  model/obs   obs - model (kR)
    0.39        1.154     0.729       0.63        +0.43
    0.78        0.913     0.485       0.53        +0.43
    1.17        0.753     0.365       0.48        +0.39
    1.56        0.530     0.291       0.55        +0.24
    1.95        0.476     0.241       0.51        +0.23
    2.34        0.367     0.205       0.56        +0.16
    2.73        0.278     0.178       0.64        +0.10
    3.12        0.244     0.157       0.65        +0.09
    3.51        0.211     0.140       0.66        +0.07
    3.90        0.211     0.126       0.60        +0.09
    4.29        0.151     0.115       0.76        +0.04

    A least-squares fit of obs = k * model + c over all 11 points gives
    k = 1.72 and c = +0.008 kR. The log-slope of model/obs against
    distance is +0.08, i.e. no trend.

  OH(0,0): model/obs = 1.25, 0.97, 0.78, 0.93 at 0.39, 0.78, 1.17, 1.56e5 km.
    The observed OH profile is slightly flatter than the model.

  H/OH brightness ratio at the four OH points:
    observed  2.15, 2.41, 2.56, 3.17
    model     1.09, 1.33, 1.58, 1.86   (0.51-0.62 of observed)

Case 2 behaves the same way at the inner points:
- fit: k = 1.75 with c = -0.23 kR;
- model/obs: 0.63 at 0.28e5 km, rising to 1.1-1.2 at the last three points;
- the observed profile drops steeply beyond 3e5 km (Section 6).


4. EXPLANATIONS EXCLUDED OR BOUNDED ON THE MODEL SIDE
-----------------------------------------------------

  candidate                          test                               effect on case-1 H/OH
  ---------------------------------  ---------------------------------  ----------------------
  finite domain (3e9 m)              collisionless reference model      x1.05 (applied above)
  collisions, incl. with CO/CO2      AMPS vs collisionless model        x1.01
  H optical depth                    sampler line-centre tau            tau <= 0.074, thin
  not converged                      three windows below 5% change      converged
  constant background offset         fit obs = k*model + c              c = +0.008 kR, absent
  rising Q(H2O) (Q x1.75 between     weight by exp(-kT), refit OH       x0.88-0.90, wrong sign
    cases 1 and 2)
  falling Q (after an outburst)      age distribution of the columns    would need Q 9x higher
                                                                         4-16 days earlier
  water from icy grains, 50-80%      extended-source test, refit OH     x0.91-1.01
    of Q, scale 1e4-1e5 km

Why time variability does not help: the H column samples water released
1-8 days earlier (34% from more than 4 days earlier), because about half of
the H comes from OH decay. OH samples water released 1-4 days earlier (14%
from more than 4 days earlier).


5. REMAINING EXPLANATIONS
-------------------------

(a) Hydrogen-bearing volatiles other than water: CH4, CH3OH, C2H6, HCN,
    NH3, H2CO, ...
    - This is the only candidate that adds H uniformly without changing OH,
      and it naturally follows the CO+CO2/H2O trend.
    - Case 1 needs about 3e28 extra H atoms/s, roughly one per water
      molecule. With about 4 H per molecule, that is about 20-25% of Q(H2O).
    - The required fraction would fall to a few percent near perihelion,
      where water dominates.

(b) Fluorescence and solar input at r_dot = -36 to -40 km/s (cases 1 and 2
    only).
    - Q(H2O), and with it the model H, scales as 1/g_OH. An OH g-factor
      overestimated by 1.5-1.7 at -40 km/s, or a Lyman-alpha g-factor or
      solar flux underestimated by that much, would explain case 1.
    - Our solar Lyman-alpha is the Earth-referenced LASP composite, not a
      measurement at the comet's solar longitude.
    - Our OH photodissociation rates have no velocity dependence.

(c) Observation-side effects: background subtraction, the FUV/MUV relative
    calibration, profile extraction (Section 7).

(d) Water photochemistry common to all cases, i.e. how the H2O -> H + OH
    excess energy is shared. With 18-20 km/s H instead of 27 km/s, H rises
    by up to 1.5x in all cases. This sets the overall level, but it cannot
    by itself make case 1 differ from case 5.

Side finding about the red model curves: our OH-matched Q divided by the
red-curve Q equals r_h^2 within +-8% in five of six cases:

  case        1     2     3     4     5     6
  AMPS / red  2.95  2.88  2.31  1.90  1.41  2.03
  r_h^2       3.06  2.76  2.37  1.85  1.85  1.88

Possible reasons: the red curves normalize fluorescence (or the rates) at
1 au, or they define Q differently, e.g. gross against net production. In
AMPS about 23% of the injected water returns to the surface. This factor
cancels in H/OH, so it does not cause the deficit, but it matters whenever
Q values are compared.


6. BACKGROUND SUBTRACTION: WHY IT MATTERS AND WHAT WE SEE
---------------------------------------------------------

Expected comet H Lyman-alpha far from the nucleus. These are our model
values beyond the observed range, scaled to match the observations. The
outer coma falls roughly as 1/b; radiation pressure is neglected.

  distance (km)   case 1   case 2   case 3   case 4   case 5   case 6   (kR)
  5e5             0.17     0.31     0.62     2.0      2.1      1.9
  1e6             0.08     0.15     0.30     0.98     1.04     0.93
  2e6             0.04     0.07     0.15     0.47     0.50     0.45
  5e6             0.016    0.03     0.06     0.18     0.19     0.17
  1e7             0.008    0.014    0.03     0.09     0.09     0.08
  2e7             0.004    0.006    0.013    0.04     0.04     0.04

Observations:

- Case 1. The deficit is multiplicative (c = +0.008 kR). A uniform error in
  the background level cannot produce it. A background or scattered-light
  error proportional to the comet signal, or a different subtraction scheme
  in the Mars-approach geometry, could.
- Cases 1-3. The observed H drops steeply at the last 1-2 points:
  - case 1: 0.21 -> 0.15 kR at 3.9-4.3e5 km;
  - case 2: 0.57 -> 0.27 kR between 2.8 and 3.9e5 km;
  - case 3: a drop near 6-7e5 km.

  The model does not drop like this. Two causes fit:
  - the background region still contains comet H (0.1-0.6 kR at 5e5 km in
    cases 1-3, see the table), so the subtraction removes real coma signal;
  - the edge of the field of view or scan.

  We therefore exclude those points from fits. Restoring them would make the
  inner deficit, if anything, larger.
- Cases 1-3 were observed while 3I/ATLAS was near Mars, with MAVEN in Mars
  orbit. The line of sight then starts inside the Martian hydrogen corona,
  and the interplanetary background and its gradients differ from the
  perihelion cases.
- In cases 4-6 the comet is 5-20 times brighter. The same absolute
  background error therefore matters much less there.


7. QUESTIONS FOR THE MAVEN/IUVS TEAM
------------------------------------

A. Background subtraction (Lyman-alpha above all)

  A1. How was the Lyman-alpha background determined for each case?
      - measured off-comet in the same scan or image, or modeled;
      - which region (angular or projected distance from the nucleus, and
        direction);
      - when, relative to the comet exposure.
  A2. What background level (kR) was subtracted in each case? Was it
      constant, or a function of position across the field? What is its
      uncertainty?
  A3. Interplanetary hydrogen: which model or measurement was used, and how
      large was the IPH gradient across the field? Were the IPH
      Doppler/velocity geometry and the comet's -40 km/s relevant?
  A4. Martian hydrogen corona (cases 1-3 especially):
      - Was the coronal contribution along each line of sight from MAVEN
        removed, and how (model, limb scans, off-comet sky)?
      - What were MAVEN's altitude and orbit phase, and how far was the
        pointing from the Mars limb, during the exposures?
      - Was data taken at low altitude or near periapsis excluded?
  A5. Could the background region contain comet H? Our model predicts
      0.1-0.6 kR at 5e5 km in cases 1-3 and about 1-2 kR at 5e5-1e6 km in
      cases 4-6. Does the steep drop of the last points in cases 1-3 come
      from the subtraction, or from the edge of the field or scan?
  A6. Was the same background procedure applied to the OH (MUV) profiles?
      What were the MUV background level and its uncertainty?
  A7. Is any scattered light or out-of-band contribution (for example
      Lyman-alpha wings, or a geocoronal-like Martian component) included
      in the Lyman-alpha signal, and is it removed?

B. Calibration

  B1. What are the absolute calibration uncertainties of the FUV channel
      (Lyman-alpha) and the MUV channel (OH 308 nm), and of their ratio?
  B2. Did the FUV or MUV sensitivity (gain, detector state, Lyman-alpha
      saturation or non-linearity) change between late September and early
      November 2025, or between observing modes or geometries?

C. Profile extraction

  C1. What does "relative distance at comet" mean exactly? Projected
      sky-plane distance, azimuthal average in annuli, or position along a
      slit or scan direction? If directional, which direction relative to
      the Sun?
  C2. Are the horizontal bars the bin widths (the "resolution" column of the
      table)? What are the PSF and the pointing and nucleus-position
      uncertainties?
  C3. What do the vertical error bars include: statistics only, or also
      background and calibration?
  C4. Could you send the numerical profiles with uncertainties? We
      currently use values digitized from the PNG figures, which have a
      0.09-kR lower plot limit and no numerical uncertainties.

D. The red "Model" curves

  D1. Which model is it (Haser, vectorial, Monte Carlo)? Which water
      outflow speed, H velocities, OH velocity, and lifetimes does it use?
      Are the lifetimes scaled to r_h?
  D2. Are the g-factors (Lyman-alpha and OH) evaluated at r_h or at 1 au?
      Which solar Lyman-alpha flux and line profile, and which OH g-factor
      table and velocity, were used? (See the r_h^2 finding in Section 5.)
  D3. Was Q(H2O) fitted to OH only, to H only, or to both? Your curves
      underpredict H in cases 1-2 as ours do. Did you consider extra H
      sources or calibration effects?
  D4. Does "Q" in the table and on the curves mean gross sublimation or net
      escaping production?

E. Solar input

  E1. Could you provide MAVEN EUVM solar Lyman-alpha irradiance at Mars for
      the six exposure windows? We currently use the Earth-referenced
      LASP/LISIRD composite; the comet was close to Mars in cases 1-3.
  E2. Which OH g-factor do you adopt at r_dot = -39.7 and -36.2 km/s, and
      with which solar spectrum?

F. Other constraints from IUVS

  F1. Are there IUVS detections of, or limits on, other emissions that
      constrain the parents? For example O I 130.4/135.6 nm, C I 156.1/
      165.7 nm, CO Fourth Positive or Cameron bands, CO2+, or CO+. These
      would pin down CO/CO2 and possibly reveal non-water H sources.
  F2. Any detections of, or limits on, H-bearing species other than water,
      from IUVS or from other observers you trust, near 1.7 au?
  F3. Was there evidence of outbursts or rapid variability in the 1-2 weeks
      before case 1, or during the exposures? The H column at the observed
      distances integrates water released 1-8 days earlier; OH, 1-4 days.


8. WHAT WE WILL DO WITH THE ANSWERS
-----------------------------------

- Background or calibration (A, B): correct the profiles, or refit with a
  background term. Restore the outer points of cases 1-3 if the background
  region is found to contain coma.
- Profile definitions (C): average the model exactly as the data were
  extracted (annulus or slit, bin widths, PSF).
- Red-curve details (D): reproduce the red curves, to separate model
  differences from data issues.
- Solar input (E): replace the LASP composite with EUVM Lyman-alpha, and test
  the OH g-factor and lifetime at -40 km/s.
- Other species (F): add measured or bounded H-bearing parents to the model,
  first in the collisionless reference model, then in AMPS.
- Archive (Section 9): reduce cases 1-3 ourselves from the archived L1B
  frames, with our own frame-by-frame background treatment, and compare
  with the digitized profiles.


9. PUBLIC ARCHIVE STATUS OF THE IUVS DATA (checked 2026-11-09)
--------------------------------------------------------------

Where the data are:

- NASA PDS, MAVEN IUVS bundles: raw L1A, calibrated L1B, processed L1C,
  derived L2, and key parameters. The latest delivery found is release 45b
  (PDS delivery 2026-06-27), with coverage through 2025-12-06. The PDS
  Atmospheres node page states coverage through 2025-12-15.
- Calibrated files:
  https://atmos.nmsu.edu/PDS/data/PDS4/MAVEN/iuvs_calibrated_bundle/l1b/
  under {disk,echelle}/2025/{09,10,11}/
- The comet observations are not labelled as such:
  - the PDS4 label Target_Identification is MARS;
  - the FITS Observation table has TARGET_NAME = MARS and
    OBSERVATION_TYPE = CDS, which ordinary Mars apoapse files also carry;
  - the release notes do not mention 3I/ATLAS.

Cases 1-3: FUV and MUV imaging found

  case  orbits        file time span (UTC)              frames per orbit
  1     24497-24500   2025-09-28 17:38 -> 09-29 05:43   45 FUV + 45 MUV
  2     24525-24528   2025-10-02 23:59 -> 10-03 12:04   45 FUV + 45 MUV
  3     24567-24570   2025-10-09 09:31 -> 10-09 21:36   45 FUV + 45 MUV

- File names: mvn_iuv_l1b_apoapse-space-orbitNNNNN-{fuv,muv}_<UTC>_v13_r01.fits
  (directory l1b/disk/2025/MM/). The time spans match the exposure windows
  of the case table.
- Related space-pointed files: echelle "apoapse-space" files on Sep 27-28
  (orbits 24491-24496) and Oct 10 (orbits 24571-24576), each with short
  4-frame FUV/MUV sets.

Verification (orbit 24497, 9 of the 45 FUV files). The pixel RA/Dec was
compared with the JPL Horizons position of C/2025 N1 seen from the Mars
centre (500@499): RA 269.35 -> 269.17 deg, Dec 7.54 -> 7.69 deg, distance
0.293 au.

- Raster across the coma. The nucleus position lies inside the pixel
  footprint in the 18:14 and 18:33 UTC frames, and the other scan steps pass
  0.2-0.5 deg from it.
- Off-comet frames. The first two frames and the last one point about
  6.1 deg away, about 4.6e6 km projected. They are probably the off-comet
  (background) reference.
- Quick-look Lyman-alpha. Integrated over 120.6-122.6 nm on the L1B data,
  median along the slit per frame:
  - the total brightness is 1.8-3.1 kR;
  - the frame containing the nucleus peaks at 3.66 kR against a 2.3 kR slit
    median, about 1.3 kR from the comet (the digitized innermost case-1 point
    is 1.15 kR);
  - at the ~6-deg pointings the background drifts from 3.05 kR (17:41) to
    1.81 kR (18:48), more than 1 kR within an hour.

What this implies:

- Comet H in the reference region is small. At 6 deg (about 4.6e6 km)
  our model gives comet H of about 0.016 kR for case 1, about 10% of the
  faintest case-1 point. Over-subtraction of coma in such a reference cannot
  produce the factor 1.7.
- The background dominates and drifts. The background is 2-20 times the
  case-1 comet signal and changes by more than 1 kR within an hour, so it
  must be removed frame by frame or with a model in time. Questions A1-A4
  matter most for cases 1-3.
- We can reduce cases 1-3 ourselves. The L1B files carry per-pixel
  RA/Dec, the calibrated spectra, and uncertainties, so the profiles can be
  rebuilt independently of the digitized figures.

Cases 4-6: imaging not found

- No FUV/MUV comet files. No space-pointed FUV/MUV files exist between
  2025-10-24 and 11-03.
- Echelle only. There are 11 single echelle "inspace" files (orbits
  24714-24734, Oct 31 -> Nov 3). The Nov 1 file (orbit 24716) points
  straight at the comet: separation 0.0 deg, with the comet 1.3 au from
  Mars. Echelle covers Lyman-alpha (H and D), not OH at 308 nm.
- Routine apoapse files are not comet pointings. A sampled FUV file
  (orbit 24700, Oct 29) points at the Mars disk, 115 deg from the comet.
- The observing plan changed. From Oct 24 the daily number of FUV
  periapse files drops from about 250-300 to about 80, but no comet-pointed
  imaging was identified by file name. A complete search needs the headers
  of about 1,000 FUV/MUV files (about 4 GB).

New questions for the MAVEN/IUVS team

  G1. Which archived files (orbits, segments, channels, frames) produced the
      six profiles? In particular, where are the case 4-6 FUV and MUV comet
      images, or were cases 4-6 built from echelle data?
  G2. Are the frames about 6 deg from the comet in each case 1-3 orbit the
      background reference? How was the background interpolated in time,
      given a drift of more than 1 kR within an hour?
  G3. Will a future PDS release identify the comet observations (target
      name, observation type)?

Sources:
- PDS release summary of 2026-06-27:
  https://pds.nasa.gov/datasearch/subscription-service/SS-20260627.shtml
- PDS Atmospheres IUVS page:
  https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/MAVEN/maven_iuvs.html
- IUVS PDS SIS:
  https://atmos.nmsu.edu/PDS/data/PDS4/MAVEN/iuvs_calibrated_bundle/document/iuvs_pds_sis.pdf
- JPL Horizons API: https://ssd.jpl.nasa.gov/api/horizons.api


APPENDIX A. SENSITIVITIES (CASE 1)
----------------------------------

- H brightness at fixed Q is inversely proportional to the H speed. The
  campaign-1 analysis found that the product energy partition changes H by
  1.5-2x.
- Doubling the OH lifetime lowers the OH-matched Q factors from
  3.95, 4.81, 2.97, 1.74, 1.25, 1.53 to 3.20, 3.65, 2.08, 0.91, 0.70, 0.80
  (cases 1-6), and lowers the model H roughly in proportion.
- The edge correction at the case-1 points is 1.014-1.094 (median 1.051).
- Collisional H enhancement in AMPS: 1.01 in case 1 and 1.08 in case 4.


APPENDIX B. INTERNAL REFERENCES (AMPS repository)
-------------------------------------------------

  srcHartley2/i3/reference/i3_cases.csv          case table and g-factors
  srcHartley2/i3/reference/i3_rate_scaling.csv   OH-matched Q factors
  srcHartley2/EDGE_CORRECTION.md                 edge correction method
  srcHartley2/STATUS.md                          campaign-2 results and the case-1 analysis
  srcHartley2/visualization/references/i3_brightness_MAVEN/
                                                 digitized IUVS profiles
  validation_runs/i3_campaign2_20261108/         run archive (not tracked):
    README.md                                    campaign-2 runs and tables
    analysis_case1_h/README.md                   case-1 diagnosis and test scripts
    analysis_case1_h/maven_archive_*.png         IUVS figures with the red model curves
    plots/i3c2_all_cases_latest_montage.png      six-case comparison
