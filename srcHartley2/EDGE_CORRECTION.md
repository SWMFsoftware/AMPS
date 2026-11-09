# Outer-boundary (edge) correction for comet brightness

All paths in this file are relative to the AMPS repository root.

A spherical AMPS comet run removes every particle that crosses the outer
boundary `mesh.outer_radius_m`. Species whose travel distance before
destruction is much larger than that radius are missing from the far coma.
The line-of-sight brightness computed by the in-run sampler is therefore too
low, and the shortfall grows with impact parameter. Atomic hydrogen is the
important case: H Lyman-alpha from a domain of a few 10^6 km misses a
measurable part of the column at every observation point.

The edge correction multiplies the AMPS brightness at each impact parameter
b by a factor C(b) > 1 from a collisionless reference model that uses the
run's own inputs. The result estimates what the same run would give with no
outer boundary. AMPS output is not changed; the correction is applied when
the output is post-processed and plotted.

| Item | Location |
| --- | --- |
| Reference model and command-line tool | `srcHartley2/visualization/comet_edge_correction.py` |
| Plot option | `srcHartley2/visualization/plot_i3_brightness_comparison.py --edge-correction` |
| Tests | `srcHartley2/tests/test_comet_edge_correction.py`, group `i3` of `srcHartley2/tests/run_all_stage_tests.sh` |
| Worked example | I3 campaign 2: `validation_runs/i3_campaign2_20261108` (run archive on the production host, not tracked) |

## 1. Why the domain edge removes brightness

AMPS integrates the particle-resolved column only over particles that are
inside the domain. Two kinds of atoms are lost:

1. **Atoms beyond the edge.** Every line of sight is cut at the sphere
   r = R. For an impact parameter b, the integral runs over
   |z| <= sqrt(R^2 - b^2) instead of over the whole line.
2. **Atoms whose history crossed the edge.** A parent (H2O, OH, or H2) that
   leaves the domain is deleted. Its daughters, which in nature are born
   outside, are never created, even when their recoil would carry them back
   inside.

The size of the loss is set by each species' scale length, speed times
lifetime, compared with R. For I3 campaign 2 (R = 3e9 m), using the run's
own rates at each heliocentric distance:

| Case | r_h (au) | H lifetime (s) | H scale length at 27 km/s (m) | OH lifetime (s) | H2O lifetime (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.75 | 2.8e7 | 7.6e11 | 2.1e5 | 2.0e5 |
| 2 | 1.66 | 2.5e7 | 6.7e11 | 1.9e5 | 1.8e5 |
| 3 | 1.54 | 2.1e7 | 5.8e11 | 1.6e5 | 1.6e5 |
| 4 | 1.36 | 1.7e7 | 4.5e11 | 1.3e5 | 1.2e5 |
| 5 | 1.36 | 1.6e7 | 4.4e11 | 1.3e5 | 1.2e5 |
| 6 | 1.37 | 1.7e7 | 4.6e11 | 1.3e5 | 1.3e5 |

27 km/s is the H speed from H2O -> H + OH with the 4.04-eV recoil energy in
`input/i3.conf`.

- **OH:** at about 1.6 km/s, OH travels about 3e8 m before destruction, a
  tenth of the radius. Its edge correction is 1.000 at every I3 observation
  point.
- **H:** H travels 100–250 times the domain radius before ionization. Inside
  the domain its density falls almost exactly as 1/r^2. The edge therefore
  acts as a geometric cutoff of a 1/r^2 coma.

## 2. The correction factor

The factor is

```
C(b) = N_full(b) / N_trunc(b),        corrected(b) = AMPS(b) * C(b)
```

- **N_full(b):** the steady-state column of an unbounded coma. Numerically it
  ends at r = 1e12 m, beyond which the I3 contribution is below 0.1%.
- **N_trunc(b):** the column counting only atoms whose whole history lies
  inside r <= R: the parent decay point, every intermediate decay point, and
  the present position. This is the population AMPS can contain.
- **Straight-line motion:** a straight segment between two points inside a
  ball stays inside it. Checking the end points of each free flight is
  therefore exact.

Only the ratio enters. The model's absolute normalization (production rate,
net-to-gross outflow, g-factor) cancels.

### Closed form for a 1/r^2 coma

For n = A/r^2 with a scale length much larger than R:

```
N_full(b)  = pi*A/b
N_trunc(b) = (2A/b) * arctan( sqrt(R^2 - b^2) / b )
C(b)       = (pi/2) / arctan( sqrt(R^2 - b^2) / b )  ~=  1 + (2/pi)(b/R)   for b << R
```

| b/R | 0.01 | 0.03 | 0.1 | 0.2 | 0.3 | 0.48 | 0.63 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C(b) | 1.006 | 1.019 | 1.068 | 1.147 | 1.241 | 1.468 | 1.766 |

This formula is a good first estimate for H. At the outermost I3
observation points it is 0.5–2.5% above the full model, because ionization
removes some of the farthest atoms. It also shows why the correction never
vanishes: even at b = 0.01 R it is 0.6%. The factor diverges as b -> R
because N_trunc -> 0. Corrected values are therefore reported only for
b <= 2R/3.

## 3. How the factor is computed

`comet_edge_correction.py` runs a collisionless steady-state test-particle
model. Its full methodology and implementation are documented in
`srcHartley2/COLLISIONLESS_REFERENCE.md`. It reads everything from the run directory `RUN`:

| Input | Source |
| --- | --- |
| Heliocentric distance, source rates, 1-au reaction rates, nucleus radius, `mesh.outer_radius_m` | `RUN/i3.in`: common settings overlaid with the selected `#begin case: N` section |
| Reaction topology (parent -> products) and recoil (excess) energy of each channel | `RUN/i3.conf` (else `input/i3.conf`): `ChemistryReaction=ID \| parent \| products \| ... \| energy_eV` |
| Radial bulk speed of each source species | `RUN/PT/plots/pic.<source>.s=<i>.out=<K>.dat` at the same output K |

The algorithm:

1. **Relevant species.** The observed species (H and OH), and every species
   that can decay into one of them through a chain. In I3 these are H2O,
   OH, H2, and H; CO and CO2 have no H/OH channels and are skipped.
2. **Parents.** For every source with a positive rate, it draws `--parents`
   decay radii from the AMPS radial speed profile v(r). The optical depth is
   tau(r) = integral of beta / v(r) dr, with beta the parent's total loss
   rate at r_h. Directions are isotropic; the velocity at decay is
   v(r) times the radial unit vector.
3. **Decays.** A channel is chosen with probability rate/beta, including
   channels with no products (ionization), which end the chain. Products
   receive the AMPS photochemical recoil law: an isotropic Gaussian draw in
   the 3N-dimensional product velocity space, centre-of-mass motion removed,
   rescaled to the channel excess energy. Momentum and energy are conserved
   exactly. Intermediate species (OH, H2) fly ballistically for an
   exponentially distributed age before their own decay.
4. **Steady state.** Each observed atom is placed at an age a drawn
   log-uniformly on [1, 1e9] s, with weight a * ln(1e9) * exp(-beta*a).
   This weight is the importance-sampled estimator of
   n(x) = integral over a of exp(-beta*a) * delta(x - x0 - v*a) da.
5. **History flag.** A particle is "inside" if the parent decay point, every
   intermediate decay point, and the present position all satisfy r <= R.
6. **Columns.** Weights are binned in 1,000 logarithmic shells: 700 from the
   nucleus to R, and 300 from R to 1e12 m. The piecewise-constant shell
   densities are projected exactly onto every impact parameter, once with
   all particles (N_full) and once with the inside particles only
   (N_trunc).

**Numerical quality, I3 case 5, 3e6 parents:**
- **Seed noise:** two random seeds differ by at most 0.14% in C at the
  observation points.
- **Run time:** about 4 s per case.
- **Independent check:** a separate implementation used during the
  campaign-2 analysis agrees within 0.7%. It lives in
  `validation_runs/i3_campaign_20261108/scripts/i3_collisionless_mc.py`,
  with the edge moved to 3e9 m.
- **Analytic limit:** in the 1/r^2 test, the full model reproduces the
  closed form within 0.6%.

## 4. Validity and limitations

The correction is valid when the material removed by the edge behaves as the
reference model assumes.

- **Collisionless exterior.** The lost atoms are in the outer coma, where the
  collision rate is negligible. AMPS itself agrees with the truncated
  reference column to within 1–7% for H and 0–2% for OH. The comparison
  covers I3 cases 1, 2, 4, 5, and 6 (campaign 2) and OH in all campaign-1
  cases. It does not use the observations. The residual 1–7% is the
  collisional enhancement of H in the inner coma, which a ratio of two
  collisionless columns leaves unchanged.
- **Same inputs as AMPS.** Rates, recoil energies, product masses, the
  recoil law, and the parent speed profile come from the run. If the product
  energy partition is wrong, AMPS and the correction share the error.
- **Optically thin.** The AMPS sampler computes thin brightness: it flags
  rays whose H line-centre optical depth exceeds 0.1 but does not attenuate
  them. The correction is a ratio of thin columns, so it is consistent with
  the sampler.
  - In campaign 2, tau exceeds 0.1 out to 2e7–6.9e8 m. That includes the
    inner observation points of cases 2–6, with up to 0.46 at the innermost
    point of case 4. Radiative transfer would lower those points; that is a
    separate issue from the edge.
  - The correction adds column mainly at large b, where tau is below 0.1.
- **Same g-factor.** The corrected brightness is AMPS(b) * C(b), which
  assumes the missing atoms have the column-weighted g-factor of the
  included atoms on that ray. The AMPS effective g-factor varies by 1–5%
  between the innermost and outermost I3 observations of a case. The error from this
  assumption is therefore at most a few percent of the added column, and less
  of the total.
- **No solar radiation pressure.** Lyman-alpha radiation pressure
  accelerates H anti-sunward by a = g * h / (m_H * lambda). With the AMPS
  g-factors this is 3.1e-3 to 4.5e-3 m/s^2 for cases 1–6, larger than solar
  gravity. A 27-km/s atom reaching radius r has been displaced by about
  a * r^2 / (2 v^2). That is 2–3% of r at 1e10 m and about 30% at 1e11 m.
  Slower H is displaced 3–15 times more. That is H at 7–15 km/s from
  H2O -> 2H + O, H2 -> 2H (0.49 eV), and OH -> O(3P) + H.
  - At b <= 2e9 m, the column beyond 1e11 m is only (2/pi)(b/1e11) <= 1.3%
    of the total.
  - Radiation pressure therefore changes C by about 1% or less at b <= 5e8 m.
  - At the outermost points of cases 5 and 6 it changes C by an estimated
    few percent.

  Section 8 describes how to include it.
- **Steady state.** The model is the steady state. Compare it with converged
  AMPS outputs, or label it.
- **Spherical symmetry.** The module assumes the 1D radial geometry of
  `geometry = spherical` and does not check it. For 2D or 3D runs, use
  Section 8.
- **What it does not fix.**
  - Model physics inside the domain (production rates, energies, lifetimes,
    collisions).
  - Statistical noise of AMPS.
  - Lack of convergence: an unconverged inner H coma is not corrected.

## 5. Example: I3 campaign 2

Campaign 2 ran the six MAVEN/IUVS cases with OH-matched production rates and
a 3e9-m domain. See `srcHartley2/STATUS.md`, sections "I3 campaign 2 launch"
and "I3 campaign 2 interim results". Cases 1, 2, and 5 converged at
iteration 70,000. The other rows use the latest output while those runs
continue.

| Case | Output | Largest observed b (m) | C at the observation points: median (range) | 1/r^2 estimate at the largest b | H model/obs, raw | H model/obs, corrected |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 1 | 6, converged | 4.3e8 | 1.051 (1.014–1.094) | 1.100 | 0.55 | 0.60 |
| 2 | 6, converged | 3.9e8 | 1.046 (1.012–1.085) | 1.091 | 0.72 | 0.75 |
| 3 | 3, running | 7.4e8 | 1.085 (1.017–1.175) | 1.189 | 0.94 | 1.02 |
| 4 | 6, running | 1.1e9 | 1.144 (1.036–1.296) | 1.321 | 1.00 | 1.14 |
| 5 | 6, converged | 1.4e9 | 1.194 (1.039–1.437) | 1.466 | 0.73 | 0.89 |
| 6 | 6, running | 1.9e9 | 1.272 (1.045–1.726) | 1.764 | 0.89 | 1.16 |

"model/obs" is the median over the H observation points. For OH the factor
is 1.000 in every case.

Point by point, case 5 (converged, `out=6`):

| b (1e5 km) | AMPS (kR) | C(b) | Corrected (kR) | MAVEN/IUVS (kR) | Raw/obs | Corrected/obs | 1/r^2 estimate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.80 | 4.743 | 1.039 | 4.928 | 4.077 | 1.163 | 1.209 | 1.040 |
| 3.60 | 2.439 | 1.076 | 2.624 | 2.853 | 0.855 | 0.920 | 1.083 |
| 5.40 | 1.569 | 1.118 | 1.755 | 1.917 | 0.818 | 0.915 | 1.130 |
| 7.20 | 1.120 | 1.166 | 1.307 | 1.467 | 0.764 | 0.891 | 1.182 |
| 9.00 | 0.850 | 1.221 | 1.039 | 1.207 | 0.704 | 0.860 | 1.240 |
| 10.80 | 0.670 | 1.283 | 0.860 | 0.970 | 0.691 | 0.886 | 1.305 |
| 12.60 | 0.541 | 1.355 | 0.733 | 0.845 | 0.640 | 0.868 | 1.380 |
| 14.40 | 0.444 | 1.437 | 0.639 | 0.760 | 0.585 | 0.840 | 1.466 |

- **Shape.** The raw ratio to the observations falls from 1.16 to 0.59
  across the profile. The corrected ratio is nearly flat, 0.84–0.92 beyond
  the innermost point.
- **Interpretation.** Most of the apparent outer-coma slope mismatch in this
  case is the domain edge, not physics.
- **Remaining deficit.** The roughly 13% that remains is a model result. It
  comes from the production rate, product energies, and lifetimes.

Figures in the run archive:
- `validation_runs/i3_campaign2_20261108/plots/i3c2_all_cases_latest_montage.png`:
  all six cases, raw and corrected.
- `validation_runs/i3_campaign2_20261108/README.md`: the interim tables.

## 6. Tools

### Plot with the correction

`plot_i3_brightness_comparison.py` draws the corrected profile, dashed orange,
when `--edge-correction` is given. It prints the factor at the observation
points:

```sh
python3 srcHartley2/visualization/plot_i3_brightness_comparison.py --no-show \
  --edge-correction --xmax 30 \
  validation_runs/i3_campaign2_20261108/case5_np128/production/PT/plots/pic.H.brightness.s=0.out=6.dat \
  -o case5_H_edge.png
```

Expected output, with the default seed 20261109 and 3e6 parents:

```text
edge correction: domain edge 3.000e+09 m; factor at the observation points median 1.1940, range 1.0390--1.4373; plotted up to 2/3 of the edge
```

Option details:
- **Run directory:** the directory three levels above the table
  (`RUN/PT/plots/...`). Override it with `--run-dir`, and the inputs with
  `--runtime-input` and `--build-input`.
- **Output index:** taken from the file name; the source profile at the same
  output is used.
- **Accuracy and speed:** `--edge-parents` and `--edge-seed` set the
  reference model's parent count and random seed.
- **Default:** without `--edge-correction` the plot is unchanged.

### Tabulate the factor

```sh
python3 srcHartley2/visualization/comet_edge_correction.py \
  validation_runs/i3_campaign2_20261108/case5_np128/production/PT/plots/pic.H.brightness.s=0.out=6.dat \
  --case 5 --csv case5_H_edge.csv
```

The command prints and writes b, AMPS kR, C(b), and corrected kR for all 200
impact parameters. C is blank at and beyond the edge. Python callers can use
`correction_for_run(run_dir, case, output, b)`, which returns `{species: C}`
and the parsed inputs.

### Tests

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest srcHartley2/tests/test_comet_edge_correction.py
```

The tests are self-contained; they synthesize a run directory. They check:
- product parsing, and case-section overlay with rejection of a missing case;
- momentum and energy conservation of the recoil law;
- the exact shell projection of a 1/r^2 density;
- that C is undefined at and beyond the edge;
- an end-to-end run that must reproduce the 1/r^2 closed form within 3%
  (achieved: 0.6%);
- that a short-lived species needs no correction.

They run in the `i3` group of `srcHartley2/tests/run_all_stage_tests.sh`.

## 7. Reproducing the effect

### 7.1 From the existing campaign-2 output

Run the two commands in Section 6. They need only Python 3, NumPy, and
matplotlib, and they read the frozen `i3.in`/`i3.conf` copies in the run
directory.

### 7.2 Build and run from scratch

Build the I3 executable once. Remove exactly `build` before every
recompilation:

```sh
./Config.pl -application=i3
rm -rf -- build
make -j32 EXTRACOMPILEROPTIONS=-O3
```

Regenerate the runtime input. This must leave `input/i3.in` byte-identical
to the tracked file:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/i3/prepare_i3_cases.py
git diff --exit-code input/i3.in
```

Run one case in its own directory. Case 5 converges in about 70,000
iterations, roughly 9.5 h on 128 ranks of one node. Keep `i3.in` and
`i3.conf` beside the output; the correction reads them from there.

```sh
RUN=validation_runs/edge_demo/case5_R3e9
mkdir -p "$RUN"
cp -p amps input/i3.in input/i3.conf "$RUN/"
(cd "$RUN" && mpiexec -n 128 ./amps -i i3.in -case 5 > run.log 2>&1)
```

On this cluster, submit it as a batch job instead. Compute-node `/tmp` is
mounted `noexec`, so point `TMPDIR` at the run directory:

```sh
cat > edge_demo.sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=edge-demo
#SBATCH --partition=superlong128cores
#SBATCH --nodes=1 --ntasks=128 --exclusive --mem=0
#SBATCH --time=2-00:00:00
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
cd "$RUN"
export TMPDIR=$PWD/tmp; mkdir -p "$TMPDIR"
mpiexec -n 128 ./amps -i i3.in -case "$CASE" > run.log 2>&1
EOF
sbatch --export=ALL,RUN="$RUN",CASE=5 edge_demo.sbatch
```

The run stops when three consecutive 10,000-iteration windows change by less
than 5%. The log then ends with `radial simulation converged after 70000
iterations`, and the last brightness table is `out=6`. Apply the correction
with the commands of Section 6, using `$RUN/PT/plots/...`.

### 7.3 The domain-radius experiment

The most direct demonstration runs the same case with two outer radii.
Changing only the two runtime keys keeps everything else identical,
including the surface resolution, which is valid for both radii:

```sh
RUN1=validation_runs/edge_demo/case5_R1e9
mkdir -p "$RUN1"
cp -p amps input/i3.conf "$RUN1/"
sed -e 's/^  mesh.outer_radius_m = .*/  mesh.outer_radius_m = 1000000000/' \
    -e 's/^  brightness.impact_parameter_grid = .*/  brightness.impact_parameter_grid = LOG | 2500.0 | 1000000000 | 200/' \
    input/i3.in > "$RUN1/i3.in"
sbatch --export=ALL,RUN="$RUN1",CASE=5 edge_demo.sbatch
```

Expected:
- **Raw brightness.** The two runs agree near the nucleus and separate with
  b. At b = 5.4e8 m the 1e9-m run should be about 1.57/1.13 = 1.4 times
  lower than the 3e9-m run, using the 1/r^2 estimates for the two edges.
- **Corrected brightness.** At b <= 6.7e8 m (2/3 of the smaller edge) the
  corrected profiles should agree within the combined Monte Carlo and AMPS
  statistical noise, a few percent.
- **Interpretation.** Agreement validates the correction. A systematic
  difference would point to physics the reference model omits near the
  smaller edge.

This experiment has not been run yet. The campaign-1 runs used a 1e9-m edge
but different production rates, so they cannot serve as the pair.

## 8. Extending the correction to 2D and 3D

In 2D (axisymmetric about the Sun–comet line) or 3D, the coma is not
spherical:
- radiation pressure pushes H anti-sunward;
- outgassing is concentrated on the dayside;
- ionization can be anisotropic;
- the domain may be a box rather than a sphere.

The factor then depends on the whole ray, the impact vector in the sky plane
and the observer direction, not on b alone:

```
C(ray) = N_full(ray) / N_trunc(ray),        corrected(ray) = AMPS(ray) * C(ray)
```

The definition of N_trunc does not change: atoms whose whole history stayed
inside the domain. There are two ways to compute it.

### 8.1 Option A: generalize the reference model (post-processing)

Change each step of Section 3 as follows:

1. **Parents.** Replace the radial speed profile by the AMPS parent field.
   Sample decay points from the local production rate n_parent(x) * beta,
   cell by cell, with the cell's bulk velocity plus a Maxwellian of the cell
   temperature. Parents that leave the domain before decaying continue
   ballistically from the outer boundary flux.
2. **Trajectories.** Use x(t) = x0 + v*t + a*t^2/2 with the radiation-pressure
   acceleration a pointing anti-sunward. The magnitude g(v_helio) * h /
   (m * lambda) uses the same g-factor model as the AMPS sampler. When g
   depends strongly on velocity, integrate in a few sub-steps.
3. **History flag.** Curved paths and non-spherical domains lose the
   convexity argument. Test the trajectory between the end points, not only
   the end points:
   - the extremum of |x(t)|^2 (a quartic in t) for a sphere;
   - the six face crossings for a box.
4. **Accumulator.** Replace the shells with a grid that resolves the
   asymmetry:
   - in 2D, (r, theta) about the Sun–comet axis;
   - in 3D, log-r times angular bins.

   Integrate each observation ray through the grid by ray marching, as the
   AMPS sampler does. Columns of points need a volume density, so a grid is
   required.
5. **Factor.** Compute C on the observation rays or on the image-plane pixel
   grid. C varies slowly, so a coarse ray set can be interpolated. About
   1e7–1e8 parents are needed for 1% noise per 3D bin.

This keeps AMPS unchanged. Its weaknesses are the separate source model and
the g-factor assumption of Section 4.

### 8.2 Option B: continue exiting particles inside AMPS (exact)

The missing population is exactly the collisionless descendants of the
particles that cross the outer boundary. AMPS knows each one when it exits:
species, statistical weight w, position, velocity, and time. The exterior
contribution to a ray in steady state is

```
N_ext(ray) = (1/T) * sum over exits during T of  w_i * integral_0^inf S_i(t) * L_ray(x_i(t)) dt
```

- **S_i(t):** the survival probability, exp(-beta*t).
- **L_ray:** the ray footprint the brightness sampler already uses (a shell
  chord in 1D, a pixel solid angle in 3D).
- **T:** the sampling-window length.

Daughters born outside are added by the same decay recursion as in Section 3.
The implementation would work as follows:

1. **Record exits.** In the outer-boundary exit handler, append (species, w,
   x, v) to an exit buffer before deleting the particle.
2. **Process the buffer.** At each sampling step, the brightness sampler
   draws K exterior ages per exit record from the log-uniform, exponentially
   weighted estimator of Section 3, step 4. It advances the record
   ballistically with radiation pressure, and deposits w * weight / (K * T)
   onto the rays.
3. **Use the real g-factor.** Each deposit uses its own Doppler g-factor.
   This removes the same-g assumption.
4. **Daughters.** For parent species, sample the exterior decay and recurse
   into the products.
5. **Report separately.** Report the interior and exterior columns as
   separate Tecplot columns, so the raw AMPS result stays visible.

Option B needs no separate source model, works for any domain shape, and
gives the exterior column with the run's own statistics. It costs one buffer
and the exterior deposit loop, which are small compared with the interior
sampler.

### 8.3 Verification in any dimension

- **Option B against Option A.** In 1D spherical runs, Option B must
  reproduce the factor of `comet_edge_correction.py` within Monte Carlo
  noise.
- **Radius independence.** In every geometry, the domain-radius experiment of
  Section 7.3 must give radius-independent corrected brightness on the rays
  inside both domains.
- **Isotropic limit.** With radiation pressure switched off and isotropic
  sources, a 3D run must reduce to the 1D factor on every ray with the same b.
