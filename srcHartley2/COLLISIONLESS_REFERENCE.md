# The collisionless reference model

All paths in this file are relative to the AMPS repository root.

The I3 comparison plots label one curve "corrected for the 3e9-m edge
(collisionless reference)". Several tables also quote an
"AMPS/collisionless" ratio. Both come from a small, independent, steady-state
Monte Carlo model of the photochemical coma, called the reference model
below.

This document covers:
- what the model takes from an AMPS run;
- where the code lives;
- how it differs from AMPS;
- its domain and grids;
- the complete mathematics;
- normalization, verification, limitations, and extensions;
- how to run and change it.

How the edge correction uses the model is in `srcHartley2/EDGE_CORRECTION.md`.

## 0. Quick answers

1. **Is it a separate model, or does it use AMPS output?** It is a separate
   model with its own particle sampling, trajectories, chemistry events, and
   tallies. It reads the same input files as the AMPS run (`i3.in`,
   `i3.conf`), so rates, recoil energies, the source rate, the nucleus
   radius, and the domain radius are identical. It also uses one AMPS output:
   the radial bulk speed of the parent gas, the H2O `Vx [m/s]` profile at the
   same output index. That profile carries the effect of the collisional
   inner coma on the parent flow. Optional comparisons use two more AMPS
   outputs (Section 6):
   - the AMPS OH production rate, for absolute normalization;
   - the AMPS effective g-factor, to convert columns into kR.

   The edge-correction factor needs only the parent speed profile.
2. **Is it part of the tree?** Yes. It is the tracked Python module
   `srcHartley2/visualization/comet_edge_correction.py` (commit 22e8d7f),
   with tests in `srcHartley2/tests/test_comet_edge_correction.py`. Nothing is
   compiled: it needs Python 3 and NumPy, plus matplotlib for plots. An
   older, I3-specific implementation lives in the untracked run archive and
   is still used by two campaign-2 scripts (Section 5.3).
3. **Does it trace particles the way AMPS does?** Free flight is the same
   kinematics: a straight 3D line, which AMPS's DIM1 spherical mover also
   follows within each step before rotating back into the radial frame.
   Everything else differs, as Section 4 tabulates:
   - no time steps and no mesh;
   - parents are not transported but sampled directly at their decay
     points;
   - each particle is followed analytically through its whole life;
   - destruction enters as an exact survival weight or an exponential decay
     age, not as a per-step probability;
   - no collisions, no particle weights or splitting, no boundary
     interactions;
   - the steady state is obtained by sampling particle ages, not by time
     averaging.
4. **How large is its domain, and does it grid?** Trajectories are
   unbounded. The AMPS radius R = 3e9 m is used only as a label ("did this
   particle stay inside R?"). The tallies reach R_big = 1e12 m, 333 times the
   AMPS radius. It uses only post-processing grids:
   - a 20,000-point log grid for the parent optical depth;
   - 1,000 log shells to bin the density (700 inside R, 300 from R to 1e12 m);
   - the 200 impact parameters of the AMPS brightness table, for output.

   No grid is used for transport.

   **How far particles go.** Each atom is placed analytically at an age of
   up to 1e9 s, so a 27-km/s H atom can sit 2.7e13 m out. Only positions
   within 1e12 m are counted. Measured in campaign 2 (cases 1 and 5):

   | Species | Median radius (m) | 99% within (m) | Beyond 3e9 m (AMPS) | Beyond 1e12 m (not counted) |
   | --- | --- | --- | --- | --- |
   | H | 2.6–4.6e11 | 2.3–3.8e12 | 99% | 10–25% |
   | OH | about 3e8 | 2e9 | 0.1% or less | none |

   Raising the counting limit from 1e12 to 1e15 m changes the H column at
   b <= 2e9 m by at most 0.02%.
5. **The mathematics** from the source boundary, through tracking and
   daughter production, to destruction is in Section 3.
6. **How is brightness calculated?** The reference model computes column
   densities only, in m^-2. It never evaluates a g-factor itself. Every
   brightness built from it borrows the AMPS effective g-factor of the same
   ray:
   - edge-corrected AMPS brightness = AMPS kR * C(b);
   - diagnostic reference brightness = eta * N_ref * g_AMPS * 1e-13 kR per
     (m^-2 s^-1).

   Section 3.13 has the details and what a self-consistent brightness would
   need.
7. **What exactly is divided by what in the edge correction?**
   C(b) = N_full(b) / N_trunc(b), both collisionless-model columns from the
   same atoms. N_full counts every atom; N_trunc counts only atoms that AMPS
   could hold. The AMPS column density is not used. Corrected brightness =
   AMPS brightness * C(b). See Section 0a.
8. **Does it have particle time stepping?** No. There is no time step, no
   time loop, and no mover. Each particle's life is evaluated in closed
   form at sampled times: decay radius, decay age, snapshot age. The result
   is the steady state directly, with Monte Carlo noise as the only error.
   Section 3.12 explains this.

## 0a. The two columns: what is divided by what

**In plain words.** The collisionless model counts the atoms on each line of
sight twice:

- once counting every atom, as if there were no domain edge (N_full);
- once counting only the atoms that AMPS could also have, because they and
  their parents never crossed the 3e9-m boundary (N_trunc).

Dividing the first count by the second gives how much bigger the column
would be without the edge. The AMPS brightness is multiplied by that number.
AMPS's own column density is not used in the correction at all. Everything
AMPS computes inside the domain (collisions, the parent flow, its brightness
factors) is kept, and only the missing far part is added: AMPS(b) * (C - 1).

The edge correction is a ratio of two column densities. Both are computed
by the collisionless model, in the same run, from the same sampled atoms.
They differ only in which atoms are counted:

| Column | Which atoms are counted | Physical meaning |
| --- | --- | --- |
| N_full(b) (numerator) | Every sampled atom on the line of sight at impact parameter b, out to 1e12 m | What a coma without a domain edge would contain |
| N_trunc(b) (denominator) | Only atoms whose whole history stayed inside the AMPS radius R = 3e9 m: parent decay point, every intermediate decay point, and present position | What AMPS can contain, since AMPS deletes particles crossing R |

```
C(b)                 = N_full(b) / N_trunc(b)          (>= 1; collisionless model only)
corrected brightness = AMPS brightness(b) * C(b)
```

The AMPS column density (the `Column Density ...` column of the AMPS
brightness table) appears nowhere in this formula. The model's absolute
normalization (gross Q, no surface return flux) also cancels, because it
multiplies numerator and denominator alike.

Worked numbers, H in case 5 of campaign 2 (`out=6`, model with 3e6 parents):

| b (m) | N_full (m^-2) | N_trunc (m^-2) | C = N_full/N_trunc | AMPS column (m^-2) | AMPS / N_trunc | AMPS kR | Corrected kR |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.8e8 | 4.078e16 | 3.924e16 | 1.039 | 3.363e16 | 0.857 | 4.743 | 4.928 |
| 7.2e8 | 1.084e16 | 9.294e15 | 1.166 | 7.840e15 | 0.844 | 1.120 | 1.307 |
| 1.44e9 | 5.253e15 | 3.653e15 | 1.438 | 3.081e15 | 0.843 | 0.444 | 0.639 |

- **Why N_trunc stands in for AMPS.** AMPS / N_trunc is nearly constant
  (0.843–0.857). It is the AMPS net-to-gross outflow (eta = 0.79), times the
  collisional increase of H in AMPS (about 1.07). Because the ratio does not
  depend on b, N_trunc has the same shape as the AMPS column. That is why it
  can represent "what AMPS contains" in the denominator.
- **Why the AMPS column is not the denominator.** N_full / N_AMPS would be
  1.21, 1.38, and 1.70 here. That would mix the edge loss with the
  gross/net normalization (1/0.79) and would remove the AMPS collisional
  increase. The result would effectively replace AMPS with the collisionless
  model.
- **What the correction keeps.** Inside the domain, everything stays AMPS's:
  collisions, the parent flow, per-particle g-factors. The model supplies
  only the fractional loss caused by the edge. The added brightness is
  AMPS(b) (C - 1).
- **Where the AMPS column is used.** Only in the separate collisional
  diagnostic ("AMPS/collisionless", Sections 3.13 and 7), never in the edge
  correction.

## 0b. Checking common readings of the method

| # | Statement | Correct? | What the model actually does |
| --- | --- | --- | --- |
| 1 | It starts from the H2O source rate and temperature | Partly | The source rate Q is used only as a weight: each test parent stands for Q/N_P molecules per second. The 175-K temperature is not used; its effect enters only through the AMPS speed profile (row 2). |
| 2 | It assumes a bulk velocity and applies the Haser model | Partly | No assumed constant speed. It uses the AMPS H2O radial bulk-speed profile v(r) at the same output: in case 1, 0.36 km/s at 3 km, 0.55 at 1e5 m, 0.63 at 1e8 m, and 1.8 at 1e9 m; in case 5, 0.43, 0.73, 1.18, and 2.6 km/s. With that speed the parent is distributed like a variable-speed Haser parent, n = Q exp(-tau)/(4 pi r^2 v) with tau = integral of beta/v dr, which reproduces the shape of the AMPS H2O density (Section 3.0). Only the parent is Haser-like; the daughters are Monte Carlo particles with real recoil velocities. |
| 3 | The profile gives the start location of a new trajectory | Yes | The start radius is sampled: tau* = -ln U, then r_d with tau(r_d) = tau*, and an isotropic direction. For H from OH or H2, the start is the sampled decay point of that OH or H2. |
| 4 | It generates a velocity vector | Yes, as follows | Parent velocity at the decay point (v(r_d) e, radial) plus a recoil from the AMPS law: Gaussian draw in 3N, centre-of-mass motion removed, rescaled to the channel excess energy. H2O -> H + OH gives H 27 km/s and OH 1.6 km/s, opposite and in random directions; H from OH gets 15–42 km/s (Section 3.4). |
| 5 | It traces the trajectory up to 1e6 s | No | Nothing is traced step by step. Each H or OH is placed once, at a random age between 1 s and 1e9 s (log-uniform), on its straight line, with its survival probability as weight (Section 3.7). The limit is 1e9 s, not 1e6 s: H lives about 2e7 s, so 1e6 s would miss most of it. Positions beyond 1e12 m are not counted. |
| 6 | It computes the column density along the line of sight | Yes | Shell binning (1,000 shells) and exact chord integration at each impact parameter, done twice from the same atoms: N_full and N_trunc (Section 0a). |
| 7 | It derives the column density from the AMPS output | No | The AMPS column is not used for the correction. Both columns in the ratio are model columns (Section 0a). AMPS supplies only the H2O speed profile (row 2). |
| 8 | It multiplies the AMPS brightness by (model column / AMPS column) | No | It multiplies by N_full / N_trunc, both from the model. Using the AMPS column instead would replace AMPS's interior physics with the collisionless model (Section 0a). |

## 1. Purpose and uses

| Use | Quantity | Where |
| --- | --- | --- |
| Edge correction | C(b) = N_full(b) / N_trunc(b) | `comet_edge_correction.py`, `--edge-correction` plots, montages |
| Collisional effect on H and OH | AMPS(b) / reference_trunc(b), normalized | "AMPS/collisionless" columns in the campaign README and `STATUS.md` |
| Fast what-if tests | Columns for modified inputs or weights | Campaign-1 Q-scaling predictions; case-1 production-history and icy-grain tests |
| Far-coma extrapolation | N_full beyond the AMPS domain | Background table in `srcHartley2/i3/COMPARISON.md` |

## 2. Assumptions

1. **Steady state.** Sources are constant in time; Section 9 gives the
   time-dependent reweighting.
2. **Spherical source.** Parents leave the nucleus (radius R_n) isotropically
   and move radially with the AMPS bulk speed v(r). Their thermal dispersion
   about v(r) is neglected.
3. **Photochemistry.** Channel c of species X has the first-order rate
   k_c = A_c(1 au) / r_h^2, from `i3.in`. The 1-au values already include the
   F10.7P solar-activity scaling. The total loss rate is beta_X = sum over c
   of k_c, including product-less ionization channels.
4. **Recoil.** It follows the AMPS law; the full excess energy from
   `i3.conf` goes into translation.
5. **Free flight.** Straight lines in the comet frame, with no collisions,
   gravity, or radiation pressure. The nucleus neither absorbs nor re-emits
   products.
6. **Observable.** The column density of the observed species (H and OH by
   default) along straight lines of sight at impact parameter b.

## 3. Mathematics, step by step

### 3.0 Step-by-step walkthrough

Line numbers refer to `srcHartley2/visualization/comet_edge_correction.py`.

**Step 0: what is read, and what is not.**

| Item | Used? | How |
| --- | --- | --- |
| H2O source rate Q (`source.H2O.rate_s-1`) | Yes | Only as the weight scale: each test parent stands for Q / N_P molecules per second |
| Surface temperature (175 K), Maxwellian injection, Knudsen layer, back-flux | No | Their effect reaches the model only through the AMPS H2O speed profile v(r) |
| AMPS H2O `Vx [m/s]` profile at output K | Yes | Outward speed v(r) of the parent gas |
| AMPS H2O density | No | Not read; the model's implied density reproduces its shape (below) |
| Rates, products, recoil energies, R_n, R | Yes | From `i3.in` and `i3.conf` |

**Step 1: H2O is not propagated through the domain.** There is no H2O
particle motion and no H2O density computation. The model uses the AMPS
speed profile only to know where a water molecule that left the nucleus
radially would be destroyed. The optical depth tau(r) = integral of
beta_H2O / v dr is computed once on a 20,000-point log grid (lines 347–349).

This implies the variable-speed Haser parent density

```
n_H2O(r) = Q exp(-tau(r)) / (4 pi r^2 v(r))
```

Compared with the AMPS H2O density, the ratio AMPS / implied is constant
with radius:

| r (m) | 3e3 | 1e4 | 1e5 | 1e6 | 1e7 | 1e8 | 1e9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| case 1 | 0.768 | 0.762 | 0.763 | 0.764 | 0.765 | 0.761 | 0.765 |
| case 5 | 0.792 | 0.790 | 0.793 | 0.795 | 0.793 | 0.791 | 0.810 |

The constant is the AMPS net/gross outflow, eta (Section 6). The parent
distribution therefore has the AMPS shape at every radius; only the gross
normalization differs, and that cancels in C.

**Step 2: start each parent at its destruction point** (lines 350–356, for
each of N_P test parents, vectorized):

```
tau* = -ln(U1)                                   (U ~ uniform(0,1); line 350)
r_d  = tau^-1(tau*)                               (linear interpolation on the grid; line 352)
mu = 2 U2 - 1,  phi = 2 pi U3                     (isotropic direction; isotropic(), line 215)
e    = (sqrt(1-mu^2) cos phi, sqrt(1-mu^2) sin phi, mu)
x_d  = r_d e,      v_P = v(r_d) e,     inside = (r_d <= R)
```

**Step 3: choose the reaction** (`decay_products`, line 322). The channel is
c with probability k_c / beta_H2O. For case 1, H2O -> H + OH has probability
0.78. Ionization channels produce nothing, and that test parent ends there.

**Step 4: give the products their start position and velocity** (lines
326–331, `recoil()` line 222). Every product starts at x_d. Its velocity is
v_P + u_i, where u_i comes from the Gaussian-in-3N, centre-of-mass-removed,
energy-rescaled recoil (Section 3.4). For H + OH: H gets 27.0 km/s and OH
1.60 km/s, opposite and in a random direction, plus v_P.

**Step 5: count observed products** (`follow` line 333, `snapshot` line 306).
For H and OH:

```
a = exp(ln 1 + U4 (ln 1e9 - ln 1))                          (log-uniform age)
x = x_start + v a
weight = (Q / N_P) * a * ln(1e9) * exp(-beta_X a)
```

The weight is added to its shell in N_full, and also in N_trunc when inside
and |x| <= R.

**Step 6: start the next generation** (`decay_products`, lines 318–321). An
OH or H2 product additionally gets a decay age a' = -ln(U5) / beta_X and the
new start point x' = x_start + v a', with inside' = inside and |x'| <= R.
Steps 3–5 repeat there with that particle's velocity. This is how H from
OH -> O + H and from H2 -> 2H gets its start position and velocity.

**Worked example (case 1, one test parent).**
1. tau* = 0.69 gives r_d of about 8.5e7 m (the median) and v_P = 0.62 km/s
   outward.
2. Channel H2O -> H + OH.
3. H starts at x_d with about 27 km/s in a random direction. Its snapshot age
   might be 2e5 s, so it is counted 5.4e9 m away. That is in N_full only,
   beyond R = 3e9 m. Its weight factor is exp(-2e5 / 2.8e7) = 0.993.
4. OH starts at x_d with about 1.7 km/s, is counted at its own snapshot age,
   and decays after a' of about 1.5e5 s, about 2.5e8 m further on.
5. There, OH -> O(3P) + H gives an H that starts at the OH decay point with
   v_OH + 15.1 km/s.

### 3.1 Frame and notation

- **Frame.** Comet-centred, non-rotating Cartesian frame; nucleus at the
  origin. x is a position, r = |x|, and e_r = x/r.
- **Species.** P is a parent (H2O). X is any species that is observed, or
  that can produce an observed species through a chain: here OH, H2, and H.
- **Rates.** Q is the gross source rate in s^-1, and beta_X the total
  destruction rate of X in s^-1.
- **Model sizes.** N_P is the number of test parents; a test particle is one
  Monte Carlo sample.

### 3.2 Boundary conditions

**Inner boundary (source).**

| | Reference model | AMPS |
| --- | --- | --- |
| Parent emission | Q parents per second leave r = R_n with isotropic directions. Their subsequent radial motion is prescribed by v(r). No Knudsen layer, back-flux, or thermal spread is modelled; their net effect on the bulk flow is inside v(r) taken from AMPS. | Surface injection is a Poisson process at rate Q / W (W is the particle weight). Velocities follow a 175-K half-Maxwellian flux. Particles that return to the surface, about 19–24% of the H2O, are re-emitted thermally, so the net outflow is 0.76–0.81 Q. |
| Daughters reaching the nucleus | They pass through it; neither absorbed nor re-emitted. | Any species reaching the inner face is re-emitted from a 175-K half-Maxwellian. |

**Outer boundary.**

| | Reference model | AMPS |
| --- | --- | --- |
| Transport | None; trajectories are unbounded. | A particle crossing r = R is deleted. |
| Use of R | Only to compute the "inside" label (Section 3.9). | Domain limit. |
| Tally limits | Positions with r > R_big = 1e12 m are not tallied. Ages are limited to [1 s, 1e9 s] (Section 3.7). | Not applicable. |

### 3.3 Parent transport and decay radius

A parent moves radially with dr/dt = v(r). Its survival probability to
radius r is

```
S_P(r) = exp(-tau(r)),    tau(r) = integral from R_n to r of  beta_P / v(r') dr'
```

The decay-radius density is

```
p(r) = beta_P / v(r) * exp(-tau(r))
```

It is normalized to 1 - exp(-tau(R_big)), which is 1 for practical purposes.

**Equivalence with AMPS.** The parent flux 4 pi r^2 n_P(r) v(r) =
Q_net exp(-tau(r)) gives p(r) = 4 pi r^2 beta_P n_P(r) / Q_net. This is
exactly the AMPS parent destruction rate per unit radius. So the reference
places parent decays where AMPS does, provided the radial flux of the AMPS
profile is conserved apart from photolysis.

**Numerical implementation.**

1. Build the grid r_j = R_n (R_big/R_n)^(j/(J-1)), with j = 0..J-1 and
   J = 20,000.
2. Set v_j by linear interpolation of the AMPS profile. Duplicate nodes at
   block faces are averaged first, and only nodes with v > 0 are used.
   Beyond the last node, v is held at the last value.
3. Integrate tau_j = sum over i < j of beta_P (r_{i+1} - r_i) / ((v_i +
   v_{i+1}) / 2) (trapezoid in 1/v).
4. Draw N_P values tau* ~ Exp(1) and discard tau* >= tau_{J-1}.
5. Invert by linear interpolation: r_d = interp(tau*, tau_j, r_j).
6. Set the direction e ~ isotropic (mu = cos(theta) ~ U(-1, 1),
   phi ~ U(0, 2 pi)), the decay position x_d = r_d e, and the parent
   velocity at decay v_P = v(r_d) e.

### 3.4 The photochemical event: channel, products, and parent removal

**Channel choice.** At a decay point the channel is chosen with probability
k_c / beta_X. The parent ceases to exist: it produces its product set once
and is never tallied or moved again. A channel with products `NONE`
(ionization) removes the particle without products, which ends the chain.
In AMPS, the parent reacts within a local step dt with probability
1 - exp(-beta dt); I3 uses AMPS's systematic-residual variant for the
product multiplicities.

**Product list.** The products are the expanded stoichiometry: `2*H + O`
becomes [H, H, O]. Each product receives

```
v_i = v_X + u_i
```

where v_X is the reacting particle's velocity and u_i its recoil velocity
in the parent frame.

**Recoil law** (identical to AMPS). For N products of masses m_i and excess
energy E:

```
g_i ~ N(0, I_3) independently, i = 1..N
G   = sum_i m_i g_i / sum_i m_i                  (centre-of-mass velocity of the draw)
u'_i = g_i - G                                   (=> sum_i m_i u'_i = 0)
s   = sqrt( E / (1/2 sum_i m_i |u'_i|^2) )
u_i = s u'_i                                     (=> 1/2 sum_i m_i |u_i|^2 = E)
```

- **Two products.** The directions are isotropic and opposite, and the
  speeds are fixed:

  ```
  |u_1| = sqrt( 2 E m_2 / (m_1 (m_1 + m_2)) ),    |u_2| = |u_1| m_1 / m_2
  ```

- **Three products.** The draw spreads the energy among the three bodies.
- **E = 0 or a single product.** u = 0.

Recoil speeds in the I3 network:

| Channel | E (eV) | Product speeds (km/s) |
| --- | ---: | --- |
| H2O -> H + OH | 4.04 | H 27.0, OH 1.60 |
| H2O -> H2 + O | 3.94 | H2 18.3, O 2.31 |
| H2O -> 2H + O | 0.697 | H mean 7.6 (5–95%: 3.7–10.5), O 0.78 |
| OH -> O(3P) + H | 1.27 | H 15.1, O 0.95 |
| OH -> O(1D) + H | 7.73 | H 37.3, O 2.35 |
| OH -> O(1S) + H | 10.0 | H 42.4, O 2.67 |
| H2 -> 2H (excited) | 0.488 | H 6.8, 6.8 |
| H2 -> 2H | 8.22 | H 28.1, 28.1 |
| H2O, OH, H2, H, O ionization | n/a | no products; the chain ends |

Masses come from the module's `MASS_AMU` table: H 1.008, H2 2.016, O 15.999,
OH 17.007, H2O 18.015, CO 28.010, CO2 44.009 u. They are not read from
`i3.conf`.

### 3.5 Daughter tracking

A product born at x_b with velocity v moves on x(t) = x_b + v t. If the
product can decay into an observed species (OH, H2), its decay age is drawn
as a ~ Exp(beta_X). The next event is at x_b + v a, and Section 3.4 repeats
there. The recursion depth equals the chain length: H2O -> OH -> H and
H2O -> H2 -> H are two levels.

### 3.5a Where each trajectory starts

No product starts at the surface. Every trajectory starts at the event that
creates the particle:

| Particle | Start position | Start velocity |
| --- | --- | --- |
| H2O (parent) | Not followed step by step. Only its decay point x_d = r_d e is sampled (Section 3.3). | v(r_d) e, the AMPS bulk speed, radial |
| H and OH from H2O -> H + OH | x_d | v(r_d) e + recoil: H 27.0 km/s, OH 1.60 km/s, opposite and isotropic |
| H from H2O -> 2H + O | x_d | v(r_d) e + three-body recoil, about 4–11 km/s |
| H2 from H2O -> H2 + O | x_d | v(r_d) e + 18.3 km/s |
| H from OH -> O + H | The OH decay point, x_d + v_OH a_OH | v_OH + 15.1, 37.3, or 42.4 km/s (O(3P), O(1D), O(1S)) |
| H from H2 -> 2H | The H2 decay point, x_d + v_H2 a_H2 | v_H2 +- 6.8 or 28.1 km/s |

The direction e of the parent is isotropic (mu ~ U(-1, 1), phi ~ U(0, 2 pi)).
Recoil directions are isotropic in the parent frame.

Typical values, campaign 2:

| Quantity | Case 1 | Case 5 |
| --- | --- | --- |
| H2O decay radius r_d, 5% / 50% / 95% | 6.0e6 / 8.5e7 / 4.2e8 m | 5.0e6 / 8.7e7 / 4.9e8 m |
| H2O decays inside 3e9 m | 99.99% | 100% |
| H2O bulk speed at the median decay radius | 0.62 km/s | 1.16 km/s |
| OH median decay age | 1.5e5 s | 8.6e4 s |
| OH median speed | 1.7 km/s | 2.0 km/s |

So first-generation H starts within about 1e7–5e8 m of the nucleus, and H
from OH starts a further few times 1e8 m out. AMPS creates daughters the
same way: at the parent's current position, with the parent's velocity plus
the same recoil. The one difference is that the AMPS parent velocity also
contains its thermal spread and collisional history.

### 3.6 Destruction ("deletion") of observed species

An observed species X (H, or OH) survives to age a with probability
exp(-beta_X a). The reference does not delete particles stochastically.
It carries this survival probability as a weight in the density estimator
(3.7), which is unbiased and has lower variance.

OH is both observed and an intermediate, so it gets two independent draws:
- a snapshot age, for its density (3.7);
- a decay age, for its products (3.5).

Both are unbiased, because every born OH eventually decays exactly once,
through one of its channels.

In AMPS, by contrast, a particle is removed when its per-step reaction
occurs (Section 3.4) or when it leaves through r = R.

### 3.7 Steady-state density

A stationary stream of species X born at x_b with velocity v at rate S
(s^-1) produces the density

```
n(x) = S * integral from 0 to infinity of  exp(-beta_X a) * delta(x - x_b - v a)  da
```

The age integral is evaluated by Monte Carlo, with one age per born test
particle, drawn log-uniformly:

```
a ~ p(a) = 1 / (a ln(a_max / a_min)),   a in [a_min, a_max] = [1 s, 1e9 s]
w = exp(-beta_X a) / p(a) = a ln(a_max / a_min) exp(-beta_X a)
```

The particle is tallied at x_b + v a with weight S w.

- **Unbiasedness.** E_a[ w delta(x - x_b - v a) ] = integral from a_min to
  a_max of exp(-beta_X a) delta(...) da.
- **Bias from the age limits.** Below 1 s a 27-km/s atom travels less than
  30 km. For H, beta_H = 3.6e-8 to 6.2e-8 s^-1, so exp(-beta_H * 1e9) is
  below 1e-15.
- **Variance.** The log-uniform choice gives a constant number of samples
  per decade of age, and so of radius, from 1e3 to 1e13 m.

**Relation to AMPS sampling.** AMPS time-averages its particle population
over a 10,000-step window. In steady state that average equals the
source-weighted age integral above (ergodic equivalence).

### 3.8 Weights, units, and tallies

Each test parent represents Q / N_P real parents per second. Every
descendant carries this scale, so a tallied weight Q/N_P * w has units of
particles: (s^-1) times (s).

**Shells.** The edges are

```
700 log shells in [R_n, R],  then  300 log shells in [R, R_big]
```

The sum of weights in shell k is W_k, and its density is n_k = W_k / V_k
with V_k = 4/3 pi (r_{k+1}^3 - r_k^3).

**Column projection.** The column of the piecewise-constant shell density
along the line of sight at impact parameter b is exact:

```
N(b) = 2 * sum over k with r_{k+1} > b of
          n_k [ sqrt(r_{k+1}^2 - b^2) - sqrt(max(r_k, b)^2 - b^2) ]
```

Discretization error of the shells:
- for n ~ 1/r^2, the model's 1,000 shells reproduce the analytic column
  within 1e-3 (at most 9e-4 at b = 1e9 m);
- the unit test, with 3,000 shells, reaches 6e-5;
- N_full and N_trunc share the shells, so most of this error cancels in C.

Material beyond R_big is not counted. For a 1/r^2 coma that is at most
(2/pi)(b/R_big), or 0.13% at b = 2e9 m.

### 3.9 The "inside" label (AMPS-equivalent truncation)

A tallied particle is "inside" when all of the following hold:
- parent decay radius r_d <= R;
- every intermediate decay point is within R;
- the tallied position is within R.

A straight segment between two points inside a ball stays inside it, so
these end-point tests guarantee that the whole history stayed inside. The
radial parent path from R_n to r_d is inside if r_d <= R.

This is the population that survives AMPS's deletion at r = R. The density
uses all particles for N_full and only inside particles for N_trunc. They
are the same particles from the same random stream, which is why C = N_full
/ N_trunc has much smaller noise than either column alone.

### 3.10 Outputs

- **Edge factor.** C(b) = N_full(b) / N_trunc(b). It is NaN where
  N_trunc = 0 or b >= R(1 - 1e-9).
- **Edge-corrected AMPS brightness.** AMPS_kR(b) * C(b). The plots show it
  to b <= 2R/3.
- **Reference brightness** (collisional diagnostics only): see Section 3.13.

### 3.11 Statistical error

- **Shell variance.** Each shell's variance comes from its weighted samples.
  The ratio C is much less noisy than either column because the columns are
  correlated.
- **Measured noise.** Two seeds for case 5 with 3e6 parents differ by at
  most 0.14% in C at the observation points.
- **How to estimate it.** Repeat with `--seed` (CLI) or `seed=` (API).

### 3.12 Time: no particle time stepping

The reference model has no time step, no time loop, and no particle mover.
Time enters only as continuous random variables of each particle's life,
and positions at those times are evaluated exactly:

| Quantity | How it is obtained | Position used |
| --- | --- | --- |
| Parent transit to its decay point | Implicit in inverting tau(r) (Section 3.3); the transit time integral of dr/v is available if needed (Section 9) | x_d = r_d e |
| Decay age of an intermediate (OH, H2) | a ~ Exp(beta_X) | x_b + v a |
| Snapshot age of an observed species (H, OH) | a ~ log-uniform on [1 s, 1e9 s], weighted (Section 3.7) | x_b + v a |

Consequences:

- There is no Courant condition and no time-discretization error. There is
  no weight adaptation, sampling window, or convergence monitoring.
- Every run is the steady state by construction. The only errors are Monte
  Carlo noise (Section 3.11) and shell discretization (Section 3.8).
- The only quadrature is the parent optical-depth integral, a spatial
  trapezoid rule on 20,000 log-spaced radii. It is not time stepping.

AMPS works the other way. It advances every particle with species- and
cell-local time steps, dt = Courant (0.8) * cell size / characteristic
speed. In case 1 of campaign 2, dt runs:

| Species | dt in the surface cell | dt in the outermost cells |
| --- | ---: | ---: |
| H | 5.7e-6 s | 3.8e2 s |
| OH | 6.8e-5 s | 4.6e3 s |
| H2O | 1.1e-4 s | 7.5e3 s |

AMPS also:
- applies photolysis with probability 1 - exp(-beta dt) per step;
- injects surface particles as a Poisson process per step;
- reaches steady state only after 70,000–80,000 iterations;
- averages over 10,000-iteration windows.

### 3.13 Brightness

The reference model's outputs are column densities N_full(b) and
N_trunc(b) (m^-2). It does not evaluate g-factors, Doppler shifts, optical
depth, or nucleus occultation. Brightness appears in three places, and all
three borrow the g-factor from AMPS.

**1. Edge-corrected AMPS brightness** (plots, `comet_edge_correction.py`
CLI, `--edge-correction`):

```
B_corr(b) = B_AMPS(b) * C(b),        C(b) = N_full(b) / N_trunc(b)
```

- The reference supplies only the ratio. The brightness itself is the AMPS
  one: per-particle Doppler g-factors, AMPS unit conversion, tau flags.
- Implied assumption: the missing atoms, beyond R or with a history that
  crossed R, emit with the same column-weighted mean g as the AMPS atoms on
  that ray. In other words, B_missing = B_AMPS (C - 1).
- The AMPS effective g varies by 1–5% across the observation points of a
  case, so the error of this assumption is a few percent of the added part.

**2. Diagnostic reference brightness** (`i3c2_mc_compare.py`; campaign-1
`i3_q_scaling_analysis.py`):

```
B_ref(b) = g_AMPS(b) * eta * N_ref(b) * K
```

- **g_AMPS(b):** the `Effective Photon g Factor` column of the AMPS
  brightness table. It is the column-weighted mean of the AMPS per-particle
  g on that ray, interpolated in log b.
- **N_ref(b):** the reference column (full or truncated), interpolated in
  log b.
- **eta:** the net-to-gross normalization of Section 6.
- **K:** the Rayleigh conversion, B[R] = g[s^-1] N[cm^-2] / 1e6, so
  K = 1e-13 kR per (m^-2 s^-1). The scripts compute it from the AMPS table:
  - `i3c2_mc_compare.py` uses mean(kR/(N g)), which is exactly K;
  - the campaign-1 script uses mean(kR/N) / mean(g) or median(kR/(g N)),
    equal to K within the 1–5% variation of g.

**3. Far-coma extrapolation** (the background table in
`srcHartley2/i3/COMPARISON.md`):

```
B(b) = B_AMPS(b0) * C(b0) * N_full(b) / N_full(b0),        b0 = 1e9 m
```

This assumes g(b) = g(b0) beyond the AMPS domain.

**Not included.** The reference particles' own Doppler g (Swings and
Greenstein effects), line optical depth, nucleus occultation, observer
geometry beyond b (1D symmetry), and radiation-pressure asymmetry.

**Adding a self-consistent brightness** would take four steps:
1. Give each tallied particle a heliocentric radial velocity
   v_h = (V_comet + v) . e_r, as in `CometBrightness.cpp`. Here e_r is the
   Sun-to-comet unit vector, so a positive v_h means receding from the Sun.
   V_comet = r_dot * e_r, from `brightness.heliocentric_radial_velocity_km_s`.
   e_r and the line of sight come from `brightness.observation_geometry` in
   `i3.in`.
2. Evaluate g(v_h) with the same models as AMPS: IKL2020 for H, and the
   Lowell OH table from -60 to +60 km/s. They are implemented in C++ in
   `src/models/resonance_fluorescence/`, so they need a Python port or a
   small helper that exports tables.
3. Tally the emission density sum(w g) in (r, mu) bins, mu = x . e_r / r.
   g depends on velocity direction, which is correlated with position
   (radial outflow), so the emission is not spherically symmetric even
   though the density is.
4. Integrate along the actual line of sight.

Validate against AMPS OH, which is collisionless in practice: agreement to
1–2% is expected.

## 4. Comparison with AMPS's spherical (DIM1) application

| Aspect | AMPS (`srcHartley2/RadialSimulation.cpp`, `src/pic/pic_mover.cpp`) | Reference model |
| --- | --- | --- |
| Time | Explicit time stepping; species- and cell-local dt = Courant * cell size / characteristic speed | No time steps; every trajectory is analytic in age |
| Mesh | Radial AMR mesh, 2.5 km to 3e9 m, about 380 blocks x 8 cells; finest cell about 0.35–0.44 m at the surface | None for transport; tally shells only |
| Particle state | Radius plus a 3-component velocity in the local radial frame. Each step moves on a straight 3D line, then the position is rotated back to (r, 0, 0) and the velocity rotated with it. | Full 3D Cartesian position and velocity in a fixed frame |
| Free-flight kinematics | Straight line within each step | Straight line for the whole flight (same kinematics) |
| Source | Poisson surface injection, 175-K half-Maxwellian flux | Parents placed directly at sampled decay radii, moving with the AMPS bulk speed |
| Parent flow | Self-consistent: collisions, Knudsen layer, back-flux | Taken from AMPS (`Vx` profile) |
| Collisions | NTC, H2O/CO/CO2 cross sections, completion law for other pairs | None |
| Photochemistry | Per-step reaction probability 1 - exp(-beta dt); channel by rate; products at the parent position with the same recoil law; systematic-residual product multiplicities; parent removed | Decay radius or age sampled exactly; channel by rate; same recoil law; one product set per test particle |
| Particle weights | Species-local weights, adapted every 10,000 steps (split/roulette); up to 20 products per reaction | One fixed scale Q/N_P for all descendants |
| Inner boundary | Re-emits every returning particle thermally at 175 K | None |
| Outer boundary | Deletes particles at r = R | None; R used only for the inside label |
| Steady state | Time average over 10,000-step sampling windows, after convergence | Direct sampling of the stationary age integral |
| Observable | Particle-resolved line-of-sight sampler with per-particle Doppler g and an optical-depth flag | Column density; brightness only through the AMPS g and unit constant |
| Cost | About 7–14 h on 128 cores per I3 case | About 4 s on one core per case |

## 5. Code, inputs, and tests

### 5.1 Location and requirements

- **Model and CLI:** `srcHartley2/visualization/comet_edge_correction.py`.
  Python 3.8 or later and NumPy. Nothing to compile.
- **Plot option:** `srcHartley2/visualization/plot_i3_brightness_comparison.py
  --edge-correction`, which also needs matplotlib.
- **Tests:** `srcHartley2/tests/test_comet_edge_correction.py`, run in the
  `i3` group of `srcHartley2/tests/run_all_stage_tests.sh`.

### 5.2 How the parameters reach the model

| Model parameter | Source | Key or column |
| --- | --- | --- |
| r_h | `RUN/i3.in`, case section | `heliocentric_distance_au` |
| k_c | `RUN/i3.in`, case section | `chemistry.rate.<ID>_1au_s-1`, divided by r_h^2 |
| Channel parent, products, E | `RUN/i3.conf` (else `input/i3.conf`) | `ChemistryReaction=ID \| parent \| products \| ... \| energy_eV` (last field) |
| Source species and Q | `RUN/i3.in` | every `source.<S>.rate_s-1` > 0 |
| R_n | `RUN/i3.in` | `nucleus.radius_m` |
| R | `RUN/i3.in` | `mesh.outer_radius_m` |
| v(r) of each source | AMPS output `RUN/PT/plots/pic.<S>.s=<i>.out=<K>.dat` | `X`, `Vx [m/s]`; duplicates averaged; only v > 0 used |
| Product masses | Module table `MASS_AMU` | Not from AMPS |
| Impact parameters b | AMPS brightness table, or user supplied | `Impact Parameter [m]` |
| eta (absolute comparisons only) | AMPS `RUN/PT/diagnostics/comet_stage0.out=<K>.csv` | `species_flux` row for OH |
| g_AMPS, K (brightness only) | AMPS brightness table | `Effective Photon g Factor ...`, `Brightness ... [kR]` / `Column Density ...` |

**Overlay rule.** Common settings in `#begin application: comet` are
overlaid by the selected `#begin case: N` section.

**Unused inputs.** Collision, mesh-resolution, cooling, brightness, and
time-step settings. Sources without H/OH chains (CO, CO2) are skipped by the
relevance analysis.

**Validation.** A missing case section, a reaction without a rate key, or
a species without a mass is an error.

### 5.3 The two implementations

| | Tracked: `comet_edge_correction.py` | Campaign: `validation_runs/i3_campaign_20261108/scripts/i3_collisionless_mc.py` |
| --- | --- | --- |
| Status | Tracked (22e8d7f), tested | Untracked run-archive script |
| Network | Generic: read from `i3.conf`, chains found by fixed-point relevance | Hard-coded I3 channels |
| R_n, R | From `i3.in` | Constants 2,500 m and `R_DOMAIN` (1e9; set to 3e9 by the campaign-2 scripts) |
| Seed | `default_rng(20261109)` | `default_rng(1 + case)` |
| Extra outputs | none | `births` (for eta), `energy_override` (energy sensitivity) |
| Used by | `i3c2_model_vs_obs_edge.png`, `--edge-correction` plots, `EDGE_CORRECTION.md` examples | `i3c2_all_cases_latest_montage.png`, "AMPS/collisionless" tables (`i3c2_mc_compare.py`), campaign-1 Q scaling |

The two agree within 0.2% in median C and within 0.7% at any I3
observation point.

### 5.4 Running it

Factor table for one AMPS brightness file:

```sh
python3 srcHartley2/visualization/comet_edge_correction.py \
  validation_runs/i3_campaign2_20261108/case5_np128/production/PT/plots/pic.H.brightness.s=0.out=6.dat \
  --case 5 --parents 3000000 --seed 20261109 --csv case5_H_edge.csv
```

- The run directory is three levels above the table; override with
  `--run-dir`, `--runtime-input`, or `--build-input`.
- The output index K is taken from the file name, and the H2O profile with
  the same K is used.

Python API:

```python
import sys; sys.path.insert(0, "srcHartley2/visualization")
import comet_edge_correction as E
inputs = E.load_run_inputs(run / "i3.in", run / "i3.conf", 5)
r, n, v = E.read_radial_profile(E.find_profile(run / "PT/plots", "H2O", 6))
cols = E.reference_columns(inputs, {"H2O": (r, v)}, b, ("H", "OH"), parents=3_000_000, seed=1)
N_full, N_trunc = cols["H"]
C = E.correction_factor(N_full, N_trunc, b, inputs.outer_radius_m)
```

`reference_columns` accepts any (radius, speed) arrays for each source. A
synthetic profile, for example a constant 800 m/s, runs the model without
AMPS output. The unit tests do exactly that.

### 5.5 Tests

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest srcHartley2/tests/test_comet_edge_correction.py
```

| Test | Checks |
| --- | --- |
| `ParsingTest.test_products` | Stoichiometry expansion; malformed terms rejected |
| `ParsingTest.test_case_overlay_and_inputs` | Case overlay, rates / r_h^2, sources, channel products; a missing case is rejected |
| `test_recoil_conserves_momentum_and_energy` | Momentum below 1e-30, energy to rtol 1e-12 |
| `test_shell_projection_of_inverse_square_density` | Exact projection against the analytic 1/r^2 column, within 2e-4 |
| `test_factor_is_undefined_at_and_beyond_the_edge` | NaN at and beyond R |
| `test_inverse_square_limit_and_short_lived_species` | End to end, synthetic run: C matches (pi/2)/arctan(...) within 3% (achieved 0.6%); 1 <= C_OH < C_H |
| `test_large_edge_needs_no_correction` | Short lifetimes give C_OH = 1 within 1e-6 and C_H within 1% |

The tests synthesize their own `i3.in`, `i3.conf`, and profile in a
temporary directory. They need neither an AMPS executable nor `validation_runs/`.

## 6. Normalization: gross against net production

The reference starts from the input Q. In AMPS, 19–24% of the injected water
returns to the surface, so its net OH production is lower.

- **Edge factor C.** The normalization cancels, so nothing is needed.
- **Absolute comparisons.** `i3c2_mc_compare.py` scales the reference by

  ```
  eta = (AMPS OH production rate) / (reference OH births inside R * Q / N_P)
  ```

  This gives 0.76–0.81 in campaign 2. With eta applied, AMPS OH agrees with
  the reference within 2% in every case.

## 7. Interpreting AMPS against the reference

Campaign 2, at the observation points:

| Case | Q(H2O) (s^-1) | OH: AMPS / reference | H: AMPS / reference |
| ---: | ---: | ---: | ---: |
| 1 | 3.2e28 | 1.016 | 1.010 |
| 2 | 6.7e28 | 1.006 | 1.019 |
| 4 | 4.5e29 | 0.998 | 1.077 |
| 5 | 3.5e29 | 0.988 | 1.066 |
| 6 | 4.0e29 | 0.999 | 1.057 |

- **OH.** Its transport is collisionless in practice at these distances.
- **H.** The 1–8% excess grows with Q. It is the collisional slowing and
  confinement of H in the inner coma, which the reference omits by
  construction.

## 8. Limitations

- **No product collisions.** Use the reference for ratios, or as the
  collisionless limit.
- **No parent thermal spread.** A few hundred m/s, against product recoils of
  1–42 km/s.
- **No radiation pressure or gravity.** This affects only the far part of
  N_full: 2–3% displacement of 27-km/s H at 1e10 m, about 30% at 1e11 m. See
  `EDGE_CORRECTION.md`, Section 4.
- **Spherical symmetry only.** For 2D or 3D, see `EDGE_CORRECTION.md`,
  Section 8.
- **Steady state.**
- **Nucleus neither absorbs nor re-emits daughters.**
- **Shared assumptions with AMPS.** Excess energy entirely translational;
  OH rates with no velocity dependence.
- **Brightness conversion.** It reuses the AMPS g on each ray.

## 9. Extending it

- **Time-dependent production.** Carry T, the time since the water left the
  nucleus: the parent transit, the integral of dr/v from R_n to r_d, plus
  every decay age, plus the snapshot age. Multiply each weight by
  Q(t - T)/Q(t). See
  `validation_runs/i3_campaign2_20261108/analysis_case1_h/production_history_test.py`.
- **Extended (icy-grain) sources.** Add parents released at
  r_0 = R_n + Exp(L) with isotropic grain-frame velocities. See
  `analysis_case1_h/icy_grain_test.py` in the same directory.
- **New parents or products** (CH4, CH3OH, ...). Add `source.<S>.rate_s-1`,
  `chemistry.rate.*`, the `ChemistryReaction` lines, and masses in
  `MASS_AMU`, and supply a speed profile for each new source.
- **Radiation pressure, non-spherical sources, other domains.** See
  `EDGE_CORRECTION.md`, Section 8.1.

## 10. Checklist for running and changing the model

- **Outputs must match.** Use the same output index K for the parent profile
  and the brightness table. `correction_for_run` does this automatically.
- **Use the frozen copies.** Use the `i3.in` and `i3.conf` frozen in the run
  directory, not the current `input/` files, which may have changed since
  the run.
- **Masses** come from `MASS_AMU`. A new species needs an entry.
- **Normalization.** C needs none; absolute comparisons need eta (Section 6).
- **Reproducibility.** It needs the seed, the parent count, and the NumPy
  version (PCG64 stream). State them with every published number.
- **Domain.** R is read from `mesh.outer_radius_m`. C is undefined at b >= R
  and is shown only for b <= 2R/3.
- **Two implementations exist.** Some campaign-2 tables come from the
  untracked one (Section 5.3). Both agree within 0.7%.
- **AMPS parent profile.** The parent speed profile is the only AMPS output
  the factor needs. For a standalone run, pass any (radius, speed) arrays to
  `reference_columns`. The CLI currently reads an AMPS run directory only.
- **Changing physics.** Physics changes (energy partition, rates, new
  channels) go into `i3.conf`/`i3.in` and so reach both AMPS and the
  reference. Changes to the reference alone (forces, extended sources,
  time dependence) must be verified against the analytic limit (Section 5.5)
  and documented as reference-only.
- **Sensitivity to the parent profile.** The edge factor is insensitive to
  it. For case 5, a constant 800 m/s parent speed gives C = 1.036 and 1.436
  at b = 1.8e8 and 1.44e9 m, against 1.039 and 1.437 with the AMPS profile.
  Absolute reference columns near the parent scale length (about 1e8 m) are
  sensitive to it.
- **Converged AMPS output.** The reference itself has no convergence, but
  comparisons with AMPS must use converged AMPS outputs, or say that they
  do not.
- **Brightness.** It always uses the AMPS g on each ray (Section 3.13). A
  reference brightness is therefore not independent of AMPS's velocity
  distribution.
- **Optically thin.** Columns are thin by construction. AMPS flags H rays
  with tau > 0.1, which in campaign 2 includes some inner observation points
  of cases 2–6. The reference does not change that.
- **Profile filtering.** Only AMPS profile nodes with Vx > 0 are used. Speed
  is interpolated linearly in r and held constant beyond the last node; it
  is not extrapolated.
- **Observed species.** The default is H and OH; pass `observed=(...)` to
  change it. A source species cannot also be observed (rejected).
  Intermediates are found automatically from the reaction network.
- **Hard limits.** The age range [1 s, 1e9 s] and R_big = 1e12 m are module
  constants (`AGE_MIN_S`, `AGE_MAX_S`, `R_BIG_M`). Raise them for species
  with lifetimes above about 1e8 s or for impact parameters approaching
  1e10 m.
- **Edge-corrected plot range.** The corrected curve is drawn only to
  b <= 2R/3, where C stays finite and well sampled.
- **Error bars.** None are produced. Estimate the Monte Carlo error by
  repeating with another seed (Section 3.11).
