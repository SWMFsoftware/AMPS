# 2013-04-11 real-event SEP test in 3-D (srcSEP3D)

This directory holds the 3-D version of the 2013-04-11 real-event test of
[`srcSEP/test/event_2013_04_11`](../../../srcSEP/test/event_2013_04_11/README.md).
The physical event is the same one:
- the shock and ambient are the shared reduced core (`event.event`, schema v1.3);
- the observations are SOHO/ERNE, ACE/EPAM and GOES-13 (`observations.csv`);
- the transport law and the source band are those of the 1-D test.

In 3-D, particles move through a 3-D mesh that covers the magnetic flux tube
connecting the observers to the shock, in the frame that corotates with the Sun.

The same three transport models are tested:

| srcSEP3D `transport` | 1-D `--particle-mover` | Equation |
|---|---|---|
| `parker` | `parker` | Parker transport, κ∥ = vλ∥/3 |
| `focused-diffusion` with `pitch_angle_diffusion_model = pitch-angle-library`, model `PA-ISO` | `fte-dmumu` | focused transport, D_μμ = v(1 − μ²)/(2λ∥) |
| `focused-scattering`, `focused_scattering_frame = alfven-wave-frame-isotropic` | `fte-mfp` | focused transport with isotropic scattering events at rate v/λ∥, split equally between the two Alfvén-wave frames (`sep_focused_transport_mfp_core.cpp` does the same for a prescribed λ∥) |

All three use λ∥ = 0.3 au (r/1 au)(R/1 GV)^(1/3), from Liu et al. (2025).

## 1. Decisions

User decisions of 2026-11-25 and 2026-11-26. `event3d_settings.json` records
each with its date.

| Topic | Decision |
|---|---|
| Shock model | The 1-D v1.3 heliocentric sphere of the shared core. It is consumed through the provider-neutral `ShockSurface` (`adapters/shock_surface.h`), so another core model or an MHD-derived front can be added later ("1: OK … able to get other shock models"). |
| Frame | Corotating with the Sun ("1: yes"). The domain is fixed along the field lines, which turn with the Sun (section 2.1). |
| Domain | The observer flux tube: the tube that the collection spheres of SOHO, ACE and Earth meet during their stay intervals ("2"). |
| Body observers | Ephemeris table or SPICE, behind an extensible observer interface ("2"). |
| Stay interval | Long tier: the whole run. Nightly tier: the first hour ("3: yes"). |
| Collection radius | Set by a convergence study ("4"). The value used now is provisional (section 6). |
| Source | The 1-D net source on every tube line, with the 1-D scale factor reused ("4: yes"; section 2.4). |
| Tiers | Nightly: first hour, coarse, all three movers, DiffNum against AMPS references. Long: full event ("5: yes"). |
| Inertial forces | The Coriolis/centrifugal option of the focused movers exists for other applications, and is off here so the equations match the 1-D movers (section 2.6). |

## 2. Model

### 2.1 Corotating frame (`core/computation_frame.h`)

`mesh.domain.coordinate_frame = HCI-corotating` selects the frame that rotates
about HCI +z at `frame_rotation_rate_rad_per_s`. That rate must equal the
event's `ambient.rotation_rate_rad_per_s` (2.86533×10⁻⁶ rad/s). The frame
equals HCI at `frame_reference_epoch_s` (0).

**Why the domain must turn with the Sun.** In HCI the solar wind is radial, so
part of its velocity is perpendicular to the spiral field. A guiding centre is
carried by that perpendicular flow at exactly the rate the spiral pattern
turns. SEPs therefore stay on their field line, and the line turns with the
Sun. AMPS cannot move its mesh, so the run is solved in the frame where the
lines are fixed:
- the reduced ambient, a rigidly rotating PFSS/Parker pattern, is static there;
- U′ = U − Ω×r is parallel to B;
- only the observers move.

The adapter rotates positions into HCI for every provider query, and rotates
B, U − Ω×r and their gradients back. The transport coefficients do not change,
because ∇·(Ω×r) = 0 and bb:∇(Ω×r) = 0. The 1-D line already turns rigidly with
the Sun (1-D README §1.3), so the 3-D tube holds the same field lines.

Because the ambient is static in this frame, the published background stays
valid to the end of the event. The decks therefore never republish it.

### 2.2 Shock

`ReducedFrontShockSurface` presents the core's front through `ShockSurface`:
- a level function;
- segment crossings;
- the local front state (support, admission, normal, speeds, upstream state,
  compression, θ_Bn);
- the release surface.

It supports the heliocentric sphere and the finite SSE cap. The mover receives
the exact radius and speed of the front at every step start.

### 2.3 Observer flux tube (`mesh/observer_flux_tube.*`)

`mesh.active_region.mode = observer-flux-tube` builds the domain in four steps:
1. Seed points lie on a latitude–longitude lattice with spacing `seed_spacing_m`
   on the reference sphere through the first tube observer.
2. Each seed's field line is traced through the event's own field (RK4 with
   arc step ≤ min(`trace_step_m`, 0.02 r)) from the inner to the outer sphere.
3. A line is kept if, at any of the observer's sampled times in its stay
   interval, the line's point at the observer's heliocentric distance lies
   within ρ + s of the observer's position. Here ρ is the collection radius and
   s is the local line spacing. The separation is measured at equal radius, not
   as the closest approach. For an offset in the plane of the line and the
   radial direction, it exceeds the closest approach by up to 1/cos ψ (ψ is the
   line's angle to the radial direction, about 45° for the Parker spiral at
   1 au). The tube is therefore narrower in that plane than ρ + s around the
   line. After the blocks are allocated, a coverage check stops the run if the
   centre or any of the six axis extremes of a sampled collection sphere falls
   outside the active blocks (`VerifyObserverTubeCoverage` in `main_lib.cpp`).
4. Each kept line is a flux element with Φ = |B_r| × (cell area). Its local
   spacing is √(Φ/|B|).

The mesh uses these lines twice:
- **mask:** the capsules around the line chords (radius = local spacing plus
  half the chord) mark active blocks; the usual buffer layers, cavity filling
  and connectivity checks follow;
- **refinement:** the requested cell size is the radial target inside the tube
  and at most twice the distance to the tube outside it.

After allocation the run checks that every observer's collection sphere stays
in active blocks at every sampled time.

How it differs from `parker-tube`:

| | `parker-tube` | `observer-flux-tube` |
|---|---|---|
| Centre | One analytic constant-speed Parker spiral | Field lines traced through the run's own field (PFSS + accelerating wind) |
| Defined by | Start longitude/colatitude and length | Observers and their stay intervals |
| Cross-section | Cone or constant radius | The real flux tube (area ∝ 1/\|B\|) |
| Width along the observers' motion | None | The union of tubes met during the stay |
| Frame | HCI | Corotating |

**Inner sphere at 2.5 R☉.** The 3-D transport inner sphere is the PFSS source
surface, not the 1-D line base at 1.1 R☉. This changes no transport, because
no particle can reach r < 2.52 R☉: the front is at R_sh ≥ 2.52 R☉ from t = 0
and absorbs every returning particle (section 2.5).

The change is needed because, in 3-D, the tube's halo blocks below 2.5 R☉
would reach open PFSS tubes where the isothermal open-tube wind has no subsonic
solution. The core rejects such points and does not clamp them.

### 2.4 Source (`adapters/tube_line_source.*`)

`source.shock.placement = release-surface-on-tube-lines` applies the 1-D srcSEP
net source to each tube line l:

    R_l = η n₁ w₁ Φ_l / |B₁·n|

Here w₁ = V_n − U₁·n is the shock-frame inflow, n₁ the upstream proton density
and B₁ the upstream field (`NetProtonRate` of the 1-D source).

For each particle:
- **release point:** where the line crosses the release sphere
  (1 + f) R_sh, with f = 0.3;
- **position:** spread over the element's footprint on that sphere;
- **spectrum:** the local DSA spectrum dN/dp ∝ p^(2−q), with q = 3r/(r − 1),
  between 20 keV and 200 MeV;
- **weights:** log-uniform importance sampling with the exact weight;
- **pitch angle:** the positive-normal-flux law for the focused movers, and
  isotropic for Parker.

Birth times are stratified, with one draw per dt/N. At each draw one line is
chosen with probability R_l/ΣR_l, so the particle represents w(ΣR) dt/N
protons. The source therefore estimates the 1-D sum over the tube's flux
elements.

Check: in the first 60 s the 3-D tube injects 4.05×10²¹ protons per weber.
The 1-D line injects 4.26×10²¹ per weber (3.34×10¹⁹ over its π m² × |B(1.1 R☉)|
flux tube).

### 2.5 Absorbing shock (`adapters/shock_return.*`)

`source.shock.returning_particles = absorb`: a particle that crosses the front
into its downstream side, where the front is an accepted shock, is removed. The
net source already counts only escaping particles (1-D README §2.2).

The test runs over each substep, with the front moving linearly:
- **ends:** a downstream start or end counts as a crossing;
- **Parker:** the 1-D Brownian bridge, with crossing probability
  exp(−g₀g₁/(κ_r dt)) and κ_r = κ∥(b·e_r)²;
- **focused movers:** the exact dip of the straight chord below the front.

The returned fraction follows the 1-D Parker ledger: 14 % at 600 s and 66 % at
one hour in 1-D.

### 2.6 Inertial forces of the focused movers (option)

`application.sep.transport.frame_inertial_forces = none | coriolis |
coriolis-and-centrifugal` is valid only in the corotating frame. It adds the
plasma-frame acceleration terms of the focused transport equation for the
frame acceleration a = 2Ω×U′ (+ Ω×(Ω×r)):

    dμ/dt   += −(1 − μ²)(b·a)/v
    d ln p/dt += −μ (b·a)/v

Together these are a uniform force −γm(b·a) along b that leaves p⊥ unchanged
(unit test E3D02). The option also adds the guiding-centre drift (F×B)/(qB²) of
the same inertial force acting on the particle.

The 2013-04-11 decks use `none`. With U′ ∥ B the Coriolis term vanishes along b,
and the 1-D movers have no inertial terms.

### 2.7 Long jumps

A 60 s step carries a fast particle across many blocks. Two mechanisms keep the
3-D movers independent of the MPI decomposition:
- **provider fallback:** a cell this rank does not hold (another rank's leaf,
  or a ghost-block centre AMPS never sends) is evaluated at its centre from the
  deterministic provider, giving exactly the value its owner published;
- **hand-off:** a particle that ends a step on another rank's active leaf is
  sent to its owner after the mover phase (`ExchangeLongJumpParticles`) instead
  of being deleted.

Particles that leave the active domain sideways are counted and reported, never
dropped silently.

## 3. Body observers (`output/observer_trajectory.*`, `output/spice_ephemeris.*`)

`kind = body` with `ephemeris = table` reads a checksummed CSV:
- columns `time_s,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s`, in the model
  HCI frame;
- cubic Hermite interpolation, with no extrapolation.

The tables in `ephemeris/` come from the pinned DE441 Horizons vectors of the
raw archive (`/data/vtenishe/sep_validation_data/event_2013_04_11/horizons`). The
model frame has +z on the IAU solar pole (Archinal et al. 2018) and +x on the
projection of the geometric Sun→Earth vector at t₀.

Check: Earth's heliographic latitude from this construction differs from the
1-D frame's Horizons sub-observer latitude by 1×10⁻⁵°.

`ephemeris = spice` asks CSPICE (`spkezr_c`, no aberration correction) for the
state of `spice_target` relative to `spice_center` in `spice_frame`, at
et = str2et(`spice_reference_utc`) + t. It rotates the result into the model
frame with `spice_to_model_rotation`. The source compiles against `SpiceUsr.h`
only when the AMPS configuration names a SPICE toolkit; a `nospice` build stops
with a typed diagnosis.

In a corotating run:
- `fixed-cartesian` is fixed in the frame, i.e. on a field line;
- `moving-cartesian` and `body` follow inertial trajectories.

**Observer units.** The `represented_particles_J-1` column of `spacecraft.csv`
is Σw/ΔE: the represented protons inside the collection sphere per joule of
kinetic energy. The srcSEP3D sampler passes the `normalization` label through
and does not apply it. The differential intensity is that value × v/(4πV),
with V = (4/3)πρ³.

**Cross-check against 1-D (nightly, t = 1 h, 2026-11-26).** Two checks:

- **Spectral shape.** Within each radius window, the fraction of particles in
  each band agrees with the 1-D energy–radius sample for every mover. Near the
  shock (0.02–0.06 au), 0.82–0.94 of the particles are in 0.02–0.2 MeV in
  1-D and 0.86–0.95 in 3-D. At 0.10–0.20 au the bands agree within a few
  percent. The 1-D table holds the weight per bin divided by Σ N ΔE, so its
  band fraction is a sum of values, not of values × ΔE.
- **Proton density.** The 3-D density is the radial–energy histogram divided
  by the tube volume ΣΦ/|B| ds in each radial bin, compared with the 1-D inner
  observers:

| r | Band | Parker | FTE-Dmumu | FTE-MFP | 3-D macroparticles |
|---|---|---|---|---|---|
| 0.1 au | 0.2–2 MeV | 0.31 | 0.57 | 0.73 | 15–22 |
| 0.1 au | 2–20 MeV | 0.44 | 0.37 | 0.48 | 77–122 |
| 0.1 au | 20–200 MeV | 1.8 | 1.5 | 2.1 | 43–63 |
| 0.3 au | 2–20 MeV | 0.94 | 0.41 | 1.6 | 10–17 |
| 0.3 au | 20–200 MeV | 1.6 | 0.93 | 1.6 | 59–105 |
| 0.5 au | 20–200 MeV | 1.7 | 1.45 | 1.26 | 51–79 |

The table gives 3-D / 1-D. The densities span four decades, and the ratios
stay within 0.3–2.1. The 1-D nightly injects 20 samples per step, and its
20–200 MeV bins at 0.1 au change by up to a factor of 5 between neighbours.
Each code uses one source realization for all its movers, so statistical
errors are shared across a row. The quantitative comparison belongs to the
long tier.

**Sampling checks (FTE-MFP, nightly).**
- A collection sphere reproduces the tube density at its radius: the
  time-summed sphere/tube density ratio is 0.96–1.22 at 0.1–0.5 au.
- The tube density does not depend on the tube width. Halving ρ (20 instead
  of 39 lines, 0.513 of the flux, source rate 0.513) changes the 0.1 au
  densities by 0.6–1.7 at 23–106 macroparticles.
- Two runs of each mover on 8 ranks give byte-identical `observers`,
  `radial_energy` and `ledger` products.

## 4. Inputs

`srcSEP/test/event_2013_04_11/prepare_event.py deck` writes every file here
except `event3d_settings.json`, `README.md` and the scripts. Its `check`
regenerates them and byte-compares.

| File | Contents |
|---|---|
| `event3d.<tier>.<mover>.in` | Run deck: `#include`s the tier part and the mover part |
| `event3d.<tier>.in` | Run, mesh (corotating frame, observer tube), background, source, observers, output |
| `event3d.convergence-r<tag>au.in` | The long tier with collection radius `<tag>` (0p005, 0p01, 0p02, 0p04 au) for the study of section 6; run as `event3d.convergence-r<tag>au.<mover>.in` |
| `event3d.mover-<mover>.in` | `application: sep` (transport and coefficients), plus the PA-ISO block for fte-dmumu |
| `ephemeris/{soho,ace,earth}.csv` | Body ephemerides in the model frame |
| `deck3d_summary.json` | Frame rows, latitude check, ephemeris digests, inner-observer positions |
| `event3d_settings.json` | 3-D settings and decisions |

The run directory also needs the shared `event.event` and
`assets/wso_cr2135_model_frame.harmonics.csv` from the 1-D directory; the Table
copies them there.

## 5. Tests

| Test | What |
|---|---|
| `run_unit_tests.sh` | E3D01–E3D06, standalone under ASan/UBSan: frame, inertial terms, ephemeris, tracing and the tube, absorbing shock and bridge statistics, tube-line source rate, weights, spectrum and pitch law |
| `run_spice_test.sh CSPICE_DIR LSK SPK` | E3D07: the SPICE body-observer path against the Horizons Earth table. On passion (CSPICE N0067, `naif0012.tls`, `de430.bsp`) Earth differs by at most 0.14 km and 3×10⁻⁵ m/s over the 67 samples. SPICE kernels for SOHO and ACE are not available here, so those bodies use the table. |
| `SEP3D--event-2013-04-11--{Parker,FTE-Dmumu,FTE-MFP}` | Nightly Table tests: first hour, DiffNum of `observers`, `radial_energy` and `ledger` (`regression_products.sh`) against `output/test_SEP3D/SEP3D--event-2013-04-11--<Mover>--<label>.ref_np<N>.gz` |
| `SEP3D--event-2013-04-11-long--{FTE-MFP,Parker,FTE-Dmumu}` | Long tier: full event, `Ranks=32` |
| `run_convergence.sh AMPS MOVER WORK_DIR [RADIUS_TAG ...]` | The collection-radius study (section 6): one long-tier run of `event3d.convergence-<tag>.<mover>.in` per radius (`MPIRUN`, default `mpirun -np 32`), then `convergence_report.py`. Finished runs are reused. |
| `convergence_report.py RADIUS_AU=RUN_DIR ...` | Pairwise comparison of consecutive radii and the selected radius (section 6) |

`make_references.sh RUN_DIR TEST_NAME REFERENCE_DIR` writes the references
from a finished run.

## 6. Provisional numerics and the collection-radius study

These values are provisional until the convergence study:

| Setting | Value |
|---|---|
| Collection radius | 0.02 au |
| Cell size | grows linearly from 10⁸ m at the Sun to 0.02 au at 1 au (2–6 % of r) |
| Seed spacing | half the cell size at the reference observer |
| Buffer blocks | 1 |
| Time step | 60 s, as in 1-D |

The study runs the long tier for collection radius 0.005, 0.01, 0.02 and
0.04 au (`event3d_settings.json` `convergence`). It keeps the largest radius
whose Earth/L1 time profiles agree with those of the next smaller radius within
their Monte Carlo uncertainty.

- **Decks:** `event3d.convergence-r{0p005,0p01,0p02,0p04}au.<mover>.in`, from
  `prepare_event.py deck`. Each is the long tier with that collection radius
  on every observer and the output prefix `event3d-convergence-r<tag>au`. The
  radius also sets the tube width (ρ + s around each observer), so every
  radius builds its own mesh.
- **Run:** `run_convergence.sh <amps> fte-mfp <work-dir>`. Each run is one
  long-tier job (about 12 h on 32 ranks, as `SEP3D--event-2013-04-11-long--*`).
  The collection radius is a sampling parameter, so one mover suffices; pass
  another mover to check that.
- **Profiles:** the tube observers `tube-earth`, `tube-soho` and `tube-ace`,
  in the bands 0.02–0.2, 0.2–2, 2–20 and 20–200 MeV, at every publication
  (10 min), as intensity j = Y v/(4πV) with V = (4/3)πρ³.
- **Agreement of a pair ρ_s < ρ_l:**
  - cells are never chosen by their own value (that selects upward
    fluctuations). An (observer, band) enters when both runs hold at least 10
    effective macroparticles in it over the run, N_eff = (Σj)²/Σσ²;
  - in an entered band, each publication where either run is nonzero is a
    cell, with z = (j_l − j_s)/√(σ_l² + σ_s²);
  - the pair agrees when Σz² does not exceed the 3σ upper quantile of χ²_ν.
  An empty sphere has σ = 0, and consecutive publications share particles.
  Both make the test reject too often rather than too rarely. `--min-count`
  and `--sigma` change the two thresholds.
- **Selection:** the largest radius that agrees with the next smaller one.

**Smoke check (2026-11-26).** The report ran on nightly-length fte-mfp runs
at 0.01 and 0.02 au, using the inner observers. Rescaling their radii by r/r⊕
cancels in z. No band reaches 10 effective macroparticles in both runs in one
hour (the thinner spheres hold at most 5–7 macroparticles), so the report
returns "not decided" (exit 1). The study needs the long tier.

## 7. Limitations and pending work

- The heliocentric sphere has no longitudinal structure. Inside the tube
  (a few degrees wide) this is the same local approximation as in 1-D.
- The background is cell-constant, so a particle follows b̂ of its cell. This
  costs accuracy at order h/r per cell.
- Long-tier references and the 3-D observation score (the 1-D
  `score_observations.pl` reads the 1-D outputs) come after the convergence
  study.
- The long-tier and convergence run times are estimates (`Time=720`, about
  12 h on 32 ranks each), not measurements. The nightly tier takes 4–13 min
  on 8 ranks for 60 steps; the long tier has 2484 steps and a growing
  population.
- Memory: AMPS replicates the whole tree on every rank. With the coarsened
  observer-tube mesh (7883 leaves, 3023 active blocks) a nightly run holds
  about 0.6 GB per rank on 8 ranks, flat over the 60 steps (measured
  2026-11-26).
