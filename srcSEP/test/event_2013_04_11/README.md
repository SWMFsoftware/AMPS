# 2013-04-11 SEP event: real-event end-to-end validation (srcSEP, reduced shock-front path)

The 3-D version of this test (srcSEP3D, the same shock, ambient and
observations in the corotating frame) is [`srcSEP3D/test/event_2013_04_11`](../../../srcSEP3D/test/event_2013_04_11/README.md).

This directory holds the inputs of the 2013-04-11 real-event tests, in two
tiers (Section 5): the nightly regression
`SEP--event-2013-04-11--{Parker,FTE-Dmumu,FTE-MFP}` (`event.nightly.in`) and
the full-event validation `SEP--event-2013-04-11--{Parker,FTE-Dmumu,FTE-MFP}--long`
(`event.in`). One field line is modelled: the line that passes through Earth
at t0. Protons are transported along it by three srcSEP movers
(`--particle-mover`):
- `parker`: the Parker equation;
- `fte-dmumu`: the focused transport equation with isotropic D_μμ;
- `fte-mfp`: the focused transport equation with mean-free-path scattering.

The shock, the solar wind and the interplanetary magnetic field come from the
shared reduced core (`src/models/sep_corona_swcme/shock_front`, schema v1.3).
Real observations are used both as input (shock heights, L1 plasma, the WSO
magnetogram) and as the comparison (SOHO/ERNE, ACE/EPAM, GOES-13).

Python (`prepare_event.py`) is used only to turn the archived observations
into the checked-in input files. The tests themselves run AMPS and Perl only.

| File | Role |
|---|---|
| `raw_manifest.json` | Every archived input with its SHA-256 and source URL. The archives are kept outside the source tree; on passion in `/data/vtenishe/sep_validation_data/event_2013_04_11` (`SHA256SUMS` there) |
| `event_settings.json` | Every physical decision and numerical choice, in one place |
| `prepare_event.py derive` | raw → `derived.json`, `observations.csv` |
| `prepare_event.py deck` | → `event.event` (core), `event.in` (srcSEP, long tier), `event.nightly.in` (srcSEP, nightly tier), `event.args`, `assets/wso_cr2135_model_frame.harmonics.csv`, `deck_summary.json` |
| `event.args` | One line of srcSEP command-line options shared by every mover and both tiers: the λ∥ law (`transport` in `event_settings.json`), log-uniform source sampling, the coefficient providers, and the provisional shock-boundary numerics. The Table recipes pass it with `$$(cat event.args)` |
| `prepare_event.py check` | Regenerates every product and byte-compares it with the checked-in file |
| `score_observations.pl` | Observation score (reported, not gating) |
| `regression_products.sh` | Assembles the DiffNum products of one run |
| `make_references.sh` | Makes a Table test's references from a finished run of that test (`gzip -n` of the products, named for the run's rank count) |

## Decisions

Every physical choice was made by the user and is recorded in
`event_settings.json` under `decisions`. The two decisions documented in detail
here were made on 2026-11-25:
- the PFSS source surface is at 2.5 R☉, and the field line is transported
  through the open PFSS corona between 1.1 and 2.5 R☉;
- the release distance is an `.in` parameter, L_ref = f R_sh(t) with
  f = 0.3 (Section 2);
- the source energy band is 20 keV – 200 MeV (Section 3);
- the tests come in two tiers: a short coarse nightly run of all three movers,
  with inner regression observers, and the full event as a separate long test
  run outside the nightly (Section 5).

## 1. Magnetic field and field-line geometry

### 1.1 Model

| Region | B | U | ρ |
|---|---|---|---|
| 1 R☉ ≤ r ≤ 2.5 R☉ (PFSS shell, open tubes) | Potential field of the WSO CR2135 radial-model coefficients, rotated into the model frame. Corotates rigidly | Ω×r + s u b̂: corotation plus flow along B | ρ = K\|B\|/u |
| r ≥ 2.5 R☉ (Parker exterior) | B_r(R_ss)(R_ss/r)²(r̂ − Ω(r − R_ss²/r)/U φ̂) | U(r) r̂ + Ω R_ss²/r φ̂ | ρ = F_m/(U r²) |

- U(r) is the isothermal Parker (1958) transonic wind. Its sound speed
  c_w = 106.4 km/s is fitted so that U(L1) equals the OMNI speed, 377.9 km/s,
  which puts the critical point at r_c = 8.43 R☉.
- F_m is the mass flux per steradian, fixed by the OMNI density at L1.
- One factor scales the magnetogram so that |B| at L1 at t0 equals the OMNI
  value. That factor is 17.48.
- The field line is the core field line through Earth at t0. srcSEP traces it
  at start-up, from Earth sunward to the 1.1 R☉ base
  (`mesh.field_line.anchor_{x,y,z}_m`, Section 1.5), so Python never
  evaluates the magnetogram.

So the line is **not** radial below 2.5 R☉. It follows the open PFSS field
line, which bends in latitude and longitude, and it becomes radial exactly at
2.5 R☉, where the PFSS boundary condition forces it. Above 2.5 R☉ it is a
Parker spiral that starts radial (B_φ = 0 at R_ss) and winds up to the 1 AU
spiral angle.

### 1.2 Why the source surface is at 2.5 R☉ and not at the 1.1 R☉ line base

The earlier choice put the source surface at the 1.1 R☉ base (decision B1),
with the kinematic Parker spiral starting there. With the slow isothermal
wind near the Sun that winds the line through 306° inside 2 R☉ (B_φ/B_r = 20 at
1.2 R☉) and through 414° in total.
- That is unphysical. Below the Alfvén surface (about 10–20 R☉; Kasper et al.
  2021) the magnetic field dominates and the plasma corotates with it. The
  field is close to potential there and is not wound up by the wind (Weber &
  Davis 1967).
- The PFSS model with R_ss = 2.5 R☉ is the standard description of that region
  (Altschuler & Newkirk 1969; Schatten, Wilcox & Ness 1969; Hoeksema 1984;
  Wang & Sheeley 1992). It is the value for which WSO publishes its source-
  surface products. WSO's notes say the source surface is "usually at 2.5R"
  (`raw/wso/pfss.txt`, X. Sun, *Notes on PFSS Extrapolation*).
- The radial-model coefficients used here are listed by WSO as "Potential
  Field Source Surface Model - Radial (Independent of Rs)"
  (`raw/wso/ghlist.rad.html`). They therefore describe the photospheric field
  for any R_ss.
- **The event paper itself uses this construction.** Lario et al. (2014,
  Sect. 2.1, after Park et al. 2013) find the L1 footpoint with "a ballistic
  projection … through a Parker spiral … up to a distance of 2.5 R⊙ … followed
  by the open field lines resulting from the potential field source surface
  (PFSS) model".

### 1.3 Physics in the open PFSS shell, and why it is a valid material region

srcSEP moves the field-line vertices with the plasma velocity U (a material
line). This is valid only where U and B satisfy ideal induction and the
density satisfies continuity. That requirement is why srcSEP refused the
SWCME1D ambient: Leblanc density with a constant wind does not conserve mass.

In the open shell:

1. **Induction.** The potential field is static in the frame that corotates
   with the Sun. With V = Ω×r, ∇·V = 0 and (B·∇)V = Ω×B, so
   ∇×(V×B) = −(V·∇)B + Ω×B = ∂B/∂t for a rigidly rotating field. The added
   field-aligned flow u b̂ contributes (u b̂)×B = 0. U = Ω×r + s u b̂ therefore
   solves ideal induction exactly, for any field-aligned speed u. A material
   line stays the same field line, rigidly rotated, while its vertices slide
   along it.

2. **Continuity.** Steady state in the corotating frame needs ∇·(ρ u b̂) = 0.
   That holds if and only if ρu/|B| is constant along each tube, since
   ∇·(K B) = B·∇K = 0. The constant, K = F_m/(R_ss²|B_r,ss|), is the
   exterior's mass flux per unit magnetic flux at the tube's source-surface
   foot, so ρ is continuous at R_ss. The radial law ρ = F_m/(U r²) used in the
   exterior would violate continuity on a non-radial tube.

3. **Field-aligned momentum.** For an isothermal flow along a tube of area
   A = Φ/|B| (Kopp & Holzer 1976, the flux-tube form of Parker 1958), the
   Bernoulli integral u²/2 + c_w² ln ρ − GM/r = const, with ρ = K|B|/u, gives

       y − ln y = 4 ln(R_ss/r_c) − 3 + 2 ln(|B_r,ss|/|B|) + 4 r_c/r,   y = u²/c_w²

   on the subsonic branch. The constant is that of the exterior transonic
   wind, so u is continuous at R_ss. Where |B| = |B_r,ss|(R_ss/r)² the relation
   reduces to Parker's radial one. In the deeply subsonic shell the result is
   hydrostatic: ρ hardly depends on the tube area, and u adjusts to the
   expansion.
   - The centrifugal potential is left out, exactly as in the exterior Parker
     wind. Below R_ss it is at most Ω²R_ss³/(2GM) ≈ 2×10⁻⁴ of the
     gravitational term.
   - The closure needs r_c > R_ss, and here r_c = 8.43 R☉. A point without a
     subsonic root is rejected.

4. **Matching at R_ss.** B is radial there in both regions, B_r is the same,
   and U = U r̂ + Ω×r on both sides. So B, U and ρ are continuous.

5. **What is not modelled.**
   - Cross-field momentum balance: PFSS is current-free, which is justified by
     low plasma β. β can approach 1 at about 1.2 R☉ above active regions
     (Gary 2001).
   - The time-dependent eruption below 2.5 R☉.
   - Closed tubes: no material vertex may enter one, and such a query is
     rejected.
   - A downstream volume: none exists in the reduced model.

The core implements this as `ambient.open_wind_density_model =
isothermal-flux-tube`. srcSEP admits it as the material profile
`reduced-pfss-open-parker-v1`: Parker exterior plus `CoronalOpen`, and nothing
else (`sep_material_background.h`).

### 1.4 Numbers for this event (t0 = 2013-04-11T07:26:44Z)

| Quantity | This model | Reference |
|---|---|---|
| Exterior winding, Earth → 2.5 R☉ | 87.6° | Lario et al. 2014 Table 1, ballistic 391 km/s spiral: 61° |
| Source-surface foot (Carrington lon, lat) | 173.0°, −5.9° | — |
| Foot at 1.1 R☉ (Carrington lon, lat) | 155.2°, −5.0° | Lario et al. 2014 Table 1, PFSS (HMI, 06:04 UT): 147°, +13°; MAS: 147°, +13° |
| n_e, u∥, \|B\| at 1.1 R☉ | 1.4×10¹⁰ cm⁻³, 0.09 km/s, 25 G | Newkirk (1961) n_e(1.1 R☉) = 3.6×10⁸ cm⁻³ |
| n_e, u∥, \|B\| at 2.46 R☉ | 3.0×10⁶ cm⁻³, 6.0 km/s, 0.34 G | Newkirk (1961): 2.4×10⁶ cm⁻³ |
| Magnetogram factor fixing \|B\|(L1) | 17.48 | — |

Notes on the table:
- The Carrington longitudes in Lario et al. Table 1 are their heliographic
  longitudes minus 40°.
- The extra 27° of exterior winding comes from the slow near-Sun acceleration
  of the isothermal wind fitted to 378 km/s (decision A1).
- In longitude, our 1.1 R☉ foot is 8° west of the PFSS and MAS footpoints of
  Lario et al. Their latitude is 18° further north. They used a different
  magnetogram (HMI rather than WSO) and a constant-speed spiral.
- The low-corona density is the hydrostatic isothermal corona at the fitted
  temperature, which carries the observed 1 AU mass flux. It is about 40×
  Newkirk at the base but matches it near 2.5 R☉.
- The density below 2.5 R☉ does not enter the source. The shock sphere starts
  at 2.52 R☉, so every source point is in the exterior.
- A factor of 17.48 is large. It reflects the known PFSS underestimate of open
  flux (Linker et al. 2017) and a footpoint near the polarity boundary.

These numbers come from the core (`REDUCED_ANCHOR` line of the run log and the
diagnostic in the Section 7 test list). `deck_summary.json` records the
exterior part.

### 1.5 The anchored field line

- `mesh.field_line.anchor_{x,y,z}_m` is the Earth position at t0. It replaces
  `initial_{x,y,z}_m`; giving both is an error.
- srcSEP calls `TraceSunwardToRadius`, which traces the core field from the
  anchor sunward to `mesh.domain.inner_radius_m`. Each leg is an arc equal to
  the remaining radial distance (|dr/ds| ≤ 1, so it cannot overshoot) and is
  traced to 1 m. The base is placed on the sphere to the line's coordinate
  tolerance.
- The line is then built outward from that base with the configured length
  and geometric spacing.
- The native constant-wind Parker builder rejects the anchor.

## 2. Release distance L_ref and the release fraction f (user decision 2026-11-25: f = 0.3)

### 2.1 Definition

f is the release distance divided by the shock radius:

  L_ref(t) = f · R_sh(t),

where R_sh(t) is the heliocentric radius of the shock sphere at the birth
epoch t. Particles are released where the field line crosses the Sun-centred
sphere of radius

  r_rel(t) = R_sh(t) + L_ref(t) = (1 + f) · R_sh(t).

At each source sample time (20 stratified samples per 60 s step),
`reduced_shock_particles.cpp` does the following:
1. It finds where the moving field line crosses the shock sphere. The core's
   local oblique MHD Rankine–Hugoniot solution there gives the upstream proton
   density n₁, the shock-frame normal inflow u₁ₙ, the upstream field B and the
   compression ratio r.
2. It sets the net rate η n₁ u₁ₙ Φ/|B·n̂| (`NetProtonRate`). Φ is the tube's
   magnetic flux, so Φ/|B·n̂| is the area the tube cuts on the shock. The
   spectrum is dN/dp ∝ p^(2−q), q = 3r/(r−1), over [E_min, E_max]
   (`MomentumSpectrum::LocalCompressionDsa`).
3. It places the particles where the same line crosses the sphere r_rel, on
   the upstream side, at the crossing nearest the shock along the line
   (`ReferenceIntersections`).
4. It draws the pitch angle: the outward-flux-weighted distribution for the
   two focused-transport movers, isotropic for Parker.

After release, r_rel is not a boundary; particles cross it freely. A particle
that comes back to the shock is absorbed, because the rate is already net. A
particle that reaches either end of the line (1.1 R☉, about 2 au) leaves the
domain.

L_ref is measured radially. Where the spiral is wound, the path along the
line between the two crossings is longer. For this event with f = 0.3
(computed with the core's `Provider::Sphere` along the Earth-anchored line):

| t | R_sh | r_rel = 1.3 R_sh | L_ref | Path along the line |
|---|---|---|---|---|
| 0 | 2.52 R☉ | 3.28 R☉ | 0.76 R☉ | 0.76 R☉ |
| 1 h | 7.4 R☉ | 9.6 R☉ | 2.2 R☉ | 2.2 R☉ |
| 6 h | 29 R☉ | 38 R☉ | 8.7 R☉ | 9.1 R☉ |
| 12 h | 0.254 au | 0.331 au | 16 R☉ | 18 R☉ |
| 24 h | 0.480 au | 0.624 au | 31 R☉ | 38 R☉ |
| 36 h | 0.682 au | 0.886 au | 44 R☉ | 60 R☉ |
| 41.4 h (run end) | 0.761 au | 0.990 au | 49 R☉ | 70 R☉ |

**Input.** `source.shock.release_distance_shock_radius_fraction = 0.3` in
`event.in`, written from `source.release_distance_shock_radius_fraction` in
`event_settings.json`. The accepted range is (0, 10]. The alternative is a
fixed length, `source.shock.release_distance_m` or
`--shock-release-distance-m`. Exactly one of the three may be given. The value
is part of the run fingerprint and of the checkpoint record, and the
`REDUCED_SOURCE` log line prints it as `reference_shock_radius_fraction`
(code: `ReducedParticles::ReleaseDistanceM`).

### 2.2 What the release surface represents

Models that resolve diffusive shock acceleration (DSA) simulate the upstream
region between the shock and a free-escape boundary (FEB). Particles there are
scattered by the turbulence and cross the shock repeatedly. A particle that
reaches the FEB is counted as escaped into the interplanetary medium. SOLPACS
is an example (Afanasiev et al. 2015), and it supplies the escaping particles
in PARASOL (Afanasiev et al. 2025).

The reduced model here does not resolve that region. It takes the FEB as
given and releases there the net escaping flux of already-accelerated protons,
with the local DSA spectrum. Consequently:
- L_ref plays the role of the FEB distance, and r_rel the role of the FEB.
- The population inside r_rel (the foreshock and its trapped particles) is
  not represented. A model-to-observer comparison is meaningful only while the
  observer is outside r_rel.
- The shock absorbs returning particles. The source rate is already net, so
  letting them cross again would count them twice.

### 2.3 Why the release distance grows with the shock radius

The lengths that decide when a particle stops being trapped near the shock all
grow with r:
- **Focusing length.** L = B/|dB/ds| = r/2 for a radial field. Particles move
  away from the shock where V + D/L exceeds V_sh (Battarbee et al. 2011).
  Vainio et al. (2014) place the FEB where λ/L = 3(u₁ − v_A)/v.
- **Diffusion length ahead of the shock.** κ/(V_sh − U) grows with r: λ∥ ∝ r
  here (Liu et al. 2025), and the shock decelerates.
- **Flux-tube expansion.** A ∝ r² gives an outward drift 2κ/r. For Parker
  diffusion with an absorbing shock, a particle released at distance L escapes
  with probability P ≈ 1 − exp[−L(2/r − ΔV/κ)], where ΔV = V_sh − U.
  - With L = f r the geometric part is 1 − e^(−2f) at every radius: 0.45 for
    f = 0.3.
  - With a fixed length it falls roughly as 1/r. With the former 10⁸ m
    placeholder, P ≈ 9% for 14 MeV protons at 2.5 R☉, 2% at 10 R☉ and 0.1% at
    1 au. That artificial decline would be absorbed by the amplitude
    calibration.
- **The published boundaries scale with distance** (Section 2.4). In absolute
  terms they differ by a factor of about 40. As fractions of the shock
  distance they lie within 0.17–0.5, and steady-state theory allows up to
  about 1.

### 2.4 What other models use

The quotations were checked in the papers' text.

| Model / study | Kind of boundary | Value | f = distance / shock radius |
|---|---|---|---|
| Afanasiev, Battarbee & Vainio 2015 (SOLPACS, coronal shock) | Fixed edge of the upstream simulation box: "The other boundary is taken to be at a distance of 1R⊙ from the shock" | 1 R☉, for conditions at 2 R☉ | ≈ 0.5 |
| Afanasiev, Wijsen & Vainio 2025 (PARASOL; SOLPACS runs with the shock at 1 au) | Fixed FEB distance along the field line, "equal to the simulation box size in SOLPACS (50 R⊙ in all simulations in this work)" | 50 R☉ | ≈ 0.23 |
| Wijsen et al. 2022 (PARADISE with EUHFORIA; energetic storm particle event of 14 July 2012) | Constant FEB length Δ_feb in the foreshock mean-free-path profile λ∥ ∝ (as + b) e^(s/Δ_feb); Δ_feb, a and b "are kept constant during the entirety of the simulation" | 0.17 au, "obtained by multiplying the modelled shock speed at Earth (∼700 km/s) by the duration of the observed ESP onset phase (∼10 h)", and "close to the presumed width of the sheath" | ≈ 0.17 near 1 au |
| Vainio et al. 2014 (semi-analytical foreshock model inside 1 AU) | Where focusing overcomes scattering: λ/L = 3(u₁ − v_A)/v, with focusing length L = R/2 (Eqs. 9, 10, 12). For the steady-state foreshock λ ∝ x + x₀ this gives x_feb ≈ R_S − 2x₀ (Eq. 11), "the boundary of the turbulent foreshock (which can be regarded as a free-escape boundary for the ions)". Eq. 11 "applies only to the late phases of the SEP event at energies below the maximum energy in the spectrum". The model itself uses an energy-dependent x_feb(E), limited by the distance v·t_max the fastest escaping particles can travel (Eq. 24). | Energy dependent; steady-state limit ≈ R_S | Up to ≈ 1 |
| Battarbee, Laitinen & Vainio 2011 (self-generated waves in coronal shocks) | Dynamic trapping boundary: "the distance where V + D/L = Vs can be regarded as the boundary of the turbulent trapping region ahead of the shock" (L the focusing length) | Energy and time dependent | — |
| PATH/iPATH (Zank, Rice & Wu 2000; Rice, Zank & Li 2003; Li, Zank & Rice 2005), as summarized by Ding et al. 2024 | Energy-dependent escape length L_esc = α L_diff, with L_diff = κ/U₁: "L_esc is typically assumed to be 2 − 4 times larger than L_diff" | α = 2–4 | — |
| Ding et al. 2024 (EUHFORIA with iPATH; two energetic storm particle events observed by Solar Orbiter) | The same law, with α fitted to the observed upstream profiles: "For the August event, α ranges from 2 to 3 … while for the September event, it ranges from 3 to 7" | α = 2–3 and 3–7 | — |
| Sokolov et al. 2004 (FLAMPA) | No imposed FEB: "the length of the tail is ∼D/ush" | — | — |

Two cautions apply to any imposed boundary, including this one:
- Lee (2005): "A 'free escape boundary' could be imposed ad hoc to allow ion
  escape, but it is unclear how to choose the location of the boundary as a
  function of v."
- Fraschetti (2021): "a momentum-independent x_FEB is in tension with the
  measurements".

No published model of 2013-04-11 states an escape distance. Liu et al.
(2025), the source of our λ∥, injects at the shock front.

### 2.5 Why f = 0.3

The user chose f = 0.3 on 2026-11-25, for these reasons.

1. **The observers stay outside the release sphere for the whole run.** r_rel
   reaches the observer radius (0.993 au) at the times below. The particle run
   ends at 41.4 h, and the last comparison window (GOES > 10 MeV) ends at
   41.2 h.

   | f | r_rel reaches 0.993 au at | r_rel at the run end (41.4 h) |
   |---|---|---|
   | 0.17 | 48.0 h | 0.891 au |
   | 0.2 | 46.3 h | 0.914 au |
   | 0.25 | 43.8 h | 0.952 au |
   | **0.3** | **41.6 h** | **0.990 au** |
   | 0.4 | 37.8 h | 1.066 au |
   | 0.5 | 34.8 h | 1.142 au |
   | 1.0 | 24.9 h | 1.523 au |

   With f = 0.3 the margin at the end of the run is small: r_rel is 0.0031 au
   (4.7×10⁸ m) inside the observer radius. With f = 0.5 the observers would
   be inside r_rel from 34.8 h, which includes the 36 h held-out windows
   (35–37 h). Computed with the core's `Provider::Sphere`.

2. **It is close to the published values at both ends of the run.**
   - At t0, L_ref = 0.76 R☉ with the shock at 2.52 R☉. Afanasiev et al. (2015)
     used 1 R☉ at 2 R☉ (f ≈ 0.5).
   - At the run end, L_ref = 49 R☉ (0.23 au) with the shock at 0.76 au. That
     matches the 50 R☉ PARASOL box for the shock at 1 au (Afanasiev et al.
     2025) and is above the 0.17 au of Wijsen et al. (2022).
   - One fraction cannot match both ends exactly: the near-Sun model gives
     0.5 and the near-1 au models 0.17–0.23. 0.3 lies between them.

3. **The geometric escape probability is the same at every radius:**
   1 − e^(−0.6) = 0.45 (Section 2.3).

### 2.6 Limitations

- f does not depend on energy. Lee (2005) and Fraschetti (2021) argue that
  the boundary should.
- The energy-dependent PATH law with this model's ambient κ would put the
  boundary far beyond the shock radius at high energy. κ/(V_sh − U) is
  0.02 r at 10 keV and 11 r at 131 MeV for this event's λ∥ (Liu et al.
  2025). That law is meant for the enhanced turbulence at the shock, which
  this model does not contain.
- The foreshock population inside r_rel is not modelled (Section 2.2).
- Suggested off-line sensitivity runs, not nightly references: f = 0.2 and
  f = 0.4.

## 3. Source energy band: 20 keV – 200 MeV (user decision 2026-11-25)

The band was 10 keV – 200 MeV. At t0 the PCHIP through the Lario et al.
heights gives a shock speed of 1459 km/s, the end slope of the first knot
interval (2.52 → 3.57 R☉ in 600 s).

Before creating any particle, srcSEP checks at the release point that every
momentum in the band has an outward shock-frame normal velocity for some
pitch angle. The check runs `BuildFocusedEscapeDistribution` at the minimum
momentum; because the condition is monotone in speed, that covers the whole
band.
- A 10 keV proton (1384 km/s) fails. Part of a 10 keV – 200 MeV net spectrum
  would have no outward escape.
- srcSEP refuses to filter or renormalize a prescribed net spectrum, because
  that would change the prescribed number and its DSA spectrum.
- The strict minimum is about 11 keV.
- 20 keV (1957 km/s) is 1.34 times the initial shock speed.

E_min is part of the source normalization. η n₁ u₁ₙ Φ/|B·n̂| is the number
released in [E_min, E_max], so the share above 1 MeV is
(E_min/1 MeV)^((q−3)/2) (non-relativistic, E_max ≫ 1 MeV). Raising E_min from
10 to 20 keV multiplies that share by 2^((q−3)/2):
- 1.41 for the strong shock (q ≈ 4);
- 5.7 for the weak phase at 4.7–5.2 R☉ (q ≈ 8).

The global amplitude fitted on the ERNE 4 h window absorbs a common factor. It
does not absorb the different weighting of the weak- and strong-shock phases.
All comparisons are at 1 MeV or above. The broadband SOHO diagnostic observer
in `event.in` spans the same 20 keV – 200 MeV band.

## 4. Numerics

- **Line grid.** Geometric spacing with chord ≈ ε r everywhere: first chord
  ε r_base and growth ratio 1 + ε give chord_k = a + ε s_k. Both λ∥ and the
  focusing length scale with r. ε = 0.02 gives 326 points, 1.5×10⁷ m at the
  base and about 3.6×10⁹ m at 1 AU. ε is provisional until the convergence
  study.
- **Open-tube foot trace.** `ambient.trace_step_m = 2e7` m (HSP12). The cost
  per open-shell query is about 0.1 ms; an exterior query costs about 1.3 µs.
- **Candidate pruning.** Line–front intersections (the source, and the
  absorbing-boundary roots) advect only the segments that
  `AmbientModel::SpeedBoundMPerS() × elapsed time` allows to reach the front.
  The roots are those of the whole line.
- **Shock-root table** (`--shock-boundary-root-table on`). A piecewise-linear
  table per AMPS step. Its nodes are exact roots, and every interval is checked
  at 1/4, 1/2 and 3/4 against exact roots to within
  `--shock-boundary-tolerance-m`. Intervals with a connectivity change are
  answered exactly.

## 5. Test tiers: the nightly and the long deck

User decisions 2026-11-25: a short coarse nightly run of all three movers,
and the full event as a separate test that runs outside the nightly; the
proposal for organizing the two tiers was approved and implemented. A Table
block with `Tier=long` and `Ranks=<N>` is written by `table_parser` to
`Makefile.test.long` (included by `Makefile.def.amps`), which the nightly
never reads; the periodic long routine of `utility/nightly`
(`runners.long.json`) runs those tests. The mechanism is generic: any
application's test can be a long test (`doc/ApplicationRegression.tex`,
section "The long tier").

`prepare_event.py deck` writes both decks from the same `event_settings.json`.
`event.nightly.in` is `event.in` with only the values in
`event_settings.json` → `nightly` changed:

| | `event.in` (long tier) | `event.nightly.in` (nightly tier) |
|---|---|---|
| Simulated time | 2484 steps of 60 s (41.4 h) | 60 steps of 60 s (1 h) |
| Grid | ε = 0.02, 326 points | ε = 0.08, 85 points |
| Source samples per step | 200 | 20 |
| Observers | ERNE, EPAM, GOES channels at L1; broadband `soho-spectrum` | the same, plus broadband `nightly-r0p1au`, `nightly-r0p3au`, `nightly-r0p5au` |
| Checked by | DiffNum (gating) and `score_observations.pl` (reported) | DiffNum (gating) |

**Why the inner observers.** In the first hour the source particles fill the
line out to about 1 AU: in an FTE-MFP run of the nightly deck the outermost
occupied `EnergyRSample` bin is at 0.44 AU after 20 min, 0.66 AU after 40 min
and 1 AU after 60 min, and the broadband SOHO observer records its first
particles (two bins near 100–200 MeV) only at the last output. The L1
observers therefore record almost nothing in the nightly window. The three
inner observers (same band, 40 logarithmic energy bins, 24 pitch-angle bins)
give the DiffNum products signal: `regression_products.sh` adds them to the
density and flux products (all observers) and puts their pitch-angle
distributions before the SOHO one. They stay outside the release sphere
1.3 R_sh, which is at 0.045 AU after one hour, and they are not compared with
observations.

### Nightly Table blocks

`MakefileTest/Table` holds `SEP--event-2013-04-11--Parker`, `--FTE-Dmumu` and
`--FTE-MFP` (`Family=sep`). They replace the former
`SEP--CME-shock-transport--*` blocks and differ only in `--particle-mover`.

- **Compile.** `Config.pl -application=test/sep_reduced_front__field_line
  -amps-test=on`, then `make amps`.
- **Run.** Copies `event.nightly.in` (as `event.in`), `event.event`,
  `event.args` and the WSO asset into `run_test_<name>`, then
  `$(MPIRUN) ./amps -i event.in --reduced-shock-event event.event
  --particle-mover <mover> $(cat event.args)`, with output in `run.log` and
  `run.err`.
- **Check.** `regression_products.sh` builds the four products; each is
  compared by `DiffNum.pl` with
  `output/test_SEP3D/SEP--event-2013-04-11--<Mover>--<label>.ref_np<N>.gz`,
  N being the run's MPI ranks. The partial results stay in
  `run_test_<name>/check/`; their concatenation is `test_<name>.diff`, empty
  when the test passes. A missing product or reference is written into it.

**References.** Made with `make_references.sh run_test_<name> <name>
output/test_SEP3D` after the run stage of the Table test itself on passion
(GNU, 8 ranks), so they carry the test's own build (`-amps-test=on`) and
command line. DiffNum uses its defaults (`-a=1e-30 -r=1e-12`), which the
byte-identical repeat runs allow; changing one value of the energy product by
0.5% makes the check fail.

Measured on passion, 8 ranks, 60 steps (run stage only; the compile stage
takes about 2.5 min with 16 jobs):

| Test | Run stage | Table `Time=` |
|---|---|---|
| `SEP--event-2013-04-11--FTE-MFP` | 100–123 s | 5 |
| `SEP--event-2013-04-11--Parker` | 407–460 s | 12 |
| `SEP--event-2013-04-11--FTE-Dmumu` | 1521–1744 s | 30 |

FTE-Dmumu costs most: its steps take 25–45 s each while the shock decelerates
(steps 20–60).

The nightly window is 60 steps (user decision 2026-11-25). Its grid and
samples per step stay provisional until the convergence study.

Two runs of the nightly deck (FTE-MFP, 8 ranks) gave byte-identical density,
flux, energy and pitch-angle products, so DiffNum against np8 references can
gate the test. Number of energy bins with particles in that run, at 12, 24,
36, 48 and 60 min:

| Observer | Bins with particles |
|---|---|
| `nightly-r0p1au` | 2, 16, 25, 30, 30 |
| `nightly-r0p3au` | 0, 5, 8, 16, 20 |
| `nightly-r0p5au` | 0, 1, 1, 8, 13 |
| `soho-spectrum` (1 AU) | 0, 0, 0, 0, 2 |

### Long Table blocks

`SEP--event-2013-04-11-long--Parker`, `--FTE-Dmumu` and `--FTE-MFP`
(`Family=sep-long`, `Tier=long`, `Ranks=16`). They run `event.in`, the full
event, with the same `event.args`, launched as
`$(TESTMPILAUNCH) $(TEST<name>RANKS) ./amps`. The check compares the same four
products with `output/test_SEP3D/SEP--event-2013-04-11-long--<Mover>--<label>.ref_np16.gz`,
fails a run on another rank count, and then runs `score_observations.pl`
(report `run_test_<name>/observation-score.{txt,csv}`, printed in the log,
never gating). Their names do not begin with a nightly test's name, because
the runner collects `test_<name>*.diff`.

By hand, from the AMPS root: `make test_SEP--event-2013-04-11-long--FTE-MFP`
(add `TESTMPILAUNCH='mpiexec -n'` where needed).

`Ranks=16` (2026-11-28): passion has 18 physical cores and Open MPI refuses
`mpirun -np 32` there; 16 ranks leave two cores for the system and the
runner. Provisional until the first long runs: the `Time=` values (1440,
1440 and 360 min), and the deck itself (the convergence study may change the
long numerics). Each `Time=` is then set to 1.5 times the measured wall time
of that test's first full 16-rank run, made in a quiet window of the machine.
The long references do not exist yet: `make_references.sh` makes them from
the same first full runs.

## 6. Pending

1. First full long runs by hand (cost at 16 ranks; references); then the
   periodic crontab line.
2. Convergence (ε, dt, samples per step, boundary tolerance) for the long
   deck; strong scaling for its rank count.

## 7. Tests

- Core: `make -C src/models/sep_corona_swcme shock-front-test`, HSP00–HSP12
  (README of `shock_front`).
- Material admission: `srcSEP/test/test_material_background.cpp`, MBC02.
- Event inputs: `python3 prepare_event.py check --raw /data/vtenishe/sep_validation_data/event_2013_04_11`
  regenerates every checked-in input from the raw archives and byte-compares
  it (`derive` needs the Python package `cdflib`). The Table tests never run
  it; the nightly and the long routine of `utility/nightly` do, as the
  runners `sep-event-inputs` (family `app-sep`) and `sep-event-inputs-long`
  (family `sep-long`), with the machine variables `sep_validation_data` and
  `sep_python_path` (on passion: cdflib 1.2.6 unpacked from its PyPI wheel in
  `/data/vtenishe/sep_validation_data/tools/python`). Any change of
  `event_settings.json` or `raw_manifest.json` requires `derive`, because
  `derived.json` records their SHA-256.

## References

- Afanasiev, A., Battarbee, M., & Vainio, R. 2015, A&A, 584, A81
- Afanasiev, A., Wijsen, N., & Vainio, R. 2025, JSWSC, 15, 3
- Altschuler, M. D., & Newkirk, G. 1969, Sol. Phys., 9, 131
- Battarbee, M., Laitinen, T., & Vainio, R. 2011, A&A, 535, A34
- Ding et al. 2024, A&A, 681, A92
- Fraschetti, F. 2021, ApJ, 909, 42
- Gary, G. A. 2001, Sol. Phys., 203, 71
- Hoeksema, J. T. 1984, PhD thesis, Stanford University
- Kasper, J. C., et al. 2021, PRL, 127, 255101
- Kopp, R. A., & Holzer, T. E. 1976, Sol. Phys., 49, 43
- Lario, D., et al. 2014, ApJ, 797, 8
- Lee, M. A. 2005, ApJS, 158, 38
- Li, G., Zank, G. P., & Rice, W. K. M. 2005, JGR, 110, A06104
- Linker, J. A., et al. 2017, ApJ, 848, 70
- Liu et al. 2025, ApJ, 985, 82 (doi:10.3847/1538-4357/adc4e3)
- Newkirk, G. 1961, ApJ, 133, 983
- Park, J., Innes, D. E., Bucik, R., & Moon, Y.-J. 2013, ApJ, 779, 184
- Parker, E. N. 1958, ApJ, 128, 664
- Rice, W. K. M., Zank, G. P., & Li, G. 2003, JGR, 108, 1369
- Schatten, K. H., Wilcox, J. M., & Ness, N. F. 1969, Sol. Phys., 6, 442
- Sokolov, I. V., et al. 2004, ApJ, 616, L171
- Sun, X., Notes on PFSS Extrapolation, Wilcox Solar Observatory
- Vainio, R., et al. 2014, JSWSC, 4, A08
- Wang, Y.-M., & Sheeley, N. R., Jr. 1992, ApJ, 392, 310
- Weber, E. J., & Davis, L., Jr. 1967, ApJ, 148, 217
- Wijsen, N., et al. 2022, A&A, 659, A187
- Zank, G. P., Rice, W. K. M., & Wu, C. C. 2000, JGR, 105, 25079
