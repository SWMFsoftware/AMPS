# Pitch-angle diffusion library: implementation status and evidence ledger

Specification: `PITCH_ANGLE_SCATTERING_MODEL.md` v1.2 (SHA-256
`1a7ecf1a599ea54a24e5e0da061a0d71c00ed75f4b74e3b561262d724ea96b92`), roadmap
Section 26 (`pitch_angle_diffusion_model_data/detailed_roadmap.md`).
Scope: Release A (user decision 2026-10-13).  Status labels follow Section
26.1: `PASS`, `FAIL`, `BLOCKED`, `SKIP` (not in Release A), `PARTIAL`.
Library evidence, srcSEP3D component evidence and native evidence are kept in
separate columns, as Section 26.8 requires.  No observational validation is
claimed.

## User decisions applied (2026-10-13)

| Item | Decision |
|---|---|
| Scope | Release A: PA-ISO, PA-QFORM, PA-EPS, PA-EPREM; EM/Milstein/rotation kernels; finite-volume reference; other inventory IDs registered and refused |
| Data archive (PS01) | not available; built-in Section 4.2 registry; PS01 BLOCKED |
| Amplitude | per particle from the pusher state: λ∥ = λ_ref (R/R_ref)^a (r/r_ref)^b, D0 from Eq. (23) |
| PA-QFORM | refused in the mover (no Itô pair) |
| Selector | new value `pitch_angle_diffusion_model = pitch-angle-library`; built-in models kept |
| Step rules | rate-fraction, van den Berg 2026 (floored everywhere by the rate fraction), Sandroos 2014; all selectable, no defaults; ν_max, a_max, b_max as sup over μ |
| Clock | host time, order 0 in U/c (Eq. 173 factor not applied) |
| srcSEP | binding deferred; the block is a D3 rejection there |

## Stage ledger

| Stage | Status | Library evidence | srcSEP3D component | Native | Open items |
|---|---|---|---|---|---|
| PS00 host audit | PASS | — | focused-diffusion mover audited: `ResolveLocalTransportImpl` → `EvaluateLocalScattering` → `AdvanceFocused` (symmetric split, Milstein/EM kick, folding); no scattering step limit existed; p, μ in plasma frame; RNG keyed per particle/substep; host owns streaming | — | — |
| PS01 data, versions, schemas | BLOCKED (archive) / PARTIAL | registry of 172 IDs from Section 4.2 with the Section 4.1 counts (PAD06, PAR10); strict schema, typed refusals (PAD05); spec digest pinned | — | — | canonical `PITCH_ANGLE_SCATTERING_MODEL_DATA` (admissibility/capability/domain per ID, 89 fixtures) not supplied |
| PS02 scalar operator, SI, normalization | PASS | PAD01, PAD02, PAD04, PAR01-PAR04, PAR06-PAR08 | PAD3D01 | — | — |
| PS02a tensor | SKIP | — | — | — | Release A has no tensor capability |
| PS03 prescribed shapes | PASS | PAD02, PAD03, PAD04, PAR03, PAR06 | — | — | printed variants (LANG, PACHECO), amplitude variants (HE11, WQ15, Z23, W19, QW15; U-7/U-9/U-10), PA-DROGE-VA (U-14), PA-KOLMO not implemented |
| PS04 scalar kernels and FV reference | PASS (Milstein at ν_max Δt ≤ 0.001; rotation as library kernel) / FAIL (EM) | PAD10-PAD12; AP-ISO-01, AP-TIME-01 (see below) | PAD3D02 | PAN01-PAN03 | EM not certifiable (endpoint-beam bias) |
| PS04a physical time and streaming | PARTIAL | AP-ISO-02 (finite-time displacement), AP-TIME-01 (rate update inside a host interval) | PRK3D08 (limiter), PAD3D01 | PAN01-PAN03 | AP-RNG-01 not applicable (no rejection); AP-CLOCK-01 not applicable (no clock factor) |
| PS05 zero pitch cosine, pair certification | PASS (PA-ISO, PA-EPREM with Milstein) / FAIL (PA-EPS with Milstein) | AP-ZERO-02 (PA-EPS), AP-ZERO-02.eps0 (barrier) | CFG3D20 (gate) | PAN04-PAN06 (refusals) | AP-ZERO-01 (PA-QFORM) not run: q-form not released |
| PS06 rotation / hemispheric / redraw | PARTIAL | rotation: PAD10 (Eq. 169 exact), AP-ISO-01 | — | — | hemispheric (U-13) and redraw (U-63) not implemented |
| PS06b host adapters, lifecycle | PARTIAL | — | PAD3D01-03, CFG3D20 | PAN01-PAN03 (1 and 4 ranks) | AP-GEOM-01..04 (equivalent 1-D/3-D, radial field, forced migration, mid-subcycle restart with the library active) not produced |
| PS07 FTE, frames, reconstruction | SKIP | — | — | — | no frame transforms, no backward pusher in Release A |
| PS08 focusing, observables | PASS (reference and composer) | PAD12 (F-FOC-01/04, F-PUSH-08 from the FV Poisson solve); AP-FOC-01 at f = 0.001 (at 0.003: ξ = +5 FAIL) | PAD3D03 (production composer: exact focusing, coth ξ − 1/ξ) | — | — |
| PS09-PS15 | SKIP | registered, refused (`UnsupportedModel`) | — | — | Release B/C |
| PS16 regression, handoff | PARTIAL | `run_all_tests.py` | `run_tests.py --all`, `run_all_test.sh` steps `pitch-angle-library`, `pitch-angle-native` | | performance profiling of a production workload not done |

## Unit and reference evidence

Command (from this directory): `python3 run_all_tests.py --skip-acceptance`.

| Stage | Result |
|---|---|
| build | PASS (C++17, `-Wall -Wextra -Wpedantic`, no warnings) |
| PAD01-PAD12 | 12/12 PASS (~9 s) |
| PAR01-PAR11 | 11/11 PASS (~4 s): I(q,H) vs mpmath 9.2e-14 (closed) / 1.6e-15 (quadrature); φ(ε) exact to roundoff; D, D′ 7.9e-17 / 1.8e-16; Eq. (23) 1.5e-15; transmissibility 1.6e-15; MFP bundle 3.0e-15; PARALLEL bundle 3.6e-16; F-KIN-01 2.0e-15; `check_physics.py` 261/261; registry digest; contract compiles |

Two reference defects were found and fixed during development and are
recorded here because they changed the evidence:

* a plain `mp.quad` of I(q,H) at q = 1.9, H = 0 is inaccurate (9.52372 vs the
  exact 2/((2−q)(4−q)) = 9.52381); the reference now uses mpmath's `hyp2f1`
  cross-checked by a substituted quadrature;
* the finite-volume transmissibility initially used a recursive integrator
  with a local relative criterion that could not converge near the t^(2/3)
  cusp of PA-QFORM with H > 0 (`Nonconverged`); both translation units now
  share one global-adaptive integrator (`pitch_angle_quadrature.h`).

## Acceptance campaign (Section 26.5)

Binary SHA-256 and results: `build/acceptance/binary.sha256`,
`build/acceptance/results.txt` (regenerate with
`python3 run_all_tests.py --rerun-acceptance`).  The pilot and the first
campaign are kept beside them (`pilot-results.txt`,
`first-campaign-results.txt`); their numbers are bit-identical to the final
run for the cases they share.

Command: `python3 run_all_tests.py --rerun-acceptance` (2026-10-13, 06:44 to
10:40, 14169 s on 16 threads of a node shared with another job).  Binary
SHA-256 `8e79b97e79e947a6361c5c9376e01fafa6a7db15cf6bebc5864f071550c4f96a`.
Result: 27 PASS, 14 FAIL, 11 SKIP (the SKIPs, with reasons, are the cases
of Section 26.5 outside Release A; none counts as passed).

Campaign history (no budget, seed or particle count of a case was changed
after it ran):

1. pilot (`pilot-results.txt`): AP-ISO-01 for the three kernels; fixed the
   scope (Euler-Maruyama kept only to document its failure);
2. first campaign, 44 cases, binary `98238c7b...`
   (`first-campaign-results.txt`, 4127 s): Milstein failed AP-FOC-01 at
   ξ = +5, f = 0.003 (101% of budget), so f = 0.003 could not be certified;
3. second campaign, the 44 cases plus the roadmap's next refinement level for
   Milstein (AP-ISO-02 f = 0.001 to νt = 50 and AP-FOC-01 at f = 0.001 for the
   seven ξ), scope frozen before it ran.  The 44 shared cases reproduce the
   first campaign bit for bit (deterministic per-batch seeds).

| Case | Status | Measured (budget use = (\|err\| + z·SE + ε_ref)/ε) |
|---|---|---|
| `AP-ISO-01.em.f0.01` | FAIL | 23 budget + 0 design failures of 60; first: P4 start=0.500000 nu*t=0.1: estimate -0.100543 +- 0.000127347 vs reference -0.10634 (ref. error 0), \|err\|+zSE+eps_ref=0.00627707 > budget 0.005 |
| `AP-ISO-01.em.f0.003` | FAIL | 9 budget + 0 design failures of 60; first: P2 start=1.000000 nu*t=0.1: estimate 0.730656 +- 9.83897e-05 vs reference 0.740818 (ref. error 0), \|err\|+zSE+eps_ref=0.0105325 > budget 0.005 |
| `AP-ISO-01.em.f0.001` | FAIL | 6 budget + 0 design failures of 60; first: P2 start=1.000000 nu*t=0.1: estimate 0.73603 +- 9.03091e-05 vs reference 0.740818 (ref. error 0), \|err\|+zSE+eps_ref=0.00512825 > budget 0.005 |
| `AP-ISO-01.milstein.f0.01` | FAIL | 2 budget + 0 design failures of 60; first: P4 start=1.000000 nu*t=0.1: estimate 0.360759 +- 0.00019049 vs reference 0.367879 (ref. error 0), \|err\|+zSE+eps_ref=0.00783756 > budget 0.005 |
| `AP-ISO-01.milstein.f0.003` | PASS | worst budget use 53.5% (P4 start=1.000000 nu*t=0.1) |
| `AP-ISO-01.milstein.f0.001` | PASS | worst budget use 32.4% (P4 start=1.000000 nu*t=0.1) |
| `AP-ISO-01.rotation.f0.01` | PASS | worst budget use 31.4% (P4 start=1.000000 nu*t=0.1) |
| `AP-ISO-01.rotation.f0.003` | PASS | worst budget use 27.7% (P1 start=0.500000 nu*t=1) |
| `AP-ISO-01.rotation.f0.001` | PASS | worst budget use 35.8% (P1 start=0.000000 nu*t=1) |
| `AP-ISO-02.milstein.f0.003.t50` | PASS | worst budget use 35.5% (kappa(t) |
| `AP-ISO-02.milstein.f0.001.t10` | PASS | worst budget use 29.4% (kappa(t) |
| `AP-ISO-02.rotation.f0.003.t10` | PASS | worst budget use 23.4% (kappa(t) |
| `AP-ISO-02.milstein.f0.001.t50` | PASS | worst budget use 30% (kappa(t) |
| `AP-TIME-01.em` | PASS | worst budget use 26.2% (P2) |
| `AP-TIME-01.milstein` | PASS | worst budget use 28.9% (P1) |
| `AP-TIME-01.rotation` | PASS | worst budget use 24.4% (P2) |
| `AP-ZERO-02.eps0.01.milstein.f0.01` | FAIL | 21 budget + 7 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0139865 +- 3.12449e-05 vs reference 0.01 (ref. error 5.39e-12), \|err\|+zSE+eps_ref=0.00410416 > budget 0.0005 |
| `AP-ZERO-02.eps0.01.milstein.f0.003` | FAIL | 13 budget + 6 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0118986 +- 2.83159e-05 vs reference 0.01 (ref. error 5.39e-12), \|err\|+zSE+eps_ref=0.00200525 > budget 0.0005 |
| `AP-ZERO-02.eps0.01.milstein.f0.001` | FAIL | 9 budget + 4 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0106754 +- 2.65704e-05 vs reference 0.01 (ref. error 5.39e-12), \|err\|+zSE+eps_ref=0.000775386 > budget 0.0005 |
| `AP-ZERO-02.eps0.048.milstein.f0.01` | FAIL | 20 budget + 6 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0124867 +- 3.09899e-05 vs reference 0.01 (ref. error 5.01948e-12), \|err\|+zSE+eps_ref=0.00260334 > budget 0.0005 |
| `AP-ZERO-02.eps0.048.milstein.f0.003` | FAIL | 9 budget + 4 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0113363 +- 2.76171e-05 vs reference 0.01 (ref. error 5.01948e-12), \|err\|+zSE+eps_ref=0.00144028 > budget 0.0005 |
| `AP-ZERO-02.eps0.048.milstein.f0.001` | FAIL | 8 budget + 4 design failures of 60; first: P(\|mu\|<0.01) uniform nu0*t=0.2: estimate 0.0106037 +- 2.87654e-05 vs reference 0.01 (ref. error 5.01948e-12), \|err\|+zSE+eps_ref=0.000711971 > budget 0.0005 |
| `AP-ZERO-02.eps1.milstein.f0.01` | FAIL | 3 budget + 4 design failures of 60; first: P(\|mu\|<0.01) zero nu0*t=0.2: z*SE=0.000143236 exceeds budget/4=0.000125 (particle count too small) |
| `AP-ZERO-02.eps1.milstein.f0.003` | FAIL | 2 budget + 3 design failures of 60; first: P(\|mu\|<0.01) zero nu0*t=0.2: z*SE=0.000131627 exceeds budget/4=0.000125 (particle count too small) |
| `AP-ZERO-02.eps1.milstein.f0.001` | FAIL | 0 budget + 4 design failures of 60; first: P(\|mu\|<0.01) zero nu0*t=0.2: z*SE=0.000157523 exceeds budget/4=0.000125 (particle count too small) |
| `AP-ZERO-02.eps0.em` | PASS | crossing fraction at nu_max*dt = 0.01/0.003/0.001: 0.076/0.0599/0.049 (reference: 0, disconnected) |
| `AP-ZERO-02.eps0.milstein` | PASS | crossing fraction at nu_max*dt = 0.01/0.003/0.001: 0.000279/0/0 (reference: 0, disconnected) |
| `AP-FOC-01.milstein.xi0.f0.003` | PASS | worst budget use 20.5% (kappa_var) |
| `AP-FOC-01.milstein.xi0.1.f0.003` | PASS | worst budget use 24.4% (kappa_var) |
| `AP-FOC-01.milstein.xi-0.1.f0.003` | PASS | worst budget use 21.7% (kappa_var) |
| `AP-FOC-01.milstein.xi1.5.f0.003` | PASS | worst budget use 32% (kappa_var) |
| `AP-FOC-01.milstein.xi-1.5.f0.003` | PASS | worst budget use 24.5% (kappa_var) |
| `AP-FOC-01.milstein.xi5.f0.003` | FAIL | 1 budget + 0 design failures of 2; first: kappa_var: estimate 0.0104103 +- 6.76537e-06 vs reference 0.0103248 (ref. error 0), \|err\|+zSE+eps_ref=0.000104505 > budget 0.000103248 |
| `AP-FOC-01.milstein.xi-5.f0.003` | PASS | worst budget use 94.7% (kappa_var) |
| `AP-FOC-01.milstein.xi0.f0.001` | PASS | worst budget use 23% (kappa_var) |
| `AP-FOC-01.milstein.xi0.1.f0.001` | PASS | worst budget use 31% (kappa_var) |
| `AP-FOC-01.milstein.xi-0.1.f0.001` | PASS | worst budget use 36.4% (kappa_var) |
| `AP-FOC-01.milstein.xi1.5.f0.001` | PASS | worst budget use 33% (kappa_var) |
| `AP-FOC-01.milstein.xi-1.5.f0.001` | PASS | worst budget use 22.5% (kappa_var) |
| `AP-FOC-01.milstein.xi5.f0.001` | PASS | worst budget use 47.8% (kappa_var) |
| `AP-FOC-01.milstein.xi-5.f0.001` | PASS | worst budget use 50.9% (kappa_var) |

Findings:

* **Milstein** passes every case at ν_max Δt = 0.001 (worst 50.9%,
  AP-FOC-01 ξ = −5); at 0.003 it fails AP-FOC-01 ξ = +5 (the explicit focusing
  drift inside the kick); at 0.01 it fails the endpoint beam of AP-ISO-01.
* **Euler-Maruyama** with mirror folding fails AP-ISO-01 at every step: the
  endpoint-beam P2 error decreases only from 0.0102 (f = 0.003) to 0.0048
  (f = 0.001), the reflection bias of Section 16.4.
* **PA-EPS with Milstein** (AP-ZERO-02) over-populates |μ| < 0.01 near the
  drift jump D′(0±) = ±D0: from the uniform (stationary) start the exact
  band probability 0.01 is measured as 0.01399 / 0.01190 / 0.01068
  (ε = 0.01) and 0.01249 / 0.01134 / 0.01060 (ε = 0.048) at f = 0.01 / 0.003 /
  0.001, against a 5e-4 budget; the bias decreases roughly as f^0.6.  An
  independent numpy Milstein of the same SDE gives 0.013991 ± 1.6e-5
  (ε = 0.01, f = 0.01) and 0.011295 ± 1.5e-5 (ε = 0.048, f = 0.003), so the
  failure is the scheme's, not the implementation's.  For ε = 1 there is no
  budget failure at f = 0.001 but design failures remain (below), so it is
  inconclusive.
* **Design defect (AP-ZERO-02).** The counts of the concentrated "zero" and
  "near-zero" starts assumed a band probability near 0.01; it is much larger
  there, so z·SE exceeds ε/4 (design failures in every AP-ZERO-02 case).
  Counts were left frozen.  The PA-EPS verdict rests on the uniform start,
  whose sizing is correct.
* **ε = 0 barrier**: the finite-volume operator is disconnected and the
  crossing fraction decreases with the step for both Itô kernels (Milstein
  crosses exactly zero at f ≤ 0.003).

## Certification matrix

`pitch_angle_certified_pairs.inc` (not linked into the acceptance binary:
its SHA-256 is unchanged by the edit, so the recorded run stays valid).
`CheckKernelPair` refuses everything else with `IncompatibleKernel`.

| Model | Kernel | Certified ν_max Δt | Evidence | Status |
|---|---|---|---|---|
| PA-ISO | Milstein | ≤ 0.001 | AP-ISO-01 f0.001, AP-ISO-02 f0.001 t10/t50, AP-TIME-01, AP-FOC-01 f0.001 (7 ξ) | certified |
| PA-EPREM | Milstein | ≤ 0.001 | same: the PA-EPREM shape is h = 1, so the kernel sees the PA-ISO D(μ)/ν_max | certified |
| PA-ISO, PA-EPREM | Euler-Maruyama | — | AP-ISO-01 fails at 0.01, 0.003, 0.001 | not certified |
| PA-EPS | Milstein | — | AP-ZERO-02 fails for ε = 0.01, 0.048 at every step; ε = 1 inconclusive | not certified |
| PA-EPS | Euler-Maruyama | — | not run (EM not certifiable) | not certified |
| PA-QFORM | Euler-Maruyama, Milstein | — | unbounded drift at μ = 0 | refused in code |
| PA-ISO | random rotation | — | AP-ISO-01/02 pass at 0.003 and 0.001, AP-TIME-01 passes; no AP-FOC-01 case run; no host scheme | not certified |
| any | finite-volume chain | — | reference operator, not a pusher | n/a |

srcSEP3D consequence: the library can be used with `reflecting-milstein`
and PA-ISO or PA-EPREM with any step rule whose worst-case ν_max Δt is at
most 0.001 (`rate_fraction ≤ 0.001`).

## srcSEP3D evidence

| Test | Result |
|---|---|
| `PAD3D01` | PASS: bridge = library for 3 rules x 5 μ with a translated solar origin |
| `PAD3D02` | PASS: `AdvanceFocused` kick = library kick (same normal, folding) |
| `PAD3D03` (Extended) | PASS: deterministic focusing order 2.009; stationary ⟨μ⟩ = 0.438232 ± 0.000411 vs 0.438125 (\|err\| + z·SE = 1.17e-3 ≤ 5e-3), 205 s |
| `PRK3D08` | PASS: nine named limiters |
| `CFG3D20` | PASS: PA-ISO/Milstein deck at ν_max Δt = 0.001 (first Milstein row); identity, installation and the refusals including the uncertified PA-EPS/Milstein pair |
| stage-1 routine (`./test/stage1 --all-tests`) | PASS 147, FAIL 0, SKIP 15 (the SWCME suspensions), 10.9 s |

## Native evidence (srcSEP3D)

Build: preflight (root confirmed, no Mars2 build/test process), `rm -rf --
build`, `./Config.pl -application=sep3d`, `make -j16` (1 min 30 s);
`build/main/mainlib.a` contains `pitch_angle_diffusion.o`,
`pitch_angle_kernels.o` and `pitch_angle_certification.o`.

Command: `python3 srcSEP3D/test/run_pitch_angle_native.py --amps ./amps
--output-dir test_output/pitch-angle-native` (reduced-front injection deck,
`transport = focused-diffusion`, PA-ISO/Milstein, `rate_fraction = 0.001`,
20 steps).

| Test | Result |
|---|---|
| `PAN01.rank1`, `PAN01.rank4` | PASS: PA-ISO installed, fingerprint `70099412c1dea376...`, 20 steps, normal end |
| `PAN02.rank1`, `PAN02.rank4` | PASS: 20088 particles injected, 0 escaped, 0 absorbed, 0 mover failures |
| `PAN03` | PASS: per-step ledger rows identical on 1 and 4 ranks |
| `PAN04` | PASS: PA-QFORM refused at start-up (exit 2, "unbounded") |
| `PAN05` | PASS: Euler-Maruyama refused at start-up (exit 2) |
| `PAN06` | PASS: PA-EPS/Milstein refused at start-up (exit 2, "no certification entry") |

This is integration evidence for the binding (installation, mover with the
library active on one and four ranks, refusals).  It is not a physics
validation of focused transport over 20 steps; the physics evidence is the
acceptance campaign and `PAD3D01`-`03`.

srcSEP (shared input layout changed): native rebuild with the preflight
(`./Config.pl -application=test/sep_parker_spiral__field_line`, `make -j16`);
`INIT01`-`INIT03` PASS through `./amps` and through
`srcSEP/test/run_initialization_tests.sh`, including the D3 rejection of the
`physics: pitch_angle_diffusion` block.

Full srcSEP3D suite: `env MAKEFLAGS=-j16 python3 srcSEP3D/test/run_tests.py
--all --amps-source . --make-config Makefile.conf --rebuild --output-dir
test_output/run-tests-all` -- PASS 170, FAIL 0, SKIP 44, ERROR 0 (SKIPs:
SWCME suspensions, evidence-root and Phase-V cases needing external data).
`BLDL3D01` passes against the fresh native build (the three library objects
are in `mainlib.a`); `CFG3D20`, `PAD3D01`-`03` and `PRK3D08` pass.
