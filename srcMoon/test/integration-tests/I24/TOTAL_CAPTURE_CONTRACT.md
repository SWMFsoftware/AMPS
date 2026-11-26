# I24 total proton–sodium capture: implementation and test contract

The user selected this case on 2026-11-25 after the
[quantitative source review](SOURCE_REVIEW_20261125.md). The same day the
user approved four implementation details: log-log interpolation of the
shell partials, a constant 0.4° center-of-mass recoil bound, radiative
decay booked at the event, and NIST ASD shell centroids. This contract
records the scope. It is supporting development evidence: formal I24 keeps
its predecessors and remains **SKIPPED / NOT QUALIFIED**.

| Item | Specification |
| --- | --- |
| ID / name | I24 supporting case `sodium-total-capture` (ensemble) and `sodium-pic-total-capture` (linked PIC). No new test ID. |
| Purpose | Verify the production provider for H⁺ + Na(3s) → Na⁺ + H(n), n = 1, 2, 3: energy-dependent total rate, sampled final shell, zero-angle two-body products, radiated-energy ledger and normal-startup coupling. |
| Production capability | `Moon::Species::SodiumTotalCapture::Provider` and the normal-startup owner `Moon::Chemistry::SodiumTotalCaptureCohort`, registered through `NativeRegistration::SetBeforeParticleBuffer` and run by the reaction controller. |
| Model identity | `hplus-na3s-total-capture-forward-1to5kev-u-v1` |
| Rate | `ν = n_p σ_tot(E) g`. `σ_tot` is the installed Morgan et al. 1985 Table 2 Na law (`ProtonNaMorgan1985LogLog`, source-knot log-log), loaded through its hash-checked publication. `g` is the relative speed and `E = ½ u g²` in eV/u. |
| Energy envelope | 1000 ≤ E ≤ 5000 eV/u, inclusive. Outside is refused with no clamping or extrapolation. |
| Final-shell split | Shingal et al. 1986 Table 4 (A34) rows at 1.0, 4.0 and 10.0 keV/u. Shell partials are σ₁ = σ1s, σ₂ = σ2s + σ2p and σ₃ = σ3s + σ3p + σ3d. Each shell partial is interpolated log-log between adjacent rows and the three are normalized to fractions. A row energy uses that row exactly. |
| Shell energies | NIST ASD H I levels. The excitation of shell n is the statistical-weight centroid `Σ g_j E_j / Σ g_j` over its nlj levels, with g = 2J + 1. NIST's separately listed n-shell rows are recorded as a cross-check only. |
| Selection | One open uniform u (purpose Kinematics, index 0) selects n = 1 if u < f₁, n = 2 if u < f₁ + f₂, otherwise n = 3. |
| Kinematics | Zero-angle exchange in the pair center-of-mass frame: H continues along the incident relative velocity and Na⁺ opposite. The available COM energy is `½ μ g² + Q₁ − E_exc(n)`, with `Q₁ = IE(H) − IE(Na)` from the NIST inventory. Momentum and energy close exactly; the longitudinal recoil follows from Q(n). |
| Radiation ledger | H enters transport as H(1s). `E_exc(n)` is booked at the event as radiated energy leaving into the untracked environment, together with the consumed proton (`−m_p`, `−m_p v_p`, `−(½ m_p v_p² + IE(H))`). Photon momentum is omitted; the H recoil it would cause, `E_exc(n)/(m_H c)`, is reported as a bound. |
| Recoil bound | Approved constant θ_max = 0.4° in the pair COM. The omitted transverse impulse bound is `2 p_f sin(θ_max/2)`, where p_f is the final COM momentum. The Na⁺ recoil-energy bound is that impulse squared over `2 m_Na⁺`. Na⁺ velocity observables are admitted only above this scale. The tail beyond the source extent is an unquantified, recorded limitation. |
| Target domain | Stationary, force-free ²³Na(3s) in a uniform monoenergetic proton beam, as in the existing I24 fixture. A moving or accelerated parent is refused on every accepted segment. |
| Configuration | Programmatic normal-startup request: provider, beam density and velocity [SI], semantic seed and cohort, explicit e = 1.602176634e-19 C and solver effort. The frame is frozen LSO. |
| Inputs | Morgan 1985 PDF (`7c412cb9…`), Shingal 1986 PDF (`58bc33f7…`), a retained NIST ASD H I response (acquired by `acquire_hydrogen_levels.py`) and the NIST mass inventory. `prepare_total_capture.py` derives two pinned CSVs whose SHA-256s are compiled into the provider. Each run regenerates both CSVs and requires byte equality. |
| Independent oracle | Python `Decimal`: the Morgan knots read from the case file and interpolated log-log independently; the Shingal rows and NIST levels parsed independently of C++; the keyed uniforms reproduced with RFC 2104 HMAC; the expected shell per event; two-body zero-angle momenta from the masses and Q(n); exponential survival. |
| Units / frame / time | SI internally; source eV/u and Å² converted explicitly. Event frame is frozen LSO, with elapsed seconds from the interval origin. |
| Stochastic controls | Seed 731993 with separate cohort names per energy, and 32,768 original parents per energy at one mean optical depth. The same histories are run at 1 and 8 MPI ranks. |
| Energies | 1.5, 3.0, 4.0 and 4.75 keV/u. Together they cover Morgan interpolation and knot points, Shingal interpolation, the exact 4.0 keV/u row and the 4–10 keV/u interval. |
| Controls | Zero density (rate 0, no event); negative or nonfinite density refused; 0.99 and 5.01 keV/u refused; a moving parent refused; a closed-endpoint uniform refused with the prior plan preserved; rescaled daughter weights preserve represented numbers. |
| Metrics | Direct rate; Morgan σ; shell fractions; per-event shell; lifetime; branch frequencies; momentum, energy and Q closure; zero angle; radiated energy; recoil bounds; MPI invariance. |
| Acceptance | Direct rate error ≤ 1%; ensemble lifetime error ≤ 3%; per-event shell identical to the oracle; shell frequencies consistent with f_n under exact binomial tests at a 1% family false-alarm probability; exact charge, nucleus and represented-number closure; momentum, energy and Q closure, transverse velocity and fractions within the arithmetic bound 5e-12; bit-identical parent histories at 1 and 8 ranks. |
| Linked PIC | One stationary Na packet at (1e7, 1001, 1003) m in frozen LSO; first- and second-order movers at 1 and 8 ranks; optional step (1, 2, 4) and weight (1, ½, ¼) refinement. Checks: event time, the bath and radiation ledger, both daughters and their ballistic continuation. Native registration refusals cover a non-LSO frame, an energy outside the envelope, negative density and zero solver effort. |
| Expected PASS | All checks above PASS at every energy and layout; formal I24 stays SKIPPED. |
| Expected failure modes | Wrong rate law or unit; wrong shell interpolation; misread table; wrong Q sign; isotropic angle leaking in; missing radiation term; MPI-dependent draws; moving target accepted; out-of-envelope energy accepted. |
| Status semantics | A missing source, compiler or MPI gives SKIPPED; a corrupt or changed source, a build failure or a launch failure gives ERROR; a measured mismatch gives FAIL. |
| Artifacts | `result.json`, `run_manifest.json`, stdout/stderr, input snapshot, source and binary hashes, per-rank parent CSVs, `comparison.csv/png/eps`, and for PIC per-layout event, product and particle CSVs. |

**Limitations.** The tail of the scattering distribution beyond 0.4° COM is
not quantified by any acquired source. Capture into n ≥ 4 is outside
Shingal's basis. Fine structure inside a shell is replaced by its centroid.
H(2s) metastability and cascade timing are not modeled. The beam is
monoenergetic, and thermal or measured proton distributions (R4.1/D06) are
not admitted.
