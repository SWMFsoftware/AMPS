# Review of AMMENDED_ROADMAP.md

Reviewed 2026-10-19 against the live CPU tree (`../LiouvilleCOPS`, edited 13:02–13:06
today), the CUDA tree and the external output root. The amendment's decisions are
incorporated into `ROADMAP.md` (revision 2). The original plan is kept unchanged as
`ROADMAP_v1.md`.

Each comment is labelled:
* **[agree]**: adopted as written;
* **[fact]**: a statement in the amendment that is stale or inaccurate, with the
  verified state;
* **[add]**: something the amendment does not cover and the plan now does;
* **[decision]**: needs the user.

## 1. The fitting scheme and the oscillations

1. **[agree] The change of primary target addresses the cause of the oscillations.**
   * **Old plan:** the original ROADMAP §3.3/C09 put the fixed-baseline harmonic
     experiments at the centre of the campaign: degree-4 control, then degrees 5–8. There,
     the 25 legacy coefficients, the three temperatures and κ are all frozen, and only
     higher-degree modes are free (degree 5 alone has 11; degrees 5–8 have 56).
   * **Why that oscillates:** with the baseline frozen, every misfit, including smooth
     long-time trends that a different baseline or temperature would explain, must be
     absorbed by short-wavelength spatial modes. Rotation and the changing illumination
     and viewing geometry project short-wavelength surface structure into time as humps
     and dips: the oscillations.
   * **New plan:** the primary families B29/N30/C30/NC31 restore the free degree-0–4
     baseline, temperatures and κ, and add at most two fixed-shape, nonnegative
     templates. The spatial freedom stays low-wavelength or physically prescribed.
   * **Status of the harmonic modes:** they become explicitly secondary, labelled
     experiments.
2. **[agree] Do not blame λ.**
   * The completed L5 scan (`fit_shap5_788k_l5_harmonic_all_available_scan`, converged
     17:24 UTC) varied λ over 1e-6 to 1e-3 at fixed parameter freedom. It has no matched
     λ = 0 degree-5 control.
   * λ = 0 is therefore the primary setting, and positive λ is a labelled extension.
3. **[add] Make "no oscillation" measurable, not a visual impression.** The amendment
   rightly forbids promising smooth curves or filtering the model, but it gives no way to
   tell an oscillating fit from a good one. The revised plan adds *reported, non-gating*
   diagnostics for every C09 fit, computed per gauge on the selected observations:
   * the number of sign changes and the lag-1 autocorrelation of the log residuals in time
     order (long runs of one sign mean the model misses a trend, and oscillation shows as
     alternating runs);
   * the second-difference ("roughness") energy of the model log pressure against that of
     the data, over common observation times, for the full 10 h and the final descent
     separately;
   * the same for the archived B29 result, as a reference.

   These are diagnostics only. **[decision]** They would become a gate only through a
   separate user-approved criterion.
4. **[add] B29 alone does not guarantee physically smooth or credible curves.**
   * The archived B29 stopped at max_nfev = 1000. Its temperatures sat at the bounds:
     Tmax = 30 K, T75 = 400 K, Tmin ≈ 49 K.
   * Tmax < T75 inverts the analytical law: temperature rises toward the terminator. In
     that regime the activity map and the temperature law can trade off, and the trade-off
     can itself produce structure in time.
   * The amendment already defers temperature plausibility to a user decision (§8).
   * The plan adds a concrete trigger: if B29 on 788k again ends at a temperature bound,
     C09 reports it, runs the already-supported ordered-temperature mode as a *labelled*
     control, and stops for the user's physical decision before any interpretation.
5. **[add] The single clip can also put kinks in time.**
   * max(0, u) with a free baseline can switch facets on and off sharply; the derivative is
     zero at and below the clip, as on the CPU.
   * C06 tests must include facets moving across the clip between Jacobian evaluations.
   * C09 reports the clipped-area fraction over the iterations, not only at the end.
6. **[agree] The N mask is applied throughout the trajectory**, not switched on at the NG
   jump. **[add]** Report the area overlap between `late_ng_not_rg` (4.9e-4 of the area)
   and `cliff_ge_60deg` (3.77% in the C03 field), and the Jacobian correlation of the two
   amplitude columns in NC31. If the overlap is large, the amplitudes are not separately
   identifiable.

## 2. Statements to correct or update

1. **[fact] §1 location.** The file is `LiouvilleCOPS_cuda/AMMENDED_ROADMAP.md`, in the
   *same* underscore directory as the implementation; there is no separate
   `LiouvilleCOPS-cuda` directory. Its `../LiouvilleCOPS_cuda/...` links resolve anyway.
   The filename spelling ("AMMENDED") is kept, since the plan links to it.
2. **[fact] §3 "C00 acceptance … pending".**
   * Policy: the numerical policy was approved on 2026-10-19 after this file was written,
     with the objective budget now derived from the pressure budget. The Jacobian and
     end-to-end fit budgets are still provisional (`VALIDATION.md`).
   * Fixtures: the amendment's *re-freeze* requirement does still hold. The CPU reference
     tree changed at 13:02–13:06 (`fit_cops.py`, `gravity_cache.py`, the builder,
     `gravity_topology.py`), after the C00 capture (06:45 UTC). The C00 fixtures are bound
     to the earlier tree digest, so a versioned **C00b** capture is required for joint-mode
     states. The original fixtures remain valid history.
3. **[fact] §3 C02 "extend to joint 29–31 modes".** The contract already supports them.
   * The column order is a(fit_j), templates, T, κ, read from the CPU `Problem`, not
     hard-coded.
   * `tests/test_backend_contract.py` now includes a synthetic NC31 configuration: 25 + 2
     + 3 + 1 = 31 columns, with the joint eigen-regularization. Residuals and Jacobians
     through the facade are bitwise identical to the CPU: 19 tests pass.
   * What remains for C02: template/prior identities in `backend.identity()`, parameter
     state serialization, and the memory plan of the live 31-column evaluator, which
     belongs with C06.
4. **[fact] CUDA STATUS statements now stale.**
   * "L5 scan running": it has completed.
   * "HEAL builder unfixed": it is fixed, and the five products are quarantined as
     `*_heal_defect`.

   `STATUS.md` is reconciled in this revision.
5. **[fact] §4.5 audit list against the CUDA audit.**
   * `cops_cuda/gravity_body.py` checks: closed, consistent orientation, Euler 2, V > 0,
     no zero-area faces, no unreferenced vertices.
   * The CPU's `gravity_topology.require_outward_closed_body` additionally checks one
     connected component, manifold vertex links and unique triangles.
   * **[add]** C04: the CUDA side calls the reference audit, so there is one definition,
     and keeps its own as a cross-check. The C03 evidence is unaffected: the 788k body
     passes both.

## 3. Gravity admission (C03/C04/C05)

1. **[agree]** HEAL artifacts stay quarantined. Cliff columns are selectable only with
   repaired provenance.
2. **[decision] Acceptance route for the production gravity field.** Recommendation: run
   the *repaired* CPU builder once for the full 788k field and compare it with the C03 GPU
   field (`cuda_c03/gravity_20261019T171225Z/full_field.npz`) facet by facet.
   * Cost: about 3.8 h on 36 threads, the measured cost of the HEAL build. That is small
     next to the campaign.
   * Benefit: it turns "sample-validated" into full-field equivalence and removes the
     question of whether a GPU import may stand in for a CPU product.
   * Interaction with CPU fits: if the CPU is busy with B/N fits, schedule it rather than
     overlap.
   * In C03 no facet lay within 1e-6° of 45° or 60°, so any classification difference
     would be a real defect.
3. **[add] Cache identity for a GPU-produced field.**
   * The amendment requires the repaired algorithm id
     (`…_v2_disable_after_topology_audit`) and a truthful record of the backend.
   * The algorithm id names the *mathematical definition*, which is identical for both
     backends. The backend, compiler and evidence belong in the cache's calculation
     metadata.
   * The CPU builder's reuse path also compares package-binary identity. A GPU-written
     cache would fail that comparison.
   * C04 must therefore extend and test the schema (a `backend` field accepted by the
     loader) rather than copy CPU metadata. That change touches the CPU reference, so it
     needs the user's go-ahead. **[decision]**, only if the GPU import route is chosen.

## 4. Evaluator (C06) consequences of the amended scheme

1. **[add] There are no fixed-temperature shortcuts in the primary modes.**
   * B29–NC31 free the temperatures, so the per-contribution temperature factors, u = T^-½
     and du/dT, are recomputed at every parameter vector. The "precompute temperature
     factors once" specialization (original C06 item 5) does not apply.
   * Performance estimates must be made for this case (amendment C10; agreed).
2. **[add] Jacobian work is dominated by the 25 baseline columns.**
   * For each contribution, the 25 Y-columns, the template columns and the three
     temperature sums must be produced.
   * Plan: one fused pass per column tile (`memory.plan_tiles`), with per-point
     fixed-order partial sums, so there is never a contribution × column array.
   * The analytic columns for templates and baseline share the masked factor
     `[u>0]·w·u`; compute it once per tile.
3. **[agree]** Return only pressures and derivative columns. Residuals, regularization and
   TRF stay in the reference CPU code; the C02 facade already enforces this.

## 5. Order of work (incorporated in ROADMAP §7)

The tracks below can run concurrently:
* **Contracts:** C00b re-freeze with B29/N30/C30/NC31 states, then C06 for B29 first.
* **Gravity:** the full repaired CPU field, if approved, then C04 cache import or build,
  then C05 repaired cliff bundle.
* **CPU fitting:** CPU fitting of B and N needs no gravity and can start any time on the
  CPU; that is the user's campaign, not a CUDA stage.

C/NC production waits for C05.

## 6. Decisions requested

1. Gravity acceptance route (§3.2): recommended, a full repaired CPU field plus a
   facet-wise comparison.
2. Whether the oscillation diagnostics (§1.3) stay reported-only (recommended) or become a
   criterion.
3. If B29 again hits temperature bounds: the physical decision of amendment §8, using the
   trigger defined in §1.4.
4. A GPU-cache schema extension in the CPU loader (§3.3), only if the import route is
   chosen.
