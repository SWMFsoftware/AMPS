# I24 sodium ground-capture development case

| Contract item | Controlled selection |
| --- | --- |
| Purpose | Verify the single-energy partial charge-exchange rate, residence-time distribution and two-heavy-product event plan using the production C++ functions. |
| Reaction and source | H⁺ + Na(3s) → H(1s) + Na⁺(2p⁶) at 4.0 keV/u projectile kinetic energy in the target-rest frame. Shingal Table 4 supplies 0.22 Å². The reviewed r14 model binds this partial to the isotropic pair-COM angular assumption. |
| Bath and units | Stationary sodium targets, constant prescribed proton density 2.5 × 10⁶ m⁻³ (the P07 development fixture), uniform monoenergetic beam. Times are elapsed SI seconds in the chosen inertial target frame. The density is a controlled input. No measured solar-wind interval is claimed. |
| Mass and energy convention | Reference masses from the retained NIST isotope/CODATA packet subtract one electron mass for a singly charged atom. Electronic formation energies are handled explicitly; Q = IE(H) − IE(Na). The isolated elemental threshold records carry no isotope-shift qualification. |
| Production path | r14 database reader, Moon ground capture provider, production parent-keyed optical-depth generator, constant-hazard clock, transport event planner and signed projectile-bath accounting. The MPI adapter supplies a homogeneous particle ensemble, explicit species slots, domain and weights. Local PIC allocation, mover trajectories and lunar fluxes are outside the adapter's scope. |
| Independent oracles | With projected area σ from the printed table and speed g = √(2 E/u), the per-target frequency is nσg and survival is exp(−nσgt). Independent parent keys select first events. Check both product momenta and kinetic energies directly from output velocities and source reference masses. |
| Random controls | 32,768 original Na targets, campaign seed 731993, one mean optical depth of prescribed exposure. Purpose-separated production streams are keyed by parent ID. Repeat each history at 1 and 8 MPI ranks by default. Compare physical records after sorting on parent ID. |
| Metrics | Reference rate, endpoint survival/lifetime, event time, realized angular CDF, represented counts, energy and vector-momentum closure including the consumed proton, and MPI invariance. Inspect accepted product velocities for every event. |
| Acceptance | Direct rate error ≤ 1%, endpoint lifetime relative error ≤ 3%, first-hit errors within floating-point tolerance, charge/nuclear identity and both daughter weights consistent per event. Isotropy is tested against its specified solid-angle distribution. No physical angular uncertainty or off-energy cross section is inferred. |
| Artifacts and status | Source and lock digests, exact command lines, compiler identity, parent-keyed records, comparison CSV/PNG/EPS and result/manifest in the selected output directory. Missing required input, compiler or MPI launcher gives SKIPPED; corruption, parser, build or launch errors give ERROR; a measured mismatch gives FAIL. A supporting PASS leaves formal I24 gated on its existing requirements. |

The direction tests use a family false-alarm probability of 1%. For each of
the polar-cosine and azimuth CDFs, the Dvoretzky–Kiefer–Wolfowitz bound is
√[log(4 / 0.01) / (2N)], where N counts captured parents. Each product has a
separate reference weight (2 and 4 represented particles); changing both by
a factor of two must preserve each daughter's represented number, momentum and
velocity for the same parent. The zero endpoint of the angular draw must be
rejected while the most recent valid event plan remains available.

Published sources and native SHA-256s are bound by r14 and the independently
replayed NIST four-source mass packet. The contract freezes the r14 lock,
inventory, seed, controls and primary numerical comparisons in
[the input case](reference/sodium_case.json). Absolute rate/lifetime thresholds
follow [the formal reference](reference/acceptance.json). The current quantity
is conditional on the specified final states and energy. The computed lifetime
is under a spatially homogeneous test bath continually maintained for one
optical depth. It is not a prediction of a lunar neutral's residence time.
[Source selection and use](../../../../utility/chemical_database/DATABASE_README.md#sodium-pickup-and-binary-charge-exchange-r14)

## Linked particle transaction contract

Use one stationary Na packet in the frozen LSO frame at (10,000,000, 1001,
1003) m. The H⁺ beam unit vector is (0.6, 0.8, 0). The packet represents
2.625 atoms. Use statistical reference weights of 7 for Na, 56 for Na⁺ and
0.875 for H. Set the trial duration to the independently calculated first
reaction time for parent 0, seed 731993 and cohort
`I24-PIC-ground-capture-v1`, plus exactly 1/16 s. This conditioning makes
both daughters observable inside the mesh. Use the independent unconditioned
ensemble above to estimate lifetime.

Configure the normal Moon initializer with the positive H reference-weight
ratio 0.125 and a characteristic H speed computed from the frozen incident
energy and the release of internal energy. Giving all available kinetic
energy to H provides an upper bound for this stationary-target fixture.
The time-step callback uses this speed for the usual cell-crossing estimate.
The one-event fixture then installs its common physical time step explicitly.
Its ledger must debit one bath proton per represented atom, remove the Na
packet, publish both daughter packets with the correct statistical corrections,
and continue each for the remaining flight. The audit compares live mesh
particles with analytical ballistic endpoints, isotope and charge budgets,
momentum, and energy. Repeat with each tracer at one and eight ranks. A
build, initialization, or corrupt receipt is an ERROR. An evaluated physical
discrepancy is FAIL. This case scores the registered reaction in the
stated incident channel and has no observational score.

### Optional native timestep and weight checks

The refined native cases share the base PIC beam, source key, stochastic
streams, position, force-free LSO frame, ET origin and 2.625 represented
parent atoms. This experiment tests whether changing numerical scales changes
the same physical realization. Its deliberately selected single event has no
independent statistical confidence interval. Use the separate 32,768-parent
cohort for the survival and angular distribution checks.

The second-order step counts are 1, 2 and 4 at the same horizon. The weight
scales are 1, ½ and ¼ at one step. Both sets run on one and eight MPI ranks.
The reference weights of Na, Na⁺ and H scale together. The parent individual
weight correction compensates before injection and the native transaction
must compute each daughter's correction. Exact packet counts, clock-step
history, represented amounts, signed bath accounting and genealogy are tested.

The independent event lies in the final step. Its exponential-law elapsed time
from injection fixes the absolute ET epoch and the duration of ballistic
continuation. The adapter checks each final position against its birth position
plus its velocity times that duration. The clock allowance is 16 binary64
spacings of the largest absolute epoch or interval. The position allowance is
the component speed times that clock bound plus eight source-coordinate
spacings. The inter-level check budgets the errors of both compared endpoints.
The direct rate, reaction energy and two-product momentum/energy checks use
the existing native commissioning relative bound of `5e-12`.
Comparison values and per-step native rank records are retained for physical
inspection. [Run controls and outputs](README.md#time-step-and-particle-weight-refinement)
