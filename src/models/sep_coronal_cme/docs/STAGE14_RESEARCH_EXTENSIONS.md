# Stage 14 research implementation and verification

This is an implementation guide for Sections 6.5, 15.8, 16.14 and 17.2 of the
specification. It does not amend their acceptance criteria. Stage 13 was added
and verified before this research layer. Stage 14 is optional and independent
of the baseline release profile. Application schema 5 and ordinary field-line
bundle major 3 retain their existing selectors and interpretation; the new
research records use schema 6 and a momentum-dependent family uses bundle
major 4. Current AMPS input cannot silently select these offline kernels.

## Commands and evidence ownership

```sh
make -C src/models/sep_coronal_cme test-stage13
make -C src/models/sep_coronal_cme test-stage14
python3 src/models/sep_coronal_cme/test/run_tests.py --test SRC3D22
python3 src/models/sep_coronal_cme/test/run_tests.py --all
```

The cumulative registries contain 209 cases through Stage 13 and 222 through
Stage 14. Thirteen new canonical Stage-14 IDs are registered in
`test/run_tests.py`; `test/test_stage14.py` invokes the public C++ kernels and
Python producers for those IDs. `make test`, individual launchers, and the
existing aggregate SEP+corona driver use the same registry. No IDs need to be
listed on the aggregate command line. With the current seven native host
initialization cases, a complete aggregate has 229 records. The native seven
remain generic Parker/SWCME host evidence, separate from shared model tests.

Shared JSON explicitly labels its results `software-verification`, with
production and observational qualification false. EVT3D01, XMD3D01 and
SLM3D01 below verify campaign **protocols on manufactured data**. They are not
reports of a completed observed event. The actual observation products and
qualified host runs must be supplied separately. An unsupported future branch
is rejected before construction; it is never replaced by a plausible zero.

## Implemented domains

| Capability / gate | Implementation and independent checks | Required before broader qualification |
|---|---|---|
| Nonradial winding / CPL3D10 | `preprocessing/nonradial.py` integrates the steady corotating spatial characteristic, differentiates its deformation and applies the Piola transform. Tests use an independently known SO(3) characteristic, inversion, orientation, divergence, area-flux identity, pushed interface normals, focusing gradients, mass-per-flux and radial reduction. | A native background adapter, event field/flow reconstruction, certified global exclusion census and convergence across complete event support. |
| Foreshock coefficient proxy / MFP3D08 | `research_extensions.cpp` provides a compact C2 factor in `(0,1]`, current-front binding, resolved-upstream support and exact ambient recovery. Line arclength and independent volume distance agree at matched plane-front points. Wave energy is untouched. | A qualified moving-front host distance authority and production mover integration. This proxy is not a solved particle-wave instability. |
| Resonant waves / TUR3D08 | `preprocessing/research.py` evolves a 1-D directional log-wavenumber spectrum with finite-volume advection, refraction, cascade, damping, focusing work and explicitly coupled particle energy/momentum reservoirs. Constant spectra give an independently integrated mean free path; periodic advection converges under refinement. Replay bindings and unsupported resonance cells are checked. | Full kinetic particle-distribution coupling, open shock boundaries, moving wave-frame qualification and multidimensional/native-MPI convergence. The implemented two-moment closure is a named small-amplitude approximation. |
| Streaming-limit comparison / SLM3D01 | Frozen response folding, registered radial power-law sensitivity, units/window/support checks, log-ratio covariance, all requested strata and explicit zero/censored/unsupported states. | Independently prepared streaming-limit observations and preregistered full runtime strata; no feedback-closure or universal agreement threshold follows from this diagnostic. |
| Coherent drift / FTE3D10–11 | Neutral `sep_common/sep_coherent_transport.*` owns separate antisymmetric Parker and full focused Hamiltonian operators in a stationary inertial frame. Independent Lorentz full-orbit references include dimensional proton/alpha examples. Charge/polarity, magnetic moment, electric work, timestep refinement, invalid-region dispatch and disabled selection are verified. | Time-dependent/plasma-frame Hamiltonian terms and integration into each production mover. The callback-disabled test verifies kernel selection, not a complete AMPS trajectory. |
| Dynamic geometry / ELL3D11 | `DynamicEllipsoid` accepts fixed-axis Rodrigues or a full proper body-to-inertial SO(3) matrix with inertial angular velocity. Translation, angular motion and axis expansion all enter surface velocity and normal speed. Independent finite differences include noncommuting rotations. | Multi-epoch observational reconstruction and covariance ensembles propagated through shocks, source rates, connectivity and observers. A supplied R/omega pair is not itself an observational fit. |
| Impulsive attribution / SRC3D20 | A normalized open-footpoint/time/momentum/pitch/species product gives physical births and separate origin/interaction number and kinetic-energy ledgers. Current-front upstream or typed pre-front births are accepted; downstream births are rejected. Overtake absorption debits once. | Native source superposition, host mapping/response ensembles and mixed-run attribution. Transfer/reacceleration needs its independently validated downstream/renewal authority. |
| Reference family / SRC3D21 | Positive normal-coordinate inflow/diffusivity is integrated and root solved under one common Pe target. A distinct geometry/joint-measure authority is required; member gaps, generation mismatches and uncertified certificates reject bundle-major-4 publication. Restart/roundtrip and legacy-reader rejection are tested. | Event geometry reach/overlap/clearance certification and rederived joint source measure, plus native family-aware export/readers. No schema-3 column insertion is allowed. |
| Conditional renewal / SRC3D22 | Conditional tabulated branches normalize and close number, kinetic/total energy, momentum and shock work. Frozen checkpoints retain first-passage identity, cohort ancestry and cycle history; manufactured splitting/restart is deterministic. Unchanged or zero-delay re-emission and original-spectrum redraw are rejected. | Independently validated physical conditional kernels, stochastic host sampling and MPI ownership/restart tests. Portable weighted-cohort splitting is not an MPI test. |
| Wind envelope / WND3D21 | Paired speed and quasi-steady advective-acceleration channels share frame/support/selector. Outward-rounded interval branch-and-bound certifies continuous polynomial extrema. An endpoint-valid but interior-invalid profile fails and reports rejected flux. | Qualified observational channel inference/covariance and extension of the production D7 product. Time-dependent material acceleration requires its own authority. |
| Event transfer / EVT3D01 | Frozen event-specific inference inputs produce a fresh run fingerprint; all seven metric families retain failed values. Compound 2020 May 27–June 2 records count as one event; 2012 May 17 is stress-labeled without a validated preceding cloud. | Real independent event runs, registered observations and complete response-derived metrics. The CLI evaluates already prepared metric values; it does not infer onset from an arbitrary time series. |
| Offline MHD / XMD3D01 | Matched sample manifests compare all required plasma/field/flow/flux/D2 variables, uncertainties, topology and connectivity. Unreported boundary mismatch, unsupported samples and residual tuning reject; permitted mismatch is explicitly combined boundary-plus-model discrepancy. | Independently supplied thermodynamic MHD samples. This is neither a runtime imported-MHD provider nor truth validation. |

These remaining gates are part of the specification. Their absence is a known
limitation, not a changed acceptance limit or an implicit release PASS.

## Numerical and lifecycle details

`CapabilityIdentity` binds algorithm/coefficient identity, background and front
generations. The JSON configuration contains a complete set of independent
boolean capability flags; each enabled flag has a version, gate, domain and
claim authority. It explicitly reports no qualified production adapter.

The characteristic winding producer solves `dF/dr=u'(F)/(u'·e_r)` from an
identity inner sphere. It checks positive radial projection and the radius
constraint; `B=(DF/J) B0` requires positive orientation. At the inner sphere,
tangential derivatives are the identity and the radial derivative comes from
the flow. A complete one-sided physical vector trace is compared against a
separate boundary authority. The inverse is bounded Newton; failed inversion
or turning rejects publication. Cartesian finite differences are used only
inside support; a separate analytic boundary trace prevents hidden centered
extrapolation. Callback field/flow authorities must be independently frozen
by their owning producer. The current implementation is offline, not an AMPS
mesh/background lifecycle provider.

For the wave step, W is total Alfvén energy per physical volume per `log|k|`;
its directional pseudomomentum is `+/-W/v_A`. Physical volume is `ds*area`, and
spectral integration includes `dlogk`. Growth takes its wave gain from particle
energy and parallel momentum. Damping heat/momentum, spectral boundary escape
and background focusing work have separate ledgers. Conservation uses an
absolute momentum inventory so counterpropagating cancellation is not an
impossible zero tolerance. Explicit positivity/CFL and reservoir positivity
are enforced. Spatial boundaries are periodic and the shock is frozen over a
coupling step; replacing either approximation requires another validated
boundary/frame authority. Resonance is broadened, interpolation is bounded,
and there is no uncovered-cell fallback. The wave asset includes all replay
bindings, including source/distribution/species/return/cadence for coupled
iteration. External prescribed waves instead type those fields inapplicable.

The focused deterministic characteristic derives streaming, mirror force,
grad-B/curvature motion and electric work from one relativistic guiding-center
Hamiltonian. It supplies the entire deterministic update. The Parker route
instead supplies its antisymmetric diffusion operator, with drift velocity
as a diagnostic. A host must select the correct ownership; adding either
velocity to an already complete focused update would double count. HCS,
separatrix, null, transition, shock, weak-field and magnetization/step-ordering
violations return typed failure/dispatch. A moving-frame snapshot is currently
unsupported. Drift reporting retains vector net displacement and noncancelling
path exposure, and labels tube-radius, one-axis and plane RMS denominators;
those widths must come from a separately registered drift-disabled diffusion
experiment. Cartesian components are signed projections in the declared frame.

For geometry, R rotates body coordinates into the inertial frame. The shape
tensor Q is `R diag(a^-2) R^T`; R and Q are not interchangeable. Angular
velocity, center velocity and axis rates enter `-F_t/|grad F|`. Matrix
orthogonality and positive orientation reject componentwise attitude mixing.
The caller must reconstruct consistent history derivatives; the supplied
manufactured histories are not an event inversion.

Renewal work uses a stable relativistic kinetic-energy expression, avoiding
subtraction of nearly equal rest energies. Work tolerance is scaled by kinetic
energy rather than rest energy. Absorbed/downstream reservoirs retain their
four-momentum. The first-passage source stays immutable across cycles; only a
conditional outgoing state creates a continuation after positive residence.
Tabulated branches are physical weighted measures rather than fresh draws
from shock `g(p)`. Event validation of that conditional authority remains
required. Impulsive origin remains distinct from shock first passage even
when a front interaction occurs.

Interval wind certificates use directed binary64 arithmetic, including an
IEEE bit implementation on Python 3.7/3.8. Monotone cells use endpoint
bounds; stationary cells are subdivided until their enclosure meets the
registered SI tolerance. Uncertifiable bounds fail instead of being sampled,
clipped or widened. The polynomial-interpolant family is explicit; an arbitrary
spline must have its own certificate. No result changes the schema-5 D7 gates.

## Offline CLI and reproducible examples

`tools/research_stage14.py` reads an explicit `sep-stage14-offline-job-v1` JSON
job and writes a checksummed result that binds the entire job. Existing output
paths are rejected. See [the examples](../examples/stage14/README.md). The
C++ providers and callback-based winding producer have public APIs; the JSON
CLI cannot invent their missing host integration. Every example is synthetic.

## Canonical campaign disposition

The three reserved campaign IDs EVT3D01, XMD3D01 and SLM3D01 report **SKIP**,
with `verification_passed=true` after their synthetic contract checks. The
actual campaign acceptance/input-loading gates remain open; the portable
contract checks cannot close them. JSON/JUnit and the aggregate preserve this
disposition instead of treating exit zero as campaign PASS. The default full
shared run selects 222 records: 219 PASS, 3 SKIP, 0 FAIL. The Stage-13 subset
is 209 PASS without these research dependencies. `--require-no-skips` returns
nonzero when actual campaign evidence is incomplete.
