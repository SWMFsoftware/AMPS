# Independent physical event-ledger reconstruction contract

This supporting extension of U14/U15/U16/U17 and linked I21–I24 does not assign
a new test ID or alter qualification order. Its purpose is independent event
bookkeeping, not source authority, daughter-distribution validation, statistical
survival/branch scoring, convergence or observational validation.

Production entries are PhysicalReactionEvents::SodiumPhotoionization::Compose,
Moon::Chemistry::ParticleTransaction, and the actual committed PIC child setter
and readback. Native Boolean conservation flags are never accepted as oracles.

Inputs: authenticated source profile JSON and its CSV/source descriptors;
explicit parent/ion/free-electron state IDs and compiled indices; actual
per-rank event scalar CSV (`original_parent_id,key,value`), with hexadecimal
binary floating values; externally frozen SI numeric bound dictionaries and
explicit event frame. Source values remain exact decimal rationals; native
hexadecimal binary numbers remain exact Fraction values without narrowing
long double to binary64. Missing actual child readback/schema/source evidence
is SKIPPED; malformed, changed, duplicate or nonfinite records are ERROR;
evaluated numerical/conservation mismatches are FAIL.

Independent reconstruction: additive isotope reference masses minus integer
charge times electron mass; source formation+excitation energies and Q;
integer nuclei/charge stoichiometry; actual reference-weight times correction;
per-event and represented mass/momentum/kinetic+internal energy; independently
signed outgoing-electron minus incoming-photon reservoir moments; proper event
rotation/velocity-offset transformation; massless photon and parent-rest
Doppler energy; actual PIC initial child state/position/velocity/weight/internal
energy/interval versus its complete plan. Exact integer and stored-value checks
have no tolerance. Every numerical comparison uses caller-frozen bounds;
there are no verifier-chosen physical/error tolerances or residual corrections.

The verifier exposes verify(paths,profile_json,mapping,bounds). It returns
criterion rows, exact rational residuals/limits, source identities and a
supporting result. A loaded partial source profile with unresolved Na+ is
SKIPPED before event reconstruction. It does not create an ion state, replace
a PHIDRATES threshold, infer isotope populations, omit electron/photon baths,
or overwrite immutable evidence. Existing ensemble and numerical-error gates
(including approved 0.2/0.3 percent bias allocations) remain in their owning
qualification scorer and are not weakened or substituted by event closure.

Artifacts: strict result/run_manifest, exact input hashes, actual CSV comparison,
PNG/EPS diagnostics when evaluated, caller's native logs and data identities.
Independent analytical toy-state fixtures may test software only and are
explicitly labelled; they are never physical profiles or observations.

An authenticated set of valid header-only rank files is a legitimate empty
event-set identity when the owning audit has reconciled all original parents
and zero reacted histories. It yields one explicit empty-set PASS criterion,
`events_evaluated=0` and no physical-event-validation claim. Absent files or no
rank file descriptors remain unavailable; unresolved source states still SKIP
BEFORE this identity. A nonzero event population is never manufactured.

## Unit scale for rotation invariants — November 19 correction

The production frame contract tests `R^T R - I` and `det(R)-1` with dimensionless
scale **one**: columns are declared unit vectors, including when their expected
mutual dot product is zero. The independent verifier must use that same declared
quantity/scale while reconstructing the Gram matrix exactly. Scaling an
off-diagonal residual by its near-zero actual dot product instead would demand
exact orthogonality of rounded SPICE coefficients and does not implement the
contract. No configured absolute/relative value is changed by correcting this
scale. The initial real-event FAIL report is retained. Analytical controls must
admit a small off-diagonal perturbation within the explicit unit-scale bound,
reject one above it, and preserve exact refusal when both bounds are zero.
