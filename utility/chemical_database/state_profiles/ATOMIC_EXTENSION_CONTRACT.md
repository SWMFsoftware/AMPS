# Opt-in source-bound Na ground-ion profile — preimplementation contract

2026-11-19; original `prepare.py` and partial-profile behavior remain unchanged.
The user approved source-derived atomic ion mass and the explicitly recorded
PHIDRATES/physical-threshold crosswalk. This extension concerns **original
research source use**, not a managed import of held direct-ASD numerical data.
No installed database/lock, canonical raw file, compile species table, application
source selection, uncertainty budget or physical acceptance threshold is changed.

Inventory: the existing preparation authenticates NIST isotope/CODATA originals,
preserves their native uncertainty strings, and produces seven available and
fourteen unavailable states. `PhysicalStateProfile::Registry::Load` accepts
explicit original-research sources with `energy` roles and calls the production
`ReducedKinematics::ResolveState`. Physical mass is separately supplied and never
replaced by the additive reference mass. This existing schema supports an
additional source-bound ion without a native-loader or application-code change.

Actual retained ASD5.12 Na ionization output explicitly supplies the neutral
ground [Ne].3s2S1/2 first ionization value and uncertainty. The retained handbook
on the SAME element page supplies NaI/NaII electronic ground labels and23Na100%
abundance. The modern isotope file provides23Na mass/uncertainty; old handbook
masses/rounded energies are NOT substituted. The ASD NaII second-ionization
energy is not the neutral threshold. No other ion/isotope/level is inferred.

New opt-in `prepare_atomic.py` authenticates these exact source files, the P02
decision register and D04 provenance. A thin C++ adapter calls the actual shared
`ReducedModels::IonPhysicalMass`; independent Decimal arithmetic verifies its
physical/reference/binding masses and converts source uncertainties. The
covariance-independent first-order sigma upper bound is explicitly distinguished
from an independent-RSS assumption, covariance estimate or accepted physical
error gate. Exact SI e/c source flags are checked. Every source record and
conversion/derivation identity is retained; generated receipt is a source-bound
derivation, not a measured ion-mass datum.

Outputs: the existing profile schema, adding `23Na+:ground` (NaII1S0) while
retaining all old unresolved IDs/reasons; exact source-bound CSV/JSON identities;
the derivation receipt; and `srcMoon-M3-Na-source-energy-crosswalk-v1`. The latter
preserves native D04 rates/means and threshold, records physical IE and the
signed energy correction for E_gamma=I_source+mean_source. No undocumented
threshold rounding interval, altered rate, inferred spectrum, numerical gate
relaxation or bounded absolute-rate bias is asserted. Existing strict source
guard behavior remains the application owner's responsibility until its explicit
crosswalk support is implemented and verified.

Preflight all source/hash/state/unit/P02/D04 checks before creating the profile
directory; failed preparation must not leave an apparently usable profile.
Optional work directory retains compiler/source/executable/log identities. Missing
sources/dependency=>SKIPPED77; corruption/parser/compiler=>ERROR2; evaluated
independent arithmetic mismatch=>FAIL1. No source file is copied into a managed
archive. No application or full AMPS build occurs.

Independent tests must compile/call the existing native profile loader and
shared mass kernel, compare source quantities and derived masses using Decimal,
retain all fourteen prior unavailable states, exercise source tamper/no-output
refusal, preserve old partial-profile identity and explicitly identify absent
formal M3 predecessors/other physical channels. A local supporting PASS is never
a linked-runtime or observational PASS. Source/data roadmap hashes remain
712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17 and
d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4.
