# Source-bound normal-compile species recipe, phase-two staging

This is a maintained supporting source-bound compile-recipe package. It does not configure
or compile an AMPS application, activate a physical process, change constants,
or qualify a missing source-state mapping. The current partial real profile
has no admissible Na+ state; requests including it must refuse BEFORE the
requested output directory exists.

The generator authenticates the complete supplied profile through the SAME
production PhysicalStateProfile::Registry and ReducedKinematics resolver.
Caller supplies an explicit compiled symbol -> resolved state mapping; names
never select an isotope/state automatically. Every selected state must resolve
and have an independently sourced physical rest mass. Compiled mass is the
resolver's additive reference mass, rounded once into the native PIC binary64
representation. Physical rest mass remains separately recorded, never used as
an unapproved compiled inertia replacement.

Use the existing Config.pl -application path with an application-owned compile
input and species include. ampsConfig.pl::ReadSpeciesBlock copies a numeric
mass token and emits a numeric charge token multiplied by legacy ElectronCharge.
The source-SI recipe therefore chooses a numeric scale token and independently
compiles its EXACT emitted expression to require equality to double(z)*source_e.
An unrepresentable expression refuses. This is a conversion from a legacy
code-unit token to coulombs, NOT fractional physical ion charge. Integer charge
continues to belong to the reaction catalog. No global constant is patched.

Acceptance: all requested profile and source hashes authenticate, exact schema
and unique symbols/states, finite positive rounded masses, finite charge token,
actual compiled expression identity, exactly one SpeciesList and species include
replacement, no output on missing/invalid source, existing output refused.
Tests use a real available neutral profile and real missing-ion refusal plus
corrupt-source refusal. No fabricated ion/profile is called physical evidence.

Phase-two Moon-only integration proposal: add a source-bound explicit charge
unit policy with legacy default to Controller::Configuration, include it in
configuration fingerprint, propagate it into every ParticleTransaction and
collector installation/impact. All compiled mass/charge checks and collector C
ledgers use that same installed unit. Reject changed/nonfinite/nonpositive unit,
empty source identity or source hash. Legacy overloads retain ElectronCharge.
Neither eV2J nor global ElectronCharge/ElectronMass nor input/species.input is
modified. Formal physical startup still refuses the unavailable Na+ state,
source-threshold crosswalk and source-error budgets.

## Explicit homogeneous reference experiment selection

The optional named reference_problem=homogeneous-zero-acceleration is an
approved reference experiment, not a physical-model fallback. Default
normal-template preserves the template force definition byte-for-byte and adds
no cohort marker. The named mode requires exactly one simple six-argument
acceleration DEFINE inside the normal Exosphere input block and no existing
cohort marker. It substitutes the EXISTING production core zero-acceleration
function and appends _MOON_COHORT_ZERO_FORCE_REFERENCE_=1 through ordinary
compile-input DEFINE directives. Ambiguous/missing definitions refuse before
output; no generated build/header/CPP source is patched and no new force law
is implemented. All native source-profile/state/charge checks still precede
output; missing Na+ remains unavailable.

Independent controls run the maintained Exosphere Perl parser against isolated
source copies, then compile its generated macro expansion and marker. These
prove normal parser/compiler binding only, not physical dynamics or source
readiness. The normal default force line and absence of a newly added marker
are separately asserted. Existing12-control results remain their old basis.
