# Oxygen five-level radiative network: contract before implementation

Scope: local supporting extension of existing U15/U16/state-cascade coverage;
no new permanent test ID. This component uses the actual production
`OxygenCascade` service and shared immutable `Catalogue`/`RadiativeBranches`.
It selects no lunar population, isotope, quenching rate or physical default.

The source is the permitted CHIANTI O I native seven-level/sixteen-entry
collection distributed with reviewed shared release r3. The selected lowest
five levels and **all nine** downward links are bound to table SHA-256
`8b0bfcbae54365ed668efa35f337e7b9ed6a60a55d0b4fbfe94693e22660c192`
and `bf4e17705e69cbaebbabd07f624e6034c5d96817b10aeb51ee8a77c96a21ed4f`.
Original level/rate file hashes are
`230bbdcd95c269fc0eea04a1fd37f534c1a42979d2e29c1b346e49dd71c688da`
and `f24efc2bd0fa5120baca6dce06ffb7831ae2a59e17cc186c233a43b5145866a8`.
Both weak excited exits and all three ground fine-structure links remain.
Complete selected-source closure does not mean a complete arbitrary O I atom.

Production entry points will load the source transactionally through the
actual shared `Catalogue`; query coherent lifetimes/branches; select a decay;
compute observed **or** theoretical state-gap metadata with explicit h/c;
and follow a supplied-random-number, radiative-only finite history. Histories
are mathematical radiative-network evolution: no unmodelled collisions,
excitation, radiation trapping or isotope interpretation is declared absent
from the actual Moon by that operation. Physical use needs an independently
admitted state/isotope crosswalk and omission/quenching-domain evidence.
No per-entry uncertainties are invented. Native signed vacuum-Angstrom
wavelengths remain diagnostics, separate from cm^-1 state energies.

A state gap is not automatically an exact emitted-photon energy. A separate
explicit-mass **parent-rest-frame** mathematical solver uses physical parent
and daughter masses, h/c, declared state-gap energy and caller-supplied
arithmetic bounds. It checks `(Mparent-Mdaughter)c²` against the state gap,
then uses `Egamma=DeltaE*(1-DeltaE/(2*Mparent*c²))`; emitted photon momentum
and opposite atomic recoil are both retained, with stable relativistic atomic
kinetic energy. Missing mass/source profiles do not fall back to infinite
mass. Inputs/outputs identify the selected source levels, energy basis,
frame and independently supplied mass/constants identities. This is generic
math verification until real isotope/state/omission evidence is admitted.

Independent oracles: Decimal arithmetic for source A sums/branches, energy
conversions and recoil/conservation; independent continuous-time generator
matrix exponential for populations and link photon counts; exact binomial
CDF bounds for the actual native supplied-uniform histories. Stochastic
protocol is fixed before execution: eight families (initial levels4/5 and
four durations1,100,1e5,1e10 s), N=20000 per family, seed1846053 with a separate
stream per family; 14 simultaneous population/link-count claims per family,
familywise alpha0.01/112. Every positive source branch additionally has a
deterministic CDF-interval query, so rare weak exits are never inferred absent
from finite sampling. Source/model counts and energies are not observations.
Source arithmetic comparisons use128 binary64 eps relative after JSON
conversion; no new physical uncertainty or observational acceptance gate.

Required dependencies: C++17 compiler, reviewed r3 shared catalogue, Python
SciPy for independent matrix/binomial controls and matplotlib for CSV/PNG/EPS.
Missing dependencies/source are SKIPPED; wrong hashes/metadata, compilation,
runtime or contradictory receipts are ERROR. An evaluated criterion mismatch
is FAIL. PASS covers only this local source/network/mathematical contract;
formal U/I/M3 qualification and physical activation remain separate.
All outputs/logs/commands/build/source/data identities are preserved under a
new caller-supplied output directory. No canonical/installed database write,
application build/configuration or reference/threshold change is permitted.

Representation guard added from independent review before final build:
Cartesian direction normalization/multiplication must also satisfy the caller's
existing momentum and energy bounds, including `E_gamma=c*|p_vector|` and
vector-derived atomic KE. No epsilon is supplied for zero bounds. A concrete
Mi=2,Mf=1,c=3 mathematical non-axis vector previously produced a represented
norm2.25000000000000000021684 against scalar2.25; the added control requires a
transactional NumericalFailure at zero bounds. This guard checks storage
arithmetic only and selects no physical O mass or angular distribution.
