Oxygen radiative cascade local service
====================================

Run from this AMPS checkout:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_oxygen_cascade.py \
  --database-root /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r3 \
  --output-dir test_output/srcMoon/oxygen-cascade-NEW
```

The output directory must not already exist. The runner compiles the maintained
`srcMoon/species/OxygenCascade.cpp` and actual shared database reader using
C++17; it does not configure or build AMPS. Python SciPy and matplotlib are
required for independent matrix-exponential/binomial checks and PNG/EPS. A
missing dependency or permitted source produces SKIPPED (exit77), a failed
criterion FAIL (1), corruption/build/runtime error ERROR (2), evaluated local
success PASS (0). `result.json`, `run_manifest.json`, native output, compiler
log, comparison CSV/PNG/EPS and exact commands/hashes are retained. The native
probe uses one process and no MPI; it is a local kernel test.

The numerical acceptance and stochastic protocol were fixed in
[CONTRACT.md](CONTRACT.md) before implementation. These results extend existing
local U15/U16/state-cascade coverage without assigning a new test ID. Formal
linked/M3 status and physical activation remain separately SKIPPED/DISABLED.

`Network::LoadFrozenSelectedChiantI` accepts an authenticated shared Catalogue
and pins the actual permitted r3 selected level/link tables and full native
originals. All nine links are retained, including weak O(1D)→3P0 and
O(1S)→3P2 exits. It uses shared `ReducedModels::RadiativeBranches` for coherent
rates, lifetimes and conditional branches. Level1 is stable *within this
selected network*. All endpoints descend, so a history has at most four
photons. `AdvanceRadiativeOnly` consumes supplied open waiting uniforms and
[0,1) branch uniforms. It selects no production RNG or population. Its isolated
radiative-only mathematical scope does not assert absence of lunar collisions,
upward excitation, quenching, trapping, or competing chemistry.

State gaps require explicit observed OR theoretical CHIANTI energy columns
and explicitly supplied SI h/c with source identities. Native wavelengths are
signed vacuum Angstrom diagnostics; they are retained independently, never
substituted as exact state energies. The source supplies no per-entry A/level
uncertainties or O isotope mapping; no values are invented. Original CHIANTI
rights do not clear independently acquired NIST ASD products.

`SolveRestFramePhotonRecoil` separately handles finite-emitter two-body
kinematics. It requires independently sourced physical parent/daughter masses,
state gap, identical explicit h/c policy, isotope/frame IDs, direction and
caller-frozen arithmetic bounds. It checks mass-gap consistency and both state
and mass energy ledgers. In the parent rest frame it emits a photon with
`Eγ=ΔE[1−ΔE/(2Mparent c²)]`, momentum `Eγ/c`, opposite atomic recoil, and stable
relativistic kinetic energy. The test uses labelled mathematical masses,
never an assumed oxygen isotope. Missing profiles do not choose infinite
mass. This generic solver supplies no angular law; its direction is explicit
input. Invalid/domain/unrepresentable inputs refuse transactionally.

A real lunar decay event additionally needs an admitted isotope/state/physical
mass crosswalk, initial state populations, source-supported omissions and
quenching-domain evidence, and linked production wiring/conservation/
convergence qualification. None is approved by a local network PASS. The
five-level subset is downward-closed for the selected native source; it is
not a complete arbitrary O I atom.

Independent references are Decimal source sums, state conversions and exact
finite-recoil conservation, a separate 14-variable continuous-time generator
for populations and accumulated per-link photons, and predeclared exact
binomial marginal acceptance bounds. The plotted comparisons are mathematics,
not experimental observations. Every positive branch also receives a
CDF-interior deterministic check, so small Monte Carlo samples cannot erase
weak published channels. All original and installed database bytes remain
read-only.

The final recoil guard also certifies the *represented Cartesian* photon norm,
massless dispersion, vector-derived atomic energy, velocity reconstruction and
subluminality against the same explicit input bounds. A zero bound has no
hidden numerical epsilon. The maintained mathematical regression exercises the
specific non-axis vector/mass combination that exposed this storage seam;
other exactly representable combinations remain admissible.
