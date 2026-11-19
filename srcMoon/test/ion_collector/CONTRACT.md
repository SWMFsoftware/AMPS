# P08 ideal ion collector implementation contract

This implements the user-approved P08 mathematical boundary. It supplies no
lunar ion neutralization, implantation, heating, sputtering, accommodation,
residence or re-emission law. It does not qualify a physical lunar surface.

Inventory: actual sphere callback registration is in main_lib.cpp; static and
moving resolved-terrain adapters call the generic Exosphere accommodation
kernel. That kernel queries neutral temperature/sticking before any ion-specific
production boundary. M2 rejects Na+/Ar+ there; legacy has neutral aliases for
those ions. New Moon dispatch must precede that entire neutral path.

Production API: MoonSurfaceEnvironment delegates ions to an optional explicitly
installed collector/state resolver. A neutral returns not-handled and follows
its original callback. An ion with the default unsupported selector fails
closed. ideal-collector requires the first/second physical-time tracer, an
actual finite impact epoch, initialized state resolver and a process-local
ledger. Legacy impact epochs cannot be inferred from a global source epoch.
The PIC-owning mover deletes a handled particle exactly once; collector code
never frees/reinserts it and never mutates a neutral inventory.

Pure production kernel: species/IonCollector validates a state with explicit
integer nuclei/charge, compiled species, positive model mass, separately
supplied formation/excitation energies and their reference/provenance. The
impact has actual represented number N=W_node,s*individualCorrection,
position/velocity/normal/surface velocity, physical epoch and geometry identity.
It records species/state material, charge, momentum and kinetic+formation+
excitation energy in the actual facet rest frame (axes labelled LSO), together
with LSO snapshot diagnostics. None is labelled deposited heat. Missing
energies/states/provenance or invalid values cause transactional refusal.

The Moon adapter checks exact compiled mass/charge/species and passes actual
collision-frame transforms. Config opt-in alone cannot invent a state provider.
The optional resolver may read preserved reaction state; initial ions need an
explicit source-state resolver. Multiple states sharing one compiled species
must never be guessed. No unresolved full physical production is admitted.

Existing IDs remain unchanged. Local analytical controls extend legitimate
U24 material-boundary bookkeeping/support; they are not observational evidence.
Independent oracle: direct N*nuclear multiplicity, N*z*e, N*m*v_relative,
and N*(m*v_relative²/2+formation+excitation). Check moving-facet invariance,
separate internal energy, exact identity retention, duplicate-state/species
refusal, strict/default/neutral paths, missing data and overflow refusal,
transactional unchanged ledger on failure. Parser tests compile actual
MoonInput and preserve all existing neutral/default controls.

Required artifacts: standalone compiler/log, result.json, source hashes,
analytical comparison CSV and parser outputs under source-tree test_output.
No external data or reaction/source assumptions are used by the analytical
fixture. Linked native callback coverage belongs to the parent's fresh matched
producer build and must not be claimed from source inspection/local tests.
