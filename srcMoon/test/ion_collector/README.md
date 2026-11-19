# P08 optional ideal ion collector

The production implementation is in `srcMoon/species/IonCollector.{h,cpp}`;
`MoonSurfaceEnvironment` adapts actual PIC impacts. This is an **artificial
absorbing boundary for verification and first-impact diagnostics**. It is not
lunar sticking, neutralization, implantation, sputtering, heating, residence or
re-emission physics. No physical daughter/source provider is enabled by it.

## Runtime configuration and production use

In the CLI `.in` application's Moon section:

```text
#begin application: moon
! Other complete Moon settings remain required.
ion_surface_boundary = unsupported
#end application
```

`unsupported` is the default when the optional setting is omitted. All actual
charged impacts must be dispatched before neutral temperature/sticking; the
default rejects them. Neutral particles follow the original callback unchanged.
To select only the mathematical boundary, explicitly request:

```text
ion_surface_boundary = ideal-collector
```

This requires `particle_tracing_mover` to be the explicitly configured first-
or second-order force-cut tracer. A source/global SPICE epoch cannot stand
in for an actual particle impact epoch. An opt-in alone is **insufficient**:
normal production must also call
`SurfaceEnvironment::InstallProductionIonCollector(states,resolver,user,ledger,error)`.
The installer is a nonfatal transactional bool API, requires an empty ledger,
and rejects missing data, repeat installation, non-MPI-only execution, wrong
compiled species/mass/charge, or mixed energy-reference zeros. Sphere mode
requires no initialized terrain service or SPICE call in this installer.

Each `StateDefinition` contains the actual integer nuclear identity, charge,
compiled species and model mass, plus explicitly supplied formation and
excitation energies [J/physical particle], a named common energy zero and
source evidence identity. NaN defaults refuse missing values. A declared zero
is allowed; no zero or ground state is inferred. `sourceEvidenceIdentity`
labels caller-supplied provenance; this mathematical API does not certify its
physical authority. Physical-source/state qualification remains separate.

The explicit resolver receives the actual particle and owning PIC node and
returns the catalog state index. It must use the actual source/product state,
including preserved reaction-state tags when relevant; it cannot infer an
excited/ground state solely from a shared compiled species. Different states
may share a compiled species only because the resolver distinguishes them.
The state catalog is copied immutably on installation. Caller-owned resolver
context and ledger must outlive all impact callbacks.

`HandleProductionIonImpact(particle,x,v,node,frame)` returns false for a neutral;
returns true only after a charged transfer is recorded; otherwise it terminates
with `exit(__LINE__,__FILE__,msg)`. The caller then returns the mover's deleted-
particle result. The collector never frees/reinserts a particle: the mover
owns deletion exactly once and normal cell-list repair. Missing provider,
source-state, frame, impact epoch, or allocation resources is not silent loss.
No callback performs an MPI collective. The parent application owns rank-wise
export/aggregation of `Ledger::records`; no implicit global sum or steady-flux
W/dt normalization is applied.

## What the ledger means

For the actual represented physical number
`N = node.block.GetLocalParticleWeight(species) * individualCorrection`, record:

```text
nuclei[a] = N * integerMultiplicity[a]
charge = N * z * elementaryChargeC
momentum_facet = N * modelMassKg * (velocityLSO - surfaceVelocityLSO)
kinetic_facet = N * modelMassKg * |velocityLSO - surfaceVelocityLSO|² / 2
formation = N * state.formationEnergyJ
excitation = N * state.excitationEnergyJ
total_facet = kinetic_facet + formation + excitation
```

Signed formation energies remain signed. All states in an installed inventory
and every appended record use the same named energy zero. These terms are an
incident transfer to an artificial collector. **Total energy is not deposited
heat** and charge transfer is not a neutralization law. Material rest energy
mc² is not added to this nonrelativistic reference-mass/state-energy ledger.
The current mechanical fixture does not bound the physical approximation of a
new isotope/reference-mass model; that remains a separate model-qualification
requirement before physical kinetics are activated.

`ImpactFrame` retains the selected outward normal, actual material surface
velocity, LSO/body state transforms, geometry identity and geometry epoch.
Normal must be nonzero, velocity non-outgoing in the material rest frame, and
all inputs finite. Transform inverse-pair validation uses an arithmetic
forward-error bound, not a tunable physical acceptance tolerance. The physical
impact epoch is the actual tracer epoch, while frozen-environment geometry
retains its separately declared frozen epoch. LSO snapshot momentum/kinetic
energy are also retained with their frame label; moving-facet conservation
uses relative velocity. Nothing labels frozen geometry as an evolving pose.

## Run local kernel and parser controls

From the repository root:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_ion_collector.py \
  --output-dir test_output/srcMoon/ion-collector-LOCAL-NEW
```

Use a new output directory for every attempt. The runner compiles maintained
`IonCollector.cpp`, `ReactionExecution.cpp`, `MoonInput.cpp` and its real rate
registry dependencies with a normal C++17 compiler. It does not configure or
recompile AMPS, change generated files, load a physical source catalog, or
write canonical data. Kernel controls need only a C++17 compiler and Python's
standard library. Parser cases use the committed bounded legacy/common-time
examples, whose required SPICE/LOLA file paths must exist. Missing example path
prerequisites give explicit parser SKIPPED results, not a false negative-case
PASS. These path checks do not qualify SPICE or LOLA data scientifically.
Matplotlib is optional for PNG/EPS; if unavailable, plotting is explicitly
SKIPPED and the numerical comparison CSV remains available.

Artifacts include `result.json`, exact compiler command and stdout/stderr,
source revision/hashes, executable hash, fixed named native control results,
`production_budget.csv`, independent `comparison.csv`, and parser input/logs.
When plotting is available, `comparison.png` and `comparison.eps` display the
analytical bookkeeping differences. They are not observations or ion-surface
validation plots. Parent-linked native callbacks require a fresh matched AMPS
producer build and their own 1/8-MPI receipts; local source inspection cannot
substitute for that execution.

## Test specification

| Required item | Contract |
|---|---|
| ID/name | Existing U24 material-boundary bookkeeping/support extension; no new permanent ID. |
| Purpose | Verify conservative accounting and fail-closed ion-boundary selection. |
| Production capability | `IonCollector::{ValidateStates,Prepare,Collect}` and actual Moon input parser. |
| Production entry points | Parent sphere, static terrain and moving terrain wrappers call `HandleProductionIonImpact` before neutral callbacks; local tests exercise the shared collector/parser kernels only. |
| Configuration | Default unsupported; explicit ideal-collector with first/second tracing selector; MPI-only native installer. |
| Inputs | Declared analytical A2B ion, N=2.5, mass4kg, charge2, explicit energy zero/formation−7J/excitation3J, moving planar surface; numerical scales are not real lunar species parameters. |
| Oracle | Independent direct number, integer multiplicity, signed charge, momentum and energy equations above, with exactly representable binary fixture values. |
| Units | SI; J per physical state versus represented totals distinguished. |
| Frame/time | LSO axes and material relative velocity; explicit impact101.25 and frozen geometry100 TDB seconds in fixture. |
| Deterministic/stochastic | No RNG, particles or observations in local kernel test. Native parent controls use actual PIC packets. |
| Procedure | Compile actual production services; fixed named positive/refusal controls; compare thirteen independent budget quantities; actual parser tests on complete inputs. |
| Metrics | Every named control, retained identity/context, thirteen exact budget comparisons, transactional unchanged ledger/output on failure, parser status and preserved previous configuration. |
| Acceptance | Exact binary fixture arithmetic and exact native count/name/summary/exit contracts; invalid cases must fail for their intended selector/data reason. No physical thresholds changed. |
| Expected PASS | All required local cases evaluated and met. This is supporting software verification only. |
| Failure modes | Lost species/state/nuclei, W/dt multiplier, incorrect moving-frame velocity, missing energy or wrong zero, neutral law for ions, malformed transform, silent publication on refusal, parser false green. |
| Required data | None for analytical science values. Parser uses example file-path prerequisites only. No authoritative physical ion-state/formation catalog is inferred. |
| Data hashes | No observations; compiler/executable/source/input hashes are retained. Physical catalogs require separate explicit source-bound provenance. |
| Predecessors | This local support does not bypass formal linked/convergence/environment/chemistry predecessor gates. |
| Artifacts | Source/build logs, named checks, actual/analytic CSV, optional PNG/EPS and per-input parser logs. |
| Status semantics | Missing compiler/path dependency SKIPPED; compile/crash/protocol corruption ERROR; evaluated missed arithmetic/guard FAIL; evaluated met software criterion PASS. Formal lunar ion-surface status remains NOT VALIDATED. |

Default neutral regression preservation and actual sphere/terrain deletion,
compiled masses/charges, real statistical weights and MPI ownership must be
verified by the parent's unified native campaign. The mathematical collector
never establishes a lunar ion-recycling capability or closes M3 by itself.
