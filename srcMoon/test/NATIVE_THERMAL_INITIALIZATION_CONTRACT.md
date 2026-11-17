# Native applied initial thermal state — I17 / I19 / I20

## Inventory and purpose

Before this change, `InitializeProductionState` allocated every predictive
thermal column at `thermal_initial_temperature_k`. I17/I20 explicitly called
the existing observation-free `SpinUpProductionState`, but I19 numerical cases
never applied a separately qualified periodic state. A proof about one state
cannot qualify cases initialized in another. This contract is written before
the dispatch/binding change and introduces no new force, material law, forcing
interval, acceptance threshold, observation, or permanent test ID.

## Explicit production configuration

`thermal_initialization_mode = declared-uniform` is the optional default and
preserves the existing uniform diagnostic initialization. It is eligible only
when an explicit qualification proof covers that very state. The other named
mode, `periodic-spinup`, requires `surface_temperature_mode = thermal-model`.
It invokes the existing production periodic solver before initial source-rate
or weight queries. Both UTC cycle endpoints and both convergence tolerances
remain required explicit physical-input fields. The cycle end and simulation
start must convert to exactly the same SPICE TDB epoch; no interval is selected
or adjusted by the implementation. No external thermal-state loader or cache
is introduced.

The complete M2 numerical verification group may optionally request
`m2_verification_record_initial_thermal_state = true` (default false). It is
not a physics switch and does not permit the runner to change physical inputs.

## Applied-state identity and reuse

Production retains the actual initialization method, epoch [TDB seconds], and
completed periodic convergence diagnostics. An explicit I17/I20 spin-up call
reuses a completed state only while that same state has not evolved. A later
positive-duration advance invalidates completed-state reuse and initial-state
export even if a caller returns to the same numerical epoch. Existing explicit
spin-up from the default uniform diagnostic state still runs and records its
actual periodic method. Initial provenance remains readable after evolution;
export requires the actual initial columns still to be present.

`GetProductionThermalInitialization` reports the applied initial state.
`WriteProductionInitialThermalStateCsv` exports the actual production columns
without modifying them. It refuses unsupported/nonthermal/uninitialized,
evolved, invalid or incomplete state and an already existing output path.
Its exact schema is:

```text
logical_surface_element,depth_index,depth_m,kinetic_temperature_K
```

Elements are ordered, each with boundary row depth_index=-1/depth_m=0 followed
by every finite-volume centre in index order. Values use 17 significant digits,
SI metres and kinetic kelvin. This is model-state evidence, never observations
or an independently generated physical truth. No file is overwritten.

## Linked receipts and independent binding

Before I19 cohort seeding or numerical dt/weight overrides, the native case
records the actual method/epoch and spin-up diagnostics. If requested, rank
zero writes `<receipt>.I19_initial_thermal_state.csv` and records its absolute
path, row count, logical-element count and layer count. I20 likewise records
the actual production initialization, not a presumed Python pilot method.
Python qualification must bind these actual values/bytes to the explicitly
reviewed proof and immutable physical-input signature. Mismatches are ERROR;
missing proof/prerequisites are SKIPPED, never automatic approval.

## Verification and status limits

The production solver and native initialization path must be rebuilt and
exercised before claiming native support PASS. Required checks distinguish
uniform default, applied periodic state, same-state reuse, evolution/export
refusal, exact full-column export and mismatch refusal. Existing U11 thermal
oracles, I17 independent spin-up tolerances and I19 predecessor/equivalence
requirements remain unchanged. Short supporting fixtures cannot establish
formal convergence, stationarity, physical feasibility or observational
validation. The current private campaign remains frozen; source changes need a
subsequent clean native build and retained build/source identity.

The native receipt names are `thermal_initialization_mode`,
`thermal_initialization_et_s`,
`thermal_initialization_spinup_completed_cycles`,
`thermal_initialization_spinup_maximum_node_temperature_change_k`,
`thermal_initialization_spinup_maximum_column_enthalpy_change_j_m2`,
`thermal_initialization_spinup_cycle_start_et_s` and
`thermal_initialization_spinup_cycle_end_et_s`. A uniform state has no spin-up
cycles (zero diagnostics, not an assertion of periodic convergence). Requested
exports add `initial_thermal_state_csv`, `initial_thermal_state_rows`,
`initial_thermal_state_element_count` and `initial_thermal_state_layer_count`.
I19 additionally refuses a case whose actual applied initialization epoch is
not its production particle start epoch before any cohort is seeded.

The exporter uses an exclusive create operation, rather than a check followed
by truncating open. An I/O failure preserves its partial artifact and asks for
a new campaign path; it never retries by overwriting evidence. Initial-state
provenance APIs leave caller outputs unchanged on failure. Export metadata is
computed from actual production profile arrays, not a separately assumed
layer count.

## U11 initialization-dispatch fixture (specified before probe changes)

The new `thermal-initialization` probe uses the existing production 40-layer
material configuration, a 2x4 logical inventory grid, a closed octahedral
terrain and the actually furnished SPICE kernels. It requires no AMR mesh or
MPI particle run. The uniform control initializes exactly 250 K and exports
all 8*(40+1) boundary/centre rows. It independently checks the exact 250 K
values, expected depth/index ordering and the 1/40 m uniform-centre geometry.
An existing export path must be refused without altering the first file.

A separate dispatch-only control explicitly declares the 1 ms UTC forcing
interval 2009-01-24T00:00:00.000 to 2009-01-24T00:00:00.001, sets simulation
start equal to its end, selects automatic `periodic-spinup`, and supplies
the unchanged approved 0.01 K/10000 J m^-2 cycle criteria. Both actual completed
metrics must meet those limits; no exact cycle count is assumed. The short
1 ms replay remains a dispatch fixture rather than an initial-state
qualification for a lunar forcing period. The probe requires
actual periodic metadata/cycles, source-BEGIN ET equal to the declared end,
and changed model bytes relative to the uniform control. An immediate explicit
spin-up request must preserve exact exported bytes and completed-cycle count.
A positive subsequent advance must prevent initial-state export. A mismatched
simulation-start/end trial must return false, leave production uninitialized,
and make both provenance retrieval and export unavailable with caller outputs
unchanged. Failure after candidate publication clears all proof/export flags;
a mishandled false return can never expose the trial uniform state as ready.

All fixtures and exports stay beneath `MOON_U11_OUTPUT`; the common runner
preserves the actual furnished-kernel hashes. A passing dispatch fixture proves
control wiring and state identity only. It does not satisfy I19 stationarity,
thermal periodic convergence at physical tolerances, M1/I13, or observations.
