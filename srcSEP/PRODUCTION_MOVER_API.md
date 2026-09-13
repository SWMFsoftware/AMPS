# Production field-line mover API

Step 4 establishes one public registry for the three supported production
movers. The registry describes physics and state capabilities independently of
implementation addresses; the PIC callback enters through one validating
adapter.

## Canonical movers

| CLI name | Public enum | Production implementation | Coefficient contract |
|---|---|---|---|
| `parker` | `ProductionMover::Parker` | `ParticleMover_Parker_Dxx` | Pitch-angle-averaged spatial diffusion `Dxx`, integrated from the configured `Dmumu` provider |
| `fte-dmumu` | `ProductionMover::FocusedTransportDiffusion` | `ParticleMover_FTE` | Explicit pitch-angle diffusion coefficient `Dmumu` |
| `fte-mfp` | `ProductionMover::FocusedTransportMeanFreePath` | `ParticleMover_FocusedTransport_EventDriven` | Parallel mean free path `lambda`, converted to event scattering |

These are one-dimensional transport models along field lines embedded in 3-D
space. The field geometry remains three-dimensional. Fully 3-D particle
trajectories, drift/Boris variants, direct-wave experimental movers, and
historical ambiguous scattering functions are not public production choices.

## Selection and discovery

```sh
../amps --list-movers
../amps --particle-mover parker
../amps --particle-mover=fte-dmumu
../amps --mover fte-mfp
```

`--help` and `--list-movers` advertise only the three canonical names. The
following unambiguous aliases are accepted for one transition release and print
a deprecation warning naming their replacement:

- `parker-dxx`, `parker-diffusion`, `dxx` → `parker`;
- `fte`, `default`, `focused-transport`, `focused-transport-equation`,
  `legacy-fte` → `fte-dmumu`;
- `focused-transport-event-driven`, `event-driven`, `event-driven-fte`,
  `fte-event-driven` → `fte-mfp`.

Names such as `coupled-fte`, `focused-transport-wave-scattering`,
`parker-mean-free-path`, `mean-free-path-scattering`, `tenishev-2005-fl`, 3-D,
drift, and Boris names are rejected. They differ physically or are ambiguous
and therefore cannot be silently mapped.

## Capabilities and callback contract

Each descriptor reports whether the mover requires field-line attachment, uses
pitch-angle state, accumulates streaming for the turbulence manager, evolves
wave state directly, and which coefficient contract it consumes. Main-loop
turbulence dispatch queries `CurrentCapabilities()`; it no longer compares raw
function addresses.

All PIC mover calls route through `DispatchProductionMover()`. Before invoking
the selected implementation, the adapter verifies:

- a finite nonnegative time interval and a valid particle record;
- enabled field-line attachment;
- a valid field-line ID and attached segment;
- a finite field-line coordinate; and
- finite parallel and normal particle velocities.

Startup output records the canonical mover, representation requirements,
turbulence streaming behavior, coefficient contract, and active provider.
Self-consistent Alfvén turbulence remains a separate production subsystem; all
three movers provide manager-consumed streaming terms and none directly owns
wave-state evolution.

## Focused tests

Run `make test-mover-api-unit` for:

- `MOVCLI01`: canonical parsing, help, and exactly three discovery entries;
- `MOVCLI02`: transition aliases and actionable warnings;
- `MOVCLI03`: capability reporting and absence of main-loop pointer comparison;
- `MOVCLI04`: rejection of unsupported or ambiguous historical movers.

The complete native AMPS regression gate must additionally exercise one
particle through each adapter mapping and confirm startup metadata on the linked
executable. That dependency is not present in a source-only `srcSEP` archive.
