# Shared AMPS Runtime Input File

This document defines the common container syntax for post-compile AMPS input.
One file may contain core settings and settings for several applications. Each
application parser reads only its named section; it must not reinterpret keys
owned by another application.

The first implementation was the `srcSEP3D` section parser in
`srcSEP3D/runtime/application_input.{h,cpp}`. It runs after
`Init_BeforeParser` and before AMPS freezes cell storage or builds the mesh.
The selected reduced provider is resolved immediately after parsing, while
the time step and statistical weights are finalized only after AMPS has
allocated the actual distributed mesh.
The shared syntax is application-independent. `srcMoon` now consumes a strict
top-level `moon` section before MPI/application initialization; it does not use
the SEP3D subsections or defaults.

### `srcMoon` section

The native lunar executable uses the same one-dash selection syntax:

```bash
./amps -input srcMoon/examples/lola_surface.in
```

Omitting `-input` deliberately retains the historical analytic-sphere
regression. The `moon` section has no implicit values; all keys below are
required even when `surface_geometry=sphere`, so a run receipt remains
complete and switching geometry cannot expose stale undeclared paths:

```text
#section begin: moon
spice_path = /home/vtenishe/SPICE
surface_geometry = lola
surface_mesh_resolution_m = 100000
lola_product_id = LDEM_4
lola_image_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
lola_label_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL
surface_cea_file = test_output/srcMoon/surface/lola_surface.cea
surface_tecplot_file = test_output/srcMoon/surface/lola_surface.dat
#section end
```

Relative paths are resolved against the file containing the section. The
SPICE root must have `cspice/include/SpiceUsr.h`, `cspice/lib/cspice.a`, and
`Kernels/`. `surface_mesh_resolution_m` is the maximum requested great-circle
edge on the pre-topography reference sphere, not a latitude/longitude step or
an AMR cell size. The CEA and Tecplot output paths must differ. Runtime path
selection does not change the build-time SPICE mode.

## File selection

The native srcSEP3D executable selects a shared file with one dash:

```bash
./amps -input path/to/run.in
```

If `-input` is absent, the shared input path is `amps.in` in the process's
current working directory:

```bash
./amps                         # reads ./amps.in
```

The existing srcSEP3D `--input FILE` interface remains available for complete
versioned INI/schema-4 decks. The two spellings are intentionally distinct
during migration: `-input` selects this shared section syntax, whereas
`--input` selects the maintained srcSEP3D-only schema documented in
`srcSEP3D/CONFIGURATION.md`.

## Sections

A section begins and ends with these directives:

```text
#section begin: application-name
key = value
#section end
```

The `#section begin`, `#section end`, and section-name comparisons are
case-insensitive. The include directive is spelled `#include` as shown below.
Sections cannot be nested, an end marker must have a preceding begin marker,
and a named application section may occur only once in the expanded input.
Application parsers ignore ordinary assignments in other well-formed
sections.

The current srcSEP3D section is:

```text
#section begin: sep3d
shock_model = reduced-shock-surface
background_plasma_model = corona-swcme-ambient
source_model = accepted-shock-incident-flux
maximum_time_steps = 4
particles_per_iteration = 1000
maximum_particle_speed_m_s = 2.0e8
time_step_margin_factor = 0.30
source_normalization_radius_m = 1.3914e10

#subsection begin: shock-particle-injection
statistical_weight_model = constant-statistical-weight
minimum_energy_j = 1.602176634e-15
maximum_energy_j = 1.602176634e-11
phase_space_power_model = compression-ratio
maximum_events_per_species_per_step = 1000000
#subsection end

#subsection begin: reduced-shock-surface
schema = shock-front-ambient-v1.1
! ...complete strict model-v1.1 assignments...
handoff.apex_radius_m = 1.3914e10
assets.harmonics_file = magnetic/pfss.harmonics.csv
assets.harmonics_sha256 = <64 lowercase hexadecimal digits>
#subsection end
#section end
```

All top-level values are mandatory; none of the physical or numerical values
below acquires a parser default:

- `shock_model=reduced-shock-surface` selects the finite prescribed front.
- `background_plasma_model=corona-swcme-ambient` selects the provider's
  maintained coronal/PFSS-to-Parker ambient plasma and IMF.
- `source_model=accepted-shock-incident-flux` selects the gross upstream
  particle flux through accepted fast-shock faces for statistical
  normalization. It is not an SEP acceleration or injection-efficiency law.
- `maximum_time_steps` is the positive standalone AMPS iteration horizon. It
  is explicit so a particle-producing run cannot silently inherit the
  library's intentionally large coupled-host default.
- `particles_per_iteration` is a positive unsigned model-particle count per
  compiled AMPS species.
- `maximum_particle_speed_m_s` is finite and in `(0,c]`.
- `time_step_margin_factor` is dimensionless and in `(0,1]`.
- `source_normalization_radius_m` is the required heliocentric apex radius at
  which the front/source rate is evaluated. It has deliberately no default.

The `shock-particle-injection` subsection is also mandatory:

- `statistical_weight_model=constant-statistical-weight` selects the currently
  implemented sampler. `log-uniform-momentum-importance` is parsed and
  fingerprinted but fails startup as an explicitly reserved feature.
- `minimum_energy_j` and `maximum_energy_j` are positive total kinetic-energy
  bounds per particle in SI, with maximum strictly greater than minimum.
- `phase_space_power_model=compression-ratio` derives the local isotropic DSA
  exponent `q=3X/(X-1)` from each accepted face. `constant` instead requires
  `phase_space_power_index=q>2`; that index must be omitted in compression
  mode.
- `maximum_events_per_species_per_step` is a positive fatal runaway guard. It
  never truncates, caps or renormalizes a Poisson realization.

The named reduced subsection is required and cannot be empty. Every assignment
inside it is passed to the strict `shock-front-ambient-v1.1` resolver. This
includes run support, ambient/PFSS normalization, plasma composition, geometry,
launch/acceleration history, handoff radius, drag continuation, surface
resolution, shock tolerances, observer endpoint, magnetic asset path and the
asset checksum. Relative asset paths resolve from the file containing the
subsection begin marker—not from the process directory. Unknown or missing
model keys and changed asset bytes fail in the shared resolver.

### Derived global time step

After mesh allocation, all owner blocks participate in the global minimum of
the actual AMPS characteristic cell length `h_min`. The one step shared by all
species is

```text
dt = time_step_margin_factor * h_min / maximum_particle_speed_m_s .
```

Ghost blocks do not add an independent restriction. Existing observer
cadences must be integer global ticks; shared-section mode moves each requested
cadence upward to the first exact tick at or after it, so output is never made
more frequent and the exact integer-cadence gate is not weakened. Startup
prints both requested and resolved cadences.

### Per-species particle weight

At `source_normalization_radius_m`, a side-effect-free production epoch is
constructed without advancing the committed background generation. For every
accepted curved triangular face,

```text
Ndot_s = sum_faces n_s [V_n - U_1 dot n]_+ A_face
W_s    = Ndot_s * dt / particles_per_iteration .
```

`A_face` is the exact curved SSE measure, not its planar visualization area.
Sub-fast, non-forward and below-support faces contribute zero and retain their
reported excluded area; numerical-unknown area makes normalization fail.
Electron, proton and alpha populations come from the upstream EOS. A compiled
species absent from that composition, or a compiled zero-abundance population,
fails rather than borrowing another species' rate. The final rank-zero receipt
prints the radius/time, accepted/excluded areas, physical rate and weight for
each compiled species.

At a live committed epoch, the same equation is evaluated on every accepted
triangle and the constant-weight macro-event rate is

```text
lambda_s(t) = Ndot_s(t) / W_s .
```

Successive waiting intervals are `-log(U)/lambda_s`. Conditional face
selection uses `Ndot_face/Ndot_total`, and the point on the chosen planar
triangle uses square-root barycentric coordinates. Every rank generates the
same keyed candidate list; only the owning rank allocates through
`PIC::BC::UserDefinedParticleInjectionFunction`. A point outside the allocated
finite domain, including the declared radial interval even if a Cartesian AMR
leaf exists there, is counted as disconnected and is not renormalized onto
another face. Momentum is anti-sunward. A birth at time `tau` is advanced for only
`dt-tau` on its first mover call.

## Comments

An exclamation mark starts a comment. Everything from `!` through the end of
that physical line is removed before directives, assignments, or continuation
markers are interpreted.

```text
particles_per_iteration = 2500  ! exact count per compiled species
```

There is currently no quoted-string escape for `!`; therefore file names and
future string values cannot contain an exclamation mark.

## Continued logical lines

A backslash that is the last non-whitespace character before a comment joins
the next physical line to the current logical line. The parser inserts one
space between nonempty fragments:

```text
particles_per_iteration = \ ! the value follows
  2500
```

A continuation at end of file is an error. Diagnostics report the first
physical line of the continued logical statement.

## Included files

`#include` performs textual inclusion at the directive location:

```text
#include common/core.in
#include "applications/sep3d.in"
#include <site/local.in>
```

Relative paths are resolved against the directory containing the including
file, not against the process working directory. Includes may be recursive.
The parser rejects cycles and nesting deeper than 64 files. A section may be
contained wholly in an included file, or an include may provide assignments at
a location inside an already-open section.

Bare paths consume the complete remainder of the logical line. Double-quoted
and angle-bracket paths support spaces. Text after the closing delimiter is an
error unless it has already been removed as an `!` comment.

## Error contract

Input is transactional: no parsed value is committed unless the complete
expanded input and the immutable application configuration are valid. Errors
terminate native initialization before mesh creation and report:

- the actual included file;
- the physical line number;
- what was missing or unrecognized; and
- the offending logical line.

For example:

```text
/case/parts/sep3d.in:4: srcSEP3D input error: unrecognized srcSEP3D setting 'particlez'
  line: particlez = 10
```

Missing files, include cycles, unmatched sections, missing required settings,
invalid unsigned integers, duplicate settings, and dangling continuations are
all fatal. The parser never substitutes a value after an error.

## Complete minimal example

`amps.in`:

```text
! Other parsers may consume this section in the future.
#section begin: core
output_directory = output
#section end

#include "input/sep3d.in"
```

`input/sep3d.in` (the complete maintained version is under
`srcSEP3D/examples/application-input/`):

```text
#section begin: sep3d
shock_model = reduced-shock-surface
background_plasma_model = corona-swcme-ambient
source_model = accepted-shock-incident-flux
particles_per_iteration = 1000
maximum_particle_speed_m_s = 2.0e8
time_step_margin_factor = 0.30
source_normalization_radius_m = 1.3914e10
#subsection begin: shock-particle-injection
statistical_weight_model = constant-statistical-weight
minimum_energy_j = 1.602176634e-15
maximum_energy_j = 1.602176634e-11
phase_space_power_model = compression-ratio
maximum_events_per_species_per_step = 1000000
#subsection end
#include "reduced-shock-surface.in"
#section end
```

At successful startup rank zero prints the resolved root file, every expanded
file, selectors and numerical values, the particle-count source line, resolved
asset directory, provider event identity, derived mesh/time/rate/weight
receipt, and final immutable configuration fingerprint. Rank zero performs
file I/O and broadcasts the canonical parsed values; every rank independently
resolves the same checksummed provider and the derived numerical values must be
bit-identical across ranks.

## Current limitations

- Only srcSEP3D consumes this syntax today.
- The shared srcSEP3D mode currently supports normal or
  `--initialization-only` execution. Allocation-free `--dry-run` and restart
  remain on the maintained complete `--input` schema path because they execute
  before the required native parser boundary.
- Only the selected reduced shock/ambient/source-normalization combination is
  implemented in shared-section mode. Other model selector values fail closed.
- Particle species remain compile-time AMPS choices. The file cannot add or
  relabel species.
- The reduced model supplies ambient volume plus a prescribed surface and
  one-sided shock limits. It supplies no sheath/ejecta/downstream volume and
  does not qualify BG3D-4.
- The accepted incident flux is only a normalization convention. Particle
  acceleration efficiency and actual reduced-front injection remain separate
  future physics choices.
