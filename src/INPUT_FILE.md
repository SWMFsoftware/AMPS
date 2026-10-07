# Shared AMPS Runtime Input File

This document defines the common container syntax for post-compile AMPS input.
One file may contain core settings and settings for several applications. Each
application parser reads only its named section; it must not reinterpret keys
owned by another application.

The first implementation is the `srcSEP3D` section parser in
`srcSEP3D/runtime/application_input.{h,cpp}`. It runs after
`Init_BeforeParser` and before AMPS freezes cell storage or builds the mesh.
The shared syntax is application-independent even though other applications do
not yet consume it.

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
particles_per_iteration = 1000
#section end
```

`particles_per_iteration` is an unsigned decimal integer. It maps to the
existing srcSEP3D exact `source.samplesPerStep` contract, which is a count per
compiled AMPS species over the complete source surface on each injection
iteration. Zero is accepted by the text parser and is meaningful for a
source-disabled/background-only run; the immutable physics configuration still
rejects zero if particle injection is enabled.

No other srcSEP3D keys are implemented in the shared-section parser yet.
Unknown or duplicate keys fail closed.

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

`input/sep3d.in`:

```text
#section begin: sep3d
particles_per_iteration = 1000
#section end
```

At successful startup rank zero prints the resolved root file, every expanded
file, the value and its source line, and the final immutable configuration
fingerprint. Rank zero performs file I/O and broadcasts the resolved integer so
all MPI ranks commit the same configuration.

## Current limitations

- Only srcSEP3D consumes this syntax today.
- The shared srcSEP3D mode currently supports normal or
  `--initialization-only` execution. Allocation-free `--dry-run` and restart
  remain on the maintained complete `--input` schema path because they execute
  before the required native parser boundary.
- The shared srcSEP3D section currently changes only the injection count. It
  does not replace the full schema-4 physics/configuration surface.
