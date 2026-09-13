# Current srcSEP CLI contract through Step 4

The authoritative parser remains `SEP::Util::CLI::ParseCommandLine` in
`util/sep_cli.cpp`.  Step 1 is additive.

## Preserved behavior

- No arguments retain production defaults: coupling/cascade/reflection on,
  integrated turbulence, the canonical `fte-dmumu` mover (the historical FTE
  implementation), 300 injected particles per
  iteration, TestManager off, and spectrum output interval 100.
- Option values are case-insensitive.  Existing Boolean values and aliases,
  including `--cascase` and `--no-cascase`, remain accepted.
- Value-taking options accept both `--option value` and `--option=value`.
- Repeated legacy options retain last-occurrence-wins behavior.
- `-h`/`--help` exits before model initialization; malformed/unknown options
  fail and cannot fall through to production.
- Parsed production options are applied after post-compile input.
- Legacy TestManager controls retain their prior continue-to-production
  behavior and do not imply component-test mode.
- SWMF library entry points do not parse arguments or execute the registry.

## Additive Step 1 behavior

- `--list-tests` lists sorted registry metadata and exits before initialization.
- `--test ID`/`--test=ID` are repeatable.
- `--test-group GROUP`/`--test-group=GROUP` are repeatable.
- `--all-tests` selects the bounded routine set only.
- IDs/groups are case-insensitive; overlaps are de-duplicated and run once in
  stable-ID order.
- Listing cannot be combined with execution.  `--all-tests` cannot be combined
  with explicit ID/group selectors.  Missing/unknown selectors fail early.
- Any execution selector is test-only and returns before TestManager or the
  production timestep loop.  FAIL/ERROR return nonzero; SKIP remains explicit.

## Step 4 mover behavior

- `--particle-mover`, `--mover`, and `--sep-mover` select only `parker`,
  `fte-dmumu`, or `fte-mfp`.
- `--list-movers` lists those three descriptors and exits before initialization.
- Unambiguous legacy names map for one transition release with a warning.
- Ambiguous, direct-wave, and 3-D historical mover names fail before
  initialization.
- Startup metadata prints the canonical name, capability contract, and active
  coefficient provider.
