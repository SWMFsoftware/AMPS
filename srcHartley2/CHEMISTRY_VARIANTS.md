# Alternative reaction definitions (chemistry variants): design and implementation plan

Status: **implemented in the master working tree, not committed** (2026-11-21); tests in section 13. The decisions
are recorded in section 12. Sections 2-11 describe the implementation as built; where it departs from the plan the
text says so. Based on the code at commit `7c05480` (`srcHartley2/CometChemistry.*`,
`src/models/chemical_solver/ReactionNetwork.h`, `input/hartley2.pl`, `input/CometReactionNetwork.pm`,
`srcHartley2/CometRuntimeInput.cpp`, `srcHartley2/RadialMain.cpp`, `srcHartley2/i3/prepare_i3_cases.py`).

## 1. Goal

Allow an alternative definition of a chemical process to be compiled into the executable and selected at run time,
by a switch in the `.in` file and/or a command-line case switch, without rebuilding. The first use is an
alternative OH loss: OH -> O + H with a lifetime that depends on the heliocentric velocity (Schleicher and A'Hearn
1982, ApJ 258, 864), as in the legacy Fortran below, instead of the present state-resolved channels.

Original legacy source as supplied by the user on 2026-11-21 (verbatim; this is the provenance of the table). The
excitation arrays `OHLON4` and `OHLON1` are out of scope (decision 12) and are kept here only for completeness.

```fortran
      BLOCKDATA OHPHOT  
      implicit double precision(a-h,o-z)
      COMMON /OHDATA/ OHTAUS(61),OHLON4(61),OHLON1(61)
C  
C     ARRAY OHTAUS(I) -- CONTAINS THE HELIOCENTRIC VELOCITY DEPENDENT
C                        LIFETIME IN 10**5 SECONDS AT R ~ .25 AU  
C  
C     ARRAY OHLON4(I) -- CONTAINS THE HELIOCENTRIC VELOCITY DEPENDENT
C                        EXCITATION RATE  
C  
C     ARRAY OHLON1(I) -- CONTAINS THE HELIOCENTRIC VELOCITY DEPENDENT
C                        EXCITATION RATE  
C  
C     GENERAL COMMENTS -- ARRAY ELEMENTS FROM 1 TO 61 CORRESPOND TO  
C                         HELIOCENTRIC VELOCITY VALUES IN KM/S FROM  
C                         -60 TO + 60 WITH INTERVALS OF 2 KM/S.
C  
C                      -- THE LIFETIME DATA ARE FROM THE PAPER BY 
C                         SCHLEICHER AND A'HEARN (1982)  ASTROPHYS. J.  
C                         VOL. 258, P 864.
C  
C                      -- THE EXCITATION DATA ARE ALSO FROM SCHLEICHER  
C                         AND A'HEARN (1984) BULL. A. A. S. VOL 15, P806
C  
      DATA OHTAUS/ 2.85, 3.15, 3.10, 2.85, 2.60, 2.45, 2.55, 2.95, 3.70,
     2 4.00, 3.80, 3.50, 3.30, 3.05, 2.80, 2.80, 3.10, 3.45, 3.45, 3.25,
     3 3.45, 3.50, 2.80, 2.45, 2.30, 2.00, 1.90, 2.00, 2.35, 2.66, 2.75,
     4 2.80, 2.90, 3.20, 3.60, 3.40, 3.10, 2.95, 3.00, 2.95, 2.85, 2.75,
     5 2.60, 2.30, 2.25, 2.25, 2.30, 2.35, 2.25, 2.20, 2.25, 2.40, 2.55,
     6 2.55, 2.45, 2.35, 2.30, 2.40, 2.55, 2.70, 2.85/
      DATA OHLON4/ 5.0, 5.0, 4.7, 4.4 ,4.0, 3.4, 2.7, 2.5, 2.4, 2.7, 
     2   3.0, 2.9, 3.0, 3.1, 3.1, 3.1, 3.0, 2.7, 2.2, 1.8, 1.4, 1.2, 
     3   1.6, 2.2, 2.5, 2.2, 1.8, 1.7, 1.7, 1.8, 1.9, 1.9, 2.2, 2.4, 
     4   2.3, 2.7, 3.5, 4.4, 4.4, 4.1, 4.3, 4.8, 5.3, 5.4, 5.1, 4.6, 
     5   3.9, 3.5, 3.0, 2.4, 2.1, 2.3, 3.7, 5.2, 5.9, 6.3, 6.4, 6.3, 
     6   5.8, 5.0, 4.5/ 
      DATA OHLON1/ 5.5, 5.5, 5.1, 4.7, 4.1, 3.4, 2.6, 2.3, 2.3, 2.6, 
     2   2.9, 2.8, 2.8, 2.8, 2.9, 2.9, 2.8, 2.4, 2.0, 1.7, 1.2, 1.1, 
     3   1.5, 2.1, 2.3, 2.0, 1.6, 1.5, 1.6, 1.7, 1.8, 1.7, 2.1, 2.3, 
     4   2.2, 2.5, 3.5, 4.5, 4.6, 4.1, 4.6, 5.3, 5.8, 5.9, 5.6, 4.8, 
     5   3.9, 3.4, 2.9, 2.3, 1.9, 2.1, 3.7, 5.2, 6.2, 6.7, 7.0, 7.0, 
     6   6.4, 5.5, 4.7/ 
      END
C
C If vr is the heliocentric velocity in km/s
C taum0 = 4.24d5
C pick ohtaus from the array above 
C rau is the heliocentric distance in au
C ivel picks the array element from OHTAUS
C tauoh is the OH lifetime from the Schleicher & A'Hearn heliocentric velocity dependent part
C taum0 is the OH lifetime from the FUV-EUV solar spectrum -- mostly Lyman alpha
C
C tohrau is the OH lifeimt at the the heliocentric distance rau and velocity vr
C
      ivel=min0(61,max0(1,int(vr+62/2.d0)))
      tauoh=ohtaus(ivel)*0.833d5
      tohrau=rau*rau/(1.d0/taum0 + 1.d0/tauoh)
```

## 2. Concept: reaction groups and variants

- A **reaction group** is one physical process with alternative definitions. Selecting a variant replaces a whole
  set of channels, not one reaction.
- A **variant** is one complete definition of the group. Exactly one variant per group is active in a run.
- The first group is `OH_LOSS`:
  - `STATE_RESOLVED` (default): the present `OH_O3P_H`, `OH_O1S_H`, `OH_O1D_H` and `OH_ION`, with rates from the
    `chemistry.rate.<ID>_1au_s-1` keys exactly as today;
  - `SA1982_VH`: `OH_O_H_SA1982` (OH -> O + H, one lumped channel, excess energy 0.36 eV) and `OH_ION_SA1982`
    (OH -> sink), with rates from the lifetime table (section 3).
- All variants are compiled in, so switching needs no rebuild. Ungrouped reactions are always active.
- **Mechanism.** At start-up every reaction of an inactive variant gets the rate coefficient A = 0. The per-parent
  channel cache in `srcHartley2/CometChemistry.cpp` already keeps only reactions with a positive rate, in input order
  ("Retain input order and the original left-to-right rate accumulation"). The active channel list, its order, the
  floating-point selection thresholds and the random-number stream are therefore those of the active variant alone.
  `ChemicalSolver::cReactionNetwork::Validate()` already accepts A = 0. `src/models/chemical_solver/ReactionNetwork.h`
  (shared with `srcMoon`) is not changed.
- The two new reactions are appended after the existing ones in `input/i3.conf`, so no existing reaction index
  changes. Reaction IDs do not appear in the physics output (`PT/plots`) or the run log. The chemistry diagnostics
  (`PT/diagnostics/comet_stage0.out=<k>.csv/.json`) have one row per reaction and therefore gain two rows with zero
  counts.
- **Default = bitwise identical:** with `STATE_RESOLVED`, every file in `PT/plots` of the nightly deck must be
  byte-identical to the build of commit `7c05480`, and the diagnostics identical apart from the two added rows
  (section 10).

## 3. Rates of `SA1982_VH`

At 1 AU, with the existing heliocentric inverse-square law (reference distance 1 AU) supplying (1 AU / r_h)^2:

```
k_d   = 1/tau_FUV + 1/(S * OHTAUS(v_h))      [s^-1]   OH -> O + H   (OH_O_H_SA1982)
k_ion = f_ion * k_d                          [s^-1]   OH -> sink    (OH_ION_SA1982)
```

- `OHTAUS(v_h)`: linear interpolation in the table file (section 4), in units of 1e5 s. Outside the table's velocity
  range the end value is used (clamped). The clamp is reported as a warning in the run log and recorded in the
  output headers.
- `v_h`: `brightness.heliocentric_radial_velocity_km_s` of the selected case, i.e. the same value as used for the
  brightness g-factor (decisions 11 and 16).
- `S = 0.833e5 s` and `f_ion = 0.05`: build-time model constants in `input/i3.conf` (section 5); changing them
  needs a rebuild.
- `tau_FUV = 4.24e5 s`: stored in the table file, because the FUV/EUV part (mostly Lyman alpha) is the part that
  depends on solar activity (proposal awaiting confirmation, section 12.3).
- No solar-activity (F10.7) scaling. Another activity level is represented by another table file.
- AMPS evaluates k_d and k_ion once at start-up (`CometRuntimeInput.cpp`), on every MPI rank. In the 1D model r_h
  and v_h are constant per case, so there is no new per-particle rate law; the two numbers enter the network through
  the existing `ConfigureReactionRate()`.

Expected values for the six I3 cases (legacy table, linear interpolation; `input/i3.in` at `7c05480`). These are
the reference values of the unit test (section 10); the last two columns compare with the present OH loss
(`OH_O3P_H + OH_O1S_H + OH_O1D_H + OH_ION`).

| Case | v_h (km/s) | OHTAUS | k_d (1 AU, s^-1) | k_ion (1 AU, s^-1) | 1/k_d (s) | 1/(k_d+k_ion) (s) | present OH loss (1 AU, s^-1) | present / new |
|---|---|---|---|---|---|---|---|---|
| 1 | -39.7 | 3.7550 | 5.555508e-06 | 2.777754e-07 | 180002 | 171430 | 1.4464e-05 | 2.48 |
| 2 | -36.2 | 3.3200 | 5.974395e-06 | 2.987197e-07 | 167381 | 159410 | 1.4642e-05 | 2.33 |
| 3 | -29.6 | 2.8600 | 6.555974e-06 | 3.277987e-07 | 152533 | 145269 | 1.4558e-05 | 2.11 |
| 4 | -5.8 | 2.0350 | 8.257656e-06 | 4.128828e-07 | 121100 | 115333 | 1.4660e-05 | 1.69 |
| 5 | +0.4 | 2.7600 | 6.708056e-06 | 3.354028e-07 | 149074 | 141976 | 1.4836e-05 | 2.11 |
| 6 | +7.1 | 3.4200 | 5.868667e-06 | 2.934333e-07 | 170396 | 162282 | 1.4515e-05 | 2.36 |

The OH lifetime is 1.7-2.5 times longer than with the present chemistry.

## 4. The table file

New file `srcHartley2/i3/reference/oh_lifetime_sa1982.csv`, next to the other per-case solar data of the case
generator:

```
# OH photodissociation lifetime: heliocentric-velocity-dependent part (array OHTAUS of BLOCKDATA OHPHOT)
# source: legacy Fortran BLOCKDATA OHPHOT (section 1 of srcHartley2/CHEMISTRY_VARIANTS.md),
#         which cites Schleicher and A'Hearn (1982), ApJ 258, 864
# ohtaus_1e5_s is multiplied by the build-time scale S (input/i3.conf) to give seconds at 1 AU
# tau_fuv_s: FUV/EUV (mostly Lyman-alpha) lifetime at 1 AU for the same solar conditions
tau_fuv_s,4.24e5
v_h_km_s,ohtaus_1e5_s
-60,2.85
-58,3.15
...
60,2.85
```

- The reader requires at least two rows, a strictly increasing velocity column, and finite positive values; it
  does not assume a uniform grid, so a future table may use another spacing. Any violation stops the run with the
  file name and line number.
- The file is named in the runtime input (section 6). A relative path is resolved against the directory of the
  `.in` file. The I3 production scripts copy the table next to `i3.in`.
- Every rank reads the file; the ranks compare a hash of its content, as they already do for the canonical runtime
  input, and the run stops if they differ. The path, the SHA-256, `tau_fuv_s`, the velocity range and the
  interpolated `OHTAUS` go to the run log and the output headers.

## 5. Build-time definition (`input/i3.conf`, backward compatible)

The existing 9-field `ChemistryReaction` record is unchanged; optional `KEY=VALUE` fields may follow it.

```
ChemistryGroup=OH_LOSS | DEFAULT=STATE_RESOLVED | VARIANTS=STATE_RESOLVED,SA1982_VH
ChemistryReaction=OH_O3P_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 1.27 | GROUP=OH_LOSS:STATE_RESOLVED
ChemistryReaction=OH_O1S_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 10   | GROUP=OH_LOSS:STATE_RESOLVED
ChemistryReaction=OH_O1D_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 7.73 | GROUP=OH_LOSS:STATE_RESOLVED
ChemistryReaction=OH_ION   | OH | NONE  | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0    | GROUP=OH_LOSS:STATE_RESOLVED
...                                                  (all other reactions unchanged, same order)
ChemistryReaction=OH_O_H_SA1982 | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0.36 | GROUP=OH_LOSS:SA1982_VH | RATE=TABLE:OH_TAU_SA1982:DISSOCIATION
ChemistryReaction=OH_ION_SA1982 | OH | NONE  | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0    | GROUP=OH_LOSS:SA1982_VH | RATE=TABLE:OH_TAU_SA1982:IONIZATION
ChemistryRateTable=OH_TAU_SA1982 | MODEL=SA1982_OH_LIFETIME | SCALE_S=0.833e5 | ION_FRACTION=0.05
```

- The tags of the four existing OH records only add metadata; their order and the order of all other records stay
  as they are.
- Preprocessor checks (`input/CometReactionNetwork.pm`, `input/hartley2.pl`): every group has exactly one `DEFAULT`
  and lists its variants; every listed variant has at least one reaction; all reactions of a group have the same
  parent species; a reaction belongs to at most one group; every `RATE=TABLE:` reference names a declared
  `ChemistryRateTable`; table-rate reactions belong to a non-default variant. Violations stop the build with the
  record and the rule.
- Generated data (`CometChemistry.generated.h`): the groups (ID, default variant, variant names), per reaction its
  group and variant index (-1 = ungrouped) and table role, and the table records (model, `SCALE_S`,
  `ION_FRACTION`). Nothing in `ReactionNetwork.h` changes.

## 6. Run-time selection

In `i3.in`, in the global block or a case block (the case block overrides the global block):

```
chemistry.variant.OH_LOSS = SA1982_VH                        ! default STATE_RESOLVED
chemistry.table.OH_TAU_SA1982.file = oh_lifetime_sa1982.csv   ! relative to the directory of the .in file
```

- `chemistry.table.<TABLE>.file` is required when a variant that uses the table is active. Otherwise it is accepted
  but the file is not read (a missing file is then not an error), so the variant can be switched on the command line
  without editing the input, and a default-variant run does not need the table next to its input.
- An input without any variant or table keys (older inputs, `input/i3-test.in`) runs the compiled default.
- Rate keys `chemistry.rate.<ID>_1au_s-1`: required for the active reactions that take their rate from the input
  (as today); optional for reactions of inactive variants, which are accepted and logged as ignored; not allowed
  for table-rate reactions (the error says that their rate comes from the table and names the table key).
- An unknown group, variant or table name is an error that lists the valid names.

On the command line (`srcHartley2/RadialMain.cpp`, next to `-i` and `-case`), repeatable:

```
./amps -i i3.in -case 3 --chemistry-variant OH_LOSS=SA1982_VH
```

- Precedence: command line > case block > global block > compiled default. `--help` lists the option, the compiled
  groups and their variants.
- Provenance: for every group the effective variant and where it came from (command line / case / global /
  default); for the active table-rate reactions v_h, OHTAUS, the clamp flag, k_d and k_ion at 1 AU, and the table
  path and SHA-256. They go to the run log and, as `DATASETAUXDATA` records, to the Tecplot output headers.

## 7. Run-time chemistry (`srcHartley2/CometChemistry.cpp`, `srcHartley2/CometRuntimeInput.cpp`)

- `ApplyCometBlock()` resolves the variant of every group, reads the table if needed, sets the table-rate reactions'
  A from section 3, sets A = 0 for every reaction of an inactive variant, and otherwise reads the rate keys as today.
- Stage 6 (daughter-packet splitting, per reaction): policies of inactive reactions are dormant, because a zero-rate
  reaction is never selected.
- Stage 7 (channel importance sampling, per parent family): `CometChemistry.cpp` already stops with "channel
  importance targets a reaction with zero runtime rate". Extend that error (decision 19) to name the group, the
  active variant and where it was selected, the configured family, and the ways forward: select the variant the
  policy was generated for, regenerate the policy for the active variant, or disable channel importance sampling
  for that parent.
- Stage 8 (total-reaction forcing, per parent species) uses the total rate of the active channels; no change.
- The validation modes (`ValidateEnsemble`, the Stage 6/7/8 validators) exercise the product kinematics of every
  compiled reaction and do not use its rate, so they need no change.
- `GetNetwork()` sets A = 0 for the reactions of inactive variants when it constructs the network, so the rule also
  holds for compile-time-only inputs (where the compiled default is active). A build input with `RATE=TABLE`
  reactions requires `RuntimeInputMode=on` (checked by `input/hartley2.pl`).

## 8. Case generator (`srcHartley2/i3/prepare_i3_cases.py`)

- New options `--oh-loss {state_resolved,sa1982_vh}` (default `state_resolved`), `--oh-loss-case CASE=VARIANT`
  (repeatable per-case override) and `--oh-lifetime-table PATH` (default
  `srcHartley2/i3/reference/oh_lifetime_sa1982.csv`); all in `--help`. The table path is written relative to the
  output directory, where AMPS resolves it.
- Departure from decision 21: the per-case choice is the option `--oh-loss-case`, not a column in `i3_cases.csv`.
  The MAVEN/IUVS reference data pin the SHA-256 of `i3_cases.csv`
  (`srcHartley2/visualization/references/i3_brightness_MAVEN/metadata.json`), so changing that file would make
  their provenance record stale.
- The generator writes `chemistry.variant.OH_LOSS`, `chemistry.table.OH_TAU_SA1982.file` and, as today, the
  state-resolved rates, so that either variant can be selected later without regenerating. It does not write k_d or
  k_ion: AMPS computes them from the named table, so the table stays the only source.
- The regenerated `input/i3.in` with the default variant differs from today's only by those two keys.

## 9. Edge correction (`srcHartley2/visualization/comet_edge_correction.py`)

- `parse_reactions()` reads the excess energy from field 9 instead of the last field and parses the `KEY=VALUE`
  extras, `ChemistryGroup` and `ChemistryRateTable` records (today a trailing `GROUP=` field breaks it).
- `load_run_inputs()` determines the active variant from `i3.in` (global/case), or from an explicit argument, or
  preferably from the run's output header when present (a command-line override is not visible in `i3.in`). It omits
  inactive reactions and computes the table-rate reactions from the table file with its own Python implementation
  of section 3.
- `srcHartley2/tests/test_comet_edge_correction.py` and `srcHartley2/i3/tests/test_prepare_i3_cases.py` are updated.

## 10. Tests and acceptance

1. **Rates (unit).** An independent Python evaluation of section 3 (not sharing code with AMPS) gives the table of
   section 3; the AMPS start-up values must agree to a relative 1e-12 for all six cases. Plus: grid points,
   mid-points, both clamps, and malformed tables (non-increasing velocity, NaN, negative value, one row, missing
   file), each with its expected error.
2. **Selection.** Precedence (command line > case > global > default); unknown group, variant or table; missing table
   key with `SA1982_VH`; a rate key given for a table-rate reaction; rate keys of inactive reactions accepted and
   logged; `--help` text.
3. **Default regression (bitwise).** Upstream's `make test_i3` (test-mode build, which is bit-for-bit reproducible,
   np = 4, cases 1 and 5 of `input/i3-test.in`) built from `7c05480` and from the new tree: all files in `PT/plots`
   byte-identical, the brightness files after removing their new `DATASETAUXDATA Chemistry.*` provenance lines
   (which must name `STATE_RESOLVED`); the chemistry diagnostics identical after removing the rows of the two new
   reactions, which must have zero selections.
4. **Policy mismatch.** A Stage-7 family configured for the state-resolved OH channels with `SA1982_VH` selected stops
   with the detailed error of section 7. Not run in this implementation: it needs a dedicated build with
   channel-importance records, which `input/i3.conf` does not have (section 13).
5. **Production path (replaces the Haser-type check, decision 25).** The nightly deck, cases 1 and
   5, 1e4 iterations, np = 8, once with each variant, read from the per-reaction chemistry diagnostics
   (`PT/diagnostics/comet_stage0.out=<k>.csv`):
   - with `SA1982_VH`, the four state-resolved OH reactions have zero selections, and with `STATE_RESOLVED` the two
     new reactions have zero selections;
   - with `SA1982_VH`, selections of `OH_O_H_SA1982` / selections of `OH_ION_SA1982` = k_d / k_ion = 20, within
     three binomial standard deviations (both channels share the same OH visits and time steps);
   - the diagnostics' own production check (`production_agreement`) holds for the O and H products of
     `OH_O_H_SA1982`.
6. **Validation runs (decision 26).** All six I3 cases with `SA1982_VH` and the campaign-3 setup; H and OH brightness
   compared with MAVEN/IUVS and with the state-resolved runs.

## 11. Implementation (2026-11-21)

New files:
- `srcHartley2/i3/reference/oh_lifetime_sa1982.csv`: the table (section 4).
- `srcHartley2/CometChemistryVariants.{h,cpp}`: table reader, linear interpolation with clamp, the SA1982_VH rates,
  SHA-256; no AMPS/PIC/MPI dependency.
- Tests: `srcHartley2/tests/test_chemistry_variants.{cpp,py}`, `srcHartley2/tests/run_chemistry_variants_tests.sh`,
  `input/test_comet_chemistry_variants.pl`.

Changed files:
- `input/CometReactionNetwork.pm`: optional `GROUP=` / `RATE=TABLE:` fields, `ChemistryGroup` and
  `ChemistryRateTable` records, `validate_variants()`, generated group/table/per-reaction tables.
- `input/hartley2.pl`: the two record keywords; `CometChemistryVariants.o` in the radial build; `RATE=TABLE` requires
  `RuntimeInputMode=on`.
- `input/i3.conf`: the `OH_LOSS` group, the tags of the four state-resolved OH reactions, the two `SA1982_VH`
  reactions (appended last), the `OH_TAU_SA1982` table record.
- `input/i3.in`: regenerated; adds `chemistry.variant.OH_LOSS = STATE_RESOLVED` and the table path.
- `srcHartley2/CometChemistry.{h,cpp}`, `srcHartley2/CometChemistry.generated.h` (placeholder): group/table types and
  accessors, active-variant state, A = 0 for inactive reactions, provenance records, the detailed Stage-7 error.
- `srcHartley2/CometRuntimeInput.{h,cpp}`: variant resolution, rate-key rules, table reading, table-derived rates,
  MPI check that all ranks read the same tables, `AddChemistryVariantOverride()`.
- `srcHartley2/RadialMain.cpp`: `--chemistry-variant`, `--help`, provenance in the run log.
- `srcHartley2/CometBrightness.cpp`: provenance records as `DATASETAUXDATA Chemistry.*` in the brightness headers.
- `srcHartley2/i3/prepare_i3_cases.py` and `srcHartley2/i3/tests/test_prepare_i3_cases.py`.
- `srcHartley2/visualization/comet_edge_correction.py` and `srcHartley2/tests/test_comet_edge_correction.py`.

`src/models/chemical_solver/ReactionNetwork.h` is unchanged (decision 22).

## 12. Decisions (2026-11-21)

Answers of the user to the 27-item decision list and to the follow-up questions of 2026-11-21. Item numbers refer to
the 27-item list.

### 12.1 Decided

| # | Decision |
|---|---|
| 1 | Table grid: -60..+60 km/s in steps of 2 km/s, as the original comment states. The legacy index `int(vr+31)` does not match this grid (it reads the table near 2 v_h). |
| 2 | The table is used as the legacy code uses it: OHTAUS x S = velocity-dependent lifetime at 1 AU, scaled by r_h^2. The expected result is the legacy table with linear interpolation (follow-up 1); the per-case lifetimes supplied by the user (12.2) are not the reference. |
| 3, 17 | S = 0.833e5 s is a tunable model constant, fixed at build time (not a run-time key). |
| 4 | No source beyond the legacy code's own citation (Schleicher & A'Hearn 1982, ApJ 258, 864); provenance is the legacy `BLOCKDATA OHPHOT` of section 1. |
| 5 | Linear interpolation; no legacy-index mode. |
| 6 | Outside the table's velocity range: clamp, and record it in the run log and the output headers. |
| 7 | One lumped channel OH -> O + H. |
| 8 | Excess energy of the lumped channel: 0.36 eV. |
| 9 | OH ionization in this variant: k_ion = 0.05 k_d, with k_d the full dissociation rate (both parts). |
| 10 | No solar-activity scaling (follow-up 2). The table is kept in a file named in the `.in`, so that another table can be used for other solar conditions. This replaces the earlier answer "scale". |
| 11, 16 | v_h is the case's `brightness.heliocentric_radial_velocity_km_s`, also used for the g-factor (option (a); proposed and not objected to). |
| 12 | `OHLON4` / `OHLON1` are out of scope. |
| 13 | Names as proposed and not objected to: group `OH_LOSS`; variants `STATE_RESOLVED`, `SA1982_VH`; reactions `OH_O_H_SA1982`, `OH_ION_SA1982`; table `OH_TAU_SA1982`. |
| 14 | Selection in `i3.in` (global and case block) and on the command line; precedence command line > case > global > default. |
| 15 | Rate records of inactive variants are optional. The default OH chemistry is the present state-resolved one. |
| 18 | Superseded by follow-up 2: the table is read at run time from the file named in the `.in` (not embedded at build time); its SHA-256 is logged. |
| 19 | A Stage-7 policy that does not match the active channels stops the run with an error that explains what happened and how to proceed. |
| 20 | Effective variant, its source, the rates and the table provenance go to the run log and the Tecplot headers. |
| 21 | Case-generator options documented in `--help`; CSV column. |
| 22 | Implementation in `srcHartley2`; `ReactionNetwork.h` (shared with `srcMoon`) is not changed. |
| 23 | `comet_edge_correction.py` and its test are updated. |
| 24 | General group/variant mechanism; `OH_LOSS` is its first use. |
| 25 | Haser-type check with a +-5% tolerance from 1e3 to 1e5 km, fixed before the run (follow-up 3). Its practical value was questioned on 2026-11-21 (third round); a replacement is proposed in 12.3. |
| 26 | All six I3 cases with the new variant; H and OH compared with MAVEN/IUVS. |
| 27 | One commit. |
| O6 | Approved (follow-up 5): the variant's rates are per-case constants, no new per-particle rate law. Modified by follow-up 2: AMPS computes them at start-up from the table file named in the `.in`, instead of the case generator writing them, so that changing the file name in the `.in` takes effect without regenerating. |

Third round (2026-11-21):

| # | Decision |
|---|---|
| tau_FUV | Stored in the table file (section 4), with the table, because it is the activity-dependent part. |
| 26 | The six validation runs start after the O-brightness runs of `validation_runs/i3_master_campaign3_oxygen_20261119` have finished. |

### 12.2 Per-case lifetimes supplied by the user (not the reference)

Dissociation lifetime at 1 AU, 1/k_d, with tau_FUV = 4.24e5 s and S = 0.833e5 s. "Implied OHTAUS" is the table value
that reproduces the user's lifetime. No standard interpolation (linear, cubic splines, PCHIP, Akima) reproduces these
values; the reference is the table with linear interpolation (decision 2).

| Case | v_h (km/s) | user (s) | linear interpolation (s) | user / linear | implied OHTAUS | linear OHTAUS |
|---|---|---|---|---|---|---|
| 1 | -39.7 | 180139 | 180002 | 1.0008 | 3.76 | 3.755 |
| 2 | -36.2 | 165846 | 167381 | 0.9908 | 3.27 | 3.32 |
| 3 | -29.6 | 150815 | 152533 | 0.9887 | 2.81 | 2.86 |
| 4 | -5.8 | 118743 | 121100 | 0.9805 | 1.98 | 2.035 |
| 5 | +0.4 | 148724 | 149074 | 0.9977 | 2.75 | 2.76 |
| 6 | +7.1 | 174501 | 170396 | 1.0241 | 3.56 | 3.42 |

### 12.3 Proposal awaiting confirmation

- Replace the Haser-type check (decision 25) by the production-path check of section 10, test 5. A comparison of OH
  profiles with a model mostly retests things that do not depend on the variant (collisions in the inner coma, the
  reference model itself), and a failure would not say what is wrong. The per-reaction diagnostics show directly
  whether the run uses the selected channels with the configured rate ratio; the size of the effect on the OH and H
  profiles is seen in the six validation runs.

## 13. Test results (2026-11-21)

Builds: GCC 11.5.0, Open MPI 4.1.7rc1, upstream `make test_i3` (test mode, `-amps-test=on`, `build/` removed), np = 4.
Scratch trees in `validation_runs/oh_loss_variants_20261121/` (not committed): `AMPS_variants` (the working tree)
and `AMPS_7c05480` (commit `7c05480`). Executables: SHA-256 `6fcb05261f70b6aa...` (new), `d11847fcb8d12a2d...`
(baseline). Scripts and logs: `scripts/`, `logs/`, `runs/` there; a first attempt with two test-script defects is
kept in `runs/first_attempt/`.

| Test | Result | Evidence |
|---|---|---|
| 1. Rates (unit) | PASS | `srcHartley2/tests/run_chemistry_variants_tests.sh`: C++ self-test (SHA-256 vectors, reader, malformed tables, interpolation, clamps); rates against the independent Python evaluation at the six case velocities and 11 further velocities, relative difference <= 1e-12; preprocessor records (24 checks). The pre-existing `input/test_comet_reaction_network.pl` passes unchanged. |
| 2. Selection | PASS, 30/30 | `scripts/selection_tests.py` with the built executable: precedence (command line > case > global > default), source reporting, unknown group/variant/table, malformed option, duplicates, missing table key/file, malformed table, table-rate key, inactive rate keys, clamp (table truncated to +-30 km/s, case 1), relative path, `--help`. All six cases x both variants through AMPS start-up (`scripts/validate_six_cases.py`): the installed SA1982_VH rates equal the table of section 3. |
| 3. Default regression | PASS | `scripts/compare_default.py`: 22 files (18 in `PT/plots`, 4 diagnostics) identical to the `7c05480` build, after removing the brightness `Chemistry.*` provenance lines and the two zero-count diagnostics entries of the new reactions. |
| 4. Stage-7 policy mismatch | NOT RUN | Needs a dedicated build with channel-importance records. The new text only extends an existing error path. |
| 5. Production path | PASS | `scripts/production_path_check.py`, nightly deck at np = 8 with `--chemistry-variant OH_LOSS=SA1982_VH`: the four state-resolved OH reactions have zero selections (and the SA1982 reactions zero in the default runs); ionization share of the OH events 7/236 (case 1, -1.30 sigma) and 46/1111 (case 5, -0.97 sigma) against 0.05/1.05; `production_agreement=true` for O and H of `OH_O_H_SA1982`. With these short runs the ratio check detects errors of the dissociation/ionization ratio of roughly 40% or more, not small ones; the ratio itself is verified exactly by tests 1 and 2. |
| 6. Validation runs | PENDING | After the O-brightness campaign (`validation_runs/i3_master_campaign3_oxygen_20261119`) has finished. |

Also run: `srcHartley2/i3/tests/test_prepare_i3_cases.py` (11 tests), `srcHartley2/tests/test_comet_edge_correction.py`
(11 tests, 4 new), `input/test_comet_brightness_config.pl`; the edge-correction factors of an existing run (O campaign
case 1, output 6) are unchanged.

Complete portable suite `srcHartley2/tests/run_all_stage_tests.sh` (default archived modes, 2026-11-21): PASS on the
new tree (46 tests in 14 groups, including the new `i3.chemistry-variants`) and on the `7c05480` tree (45 tests),
0 skipped, 0 failed (`validation_runs/oh_loss_variants_20261121/logs/suite_*_13777[89]*`). In archived mode the Stage
4-9 integrated tests check frozen evidence logs; they do not build or run the new executable. Two environment notes:
master lacks the suite's data, `srcHartley2/tests/data/` and `input/halley1d.input` (master's `.gitignore` rules `data`
and `*.input` dropped them in the merge), so both were copied from the dev clone into the scratch trees; and `/tmp` on
the compute nodes does not allow executing files, so the jobs set `TMPDIR` to a directory on the shared file system.
