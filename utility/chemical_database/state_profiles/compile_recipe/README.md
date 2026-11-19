# Source-bound normal compilation species recipe

This maintained helper stages an application-owned species table through the
SAME PhysicalStateProfile/ReducedKinematics production resolver. It does not
run Config.pl, rebuild an application, edit generated headers/global constants,
modify input/species.input, or qualify a physical provider. Missing actual Na+
state/physical rest mass refuses BEFORE creating the requested output folder.

```bash
python3 -B utility/chemical_database/state_profiles/compile_recipe/prepare.py \
  --amps-root /absolute/AMPS --profile /absolute/profile.json \
  --mapping /absolute/explicit-state-map.json --template /absolute/AMPS/input/moon.input \
  --output-dir /absolute/new-package
```

Mapping is a JSON object from unique explicit compiled symbols to unique
profile state IDs, for example `{"Na":"23Na:ground","ELECTRON":"electron:free"}`
for an available-state SOFTWARE table demonstration. It is not a complete
physical sodium reaction configuration. No species name infers a state/isotope.
Every selected state and independent physical rest mass must exist. All source
and CSV hashes authenticate through the real native loader before publication.

The compiled mass is the resolved additive reference mass rounded once to the
actual binary64 PIC table. Independently sourced physical rest masses remain
separate. Source-SI elementary charge is read from authenticated CODATA. Native
ReadSpeciesBlock emits charge-token times legacy ElectronCharge, so a numeric
conversion token is independently compiled/checked to reproduce the requested
source C value exactly in binary64. This is NOT fractional physical ion charge;
integer charge remains in the material catalogue. Global ElectronCharge/eV2J
and baseline modes are preserved. Source authority, isotope/state populations,
threshold crosswalks and model-error bounds remain separate admission gates.

Output contains `source_species.input`, `moon_source_profile.input`, its
byte-identical normal compiler `.input.pl` and `.input.spicekernels` companions,
and a detailed manifest
with actual source, profile, template, code/header identities and native charge
expression checks. Native include paths resolve against process cwd, so the
species include is absolute; paths containing whitespace or native comment
characters refuse. Regenerate at the intended persistent destination if moved.

In a DISPOSABLE checkout, copy all moon_source_profile.input* files into its
input directory; retain the generated species include package. Then use normal
configuration, removing build before rebuilding as required:

```bash
rm -rf build
./Config.pl -application=moon_source_profile \
  -spice-path=/path/to/SPICE/cspice -spice-kernels=/path/to/SPICE/Kernels
make -j
```

The runtime .in remains a separate layer. A changed source mass table defines
a separate application build and does not authorize refreshing ordinary Moon's
reviewed I04 reference. The explicit Moon-only charge policy must be shared by
controller checks/fingerprint, transaction parent/products, collector install
and actual impact charge ledger. Legacy overloads retain their original unit.
The actual normal build/startup/Lorentz/table binding needs native verification;
passing this helper alone is not that evidence.

Run real-source controls (no AMPS application build):

```bash
python3 -B utility/chemical_database/state_profiles/compile_recipe/run.py --output-dir /absolute/new-test
```

`--raw-dir`, `--decision-register` and `--canonical-root` are adjustable. Missing
retained sources/compiler reports SKIPPED. Standard result/log/manifest and
0/1/2/77 exits retain PASS/FAIL/ERROR/SKIPPED. Current real missing Na+ is tested
as an unavailable state, never replaced by a fabricated physical profile.

The selected application basename is used in both parser `#block` and SPICE
kernel include references; no stale canonical moon.input.pl/kernel copies from
an earlier configuration are relied upon. Both companions are hashed and
copied byte-identically.

## Explicit stationary chemistry reference compile directives

For the separately declared `homogeneous-fixed-position` I21 reference problem,
copy the staged files into an ISOLATED checkout as above. In the generated
`input/moon_source_profile.input`, **replace** its existing total-acceleration
`define` line inside the existing `#block ../exosphere.pl` block with these
lines (do not retain a contradictory Moon acceleration definition):

```text
define _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode) TotalParticleAcceleration_default(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode)
define _MOON_COHORT_ZERO_FORCE_REFERENCE_ 1
```

Use the unchanged normal Config.pl/make commands shown above after removing
`build`. This selects the EXISTING core zero-acceleration function; it introduces
no duplicate test-only force law. SPICE/geometry remain enabled. Keep the exact
derived compile input, both companions, generated header, executable and command
hashes. The native owner admits only the complete bare or qualified default call
with at most one trailing semicolon; a marker, bare function name or extra
statement/suffix cannot pretend to establish the zero-force reference.

These directives define a controlled stationary-parent chemistry experiment,
not a full lunar-force run. Physical state/mass/charge/crosswalk and model-error
admission remain mandatory. Current unavailable Na+ still refuses before source
recipe output. Spatially varying trajectory runs retain actual Moon forces and
require an independent integrated-hazard reference, not this homogeneous fit.

## Explicit homogeneous zero-acceleration reference problem

`--reference-problem homogeneous-zero-acceleration` selects the approved
homogeneous reference experiment through normal compile-input DEFINE syntax.
It replaces the ONE simple acceleration definition inside the Exosphere block
with the existing `PIC::Mover::TotalParticleAcceleration_default` call and adds
`_MOON_COHORT_ZERO_FORCE_REFERENCE_=1`. It defines no force law and patches no
build/header/CPP output. The default `--reference-problem normal-template`
preserves the supplied force definition byte-for-byte and adds no marker.
Duplicate/missing definitions, existing marker or unsupported signatures refuse.
The named experiment retains ALL source/state/rest-mass/hash and exact SI-charge
checks. A marker never overrides a failed source or actual acceleration-binding
check in the conditional native owner.

From the SAME disposable checkout whose source and build will be recorded:

```bash
cd /ABSOLUTE/PATH/TO/DISPOSABLE-AMPS
MOON_PROFILE=/ABSOLUTE/PATH/TO/REVIEWED-SOURCE-PROFILE/profile.json
MOON_STATE_MAP=/ABSOLUTE/PATH/TO/EXPLICIT-SOURCE-STATE-MAPPING.json
MOON_COMPILE_PACKAGE=/ABSOLUTE/PATH/TO/NEW-COMPILE-PACKAGE

python3 -B utility/chemical_database/state_profiles/compile_recipe/prepare.py \
  --amps-root "$PWD" --profile "$MOON_PROFILE" --mapping "$MOON_STATE_MAP" \
  --template "$PWD/input/moon.input" --output-dir "$MOON_COMPILE_PACKAGE" \
  --reference-problem homogeneous-zero-acceleration
```

For a physical Na cohort the map must explicitly contain its admitted neutral
and ion states, for example the actual IDs `23Na:ground` and
`23Na+:unresolved` only after the latter has authoritative usable data. With the
CURRENT partial real profile that request exits77/SKIPPED and creates no compile
package. Stop there; do not continue with missing files or substitute toy ions.
Available-state Na/electron maps can exercise the compile-table tooling only;
they are not a complete physical Na reaction configuration.

Only after successful source admission, stage and build the emitted recipe in
that same disposable checkout using the normal configuration path:

```bash
cp "$MOON_COMPILE_PACKAGE"/moon_source_profile.input* input/
rm -rf build
./Config.pl -application=moon_source_profile \
  -spice-path=/ABSOLUTE/PATH/TO/SPICE/cspice \
  -spice-kernels=/ABSOLUTE/PATH/TO/SPICE/Kernels
make -j
```

Retain the normal producer's actual source/compiler/configuration, library and
generated-ABI evidence from THIS checkout when invoking the conditional cohort
wrapper there. This recipe does not fabricate a producer manifest, runtime .in,
model-error budget, admitted Na+ state or completed qualification. It selects a
reference problem; actual lunar-force transport uses normal-template and its
separate trajectory/reference gates.

The controls additionally invoke the maintained `exosphere.pl` DEFINE parser on
isolated source copies, then compile/stringify its generated macro and check
the actual native owner's exact-call guard plus marker. That is parser/compiler
binding evidence, not another implementation or execution of a zero-force law.

# Complete canonical species declarations

`source_species.input` preserves `input/species.input` component order and
unselected definitions, replacing only mass and charge for explicitly mapped
states. The normal parser needs unselected components to emit inactive negative
species macros used by shared modules, including the three CEX ENA names. This
does not activate background species. The recipe manifest pins the canonical
table. The local runner executes the maintained `ReadSpeciesBlock` and compiles
its generated inactive CEX identifiers to verify that interface directly.

Conservative reaction products additionally need the ordinary per-particle
weight-correction field. Use `--packet-layout conservative-products` (Python:
`prepare(..., packet_layout='conservative-products')`) and clean-build the
result. This explicitly changes the unique normal
`ParticleWeightCorrectionMode=off` to `on`; the default `template` preserves
all original bytes. It does not change `TimeStepMode` or species timestep
policy. The maintained M3 source-qualification build helper always requests
this storage layout and retains a preprocessor proof of the actual generated
macro alongside the normal library/build receipt.
