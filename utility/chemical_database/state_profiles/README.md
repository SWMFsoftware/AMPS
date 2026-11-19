# Source-bound M3 reference-state profiles

This component loads explicit isotope/state/energy mappings through the SAME
production `ReducedKinematics::ResolveState` used by event kernels. Read the
preimplementation [contract](CONTRACT.md) for schemas, source roles, units,
refusal behavior and test criteria. It changes no normal AMPS configuration,
compiled species mass, reaction rate, immutable dataset or installed database.

The current source-supported preparation is deliberately partial. It provides
six specified neutral-isotope ground reference cohorts, a free electron,
their measured nominal neutral/electron masses and four CODATA constants.
Formation zeros relative to separated neutral ground isotope atoms/free
electrons are definitions approved in P02. They are not fitted energy values.
Fourteen ion, excited and molecular states remain explicitly unavailable.
No complete physical photon/electron/CEX channel is activated by these states.
Native mass/constant uncertainties are preserved as source strings. No
uncertainty, covariance, classical-speed guard or physical-error budget is
invented from the printed precision.

Run preparation from the AMPS checkout, using the retained native sources:

```sh
python3 -B utility/chemical_database/state_profiles/prepare.py \
  --raw-dir test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/raw \
  --decision-register srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json \
  --canonical-root /data/vtenishe/moon_validation_data/photochemistry \
  --output-dir test_output/srcMoon/NEW-state-profile
```

The output must not exist. `profile.json` pins three CSVs and every original
source; source bytes are not copied. `canonical_gap_audit.json` verifies the
actual D04/D05 products against native provenance and retains missing
molecular thresholds and positive subthreshold atomic samples. Audit success
is not physical state qualification. In particular, the D04 Na threshold
differs from the acquired ASD ground energy; the profile neither replaces
that threshold nor creates an ion state from an unreviewed effective value.
`OH_exp` denotes a compilation, not an excited/ground OH population. Unlabelled
oxygen products do not become O(3P) by inference.

The published `catalogue-v1` database holds several rate/state sources because
managed-import permission is unresolved. This utility reads the two original
permission-cleared NIST mass/constant responses directly and records the AMPS
P02 reference definition. It does not import held ASD thresholds or the mixed
mass/threshold inventory into that archive, mutate the installed release, or
override any source licence. The independent test reads only permitted mass
and constant quantities from the earlier research inventory.

Run all component controls without configuring/rebuilding AMPS:

```sh
python3 -B utility/chemical_database/state_profiles/run.py \
  --output-dir test_output/srcMoon/NEW-state-profile-controls
```

Optional explicit `--raw-dir`, `--reference-inventory`, `--decision-register`
and `--canonical-root` arguments relocate the retained inputs. Missing input
or compiler gives SKIPPED77; corruption/build/runtime gives ERROR2; evaluated
numerical/refusal failure gives FAIL1; evaluated supporting checks met gives
PASS0. The formal M3 field stays SKIPPED/NOT QUALIFIED. Outputs preserve
commands/compiler/source/data/executable hashes, native logs, `result.json`,
`run_manifest.json` and `comparison.csv`. Native kg values are compared with
the independently prepared reference using Decimal and a fixed limit of two
reported long-double machine epsilons. No observational comparison is made.

For C++ callers, construct `PhysicalStateProfile::Files` from preserved
configuration and call `Registry::Load`. `GetState` returns a source-bound
`ReducedKinematics::ResolvedState`; it refuses unavailable entries without
changing output. `GetPhysicalRestMass` is separate and never substitutes the
additive reference mass. `GetConstant` requires the expected native unit.
Loaded identities are byte/source provenance, not approval signatures.
Physical providers must additionally validate each channel threshold and
product-state crosswalk, physical/reference-mass errors, reservoir/PDF,
frame/epoch and the unchanged numerical/qualification gates.

The accompanying source audit found a useful correction to the earlier
missing-metadata inventory: [NIST ASD with its uncertainty column enabled](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&format=1&level_out=on&seq_out=on&shells_out=on&spectra=Na&unc_out=on&units=1)
does supply the Na ground ionization uncertainty. The previous elemental
landing-page extraction omitted that column. This resolves an acquisition
gap, not the D04 threshold/state compatibility question. The source audit
retains both distinct numerical identities and does not broaden a tolerance
or infer a rounding interval. Its native H output also explicitly preserves
the theoretical-parentheses flag and supplied uncertainty.

Fresh research sources, failures, exact queries/hashes and the maintained
artifact extraction are in
`test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/`.
`source_audit.json` is documentary evidence; it is not an admitted state
profile. Seven direct ATcT requests returned HTTP403 challenge bodies; these
are retained as unavailable acquisitions and never parsed as thermochemical
tables. Indexed [ATcT API documentation](https://atct.anl.gov/API%20Calls/)
identifies a useful future native JSON route, while its most-frequent-isotomer
statement concerns the displayed geometry and must not be silently transferred
to every enthalpy or mixture-mass entry. The publicly indexed latest ATcT
network is1.228, distinct from the previously proposed1.220. A new version,
isotopologue/ortho-para interpretation and covariance/energy-zero crosswalk
require review before use. No author/provider was contacted and no access
control was bypassed.

The maintained [compile_recipe](compile_recipe/README.md) helper stages an
application-owned source-bound species table through the same native resolver,
without changing canonical species/global constants. Its own real-source tests
are registered as `source-compile-recipe`; missing Na+ states refuse before output.

## Explicit Na I/Na II research extension

`prepare.py` remains the original partial/default preparation described above.
The opt-in `prepare_atomic.py` additionally prepares `23Na+:ground` while keeping
every old `*:unresolved` declaration unavailable, including the old Na ion ID.
Its [preimplementation contract](ATOMIC_EXTENSION_CONTRACT.md) records the
source and numerical scope. This is original research use of authenticated
retained ASD/handbook sources; it is **not** a managed ASD import, new database
release, changed installed lock or automatic Moon source activation.

The same retained handbook explicitly identifies23Na100% abundance and the
NaI3s2S1/2/NaII1S0 ground labels. Modern NIST isotope and CODATA originals supply
the masses/constants; the old handbook mass/rounded energy are not used.
The actual ASD5.12 neutral row supplies its first IE and reported uncertainty;
the next row's NaII ionization energy is not substituted. A thin native adapter
calls shared `ReducedModels::IonPhysicalMass`, and independent Decimal checks
the reference/binding/physical masses within the existing two-native-epsilon
criterion. The physical ion mass is separately stored from the additive
reference mass. Original uncertainty strings and first-order contributions
are retained; no unspecified ASD confidence level, covariance or physical
acceptance tolerance is invented.

From the checkout, prepare a NEW profile and separate NEW work directory:

```sh
python3 -B utility/chemical_database/state_profiles/prepare_atomic.py \
  --raw-dir test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/raw \
  --decision-register srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json \
  --energy-file test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/raw/nist_na_ionization_uncertainty.html \
  --label-file test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/raw/nist_na_handbook.html \
  --help-file test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/research/raw/nist_ie_help.html \
  --canonical-root /data/vtenishe/moon_validation_data/photochemistry \
  --output-dir test_output/srcMoon/NEW-atomic-profile \
  --work-dir test_output/srcMoon/NEW-atomic-profile-work
```

The work directory preserves original default preparation, compiler commands,
source snapshots, executable and native logs. The profile contains the same
three CSVs plus `atomic_mass_derivation.json`, `na_source_energy_crosswalk.json`
and its exact eight-column CSV. All source paths/hashes and native state/record
anchors are retained. No original source bytes are copied into a managed archive.
`profile.json` is written last after source/state/independent arithmetic checks.
Missing/corrupt inputs never produce a usable profile.

The crosswalk's `source-effective-photon` policy defines
`Egamma=I_D04+mean_excess_D04`, retains actual native threshold/rate strings and
records `I_physical-I_D04` and its opposite signed translational correction.
The exact `physical_threshold_source_identity` equals the formation-source
identity returned by the production native loader. No legacy threshold is
overwritten, no source rounding interval guessed and no source-consistency
tolerance broadened. A small threshold shift is not an absolute-rate error
bound; that remains unavailable without the actual spectral convolution.
The application must explicitly implement/admit this crosswalk and its other
physical controls before a linked reaction can run. This preparation alone
does not qualify a physical event, source model, convergence or I21.

Run all new supporting controls without an AMPS application build:

```sh
python3 -B utility/chemical_database/state_profiles/run_atomic.py \
  --output-dir test_output/srcMoon/NEW-atomic-profile-controls
```

Optional source arguments match the preparation CLI and relocate actual retained
inputs; unavailable originals give SKIPPED77, never substituted science data.
The runner exercises the actual profile loader, `ResolveState`, `FormationQ`,
shared mass kernel, exact CSV/source identities, signed crosswalk, independent
Decimal arithmetic, old unavailable states, corruption/no-output refusal and
read-only source invariants. Result/status/count/log/source/resource manifests
are retained. A supporting PASS leaves formal M3 SKIPPED/NOT QUALIFIED.
