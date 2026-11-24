# I22 — Photoionization for He, Ne, Ar, and future H2O/OH

The [atomic-family preparer](../../../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
emits six explicit source-effective photon candidates, one for each
isotope and source endpoint. Each retains the D04 printed threshold and
mean, physical ground-ion energy and native row identity. The worksheet
supports source/state review before recoil or ground-ion yields are selected.
The U15/U24 and full linked I22 gates remain applicable.

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Ensemble rate error <=2%; linked branching statistically consistent at95% confidence; complete declared family sums to1 within numerical precision. U15 separately retains99% intervals and source-precision transcription.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
fixes total-removal endpoint rate `-log(S/N)/T` as primary, relative0.002
numerical bias, simultaneous95% exact branch intervals and10% relative
half-width for every positive source branch. Partial rates/event-time MLE are
separately named diagnostics. The complete family and full U24 dependency stay
required. Protocol approval does not supply missing physical/native evidence.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I22/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I22/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.

The r11 shared source packet supplies 1 au solar-cycle frequency fits for
all three charged water branches and laboratory checks of argon
fine-structure splitting and OH uptake on ice. The rates share PHIDRATES
cross-section lineage with D04. I22 still requires solar-weighted populated
ionic states and joint product kinematics. Keep the 800 argon posterior
samples distinct from independent measurements.
[Photon sources](../../../../utility/chemical_database/NATIVE_M3_PHOTON_BENCHMARK_SOURCES.md)
[OH source](../../../../utility/chemical_database/NATIVE_M3_OH_ICE_SOURCES.md)

Kronebusch and Berkowitz supply a 21.2 eV laboratory benchmark for the
three water fragment-ion yields, relative to total ionization. The
instrument discriminates against ions above 1 eV. The experiment gives a
finite-energy check of the relative signals. Electron, ion and neutral
coincidences and the solar-spectrum average require additional evidence.
[ASMS proceedings, printed page 114](https://www.asms.org/docs/default-source/proceedings-archive/1975_asms_23rd_conference.pdf#page=117)

### Conditional products for charged water breakup

`ReducedKinematics::PhotonChargedBreakup` constructs one correlated electron,
ion and neutral from a specified parent at rest. It checks the selected
states' isotope counts, integer charge and common mass reference, computes
formation/excitation Q, includes the incoming photon momentum and removes the
center-of-mass kinetic energy. The caller supplies the fraction of the
remaining energy in electron-versus-heavy-pair motion, the electron direction
relative to that pair, and the ion direction relative to the neutral. Both
directions use the same inertial axes. Every positive-energy mode is checked
after reconstructing the delivered velocities. The three output occurrences
retain their identities and the source of the supplied partition.

The deterministic test sets those inputs explicitly and verifies all three
velocities and both relative energies independently. Its nonzero center-of-mass
and threshold cases exercise drift, endpoint sharing and photon recoil. From
the AMPS root run the following with a fresh output directory.

```sh
python3 -B srcMoon/test/test_m3_reduced_kinematics.py --output-dir test_output/srcMoon/i22-kinematics-NEW
```

The output includes `comparison.csv`, `result.json`, the build and run manifest
and each product's native momentum. Physical activation of H₂O04–06 requires
the channel-conditioned states, energy sharing and joint directions at each
solar endpoint, or a source-supported approximation with its declared error.
Those inputs are external to the kinematic construction. The same-campaign
predecessors and U24 surface requirements continue to apply.

### Scope of the r14 modeled inputs

The authenticated r14 catalogue contains three EXOSPHID noble-gas records.
Each produces one singly charged ion. The ion has no resolved electronic
state. This model can benchmark rates and product counts. Broadband
ground-state ion recoil still needs a common spectral weighting and a
state-energy or omitted-excitation bound at both D04 endpoints. The NORAD
ground and total He cross sections give a possible finite-band comparison.
For charged H₂O breakup, EXOSPHID supplies band-limited product fractions
and heavy-product velocities. The Europa/Callisto tables supply per-channel
mean energies. Neither source selects the electron–heavy-pair energy sharing
and joint directions for a **single** conservative event. The implemented
Jacobi kernel accepts an independently specified partition for that purpose.
[Cross-model evidence](../../M3_LUNAR_MODEL_REFERENCE_CROSSWALK.md)

U24 now evaluates the fully specified Laferriere residence scenario for
already-bound OH/H₂O against the lunar rate and finite-time release kernels.
Jones' returning-radical chemisorption, Smolka's thermal flux and Laferriere's
mobile hopping represent distinct surface processes. A full I22 case needs
a consistent selection for the impacted material, its scope or discrepancy
bound and native species, reservoir and gas chemistry dispatch. The r14
source-bound configurations for **Na pickup and proton–Na capture** make
no selection for these water processes. [U24 residence verification](../../U24_water_chemistry/README.md)
