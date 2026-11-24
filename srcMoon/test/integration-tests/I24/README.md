# I24 — Charge-exchange reactions with solar-wind ions

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Direct rate error <=1%; ensemble lifetime error <=3%; correct products and exact charge/mass conservation under the reviewed physical convention; unsupported channels OFF. U17 retains >=100000 probabilistic-branch events.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.003 numerical bias and fixes endpoint-survival lifetime
as primary. Complete source/state products, U17's existing event minimum and
actual native/numerical/physical-ledger/predecessor evidence remain required.
An experimental partial channel cannot replace this complete formal scope.

The [angular-source follow-up](ANGULAR_SOURCE_FOLLOWUP_20261124.md)
records published ground-product scattering at higher energies,
the 4 keV/u rate discrepancy, and the scope of the conditional
isotropic center-of-mass sensitivity case.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I24/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I24/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.

The [first-channel source review](SOURCE_ADMISSION_20261124.md) identifies
the retained D119990 gaps, the state-resolved Na calculation candidates
and a higher-energy proton–He alternative with documented lab-frame angle
measurements. A source-matched angular law or independently bounded
approximation is required for linked transport.

The r11 review traces D119990 to Knoop et al. Fig. 2 using the IAEA v2024.1
release and NIFS CHART record 38917. That figure groups capture by target
shell. Reject its twelve knots as an H(1s)-resolved measurement, while
preserving the original record and its label. Shingal et al. provide
separate theoretical ground-state cross sections near the controlled
energy. Their basis dependence and missing ground-product angular data
still limit transport use. [Source evidence](../../../../utility/chemical_database/NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md)

### Selectable database reference for modeled recoil

Use `hplus-na3s-ground-capture-4kev-u-isotropic-v1` in the r14 shared
catalogue for the 4.0 keV/u single-channel development case. The integral
cross section is atomic theory. The solid-angle distribution is the approved
sensitivity assumption. Use the independent two-body energy/momentum oracle
for both products. A separate model, `poppe2022-na-pickup-median-v1`, supplies
total Na⁺ production and zero lunar-frame birth speed for a bulk pickup-flux
comparison. The two represent different state and energy selections. Publication
of the selections is a database result. Formal I24 requires the linked
rate, ensemble, conservation and predecessor evidence in this test contract.
[Physics, units and source queries](../../../../utility/chemical_database/DATABASE_README.md#sodium-pickup-and-binary-charge-exchange-r14)

### Running the 4.0 keV/u MPI provider comparison

From the AMPS root, supply the release lock and request the supporting case
on its existing I24 owner. The runner locates the NIST reference mass inventory
in the retained acquisition; use `--i24-mass-inventory /path/to/reference_mass_inventory.json`
to point to a relocated copy together with its four `raw/` files.

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
out=$(mktemp -d)
python3 -B srcMoon/test/run_all_test.py --only I24 --skip-build \
  --campaign-mode regression --mpi-ranks 8 --i24-sodium-lock "$db/database-lock.json" \
  --test-runtime-seconds 3600 --output-dir "$out"
```

In `$out/I24/sodium-ground-capture/`, inspect `result.json` for the supporting
score, `comparison.csv` for both rank layouts, and `comparison.png` or
`comparison.eps` for the survival and recoil-direction CDFs. Per-parent native
histories are in `mpi-1/` and `mpi-8/`. A linked MPI adapter samples the
production event clock and calls the provider for each actual hit. It checks the
two transported product packets and the consumed incoming projectile. This
case covers the event plan before PIC allocation or transport. The I24 owner
reports the supporting outcome separately from its formal predecessor and
physical-coverage status. [Physics, source assumptions and acceptance](SODIUM_CASE.md)

### Running both the MPI cohort and the linked PIC movers

Set `producer` to the `normal_build_evidence.json` written by the
source-qualification builder for the force-free Na, Na⁺ and H build.
The receipt must still point to its original library, generated headers,
executable, build log and state mapping. The runner checks their hashes.
An unrelated Moon executable cannot supply this species layout.
If a reboot removed the disposable configured library or generated build
directory, the PIC comparison reports SKIPPED and requests a clean producer.
Changed bytes in the retained evidence or any available build member give
ERROR, including when another disposable build member is missing.
[Build the Na, Na⁺ and H producer](../../source_qualification_build/README.md#build-the-i24-sodium-capture-case)

```bash
db=/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261124-r14
producer=/path/to/source-qualification-build/normal_build_evidence.json
out=$(mktemp -d)
python3 -B srcMoon/test/run_all_test.py --only I24 --skip-build \
  --campaign-mode regression --mpi-ranks 8 --test-runtime-seconds 3600 \
  --i24-sodium-lock "$db/database-lock.json" --i24-pic-producer "$producer" \
  --output-dir "$out"
```

The added PIC matrix runs one controlled event with the first-order and
second-order production tracing movers, each at one and eight MPI ranks.
The matrix checks the native incident-proton ledger, two allocated daughter
particles, their species weights, and their final clocks, velocities and
positions after continuing the unused part of the parent step. The
independent 32,768-parent MPI cohort supplies the rate, survival and angular
checks. Both children share one total wall allowance. The PIC matrix always
uses 1/8 ranks. `--mpi-ranks` selects the larger layout of the cohort.

The PIC records are in `$out/I24/sodium-pic-ground-capture/`. Its
`result.json` lists each mover/layout and its named checks. Directories
`first-1mpi/`, `first-8mpi/`, `second-1mpi/` and `second-8mpi/` hold the
individual input decks, event, product and final-particle CSVs by rank.
The I24 owner binds the hashes of both child result/manifest pairs and
retains a separate CLI log for each. A missing build or release is a
supporting SKIPPED, altered input or a launch error is ERROR, and a
measured mismatch is FAIL. An ERROR or FAIL propagates to the aggregate
even when the other child succeeds. Both supporting checks can pass while
formal I24 still awaits its predecessors and remaining physical coverage.

### Time-step and particle-weight refinement

Add `--i24-pic-refinement` to the aggregate command above. For the standalone
PIC runner the switch is `--refinement`. The four original mover/rank cases
remain the default. The optional matrix reuses its second-order one-step
results as anchors and runs eight more cases.

| Series | Production steps from the same source epoch | Reference weights (Na, Na⁺, H) | MPI ranks |
|---|---|---|---|
| Time step | 1, 2, 4 equal calls to the normal mover | (7, 56, 0.875) | 1, 8 |
| Reference weight | 1 | (7, 56, 0.875) multiplied by 1, ½ and ¼ | 1, 8 |

All cases use the same parent key, original beam density, 4.0 keV/u source
point, represented parent population of 2.625, and final time. Halving the
reference weights doubles the individual statistical corrections. The
independent exponential clock fixes the absolute event epoch. The parent
survives each earlier complete step and reacts during the final step. Both
daughters then coast for approximately 1/16 s. The clock oracle separately
checks elapsed time within the reacting step and elapsed time since injection.
The trajectory check uses the independently calculated remaining flight and
the native birth velocities. Component-wise comparisons across refinements
allow the rounding error of storing absolute ET in binary64.

`refinement.csv` records each tested level, native and reference clock
metrics, maximum ballistic displacement error, and energy/momentum
residuals. Each `steps-rank-*.csv` records actual post-exchange particle
counts, represented populations, committed event counts and physical end
time after every step. This comparison verifies one conditioned,
force-free collision and its two descendants. It does not estimate an
ensemble lifetime or qualify I24 across the solar-wind energy distribution.

Each completed PIC layout also writes `continuation.csv`, `continuation.png`
and `continuation.eps`. These compare the two live products' displacements in
the frozen LSO frame with analytical force-free flight over the unconsumed
part of the step. The native per-reaction charge changes are integers in
units of the elementary charge. Multiply each by its represented population
to obtain packet charge. The run manifest records the actual loaded SPICE
kernel paths and hashes. [Source-bound build command](../../source_qualification_build/README.md#build-the-i24-sodium-capture-case)
