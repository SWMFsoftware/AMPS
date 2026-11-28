# Automatic fitted-source export and 2026-11-11 campaign

Status: requested implementation, 124-fit execution, native integration checks
and Word/PDF delivery are complete. Physical/identifiability limits and the
one marginal numerical stationarity diagnostic remain explicit. Git remains
read-only.

## Final delivery

Campaign root:
`/data/vtenishe/cassini_validation_data/studies/fitting_compain_111126`.

- `report/Enceladus_fitting_campaign_111126.docx`: 12-page Word report, eight
  embedded comparison figures, all target-rate tables, methods, input sample
  counts/kinematics, calibration/geometry/identifiability limits and validation.
- Adjacent `.pdf` and self-contained `.html`; `report/figures/` contains all
  eight PNG/PDF pairs comparing every strategy with native observed profiles.
- `report/all_target_rates.csv`: 124 targets (64 H2O/UVIS, 60 CO2).
- `report/all_component_rates.csv`: 7843 canonical component rows, including zeros.
- `report/export_inventory.csv`: 224 source-file paths/hashes: 124 conditional
  shape files and 100 absolute files; 24 E3/E5 CO2 absolute files are deferred.
- `report/actual_pooled_training_losses.csv`: 20 original balanced losses and
  KKT diagnostics. The 99-source H2O INMS-only result has normalized residuals
  about 8.6e-10, above the 1e-12 audit threshold; it is explicitly flagged,
  not silently replaced or claimed to be fully stationary.
- `report/candidate_objective_input_audit.csv`: actual 98/99 target/uncertainty
  comparisons. UVIS proxy values/weights differ with independent baseline
  calibration; their pooled losses do not constitute one nested objective.
- `report/manuscript_rate_comparison.csv`: 93 comparisons with the rounded
  manuscript Tables 2–5. These are historical values, not test references.
- `execution/`, `validation/`, `provenance/`: exact commands, successful and
  failed-attempt logs, receipts, paper, input products and software snapshots.

Final checks: **62 standalone PASS; 56 global cases / 2 suites PASS; six
native startup gates PASS; 224 production source validations PASS; eight
one-rank native runs PASS**. Every evolving smoke run completed 20 steps,
injected exactly 200 particles and closed its ledger at zero. The final
global receipt is
`test_output/enceladus-fitting-111126/global/run-20261111T073812Z-235/amps-core-tests.json`.
The generated global/manual catalogs pass their checks (56 suites / 1193
cases). The final clean native build repeated the mandatory preflight after
the subnormal parser fix; full logs are `native-final-build/` in the local
evidence root and copied under the campaign's `validation/native_build/`.

LibreOffice conversion uses local temporary staging/profile directories;
some installations cannot directly use the observation mount for that state.
Figures are embedded as data URIs in HTML and inside the DOCX, not external
links. The DOCX ZIP integrity/media count and PDF layout were checked; the
final document has no cropped figures or isolated component-table rows.
Prior failed conversions/layout drafts are retained under `provenance/`.

Reproduction from the repository root:

```bash
python3 srcEnceladus/standalone/run_fitting_campaign.py --validation-data-root /data/vtenishe/cassini_validation_data --paper /home/vtenishe/Downloads/Enceladus.pdf --output-dir /absolute/path/to/NEW_campaign --jobs 3
python3 srcEnceladus/standalone/build_campaign_report.py --campaign /absolute/path/to/NEW_campaign --validation-receipt /absolute/path/to/actual_validation.json
```

The source reader accepts all four cases; choosing a forward-run source
file/strategy remains explicit. No checked-in equal-rate source fixture or
default run deck was replaced, and no new long AMPS flyby comparison is
claimed by these short native checks. Scientific limits listed below remain
part of every interpretation of the campaign.

### Evidence checkpoint (2026-11-11)

- The three base geometries completed all 93 requested target fits. The E17
  search completed and selected the same radial candidate coordinates as the
  packaged input; the full 99-source H2O/CO2 refits are running.
- Automatic exports, exact fit/export stripe coordinates, all canonical
  memberships and negative speed tails are implemented. Published jet
  `porco_085` uses azimuth 360 degrees; export records the identical 0-degree
  representation required by the half-open runtime grammar. Failed first
  exports and their logs are retained; saved fit amplitudes were not changed.
- The 98-source E3/E5-pooled CO2 TRF solve stagnated. Active-set NNLS recovery
  on the identical weighted matrix passes objective nonincrease and normalized
  KKT gates <=1e-12. Its original failure, retry directory and solver-code hashes
  are retained. No objective, data weight, source constraint or rate prior changed.
- Standalone checks: 61 PASS before the final single-profile CO2 deferral
  regression was added. Selected global suites: 56/56 PASS, 2/2 suites PASS,
  no SKIP/FAIL/ERROR. Native clean/configure/-j16 build: PASS. All six native
  startup gates: PASS. Six actual one-rank base-geometry H2O/CO2 smoke runs:
  PASS; each completed 20 steps, injected exactly 200 particles and closed its
  particle ledger at zero. These are native integration checks, not plume
  convergence/observation validation.
- Final campaign audit, all comparison figures, Word/PDF report and final
  evidence consolidation remain pending at this checkpoint.

Evidence root: `/home/vtenishe/Mars2/AMPS/test_output/enceladus-fitting-111126/`.
Selected global receipt:
`global/run-20261111T072701Z-3/amps-core-tests.json`.
Native receipts: `native-startup/test_receipt.json` and
`native-catalog/<case>/<species>/receipt.json`.

## Requested deliverables and scientific decisions

- Implement automatic fitted-result to schema-2 AMPS source export.
- Extend the AMPS reader/exporter to hybrid 8+7, Porco 98, Portyankina 41 and
  the explicitly selected Porco 99 E17-candidate variant.
- Repeat the complete fit matrix in
  `/home/vtenishe/Downloads/Enceladus.pdf` (38 pages), preserving its numerical
  method and reporting physical/identifiability limitations.
- User explicitly selected **projected 252100 m stripe geometry in both the
  new fits and exports**. Old results and default legacy configurations remain
  intact; their rate differences are scientific evidence, not test failures.
- User explicitly selected **defer absolute E3/E5 CO2 export** until the
  fragmentation correction is chosen. Preserve apparent-rate diagnostics and
  normalized conditional shapes. Never silently select a factor in [2,5].
- Save the complete campaign, comparisons, rates, source files, provenance and
  a Word report at the exact requested directory:
  `/data/vtenishe/cassini_validation_data/studies/fitting_compain_111126`.

The draft contains revision markup and interleaved removed/inserted text.
Its original PDF/hash and extracted text are retained under `provenance/`.
The executed case manifest, numerical outputs and declared caveats determine
what was actually repeated; an old paper table is not a generated reference.

## Fit matrix

| Family | Geometry cases | Strategies | Targets | Target fits |
| --- | --- | --- | --- | ---: |
| H2O | hybrid, Porco 98, Portyankina 41 | independent, all-data shape, INMS-only shape, UVIS-only shape | E3, E5, UVIS 2007, UVIS 2010 | 48 |
| CO2 | same three bases | independent, E3+E5 shape, all-five shape | E3, E5, E14, E17, E18 | 45 |
| Explicit E17 candidate | Porco 99 | same four H2O and three CO2 strategies | same targets | 31 |

Also retain instrument-combined fits, component ablations, cross-predictions,
quadrature checks, orbital-proxy diagnostics and E17 location/window/direction
sensitivity. CH4 is not a fitting family in this H2O/CO2 manuscript.

## Implementation and validation work

1. Generalize the explicit exporter and C++ receipt reader without weakening
   canonical membership, duplicate/unknown-key, shape-sum or surface guards.
2. Add a selected projected-geometry option to standalone configuration;
   export the same resolved geometry/kinematics that produced each fit.
3. Integrate automatic shape/absolute source files and an index/provenance
   receipt into single-profile and complete-study fitting entry points.
   Exclude empirical backgrounds and do not select a preferred strategy.
4. Qualify negative shifted-Maxwellian speed tails needed by tilted catalog
   sources. Preserve the old recurrence in its existing range and independently
   compare the extension with quadrature and inverse-CDF references.
5. Execute application-independent tests and native reader/sampler smoke gates
   through required clean/configure/-j16 preflight. Do not overlap native
   builds with other make/test activity or call a source-only compile native.
6. Run all declared fits into new directories, compare outcomes and construct
   a self-contained Word report with plots, totals, component tables,
   diagnostics, exact commands, hashes and limitations.

Known inherited limits remain explicit: steady collisionless/gravity-free
fitting, E3/E5 digitized outbound water profiles, correlated calibration,
ill-conditioned source partitions, pooled-template versus free-amplitude
shared-shape distinction, and the approximate 2007 ray / 308 outward 2010
rays. Do not claim that this campaign repairs those observation operators.

Local preparation and source backups:
`test_output/enceladus-fitting-111126/preflight/`.
Final execution receipts and test outcomes will be recorded here and in
`STATUS.md` as implementation progresses.
