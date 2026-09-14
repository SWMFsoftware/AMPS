# Step 15 change manifest: scientific validation

Step 15 introduces numerical, cross-mover, cross-model, coupled-background,
and observational-readiness evidence without weakening the three-mover public
scope established in Steps 5 and 14.

## Production validation descriptors

- `util/sep_scientific_validation.{h,cpp}` registers `VAL01`–`VAL03` with the
  existing component-test result contract.
- `component_tests.cpp` exposes those stable IDs through the production CLI.
- `makefile` archives the descriptors into `mainlib.a`; no fourth mover or new
  coefficient authority is added.

## Focused runners and independent references

- `VAL01` checks Parker diffusion and cooling against closed-form solutions.
- `VAL02` compares `fte-dmumu` and `fte-mfp` under a matched mean free path.
- `VAL03` contains a separately coded conservative pitch-angle solver.
- `VAL04-SWCME` feeds actual SWCME adapter records through the production
  srcSEP background-snapshot and Parker-core boundaries.
- `test/run_step15_tests.sh` runs these cases with strict warnings, ASan, and
  UBSan, then checks the release gate.

## Campaign and external evidence

- `validation/campaign.schema.json` defines the release report.
- `validation/run_campaign.py` authenticates reports and external input bytes,
  captures source/compiler/schema provenance, separates evidence classes, and
  writes JSON/Markdown transactionally.
- `validation/manifests/` provides incomplete SWMF and spacecraft templates;
  templates are instructions, never passing evidence.
- `validation/cases/` records equations, exact configurations, metrics, and
  scientific limitations for `VAL01`–`VAL04`.

Real SWMF replay and held-out observations are intentionally not bundled or
claimed. `--release` remains nonzero until those checksum-verified records are
supplied and all required observational metric families close.

The validation-gate enhancement adds independent Make targets for controlled
analytical, native AMPS, real SWMF, and held-out observational evidence.
`validation/check_external_gate.py` validates SWMF or observational manifests
without merging their status; the observational command additionally requires
complete metric-family coverage across the supplied held-out events.
