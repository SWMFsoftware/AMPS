# Synthetic Stage 14 offline jobs

These JSON jobs are manufactured verification inputs. They use dimensionally
explicit toy values and frozen synthetic identities, not observations or an
AMPS production deck. Run from the AMPS root with a fresh output path:

```sh
python3 src/models/sep_coronal_cme/tools/research_stage14.py --job src/models/sep_coronal_cme/examples/stage14/wave_step.json --output test_output/stage14/wave_step.json
```

The same command accepts `wind_pass.json`, `reference_family.json`,
`renewal.json`, `impulsive_source.json`, `transfer_synthetic.json` and
`mhd_synthetic.json` by changing the explicit job/output path. Results retain
the full input-job identity and state that production-adapter and observed
campaign qualification are false. JSON identities are content hashes; after
changing a frozen inner asset, re-register it with `freeze_record` instead of
editing its old identity string.

`wind_hidden_overshoot.json` intentionally returns exit **1** with a finite
failed wind product. Both sampled endpoints satisfy the nominal bounds, while
the continuous speed and advective acceleration violate them internally.
Invalid provenance, unsupported inputs and existing output paths return **2**.
Successful supported processing returns **0**. A failed metric remains in the
result file; it is not discarded or relabeled PASS.

The [Stage-14 guide](../../docs/STAGE14_RESEARCH_EXTENSIONS.md) records each
implemented domain and the additional host, covariance, MPI and observational
qualification requirements. The examples do not enable schema-5 input
selectors or overwrite the shock-only production source.
