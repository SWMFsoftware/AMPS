# Step-8 standalone campaign example

This three-epoch analytic-dipole example demonstrates the campaign manifest,
content-addressed references, input rendering, restart layout, AMPS product extraction,
cadence averaging, and machine-readable summaries. The two synthetic GOES P5 response
channels are identical, so their expected east/west ratio is exactly one. This is a
workflow/reference check, not an O1/O2 observation claim.

The manifest fixes the complete 18-column `gridless_points_flux.dat` schema and one
expected output row before execution. The extractor reads only the declared east and
west response columns, records dimensionless units, and fails if the producer changes
the variable order, row count, or values needed for the ratio. Every resource has a
provenance statement and digest; editing a resource requires deliberately updating its
manifest digest after review.

Run preflight and input rendering first:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --output-dir test_output/standalone_step8_example \
  --validate-only

python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --amps ./amps \
  --output-dir test_output/standalone_step8_example \
  --dry-run
```

For a linked standalone executable, run:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --amps ./amps \
  --output-dir test_output/standalone_step8_example \
  --restart
```

The example keeps the Step-7 unresolved gate enabled. Do not turn it off to make the
campaign pass. Increase the energy/angular resolution and add the roadmap convergence
JSON gates before using the pattern for a scientific event. Production O1/O2/O4
manifests must use role `VALIDATION`/`HOLDOUT`, explicit event IDs, frozen exclusions,
and the required 0.02 energy and angular gates; the loader rejects an incomplete or
weakened contract.
