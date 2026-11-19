# Isotope/constants evidence for the approved M3 reference-mass convention

This preparation script acquires public NIST source responses into a new
research output directory and reproducibly extracts the explicitly selected
reference isotopes and SI constants. It never writes canonical data, edits a
raw response, changes compiled PIC masses, or supplies missing molecular states.
The contract is [the approved implementation contract](../M3_APPROVED_MODELS_IMPLEMENTATION_CONTRACT.md).

```sh
python3 -B srcMoon/test/physical_state_sources/prepare.py \
  --output-dir test_output/srcMoon/m3-state-sources-NEW --fetch
```

For an offline reproduction, omit `--fetch` and use `--raw-dir` pointing to
the earlier output's `raw` directory. Use a new output directory in both cases.
Failed responses and source-interpretation failures are preserved; no manual
replacement value or uncertainty is allowed. Exits are0 for successful source
extraction,2 for malformed/corrupt infrastructure,77 for unavailable acquisition.
This is documentary/source QA, not an I-series qualification or observation.

`reference_mass_inventory.json` preserves the source strings, uncertainties,
isotope identities, neutral masses and constants. Conversion to kg uses the
same acquired atomic mass constant, not a mixture-average mass. This file
supplies reference-mass input only; a complete state/formation/threshold
crosswalk and independently reviewed physical error bounds are still required
before event transport. Standard isotopic abundances in the source are not
adopted as lunar populations. The script records the native source hash and
retrieval URL/time, including unsuccessful requests, in `run_manifest.json`.
