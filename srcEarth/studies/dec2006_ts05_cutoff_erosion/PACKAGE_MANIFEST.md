# Package manifest

## Intended Git location

Place this directory at:

```text
srcEarth/studies/dec2006_ts05_cutoff_erosion/
```

All runners share the repository-level default output root:

```text
test_output/dec2006_ts05_cutoff_erosion/
```

The location is resolved by finding the nearest ancestor that contains the
`srcEarth` directory. When the package is used outside an AMPS checkout, the
launch directory is used as the fallback codebase root. Passing
`--output-root PATH` overrides the default without changing any scientific
configuration.

## Included scientific inputs and observations

| File | Records | Purpose |
|---|---:|---|
| `data/drivers/ts05_dec2006_5min.txt` | 865 | Five-minute IGRF+TS05 event driver |
| `data/observations/pamela_table_s1.csv` | 259 | PAMELA orbit-scale cutoff latitude |
| `data/observations/poes_metop_meped_boundaries.csv.gz` | 3,904 | POES/MetOp MLT-resolved boundary cells |
| `inputs/AMPS_PARAM_DEC2006_475km.in` | 1 | Common-rigidity morphology at PAMELA altitude |
| `inputs/AMPS_PARAM_DEC2006_850km.in` | 1 | Common-rigidity morphology at POES/MetOp altitude |

Expected hashes, row counts, temporal coverage, and source descriptions are in
`data/provenance.json`. The `vendor/C9` and `vendor/C10` directories contain
the observation-equivalent runners and source-preserving reconstruction tools.

## Runner and analysis programs

| Program | Output subdirectory |
|---|---|
| `scripts/run_study.py` | Entire shared run tree |
| `scripts/run_morphology.py` | `morphology/` |
| `scripts/compare_observations.py` | `comparison/` |
| `scripts/analyze_dynamics.py` | `dynamics/` |
| `scripts/make_figures.py` | `figures/` |
| `scripts/run_sensitivity_suite.py` | `ts05_sensitivity/` |

## Verified invocation

From the AMPS repository root:

```bash
python3 srcEarth/studies/dec2006_ts05_cutoff_erosion/scripts/run_study.py \
  --profile SMOKE --prepare-only --amps ./amps -np 4 -nt 16
```

This creates generated inputs and command inventories beneath
`test_output/dec2006_ts05_cutoff_erosion/` without launching AMPS. The same
command without `--prepare-only` executes the SMOKE calculation.
