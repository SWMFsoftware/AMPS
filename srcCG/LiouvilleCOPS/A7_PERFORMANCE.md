# Stage A7: fitter performance characterization

This document records the Stage-A7 compute envelope for the degree-5--8
activity-enhancement fitter.  It describes measured implementation cost, not
fit convergence or scientific validity.  Timings and memory are specific to
the host, libraries, filesystem cache, selected data, and static base named in
the report; rerun the benchmark before a large campaign on another machine or
with another surface mesh.

## What the benchmark measures

`tools/benchmark_fit_performance.py` uses the production readers and the
production `Model`, `Problem.residuals`, and analytic `Problem.jacobian` from
`tools/fit_cops.py`.  Each harmonic degree is measured in an isolated process
so Linux maximum-resident-set-size high-water marks do not leak from one
degree to the next.  It reports:

- cache validation and model-startup time;
- construction of the facet harmonic basis;
- median residual and Jacobian evaluation time;
- process peak RSS and a live NumPy-array accounting estimate;
- an equal-work, fork-based scaling test over candidate independent-fit job
  counts;
- the size of a complete Stage-A5 product package made by a deliberately
  truncated one-evaluation production fit.

The temporary output probes are deleted.  They are not converged fits and
must not be used as scientific results.  JSON or Markdown is written only
when explicitly requested; relative report paths are resolved under the
normal `<AMPS>/cops_output` root, and existing files are refused.

The scaling test represents independent lambda fits.  It does not make one
SciPy fit parallel: a single lambda is serial at the fitter-process level.
Every worker sets OpenMP, OpenBLAS, and MKL to the requested numerical thread
count, normally one, so process and numerical-library parallelism are not
multiplied accidentally.

## Reproduce the measurement

Run from the LiouvilleCOPS directory, using a current manifested base and its
matching housekeeping table.  Absolute input paths avoid ambiguity and do
not alter the configured output root:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 tools/benchmark_fit_performance.py \
  /absolute/path/to/output_fit_full_default_base \
  --hk-table /absolute/path/to/COPS_HK_20160930.txt \
  --profile all_available \
  --degrees 5 6 7 8 \
  --repeat 3 \
  --jobs 1 2 4 8 14 \
  --threads-per-process 1 \
  --campaign-profiles 3 \
  --campaign-lambdas 5 \
  --campaign-max-nfev 1000 \
  --campaign-id a7measured \
  --json a7_performance/performance_YYYYMMDD.json \
  --markdown a7_performance/performance_YYYYMMDD.md
```

Use new report names for every measurement.  `--skip-output-probe` avoids the
four one-evaluation fits when only kernels are being timed.

## Measured 100,000-facet result

The following run used the full 1,150-point manifested SHAP5 base, the
`all_available` profile, one numerical thread per process, three timing
repetitions, and the candidate job counts 1, 2, 4, 8, and 14.  The host
reported 36 schedulable CPUs and about 242.5 GiB available memory when the
recommendation was computed.

| Maximum degree | Free residual coefficients | Ready/startup | Basis | Residual | Jacobian | Peak RSS | Output package |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 11 | 90.083 s | 0.0265 s | 1.9233 s | 5.8742 s | 11275.5 MiB | 32.4 MiB |
| 6 | 24 | 89.211 s | 0.0408 s | 1.9048 s | 6.4069 s | 11275.4 MiB | 32.4 MiB |
| 7 | 39 | 89.103 s | 0.0609 s | 1.9081 s | 7.1977 s | 11275.7 MiB | 32.4 MiB |
| 8 | 56 | 89.619 s | 0.0819 s | 1.9303 s | 8.1912 s | 11275.7 MiB | 32.4 MiB |

`Ready/startup` includes Python/import overhead, strict base validation,
element and kernel loading, the base-consistency calculation, selection, and
problem construction.  It is not merely basis time.  Peak RSS is the
per-process high-water mark after preparing the full fitter; it must not be
added to the dense-basis byte count as though both were independent.

| Independent jobs | Evaluations/s | Speedup | Parallel efficiency | Maximum worker RSS |
|---:|---:|---:|---:|---:|
| 1 | 0.099 | 1.000 | 1.000 | 11275.7 MiB |
| 2 | 0.184 | 1.862 | 0.931 | 11255.7 MiB |
| 4 | 0.306 | 3.098 | 0.774 | 11255.6 MiB |
| 8 | 0.527 | 5.334 | 0.667 | 11255.6 MiB |
| 14 | 0.627 | 6.348 | 0.453 | 11255.6 MiB |

The measured recommendation is `--jobs 8 --threads-per-process 1` for a scan
with at least eight independent lambdas.  Four jobs are appropriate for four
lambdas; one lambda cannot use fitter process parallelism.  Fourteen jobs had
higher aggregate throughput but only 45.3% efficiency and therefore failed
the declared 60% efficiency gate.  This host-specific finding does not change
the portable runner default selected earlier; pass the measured setting
explicitly for this host.

The benchmark also applies a conservative memory gate:

```text
jobs * maximum_worker_peak_RSS <= 0.70 * MemAvailable
```

Fork copy-on-write can make actual aggregate RSS smaller than this bound, but
the conservative check avoids promising memory sharing that a later write or
library change could remove.

For the declared planning case of three profiles, four degrees, five lambdas,
and 1,000 optimizer evaluations, the script reported 147.53 CPU-hours of
residual-plus-Jacobian capacity, 27.90 hours at the recommended measured
throughput, and 1.9 GiB of product storage.  This is deliberately an upper
capacity scenario: SciPy call counts, convergence, retries, C++ verification,
plot aggregation, and filesystem contention are not predicted by it.

## Output naming and safety

Campaign names have the bounded form
`a8_<profile>_l<degree>_<campaign-id>`.  The tool rejects unsafe campaign IDs,
and the profile runner and fitter refuse an existing destination rather than
merge or overwrite it.  Keep the benchmark JSON with the campaign so the
resource decision is auditable.

## Shape resolution and performance

The current external `SHAP5_Volume.bdf` has exactly 100,000 `CTRIA3` faces.
LiouvilleCOPS can directly read either a NASTRAN `CTRIA3` surface or a
single-segment SPICE type-2 DSK whose frame equals the configured `BodyFrame`;
DSK vertices are converted from kilometres to metres.  PLY, OBJ, WRL, STL,
and VTK are not supported by the current reader.

NASA PDS Small Bodies Node collection
[`RO-C-MULTI-5-67P-SHAPE-V2.0`](https://pdssbn.astro.umd.edu/holdings/ro-c-multi-5-67p-shape-v2.0/dataset.shtml)
is the preferred downloadable archive. Its product inventory includes these
useful DSK levels in a container format supported by LiouvilleCOPS. Format
compatibility alone is insufficient: each candidate remains subject to the
single-segment, exact-`BodyFrame`, closed-topology, and consistent-winding
checks described below:

| Family | Archived DSK plate counts above the current mesh | Practical role |
|---|---|---|
| SPC SHAP5 | 199k, 384k, 788k, 1M, 3M | Best controlled resolution study because it keeps the SHAP5 family; start with 199k, then 384k, after verifying the embedded frame. |
| SPG SHAP4S | 200k, 1M, 4M | Independent shape reconstruction and resolution sensitivity, not just refinement of the active SHAP5 file. |
| MSPCD SHAP2 | 191k, 391k, 760k, 1M | Independent-method sensitivity test; its geometry/version difference is scientifically confounded with plate count. |
| ESA SPC MTP019 | archived high-resolution DSK | Mission geometry alternative; read its PDS label and verify plate count/frame before use. |

The archive also supplies lower resolutions and non-DSK variants.  Each
product has its own citation in the adjacent PDS label; preserve that label
and file hash in the run record.

The five SPC SHAP5 candidates were downloaded on 2026-10-18 to
`/data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2`.
Each binary has its PDS label. CSPICE N0067 `dskbrief -full`, with the mission
frame kernel, verified one segment, data type 2, and reference frame
`67P/C-G_CK` in every binary.

They were then converted with `tools/convert_dsk_to_bdf.py` into the
metre-unit `GRID*`/`CTRIA3` surface convention consumed from the original
SHAP5 BDF. LiouvilleCOPS ignores the historical BDF's `CTETRA` volume cells,
so no volume elements are invented. The raw lower four DSK reductions contain
small numbers of exact reversed duplicate triangle pairs: zero-thickness fins
that make an attachment edge non-manifold. The converter removes both copies
of only those exact pairs and rejects any remaining open, non-manifold, or
inconsistently wound edge.

| Resolution/file stem | DSK vertices | DSK plates | Removed fin pairs/plates | Retained BDF triangles | BDF bytes |
|---|---:|---:|---:|---:|---:|
| `cg_spc_shap5_199k_cart` | 99,729 | 199,452 | 2 / 4 | 199,448 | 20,344,531 |
| `cg_spc_shap5_384k_cart` | 191,997 | 383,988 | 13 / 26 | 383,962 | 39,166,127 |
| `cg_spc_shap5_788k_cart` | 393,963 | 787,920 | 83 / 166 | 787,754 | 80,360,332 |
| `cg_spc_shap5_001m_cart` | 597,251 | 1,194,496 | 123 / 246 | 1,194,250 | 121,827,165 |
| `cg_spc_shap5_003m_cart` | 1,572,866 | 3,145,728 | 0 / 0 | 3,145,728 | 320,864,720 |

The model's actual NASTRAN reader and strict `FinalizeGeometry` check loaded
all five converted files as closed, consistently oriented surfaces with zero
open edges. Every `.bdf.conversion.json` records source and output SHA-256,
counts, unit scale, NAIF commands, PDS label identity, and exact removed input
face IDs. The directory `SHA256SUMS` covers DSKs, labels, BDFs, and manifests.

Select a converted file in a copied input with, for example:

```ini
ShapeModel = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bdf
ShapeModelFormat = nastran
BodyFrame = 67P/C-G_CK
```

Use a new `OutputDirectory` and regenerate the complete static fitting base;
do not pair old geometry tables or kernels with the new mesh.

Much finer models exist but require conversion or access coordination:

- DLR SPG SHAP7 has 44 million facets at roughly 1--1.5 m point spacing and
  is available by request; DLR also describes 500k textured PLY reductions.
- MiARD distributes binary-PLY SHAP8 variants with 12, 20, and 44 million
  facets.
- CNES describes a 132-million-facet binary-PLY model available by request.

Those products are not drop-in LiouvilleCOPS inputs.  A converter would have
to preserve units, orientation, topology, body-fixed frame, and provenance,
and the converted mesh would need the same geometry and zero-enhancement
validation gates as a native DSK.

Do not infer high-resolution feasibility by plate count alone.  Static
point--facet transport caches and startup work grow with the number of
contributing facets, and the measured 100k base already uses about 11 GiB per
prepared fitter process.  The code has not measured a 1M, 3M, or larger
fitting base, so no memory or runtime number for those meshes is asserted
here.  Benchmark a 199k base first, with `--jobs 1`, and inspect actual cache
sizes and RSS before increasing either resolution or worker count.

Finally, a finer geometric mesh changes facet positions, normals, shadowing,
visibility, altitude, and near-surface transport.  It does **not** make the
degree-5--8 activity parameterization more local: spatial source-map
resolution is controlled separately by the harmonic degree/basis.  Compare a
matched-resolution sequence with the same data profile, physics, fitting
controls, and activity basis before attributing a pressure mismatch to local
topography.

## Interpretation limits

- Timing results do not establish optimizer convergence.
- A one-evaluation output probe is not a scientific fit.
- Fork scaling is not OpenMP speedup inside one fit.
- Peak RSS is host- and allocator-dependent.
- The campaign estimate excludes any work not named by its formula.
- A finer mesh can change geometry and visibility but cannot by itself prove
  that a near-surface pressure jump is topographic.
