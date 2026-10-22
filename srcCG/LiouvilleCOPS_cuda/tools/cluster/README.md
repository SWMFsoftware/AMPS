# Cluster tools (C15/P0: Voyager)

| File | Role |
|---|---|
| `cops_pipeline.py` | The campaign runner: one surface model and emission law, from the fitting base to the fits and the fine tuning, resumable from a restart directory |
| `cops_pipeline.slurm` | Template of the runner's batch scripts (`cops_pipeline.py submit` fills it) |
| `fit_mpi.slurm` | One fit over MPI ranks, for direct use (`sbatch --nodes=N --ntasks-per-node=<GPUs per node> ...`) |

Status: CPU-tested (fake tools, fake clock). A real run on GPUs and on Voyager is pending (STATUS.md C15, P0).

## The runner

```tcsh
# what would run, and what exists already (no GPU, writes nothing)
.venv/bin/python tools/cluster/cops_pipeline.py plan --restart-dir $SCRATCH/shap8_020m_flat \
    --shape .../cg-spg_mspcd-shap8-v0.4cmb-20Mfacets.bdf --surface-bc flat --max-nfev 2000 --fine-tune \
    --allow-genus 4
# inside an allocation: every incomplete stage, in order
.venv/bin/python tools/cluster/cops_pipeline.py run --restart-dir $SCRATCH/shap8_020m_flat [same options]
# progress; --check exits 0 (complete), 75 (resumable) or 1 (failed)
.venv/bin/python tools/cluster/cops_pipeline.py status --restart-dir $SCRATCH/shap8_020m_flat
```

The five campaign settings: `--shape` (surface model), `--surface-bc flat|rough` (emission law of the base),
`--restart-dir` (all products, stamps, logs and restart files), `--max-nfev` (evaluations per fit, default 2000),
`--fine-tune`/`--no-fine-tune` (default no). After the first call, `run --restart-dir DIR` alone repeats the stored
request; changing the mesh, the emission law, the base input, the gravity body or the precision request is refused
(use a new restart directory), and other fit settings need `--tag`.

| Stage | Group | Tool | Runs when |
|---|---|---|---|
| `base` | prepare | `../LiouvilleCOPS/tools/rebuild_cops_boundary.py` with `native/liouville_cops_cuda` (one GPU) | no `--base` |
| `store` | prepare | `tools/c15_record_store.py` | precision mixed, no `--store` |
| `visibility` | prepare | `tools/c05_write_bundle.py --visibility-only` | N30 requested |
| `gravity`, `cache`, `bundle` | prepare | `c04_gravity_equivalence.py gpu` (all GPUs of the node), `c04_write_gpu_cache.py`, `c05_write_bundle.py` | C30 or NC31 requested, no `--bundle` |
| `translate:<FIT>` | prepare | `tools/restart_tool.py translate` | `--restart-from` of another mesh with `--translate-restart` |
| `fit:<FIT>_n<N>` | fit | `tools/gpu_fit.py`: one node, local workers; several nodes, one MPI rank per GPU | always |
| `finetune:<FIT>_n<N>` | finetune | `tools/fine_tune_gpu.py` (GPUs of one node) | `--fine-tune` |

Restart directory:

```
pipeline.json         request, identity, hardware, precision decisions and policies, fit segments
stamps/<stage>.json   inputs, command, code digests, hardware, start/end, exit status (fits: segments)
logs/<stage>.log      timestamped tool output, every attempt
base/ store/ visibility/ gravity/ cache/ bundle/ translate/
fits/<FIT>_n<N>/segNN/   one gpu_fit output per segment; fit_restart.npz in each
finetune/<FIT>_n<N>_m<min>_i<n>_<method>/
slurm/                batch scripts written by submit
*.incomplete_<UTC>    earlier failed or killed attempts (moved aside, never deleted)
```

### Bases and stores built elsewhere

`--base DIR` and `--store DIR` (and `--visibility-bundle`, `--bundle`) use products copied from another machine.
They are checked, not rebuilt: `rebuild.json` complete with the same SurfaceBC and mesh (path, else SHA-256 of the
recorded mesh file, else its name), `run_manifest.json` process complete and cache verified and equal to the
digest `rebuild.json` recorded, the `elements.bin` header; the store's `store_manifest.json` must name this base's
`run_manifest.json` digest, FP32 records and the same facet count; a bundle's sidecar must list the needed columns
and the store's mesh identity (or the base path). The kernel files are never read for this.

To build on the cluster instead, the cluster needs the `../LiouvilleCOPS` tree, the SPICE kernels and the
observation file named in the base input (`--fit-input`, default `../LiouvilleCOPS/fit_base.in`; `plan` lists the
missing ones) and `native/liouville_cops_cuda` built for the GPUs (`make -C native ARCH=sm_80` on A100).

### Hardware and precision

Nothing about the hardware is fixed in the code. `run` detects the nodes, GPUs per node and ranks (Slurm),
the GPU memory (nvidia-smi, else the CUDA device properties) and the host memory per node (`SLURM_MEM_PER_NODE`,
else `/proc/meminfo`); `--nodes`, `--gpus-per-node`, `--gpu-memory-gb` and `--host-memory-gb` override them.
`plan` and `status` never touch a GPU, so without these options they report the GPU memory as unknown.

Record precision `--precision auto` (default; user decision 2026-10-22): FP64 records (94 B each on the device) when
they fit resident on all fit GPUs with the headroom; otherwise mixed (FP32 records of 20 B, FP64 arithmetic), with
as many records resident as fit and the rest streamed from pinned host memory, copies overlapped with computation
(P2; on several nodes every MPI rank streams its own rows on its own node). FP64 has no streaming path. The pinned
buffers (28 B per streamed record, measured on 20M), 1.5 GB per GPU worker and 500 B per facet for the rank-0 driver
must fit 85% of a node's memory, otherwise the fits are refused with the numbers. The decision is made once the record
count is known (after `base`) and recorded with the hardware in `pipeline.json`; it is checked again with the fit
job's hardware before every fit. On the local 20M fit the model predicts 3.56e9 records resident, 1.50e9 streamed
and 54.9 GB of host memory; measured were 3.54e9, 1.52e9 and about 54 GB.

Gravity precision `--gravity-precision mixed|fp64` is a separate policy (default mixed, user decision 2026-10-22):
mixed = `mixed_near_fp64_exact_far_fp32_quad3_v1` (exact FP64 near faces, FP32 quadrature far), fp64 =
`fp64_reference_v1`. The record precision's mixed is `fp32_records_fp64_arithmetic_v1` (FP32 storage, FP64
arithmetic). `plan` estimates the gravity time per precision from recorded runs (`gpu_manifest.json` under the restart
directory and `$COPS_OUTPUT_ROOT/cuda_c04`), scaled by facet pairs and GPU count from the recorded run of the nearest
size; it names that run's GPU model, and prints "unknown" when there is none.

### Wall time and segments

Default: no limit. `--wall-time HH:MM:SS` (from the start of `run`) or `--wall-time slurm` (the job's end:
`SLURM_JOB_END_TIME`, else `squeue -h -j $SLURM_JOB_ID -o %e`), `--wall-margin` (default 300 s) and
`--save-seconds` (default 600 s, the time a fit needs to write its results after stopping).

* Fits receive the deadline (`gpu_fit.py --deadline`). The driver times every residual and Jacobian evaluation and
  stops the optimizer at the best evaluated state when the next evaluation plus the reserve would pass the deadline.
  The fit is written normally (optimizer status -2), the restart file gets `status = "interrupted"` and a `segments`
  record, and `gpu_fit.py` exits 75. The next `run` continues it as segment k+1 from that restart file with the
  remaining evaluations (`--max-nfev` minus those used). **A segmented fit is not identical to one continuous run**:
  every segment is a new `least_squares` run, so the trust region and the Jacobian scaling start again.
* Other stages have no deadline option: a stage starts only if its recorded duration (an earlier attempt) fits the
  time left; otherwise `run` stops with exit 75. A killed stage runs again from its start.
* The segments (index, status, evaluations, cost) are in the fit's stamp, in `pipeline.json` (`fits`) and in every
  segment's `fit_manifest.json` (`compute_backend.wall_time.segments`).

### Slurm

```tcsh
# nodes of the fit job; GPUs per node from the partition's GRES (sinfo -h -o %G), or --gpus-per-node N
.venv/bin/python tools/cluster/cops_pipeline.py submit --restart-dir $SCRATCH/shap8_020m_flat \
    --nodes 3 --slurm-partition <p> --slurm-account <a> --slurm-module cuda/12.4 --slurm-module openmpi/4.1
#   -> slurm/prepare.sbatch, fit.sbatch, finetune.sbatch and the sbatch commands chained with afterok
# with a wall time, each job continues itself (afterany) while its group is resumable, then submits the next
.venv/bin/python tools/cluster/cops_pipeline.py submit --restart-dir $SCRATCH/shap8_020m_flat \
    --nodes 3 --wall-time slurm --resubmit 5
```

`submit` never runs `sbatch`. The prepare and fine-tuning jobs use one node (all its GPUs, `--exclusive`); the fit job
uses `--nodes` with one rank per GPU (`--ntasks-per-node` = GPUs per node, `--gpus-per-task=1`, `--gpu-bind=closest`,
`--cpus-per-task` from `sinfo` when known). Set the SITE lines (partition, account, modules) before submitting.
Time limits: `--time-prepare` (12 h), `--time-fit` (8 h), `--time-finetune` (4 h) until P0 measures the A100 rates.
