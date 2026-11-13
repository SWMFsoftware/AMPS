# Tenishev et al. (2013) Case-1 sodium demonstration

This runner executes the normal production `srcMoon` particle model for the
recommended Case 1 of Tenishev et al., *Icarus* 226 (2013), 1538--1549,
doi:10.1016/j.icarus.2013.08.021. It is a bounded meeting demonstration outside
the frozen U01--U26/I01--I46 registry.

The physical receipt uses a spherical Moon, the paper's fixed uniform sodium
reservoir of `2.65e16 atoms m^-2`, complete sticking with permanent removal of
returning particles, the `100 + 300 cos(chi)^0.25 K` temperature law, MIV at
`1.67e22 atoms s^-1` and `6000 K`, PSD and TD production kernels, no direct
sputtering, the Combi et al. radiation-pressure kernel, and a 36-hour
illuminated photoionization lifetime. The production adapter checks the
generated PSD/TD constants and source switches before running.

The paper specifies 30 September 1991 for Figs. 8--9 but gives no exact UTC.
Consequently `--simulation-start-utc` is mandatory and must be an independently
justified UTC on that date. The committed midnight value is only a parser-valid
template placeholder and is always replaced by the runner.

The run stops at 5000 particle steps and uses fixed 500-step sampling windows.
These bounds do not establish equilibration, independent samples, statistical
precision, or mesh convergence. The Flynn and Mendillo (1993) image/profile,
its uncertainties, and the exact observation time are not available as a
qualified machine-readable package, so no observational residual or validation
PASS is produced.

Run from the AMPS repository root after a clean Moon/SPICE build, or pass
`--rebuild` to request that exact clean build sequence:

```bash
python3 -B srcMoon/test/sodium_exosphere_demo/run.py \
  --simulation-start-utc 1991-09-30THH:MM:SSZ \
  --mpi-ranks 8
```

Replace the timestamp with an independently supported time on `1991-09-30`.
Outputs are written beneath
`test_output/srcMoon/sodium-exosphere-demo/<UTC-run-id>/`: the derived
`application.in`, full stdout/stderr, `result.json`, `run_manifest.json`, local
kernel-gate results, and the untouched production files under
`runtime/PT/plots/`. `execution_status=PASS` means only that the bounded run
completed and produced nonempty volumetric Na density plus the per-stage limb,
spherical-map, and anti-sunward-map products. Each of those observation files
contains production Na column density and both 5891.58/5897.56-A brightness
variables; the limb file also carries the tangent altitude in metres. The
capability remains `EXPERIMENTAL`.

The month-scale `pic.Moon.LimbIntegrals.*` and
`pic.Moon.AltitudeDistribution.*` files are deliberately not required. Their
production sampler writes only after the Moon completes an orbital-direction
cycle, whereas this resource-bounded run spans about 19 model hours. Requiring
those files would incorrectly turn successful bounded production into a runner
error; it would not add an independent physics gate.

If AMPS completed but an older runner rejected the preserved files during
post-processing, recover only the deterministic receipt generation with:

```bash
python3 -B srcMoon/test/sodium_exosphere_demo/run.py \
  --simulation-start-utc 1991-09-30THH:MM:SSZ \
  --mpi-ranks 8 \
  --output-dir <existing-run-directory> \
  --finalize-existing
```

Recovery does not rerun AMPS. It fails unless the derived input has the exact
requested UTC, stdout has the normal AMPS completion trailer, all six local
kernel results say `PASS`, and all four production product families are
nonempty. The manifest records `existing-completed-run-recovery` explicitly.
