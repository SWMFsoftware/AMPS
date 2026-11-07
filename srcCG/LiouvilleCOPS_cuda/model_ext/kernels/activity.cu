// ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/model_ext/kernels/activity.cu (on Voyager: /data001/highresolu_vt/cops_output/model_ext/kernels/activity.cu)
// SHA256 of the run copy = SHA256 of this file without these header lines: 4d53c184e697967effa6d8311577159666d0f5d623ad9024ad43bf7ee26357a4
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: CUDA source. activity_linear_clip with the exp option (af = exp(act_lin), active = 2) for
//   positive_activity.py; user decisions 2026-11-07. Compiled at run time.
// USED BY: positive_activity.py.
// ==== end of header ====
// =============================================================================================
// model_ext COPY (positive_activity.py, user decisions 2026-11-07) of
//     srcCG/LiouvilleCOPS_cuda/cops_cuda/kernels/activity.cu
//     sha256 61dc45410249f76d42c7b27fa773a08840f6ae7b759ae911c9335322de913048
//     (AMPS master 7ea4885; the file is unchanged since 47fccd5)
// Loaded instead of the original (cuda_sources.source_path patched for the module name
// "activity") only in processes where model_ext/positive_activity.py install() ran with
// COPS_POS_ENABLE=1.
//
// What changed, and why: activity_linear_clip is byte-for-byte the original. ADDED
// activity_linear_exp (same signature, same act_lin with the same summation order) for
// COPS_POS_ACTIVITY=exp (user 2026-11-07: positive activity by construction):
//     af[f] = exp(act_lin[f])  (no clip, no reference constant; CUDA exp, <= 1 ulp),
//     active[f] = 2            (pressure.cu copy: column factor dA/dc = af * basis).
// positive_activity.py calls it from CudaBackend._device_activity (activity="device").
// =============================================================================================
//
// Activity field on the device (SCALING_PLAN.md 3.5, user decision Q3: each GPU forms the
// activity of its own facets; rank 0 combines only the per-point sums).
//
// act_lin[f] = Y[f,:] a (+ A_scale * (B[f,:] b)), the expression of fit_cops Model.evaluate,
// with one fixed summation order per facet: terms in column order, left to right, the template
// sum formed first and scaled once. Compiled with --fmad=false, so every product and sum is
// rounded on its own: the result is deterministic and the same for a facet on every GPU and in
// every partition. It is not the BLAS rounding of the host matvec (a few ulp of the terms; the
// tests bound it), which is why the exact host path stays the default on one node.
//
// The clip is NumPy's: af = np.maximum(act_lin, 0.0) (NaN propagates; -0.0 gives +0.0) and
// active = act_lin > 0.0.
//
// basis: term-major, basis[j * nf + f]; the first ny terms are Y's columns, the next nb are the
// UNSCALED template columns.

extern "C" __global__ void activity_linear_clip(
    const double* __restrict__ basis, long long nf, int ny, int nb,
    const double* __restrict__ coef,      // a[0..ny) then b[0..nb)
    double scale,
    double* __restrict__ af, unsigned char* __restrict__ active,
    double* act_out, int write_act)       // act_lin itself when write_act (checks)
{
    for (long long f = (long long)blockIdx.x * blockDim.x + threadIdx.x; f < nf;
         f += (long long)gridDim.x * blockDim.x) {
        double s = basis[f] * coef[0];
        for (int j = 1; j < ny; ++j)
            s = s + basis[(long long)j * nf + f] * coef[j];
        if (nb > 0) {
            const double* B = basis + (long long)ny * nf;
            double t = B[f] * coef[ny];
            for (int k = 1; k < nb; ++k)
                t = t + B[(long long)k * nf + f] * coef[ny + k];
            s = s + scale * t;
        }
        af[f] = (s != s) ? s : (s > 0.0 ? s : 0.0);
        active[f] = s > 0.0 ? 1 : 0;
        if (write_act)
            act_out[f] = s;
    }
}

// positive_activity: exponential activity (see the header).
extern "C" __global__ void activity_linear_exp(
    const double* __restrict__ basis, long long nf, int ny, int nb,
    const double* __restrict__ coef,      // a[0..ny) then b[0..nb)
    double scale,
    double* __restrict__ af, unsigned char* __restrict__ active,
    double* act_out, int write_act)       // act_lin itself when write_act (checks)
{
    for (long long f = (long long)blockIdx.x * blockDim.x + threadIdx.x; f < nf;
         f += (long long)gridDim.x * blockDim.x) {
        double s = basis[f] * coef[0];
        for (int j = 1; j < ny; ++j)
            s = s + basis[(long long)j * nf + f] * coef[j];
        if (nb > 0) {
            const double* B = basis + (long long)ny * nf;
            double t = B[f] * coef[ny];
            for (int k = 1; k < nb; ++k)
                t = t + B[(long long)k * nf + f] * coef[ny + k];
            s = s + scale * t;
        }
        af[f] = exp(s);
        active[f] = 2;
        if (write_act)
            act_out[f] = s;
    }
}
