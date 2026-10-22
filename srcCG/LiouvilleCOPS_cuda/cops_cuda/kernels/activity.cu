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
