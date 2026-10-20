// Per-facet influence of an activity change on the pressures of selected points (ROADMAP C12).
//
// For a change dAF_f of the clipped H2O activity of facet f, the pressure changes of point p are
// (fit_cops.Model.evaluate, derivative of the H2O terms with respect to af):
//   dpNG = (cNG_H bm_H) * wn_H * u * dAF_f
//   dpRG = cRG_H[p] * (wp_H + bm_H * wv_H * u) * dAF_f
// with the per-record static factors of the evaluator (wn_H, wv_H, wp_H, xlit, lit; pressure.cu
// precompute_static) and u = 1/sqrt(Te) of the state's temperatures (record_factors).
// For every record of a flagged point the relative changes |dpNG|/pNG[p] and |dpRG|/|pRG[p]| are
// max-reduced per facet.  The maxima are non-negative doubles, so atomicMax on their IEEE bit
// patterns (as unsigned 64-bit integers) is exact and order independent.

extern "C" __global__ void facet_point_shares(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ dAF, const unsigned char* __restrict__ flagged,
    const double* __restrict__ pNG, const double* __restrict__ pRG, const double* __restrict__ cRGH,
    double cNbmH, double bmH, double Tmax, double T75, double Tmin, double KT,
    unsigned long long* __restrict__ maxNG, unsigned long long* __restrict__ maxRG)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    if (!flagged[p]) return;
    const double png = pNG[p], prg = fabs(pRG[p]), crg = cRGH[p];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int f = face[r];
        const double d = dAF[f];
        if (d == 0.0) continue;
        const double x = xlit[r];
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = lit[r] && (tl > Tmin);
        const double te = use_lit ? tl : Tmin;
        const double u = 1.0 / sqrt(te);
        const double dng = cNbmH * wnH[r] * u * d;
        const double drg = crg * (wpH[r] + bmH * wvH[r] * u) * d;
        const double sng = (png > 0.0) ? fabs(dng) / png : (dng != 0.0 ? 1.0 / 0.0 : 0.0);
        const double srg = (prg > 0.0) ? fabs(drg) / prg : (drg != 0.0 ? 1.0 / 0.0 : 0.0);
        atomicMax(&maxNG[f], (unsigned long long)__double_as_longlong(sng));
        atomicMax(&maxRG[f], (unsigned long long)__double_as_longlong(srg));
    }
}

// Per-block partial sums of dpNG and dpRG (same terms as above) for the records of each block;
// the host adds the block partials of a point in block order (deterministic).  The pressures are
// linear in each facet's clipped activity, so the sum is the exact first-order AND total change.
#define RS_BLOCK 256
extern "C" __global__ void __launch_bounds__(RS_BLOCK) point_change_partials(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ dAF, const double* __restrict__ cRGH,
    double cNbmH, double bmH, double Tmax, double T75, double Tmin, double KT,
    double* __restrict__ partial)
{
    __shared__ double sn[RS_BLOCK], sr[RS_BLOCK];
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    const double crg = cRGH[p];
    double an = 0.0, ar = 0.0;
    for (long long r = begin + threadIdx.x; r < end; r += RS_BLOCK) {
        const double d = dAF[face[r]];
        if (d == 0.0) continue;
        const double x = xlit[r];
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = lit[r] && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        an = an + cNbmH * wnH[r] * u * d;
        ar = ar + crg * (wpH[r] + bmH * wvH[r] * u) * d;
    }
    sn[threadIdx.x] = an; sr[threadIdx.x] = ar;
    __syncthreads();
    for (int w = RS_BLOCK / 2; w > 0; w >>= 1) {
        if (threadIdx.x < w) { sn[threadIdx.x] += sn[threadIdx.x + w]; sr[threadIdx.x] += sr[threadIdx.x + w]; }
        __syncthreads();
    }
    if (threadIdx.x == 0) { partial[2 * b] = sn[0]; partial[2 * b + 1] = sr[0]; }
}
