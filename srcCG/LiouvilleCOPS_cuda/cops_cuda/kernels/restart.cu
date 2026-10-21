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

// ---------------------------------------------------------------------------------------------
// Late-interval fine tuning (cops_cuda/fine_tuning.py; definitions of ../LiouvilleCOPS
// FINE_TUNING_STEP.md sections 2-3).
//
// first_structural_row: the earliest local row p at which facet face[r] has a structural response,
// Kn != 0 or KPhi != 0 or KV != 0 (the CPU predicate fine_tuning.row_response "structural").
// atomicMin of an integer is order independent.  ``first`` starts at 0xFFFFFFFF (never seen).
extern "C" __global__ void first_structural_row(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ kn, const double* __restrict__ kphi, const double* __restrict__ kv,
    unsigned int* __restrict__ first)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        if (kn[r] != 0.0 || kphi[r] != 0.0 || kv[r] != 0.0) atomicMin(&first[face[r]], (unsigned int)p);
    }
}

// first_interval_scores: for ONE row (the blocks [block0, block0 + gridDim.x) of local point p),
// add weight * |H_NG| and weight * |H_RG| to the facet's scores when the row's interval ``bin``
// is the facet's first interval.  H per unit clipped activity at the frozen temperatures:
//   H_NG = (cNG_H bm_H) * wn_H * u,   H_RG = cRG_H[p] * (wp_H + bm_H * wv_H * u)
// (the terms of point_change_partials).  The loader rejects duplicate facets within a point,
// so within one launch every facet is written at most once: plain read-add-write, no atomics.
// Launching the rows in chronological order reproduces the CPU's sequential np.add.at order.
// ``stamp`` re-checks the uniqueness (flag bit 0); a negative response sets flag bit 1
// (the CPU tool rejects negative gauge responses).
extern "C" __global__ void first_interval_scores(
    const long long* __restrict__ blocks, long long block0, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const int* __restrict__ first_interval, int bin, double weight, double crg,
    double cNbmH, double bmH, double Tmax, double T75, double Tmin, double KT,
    double* __restrict__ score_ng, double* __restrict__ score_rg,
    unsigned int* __restrict__ stamp, unsigned int* __restrict__ flags)
{
    const long long b = block0 + blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int f = face[r];
        if (atomicExch(&stamp[f], (unsigned int)(p + 1)) == (unsigned int)(p + 1)) atomicOr(&flags[0], 1u);
        if (first_interval[f] != bin) continue;
        const double x = xlit[r];
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = lit[r] && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        const double ng = cNbmH * wnH[r] * u;
        const double rg = crg * (wpH[r] + bmH * wvH[r] * u);
        if (ng < 0.0 || rg < 0.0) atomicOr(&flags[0], 2u);
        score_ng[f] = score_ng[f] + weight * fabs(ng);
        score_rg[f] = score_rg[f] + weight * fabs(rg);
    }
}

// first_significant_row: the earliest local row p with use_row[p] at which the change dAF_f of
// facet f alone changes the point's NG or RG pressure by more than tau (relative; the terms and
// the share of facet_point_shares).  Fine tuning with --confine total: tau = 0 over the window
// rows gives each facet's first window row with a nonzero response.  atomicMin of an integer is
// order independent.
extern "C" __global__ void first_significant_row(
    const long long* __restrict__ blocks, const unsigned char* __restrict__ use_row, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ dAF, const double* __restrict__ pNG, const double* __restrict__ pRG,
    const double* __restrict__ cRGH, double cNbmH, double bmH, double Tmax, double T75, double Tmin,
    double KT, double tau, unsigned int* __restrict__ first)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    if (!use_row[p]) return;
    const double png = pNG[p], prg = fabs(pRG[p]), crg = cRGH[p];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int f = face[r];
        const double d = dAF[f];
        if (d == 0.0) continue;
        const double x = xlit[r];
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = lit[r] && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        const double dng = cNbmH * wnH[r] * u * d;
        const double drg = crg * (wpH[r] + bmH * wvH[r] * u) * d;
        const double sng = (png > 0.0) ? fabs(dng) / png : (dng != 0.0 ? 1.0 / 0.0 : 0.0);
        const double srg = (prg > 0.0) ? fabs(drg) / prg : (drg != 0.0 ? 1.0 / 0.0 : 0.0);
        if (sng > tau || srg > tau) atomicMin(&first[f], (unsigned int)p);
    }
}

// ---------------------------------------------------------------------------------------------
// FP32 kernel records (C14 mixed-precision experiment; device_store record_precision "fp32").
// No stored per-record factors: the H2O factors are recomputed in FP64 from the FP32 records
// exactly as pressure.cu static_factors_m does (H2O part; shadow folded into cos <= 0 -> night),
// so these kernels see the factors of the mixed evaluator.  Per-point constants
// pt[0..4][p] = qmin_H, qmax_H, qmin_C, qmax_C, VdotB (component-major, evaluator array "pt").
// Each _m kernel is the FP64 kernel above with (wnH, wvH, wpH, xlit, lit) replaced by RS_FACTORS_M.

__device__ __forceinline__ void factors_h_m(
    long long r, long long p,
    const float* __restrict__ kn32, const float* __restrict__ kphi32,
    const float* __restrict__ kv32, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt, double cos75,
    double& wn_h, double& wv_h, double& wp_h, double& x, bool& is_lit)
{
    const double qmin_h = pt[p], qmax_h = pt[npt + p], vb = pt[4 * npt + p];
    const double c0 = (double)cos32[r];
    const bool night = c0 <= 0.0;
    const double cz = c0 < 1.0 ? c0 : 1.0;
    const double qh = night ? qmin_h : qmin_h + (qmax_h - qmin_h) * cz;
    wn_h = qh * (double)kn32[r];
    wv_h = qh * vb * (double)kv32[r];
    wp_h = qh * (double)kphi32[r];
    const double clipped = c0 < -1.0 ? -1.0 : (c0 > 1.0 ? 1.0 : c0);
    const double c = night ? -1.0 : clipped;
    is_lit = c >= cos75;
    x = is_lit ? 1.0 / c - 1.0 : 0.0;
}

#define RS_FACTORS_M_ARGS \
    const float* __restrict__ kn32, const float* __restrict__ kphi32, \
    const float* __restrict__ kv32, const float* __restrict__ cos32, \
    const double* __restrict__ pt, long long npt, double cos75
#define RS_FACTORS_M(r, p) \
    double wn, wv, wp, x; bool islit; \
    factors_h_m(r, p, kn32, kphi32, kv32, cos32, pt, npt, cos75, wn, wv, wp, x, islit)

extern "C" __global__ void facet_point_shares_m(
    const long long* __restrict__ blocks, const int* __restrict__ face, RS_FACTORS_M_ARGS,
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
        RS_FACTORS_M(r, p);
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = islit && (tl > Tmin);
        const double te = use_lit ? tl : Tmin;
        const double u = 1.0 / sqrt(te);
        const double dng = cNbmH * wn * u * d;
        const double drg = crg * (wp + bmH * wv * u) * d;
        const double sng = (png > 0.0) ? fabs(dng) / png : (dng != 0.0 ? 1.0 / 0.0 : 0.0);
        const double srg = (prg > 0.0) ? fabs(drg) / prg : (drg != 0.0 ? 1.0 / 0.0 : 0.0);
        atomicMax(&maxNG[f], (unsigned long long)__double_as_longlong(sng));
        atomicMax(&maxRG[f], (unsigned long long)__double_as_longlong(srg));
    }
}

extern "C" __global__ void __launch_bounds__(RS_BLOCK) point_change_partials_m(
    const long long* __restrict__ blocks, const int* __restrict__ face, RS_FACTORS_M_ARGS,
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
        RS_FACTORS_M(r, p);
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = islit && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        an = an + cNbmH * wn * u * d;
        ar = ar + crg * (wp + bmH * wv * u) * d;
    }
    sn[threadIdx.x] = an; sr[threadIdx.x] = ar;
    __syncthreads();
    for (int w = RS_BLOCK / 2; w > 0; w >>= 1) {
        if (threadIdx.x < w) { sn[threadIdx.x] += sn[threadIdx.x + w]; sr[threadIdx.x] += sr[threadIdx.x + w]; }
        __syncthreads();
    }
    if (threadIdx.x == 0) { partial[2 * b] = sn[0]; partial[2 * b + 1] = sr[0]; }
}

// The loader keeps zero-ness in FP32, so the structural predicate is unchanged.
extern "C" __global__ void first_structural_row_f32(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const float* __restrict__ kn, const float* __restrict__ kphi, const float* __restrict__ kv,
    unsigned int* __restrict__ first)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        if (kn[r] != 0.0f || kphi[r] != 0.0f || kv[r] != 0.0f) atomicMin(&first[face[r]], (unsigned int)p);
    }
}

extern "C" __global__ void first_interval_scores_m(
    const long long* __restrict__ blocks, long long block0, const int* __restrict__ face, RS_FACTORS_M_ARGS,
    const int* __restrict__ first_interval, int bin, double weight, double crg,
    double cNbmH, double bmH, double Tmax, double T75, double Tmin, double KT,
    double* __restrict__ score_ng, double* __restrict__ score_rg,
    unsigned int* __restrict__ stamp, unsigned int* __restrict__ flags)
{
    const long long b = block0 + blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int f = face[r];
        if (atomicExch(&stamp[f], (unsigned int)(p + 1)) == (unsigned int)(p + 1)) atomicOr(&flags[0], 1u);
        if (first_interval[f] != bin) continue;
        RS_FACTORS_M(r, p);
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = islit && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        const double ng = cNbmH * wn * u;
        const double rg = crg * (wp + bmH * wv * u);
        if (ng < 0.0 || rg < 0.0) atomicOr(&flags[0], 2u);
        score_ng[f] = score_ng[f] + weight * fabs(ng);
        score_rg[f] = score_rg[f] + weight * fabs(rg);
    }
}

extern "C" __global__ void first_significant_row_m(
    const long long* __restrict__ blocks, const unsigned char* __restrict__ use_row, const int* __restrict__ face,
    RS_FACTORS_M_ARGS,
    const double* __restrict__ dAF, const double* __restrict__ pNG, const double* __restrict__ pRG,
    const double* __restrict__ cRGH, double cNbmH, double bmH, double Tmax, double T75, double Tmin,
    double KT, double tau, unsigned int* __restrict__ first)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    if (!use_row[p]) return;
    const double png = pNG[p], prg = fabs(pRG[p]), crg = cRGH[p];
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int f = face[r];
        const double d = dAF[f];
        if (d == 0.0) continue;
        RS_FACTORS_M(r, p);
        const double tl = Tmax + KT * (Tmax - T75) * x;
        const bool use_lit = islit && (tl > Tmin);
        const double u = 1.0 / sqrt(use_lit ? tl : Tmin);
        const double dng = cNbmH * wn * u * d;
        const double drg = crg * (wp + bmH * wv * u) * d;
        const double sng = (png > 0.0) ? fabs(dng) / png : (dng != 0.0 ? 1.0 / 0.0 : 0.0);
        const double srg = (prg > 0.0) ? fabs(drg) / prg : (drg != 0.0 ? 1.0 / 0.0 : 0.0);
        if (sng > tau || srg > tau) atomicMin(&first[f], (unsigned int)p);
    }
}
