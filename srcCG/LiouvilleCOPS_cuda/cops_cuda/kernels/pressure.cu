// FP64 COPS pressure and analytic-Jacobian sums (ROADMAP C06).
//
// Reference: ../LiouvilleCOPS/tools/fit_cops.py, Model.load_kernels / prepare /
// evaluate.  The host (cops_cuda/evaluator.py) assembles pressures and Jacobian
// columns from the per-point SUMS produced here with the reference's own
// NumPy expressions; these kernels compute only the per-contribution terms and
// their per-point sums.
//
// Exactness contract
// ------------------
// * Every per-contribution expression is written in the evaluation order of
//   the NumPy reference (left to right, e.g. q * VB * KV == (q*VB)*KV) and the
//   module is compiled with --fmad=false, so each term is the same correctly
//   rounded IEEE result as on the CPU (sqrt and division are IEEE on both).
// * The per-face activity field act_lin = Y a + A_scale B c, its clip
//   max(act_lin, 0) and the active set act_lin > 0 are computed on the host
//   with the reference's NumPy expression and uploaded, so clipping decisions
//   are bitwise those of the CPU (a BLAS matvec order cannot be reproduced on
//   the device; a per-face upload of 2 x 8 bytes is negligible).
// * eval_sequential: one thread per point folds its records in file order,
//   which is the order of np.bincount and of SciPy's CSR x dense product.
//   It therefore reproduces the CPU sums bitwise and serves as the
//   verification mode (slow: one thread per point).
// * eval_partial + eval_finish: the production mode.  Records of a point are
//   split into fixed blocks (host block table); inside a block thread t sums
//   records t, t+B, ...; a fixed shared-memory tree combines the threads; the
//   block partials are added in block order.  Deterministic for a pinned
//   configuration and within the approved budget (<= 64 eps M vs the exact
//   sum, <= 512 eps M vs the CPU).
//
// Sum codes (one per requested sum; the host decodes them):
//   0..4   S1 = sum wn_H*af*u, S2 = sum wp_H*af, S3 = sum wv_H*af*u,
//          S4 = sum wn_C*u,     S5 = sum wv_C*u
//   5..13  temperature sums: weight w in (gN, gRH, gRC) x kind in (dl, dlx, dn):
//          code = 5 + 3*w + kind
//   14     PhiCO2 = sum qC*KPhi (parameter independent, computed once)
//   1000+c NG column c: sum (wn_H*u) * mb_c(face)
//   2000+c RG column c: sum (wp_H + bm_H*wv_H*u) * mb_c(face)
//   where mb_c(f) = active(f) ? colbasis[c][f] : 0.0*colbasis[c][f]
//   (numpy's bool * float keeps the sign of zero).

#define MAX_SUMS 32
#define BLOCK 256

struct Params {
    double Tmax, T75, Tmin, KT;      // temperature law
    double cNbmH;                    // cNG_H * bm_H
    double kcNbmC;                   // (kappa * cNG_C) * bm_C
    double bmH, bmC;                 // bm_H, bm_C
};

// The host passes the eight parameters as a device array in this order (a
// struct argument is not portable through CuPy); values are copied unchanged.
__device__ __forceinline__ Params load_params(const double* __restrict__ pp)
{
    Params P;
    P.Tmax = pp[0]; P.T75 = pp[1]; P.Tmin = pp[2]; P.KT = pp[3];
    P.cNbmH = pp[4]; P.kcNbmC = pp[5]; P.bmH = pp[6]; P.bmC = pp[7];
    return P;
}

// Static per-contribution factors (Model.load_kernels + prepare), once per bind.
// One block row (blockIdx.y) per point; records of the point in grid-stride.
extern "C" __global__ void precompute_static(
    const long long* __restrict__ offsets, const int* __restrict__ face,
    const double* __restrict__ kn, const double* __restrict__ kphi,
    const double* __restrict__ kv, const double* __restrict__ cosine,
    const unsigned char* __restrict__ shadow,
    const double* __restrict__ qminH, const double* __restrict__ qmaxH,
    const double* __restrict__ qminC, const double* __restrict__ qmaxC,
    const double* __restrict__ vdotb, const double* __restrict__ aco2,
    double cos75,
    double* __restrict__ wnH, double* __restrict__ wvH, double* __restrict__ wpH,
    double* __restrict__ wnC, double* __restrict__ wvC, double* __restrict__ phiC,
    double* __restrict__ xlit, unsigned char* __restrict__ lit)
{
    const long long p = blockIdx.y;
    const long long begin = offsets[p], end = offsets[p + 1];
    const double qmin_h = qminH[p], qmax_h = qmaxH[p], qmin_c = qminC[p], qmax_c = qmaxC[p];
    const double vb = vdotb[p];
    for (long long r = begin + (long long)blockIdx.x * blockDim.x + threadIdx.x; r < end;
         r += (long long)gridDim.x * blockDim.x) {
        const double c0 = cosine[r];
        // night = (cos <= 0) | (shadow != 0)
        const bool night = (c0 <= 0.0) || (shadow[r] != 0);
        // SZAMode continuous: cz = np.minimum(cos, 1.0)
        const double cz = c0 < 1.0 ? c0 : 1.0;
        // Q = where(night, qmin, qmin + (qmax - qmin) * cz)
        const double qh = night ? qmin_h : qmin_h + (qmax_h - qmin_h) * cz;
        const double qc0 = night ? qmin_c : qmin_c + (qmax_c - qmin_c) * cz;
        // q["CO2"] = Q["CO2"] * ACO2[face]
        const double qc = qc0 * aco2[face[r]];
        wnH[r] = qh * kn[r];
        wvH[r] = qh * vb * kv[r];          // (q * VB) * KV
        wpH[r] = qh * kphi[r];
        wnC[r] = qc * kn[r];
        wvC[r] = qc * vb * kv[r];
        phiC[r] = qc * kphi[r];
        // prepare(): c = where(night, -1, clip(cos, -1, 1)); lit = c >= COS75;
        // xlit = where(lit, 1/where(lit, c, 1) - 1, 0)
        const double clipped = c0 < -1.0 ? -1.0 : (c0 > 1.0 ? 1.0 : c0);
        const double c = night ? -1.0 : clipped;
        const bool is_lit = c >= cos75;
        lit[r] = is_lit ? 1 : 0;
        xlit[r] = is_lit ? 1.0 / c - 1.0 : 0.0;
    }
}

// Value of one sum code for record r (see the header).  All factors that are
// shared between codes are passed in, computed once per record.
__device__ __forceinline__ double term(int code, long long r, int f,
    const double af, const double u, const double dl, const double dlx, const double dn,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ wnC, const double* __restrict__ wvC, const double* __restrict__ phiC,
    const unsigned char* __restrict__ active, const double* __restrict__ colbasis, long long nfaces,
    const Params& P)
{
    if (code < 5) {
        switch (code) {
        case 0: return wnH[r] * af * u;
        case 1: return wpH[r] * af;
        case 2: return wvH[r] * af * u;
        case 3: return wnC[r] * u;
        default: return wvC[r] * u;
        }
    }
    if (code < 14) {
        const int w = (code - 5) / 3, kind = (code - 5) % 3;
        double weight;
        if (w == 0) {
            // gN = cNG_H*bm_H*wn_H*af + kappa*cNG_C*bm_C*wn_C
            weight = P.cNbmH * wnH[r] * af + P.kcNbmC * wnC[r];
        } else if (w == 1) {
            weight = P.bmH * wvH[r] * af;            // (bm_H * wv_H) * af
        } else {
            weight = P.bmC * wvC[r];
        }
        return weight * (kind == 0 ? dl : (kind == 1 ? dlx : dn));
    }
    if (code == 14) return phiC[r];
    const int column = code >= 2000 ? code - 2000 : code - 1000;
    const double b = colbasis[(long long)column * nfaces + f];
    const double mb = active[f] ? b : 0.0 * b;
    if (code < 2000) return wnH[r] * u * mb;
    return (wpH[r] + P.bmH * wvH[r] * u) * mb;
}

// Per-record parameter-dependent factors (Model.evaluate).
__device__ __forceinline__ void record_factors(long long r, int f,
    const double* __restrict__ af_face, const double* __restrict__ xlit,
    const unsigned char* __restrict__ lit, const Params& P,
    double& af, double& u, double& dl, double& dlx, double& dn)
{
    af = af_face[f];                                  // np.maximum(act_lin, 0)[face]
    const double x = xlit[r];
    // Tl = Tmax + KT * (Tmax - T75) * xlit
    const double tl = P.Tmax + P.KT * (P.Tmax - P.T75) * x;
    const bool use_lit = lit[r] && (tl > P.Tmin);
    const double te = use_lit ? tl : P.Tmin;
    u = 1.0 / sqrt(te);
    const double du = -0.5 * u / te;                  // (-0.5 * u) / Te
    dl = use_lit ? du : 0.0;
    dn = use_lit ? 0.0 : du;
    dlx = dl * x;
}

// Verification mode: one thread per point, records in file order (bitwise CPU).
extern "C" __global__ void eval_sequential(
    const long long* __restrict__ offsets, long long npoints, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ wnC, const double* __restrict__ wvC, const double* __restrict__ phiC,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,
    const double* __restrict__ colbasis, long long nfaces,
    const int* __restrict__ codes, int nsums, const double* __restrict__ pp,
    double* __restrict__ out)
{
    const long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= npoints) return;
    const Params P = load_params(pp);
    double acc[MAX_SUMS];
    for (int k = 0; k < nsums; ++k) acc[k] = 0.0;      // np.bincount starts from 0
    for (long long r = offsets[p]; r < offsets[p + 1]; ++r) {
        const int f = face[r];
        double af, u, dl, dlx, dn;
        record_factors(r, f, af_face, xlit, lit, P, af, u, dl, dlx, dn);
        for (int k = 0; k < nsums; ++k)
            acc[k] = acc[k] + term(codes[k], r, f, af, u, dl, dlx, dn, wnH, wvH, wpH, wnC, wvC,
                                   phiC, active, colbasis, nfaces, P);
    }
    for (int k = 0; k < nsums; ++k) out[(long long)k * npoints + p] = acc[k];
}

// Production mode, level 1.  blocks[3*b + 0..2] = (point, first record, end).
extern "C" __global__ void __launch_bounds__(BLOCK) eval_partial(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ wnC, const double* __restrict__ wvC, const double* __restrict__ phiC,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,
    const double* __restrict__ colbasis, long long nfaces,
    const int* __restrict__ codes, int nsums, const double* __restrict__ pp,
    double* __restrict__ partial)
{
    __shared__ double tree[BLOCK];
    const Params P = load_params(pp);
    const long long b = blockIdx.x;
    const long long begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    double acc[MAX_SUMS];
    for (int k = 0; k < nsums; ++k) acc[k] = 0.0;
    for (long long r = begin + threadIdx.x; r < end; r += BLOCK) {
        const int f = face[r];
        double af, u, dl, dlx, dn;
        record_factors(r, f, af_face, xlit, lit, P, af, u, dl, dlx, dn);
        for (int k = 0; k < nsums; ++k)
            acc[k] = acc[k] + term(codes[k], r, f, af, u, dl, dlx, dn, wnH, wvH, wpH, wnC, wvC,
                                   phiC, active, colbasis, nfaces, P);
    }
    // Fixed-shape pairwise tree over the block, one sum at a time.
    for (int k = 0; k < nsums; ++k) {
        tree[threadIdx.x] = acc[k];
        __syncthreads();
        for (int width = BLOCK / 2; width > 0; width >>= 1) {
            if (threadIdx.x < width) tree[threadIdx.x] = tree[threadIdx.x] + tree[threadIdx.x + width];
            __syncthreads();
        }
        if (threadIdx.x == 0) partial[b * nsums + k] = tree[0];
        __syncthreads();
    }
}

// Production mode, level 2: per point, add its block partials in block order.
extern "C" __global__ void eval_finish(
    const long long* __restrict__ first_block, long long npoints,
    const double* __restrict__ partial, int nsums, double* __restrict__ out)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;   // i = k*npoints + p
    if (i >= (long long)nsums * npoints) return;
    const long long k = i / npoints, p = i % npoints;
    double s = 0.0;
    for (long long b = first_block[p]; b < first_block[p + 1]; ++b) s = s + partial[b * nsums + k];
    out[i] = s;
}

// ---------------------------------------------------------------------------
// C08: FP64 restructure (fused, compile-time specialized; production tree mode)
// ---------------------------------------------------------------------------
// eval_partial keeps its accumulators in a runtime-indexed array, which ptxas
// places in local memory (256-byte stack frame), and dispatches every term on a
// runtime code.  The kernels below compute the SAME sums for the fixed code sets
// the evaluator requests, with a fixed number N of accumulators known at compile
// time (fully unrolled, held in registers):
//   NP = 5   : pressure codes 0..4
//   NP = 14  : pressure codes 0..4, then temperature codes 5..13
//   G  > 0   : activity/template columns col0 .. col0+G-1: first the G NG sums
//              (codes 1000+c), then the G RG sums (codes 2000+c)
// Exactness: every term is written exactly as in term() (same operands, same
// left-to-right order, --fmad=false); each accumulator folds the thread's
// records r = begin+t, begin+t+BLOCK, ... in the same order; the shared-memory
// tree and the partial layout partial[b*N + k] are those of eval_partial, and
// eval_finish adds the block partials in the same order.  The per-sum results
// are therefore bitwise identical to eval_partial for the same block table.
template <int NP, int G>
__device__ __forceinline__ void fused_body(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH,
    const double* __restrict__ wnC, const double* __restrict__ wvC,
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,
    const double* __restrict__ colbasis, long long nfaces, int col0,
    const double* __restrict__ pp, double* __restrict__ partial)
{
    constexpr int N = NP + 2 * G;
    __shared__ double tree[BLOCK];
    const Params P = load_params(pp);
    const long long b = blockIdx.x;
    const long long begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    double acc[N];
#pragma unroll
    for (int k = 0; k < N; ++k) acc[k] = 0.0;
    for (long long r = begin + threadIdx.x; r < end; r += BLOCK) {
        const int f = face[r];
        double af, u, dl, dlx, dn;
        record_factors(r, f, af_face, xlit, lit, P, af, u, dl, dlx, dn);
        if (NP >= 5) {
            const double wn_h = wnH[r], wv_h = wvH[r], wp_h = wpH[r], wn_c = wnC[r], wv_c = wvC[r];
            acc[0] = acc[0] + wn_h * af * u;                       // code 0
            acc[1] = acc[1] + wp_h * af;                           // code 1
            acc[2] = acc[2] + wv_h * af * u;                       // code 2
            acc[3] = acc[3] + wn_c * u;                            // code 3
            acc[4] = acc[4] + wv_c * u;                            // code 4
            if (NP >= 14) {
                const double w0 = P.cNbmH * wn_h * af + P.kcNbmC * wn_c;   // gN
                const double w1 = P.bmH * wv_h * af;                     // (bm_H * wv_H) * af
                const double w2 = P.bmC * wv_c;
                acc[5] = acc[5] + w0 * dl;  acc[6] = acc[6] + w0 * dlx;  acc[7] = acc[7] + w0 * dn;
                acc[8] = acc[8] + w1 * dl;  acc[9] = acc[9] + w1 * dlx;  acc[10] = acc[10] + w1 * dn;
                acc[11] = acc[11] + w2 * dl; acc[12] = acc[12] + w2 * dlx; acc[13] = acc[13] + w2 * dn;
            }
        }
        if (G > 0) {
            const bool on = active[f] != 0;
            const double wn_u = wnH[r] * u;                          // (wn_H * u), then * mb
            const double rg = wpH[r] + P.bmH * wvH[r] * u;           // wp_H + (bm_H * wv_H) * u
#pragma unroll
            for (int g = 0; g < G; ++g) {
                const double bv = colbasis[(long long)(col0 + g) * nfaces + f];
                const double mb = on ? bv : 0.0 * bv;
                acc[NP + g] = acc[NP + g] + wn_u * mb;
                acc[NP + G + g] = acc[NP + G + g] + rg * mb;
            }
        }
    }
#pragma unroll
    for (int k = 0; k < N; ++k) {
        tree[threadIdx.x] = acc[k];
        __syncthreads();
        for (int width = BLOCK / 2; width > 0; width >>= 1) {
            if (threadIdx.x < width) tree[threadIdx.x] = tree[threadIdx.x] + tree[threadIdx.x + width];
            __syncthreads();
        }
        if (threadIdx.x == 0) partial[b * N + k] = tree[0];
        __syncthreads();
    }
}

#define COPS_FUSED_KERNEL(NAME, NP, G)                                                         \
extern "C" __global__ void __launch_bounds__(BLOCK) NAME(                                       \
    const long long* __restrict__ blocks, const int* __restrict__ face,                       \
    const double* __restrict__ wnH, const double* __restrict__ wvH, const double* __restrict__ wpH, \
    const double* __restrict__ wnC, const double* __restrict__ wvC,                           \
    const double* __restrict__ xlit, const unsigned char* __restrict__ lit,                   \
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,             \
    const double* __restrict__ colbasis, long long nfaces, int col0,                          \
    const double* __restrict__ pp, double* __restrict__ partial)                              \
{                                                                                              \
    fused_body<NP, G>(blocks, face, wnH, wvH, wpH, wnC, wvC, xlit, lit, af_face, active,       \
                      colbasis, nfaces, col0, pp, partial);                                    \
}

COPS_FUSED_KERNEL(eval_fused_p5, 5, 0)
COPS_FUSED_KERNEL(eval_fused_p14, 14, 0)
COPS_FUSED_KERNEL(eval_fused_c1, 0, 1)
COPS_FUSED_KERNEL(eval_fused_c2, 0, 2)
COPS_FUSED_KERNEL(eval_fused_c3, 0, 3)
COPS_FUSED_KERNEL(eval_fused_c4, 0, 4)
COPS_FUSED_KERNEL(eval_fused_c5, 0, 5)
COPS_FUSED_KERNEL(eval_fused_c6, 0, 6)
COPS_FUSED_KERNEL(eval_fused_c7, 0, 7)
COPS_FUSED_KERNEL(eval_fused_c8, 0, 8)
COPS_FUSED_KERNEL(eval_fused_c9, 0, 9)
COPS_FUSED_KERNEL(eval_fused_c10, 0, 10)
COPS_FUSED_KERNEL(eval_fused_c11, 0, 11)
COPS_FUSED_KERNEL(eval_fused_c12, 0, 12)
COPS_FUSED_KERNEL(eval_fused_c13, 0, 13)
COPS_FUSED_KERNEL(eval_fused_c14, 0, 14)
COPS_FUSED_KERNEL(eval_fused_c15, 0, 15)
COPS_FUSED_KERNEL(eval_fused_c16, 0, 16)
