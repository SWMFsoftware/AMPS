// ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/model_ext/kernels/pressure.cu (on Voyager: /data001/highresolu_vt/cops_output/model_ext/kernels/pressure.cu)
// SHA256 of the run copy = SHA256 of this file without these header lines: 85842e95149817a73664f6e1afe543fd2dff5d884063124398f5887bdb08faf7
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: CUDA source. The mixed-precision static_factors_m / fused_body_m kernels of the GPU forward model with the
//   positive-activity options of positive_activity.py (night/shadow ratio, exp activity, H2O floor block indexed
//   aco2[nfaces + f]); user decisions 2026-11-07. Compiled at run time by positive_activity.py (NVRTC via CuPy).
// USED BY: positive_activity.py.
// 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
// ==== end of header ====
// =============================================================================================
// model_ext COPY (positive_activity.py, user decisions 2026-11-07) of
//     srcCG/LiouvilleCOPS_cuda/cops_cuda/kernels/pressure.cu
//     sha256 dc4b6f4ceeba96bcd2df335c6dff607d6e545c2c363b8bf897836652010798f2
//     (AMPS master 7ea4885; the file is unchanged since 47fccd5)
// Loaded instead of the original (cuda_sources.source_path patched for the module name
// "pressure") only in processes where model_ext/positive_activity.py install() ran with
// COPS_POS_ENABLE=1. Kernel names and kernel signatures are those of the original.
//
// What changed, and why (only the MIXED PRECISION section, the path of the --precision mixed
// fits; the FP64 kernels below are byte-for-byte the original, and positive_activity.py refuses
// FP64 records whenever a modification is active):
//   1. static_factors_m returns qh, the record's H2O production per unit activity
//      (q_H2O = qmin_H + (qmax_H - qmin_H) cz, or qmin_H at night), which it already computed.
//   2. fused_body_m: H2O flux floor ("floor form (a)"): with F = q_H2O * af the record's H2O
//      flux, F_min = qmax_H(p) * aco2[nfaces + f] (the second block of the per-face aco2 array is
//      floor_ratio * A_orig(f), A_orig = elements.bin column 32 "AH2O" of the base, uploaded by
//      positive_activity.py); when F < F_min the record uses F_min: wn_h, wv_h, wp_h become
//      F_min*Kn, (F_min*VB)*KV, F_min*KPhi and its activity factor becomes 1.0 (exact), so every
//      pressure and temperature sum uses F_min in place of q_H2O*af with the record's own u.
//      Records at the floor contribute 0 to the activity/template columns. CO2 has no floor.
//      With floor_ratio = 0 the condition F < 0 is never true (q > 0, af >= 0; NaN compares
//      false), so the original expressions run unchanged: bitwise the original.
//   3. fused_body_m: active[f] == 2 marks the exponential activity (af = exp(act_lin), written by
//      positive_activity.py / activity.cu activity_linear_exp): the column factor becomes
//      mb = af * basis (dA/dc = A * basis) instead of basis. active[f] in {0, 1} (clip) keeps
//      mb = active ? basis : 0 * basis, bitwise the original.
//   4. Layout: aco2 holds 2 * nfaces values: [0, nfaces) the run's CO2 map ACO2 (unchanged, read
//      by every kernel as before), [nfaces, 2 nfaces) floor_ratio * AH2O (read only by
//      fused_body_m). Kernels that index aco2[f] only (precompute_static, eval_fused_m_phic)
//      are unaffected.
//   5. (2026-11-08, user: "use the same floor for CO2"; positive_activity.py COPS_POS_FLOOR_CO2=1)
//      CO2 flux floor in the NEW kernels eval_fused_m_p5f / eval_fused_m_p14f only
//      (fused_body_m<NP, G, CF = true>): the record's CO2 flux is
//          F_C = max(kappa * qc, F_min,C),  qc = Q_CO2 * ACO2(f)  (the reference's q["CO2"], the
//          production per unit kappa, Q_CO2 with the night/shadow ratio),
//          F_min,C = qmax_C(p) * (r * ACO2(f))  (r = floor_ratio, the H2O floor's ratio; ACO2 = the
//          run's CO2 map = the CO2 activity at kappa = 1),
//      i.e. CO2 emission never drops below r times the run's maximum-production CO2 emission, even
//      for kappa -> 0. A record is at the CO2 floor when kappa * qc < F_min,C (strict; NaN false).
//      Records NOT at the floor enter the CO2 sums of codes 3, 4 and the temperature weights
//      exactly as before (the host multiplies them by kappa); their PhiCO2 is accumulated per
//      evaluation (code 15; the bind-time PhiCO2 of all records is not used). Records AT the floor
//      are left out of codes 3, 4, 15 and of the CO2 terms of the temperature weights (gN keeps its
//      H2O term, w2 none) and go into kappa-independent sums with F_min,C in place of q:
//          16  sum F_min,C Kn u             (NG)
//          17  sum (F_min,C VB) KV u        (RG, VdotB part)
//          18  sum F_min,C KPhi             (RG, Phi part)
//          19..21  sum F_min,C Kn x (dl, dlx, dn)                 (p14f only: NG temperature)
//          22..24  sum (bm_C (F_min,C VB) KV) x (dl, dlx, dn)     (p14f only: RG temperature)
//      positive_activity.py assembles p = p_H + kappa p_C,active + p_C,floor and the kappa column
//      from the active sums only: records at the CO2 floor contribute 0 to d p / d kappa. The
//      H2O part (static_factors_m, the H2O floor, the activity/template columns) is unchanged
//      (finetune.cu copies it). Layout with COPS_POS_FLOOR_CO2=1: aco2 holds 2 nfaces + 2 values,
//      aco2[2 nfaces] = kappa (written by positive_activity.py before every evaluation) and
//      aco2[2 nfaces + 1] = r; the CF kernels read them. The existing kernels (CF = false) do not
//      read them and are the 2026-11-07 code: the CF additions are compile-time dead there.
// =============================================================================================
//
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

// ---------------------------------------------------------------------------
// MIXED PRECISION (experiment, 2026-10-21; precision.py "fp32" record storage)
// ---------------------------------------------------------------------------
// Kn, KPhi, KV and cos are stored in FP32, and the shadow flag is folded into cos (cos = -2
// for shadowed records, device_store.SHADOW_FOLDED_COS; a night record's cos is never read):
// 20 bytes per record instead of 37 plus the 57 bytes of static factors, so meshes of more
// than 10M facets fit the devices. Each FP32
// value is widened to FP64, and the static factors of precompute_static are recomputed per
// record with exactly its expressions; every term and sum below is the FP64 code of
// fused_body. The only difference from the FP64 path is the FP32 rounding of the four
// inputs (relative 2^-24), so the results are not bitwise those of the reference.
// Per-point constants: pt[0..4][p] = qmin_H, qmax_H, qmin_C, qmax_C, VdotB (component-major).

__device__ __forceinline__ void static_factors_m(
    long long r, long long p, int f,
    const float* __restrict__ kn32, const float* __restrict__ kphi32,
    const float* __restrict__ kv32, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt,
    const double* __restrict__ aco2, double cos75,
    double& wn_h, double& wv_h, double& wp_h, double& wn_c, double& wv_c, double& phi_c,
    double& x, bool& is_lit, double& qh_out)      // positive_activity: + qh_out
{
    const double qmin_h = pt[p], qmax_h = pt[npt + p], qmin_c = pt[2 * npt + p],
                 qmax_c = pt[3 * npt + p], vb = pt[4 * npt + p];
    const double c0 = (double)cos32[r];
    const double kn = (double)kn32[r], kphi = (double)kphi32[r], kv = (double)kv32[r];
    const bool night = c0 <= 0.0;                  // shadow folded into cos (-2)
    const double cz = c0 < 1.0 ? c0 : 1.0;
    const double qh = night ? qmin_h : qmin_h + (qmax_h - qmin_h) * cz;
    const double qc0 = night ? qmin_c : qmin_c + (qmax_c - qmin_c) * cz;
    const double qc = qc0 * aco2[f];
    qh_out = qh;                                   // positive_activity: q_H2O of the record
    wn_h = qh * kn;
    wv_h = qh * vb * kv;
    wp_h = qh * kphi;
    wn_c = qc * kn;
    wv_c = qc * vb * kv;
    phi_c = qc * kphi;
    const double clipped = c0 < -1.0 ? -1.0 : (c0 > 1.0 ? 1.0 : c0);
    const double c = night ? -1.0 : clipped;
    is_lit = c >= cos75;
    x = is_lit ? 1.0 / c - 1.0 : 0.0;
}

// positive_activity (2026-11-08, CO2 flux floor, CF kernels only): qc = Q_CO2 * ACO2(f) of record r,
// with the expressions of static_factors_m (left unchanged: finetune.cu copies it verbatim), so the
// value is bitwise the qc that static_factors_m formed wn_c, wv_c, phi_c from.
__device__ __forceinline__ double co2_flux_m(
    long long r, long long p, int f, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt, const double* __restrict__ aco2)
{
    const double qmin_c = pt[2 * npt + p], qmax_c = pt[3 * npt + p];
    const double c0 = (double)cos32[r];
    const bool night = c0 <= 0.0;                  // shadow folded into cos (-2)
    const double cz = c0 < 1.0 ? c0 : 1.0;
    const double qc0 = night ? qmin_c : qmin_c + (qmax_c - qmin_c) * cz;
    return qc0 * aco2[f];
}

// CF (positive_activity 2026-11-08): CO2 flux floor (header item 5); CF = false (the default, every
// kernel that existed before) is the 2026-11-07 code: each CF addition is compile-time dead there.
template <int NP, int G, bool CF = false>
__device__ __forceinline__ void fused_body_m(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const float* __restrict__ kn32, const float* __restrict__ kphi32,
    const float* __restrict__ kv32, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt,
    const double* __restrict__ aco2, double cos75,
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,
    const double* __restrict__ colbasis, long long nfaces, int col0,
    const double* __restrict__ pp, double* __restrict__ partial)
{
    // CF: the CO2-floor sums (codes 15..18, and 19..24 with the temperature sums) after the others
    constexpr int NF = CF ? (NP >= 14 ? 10 : 4) : 0;
    constexpr int N = (NP > 0 ? NP : 0) + 2 * G + (NP < 0 ? 1 : 0) + NF;
    __shared__ double tree[BLOCK];
    const Params P = load_params(pp);
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    // CF: kappa and the floor ratio r, the two values after the H2O floor block (positive_activity.py)
    const double kappa_c = CF ? aco2[2 * nfaces] : 0.0;
    const double ratio_c = CF ? aco2[2 * nfaces + 1] : 0.0;
    const int bf = (NP > 0 ? NP : 0) + 2 * G;          // CF: first CO2-floor accumulator
    double acc[N];
#pragma unroll
    for (int k = 0; k < N; ++k) acc[k] = 0.0;
    for (long long r = begin + threadIdx.x; r < end; r += BLOCK) {
        const int f = face[r];
        double wn_h, wv_h, wp_h, wn_c, wv_c, phi_c, x, qh;
        bool is_lit;
        static_factors_m(r, p, f, kn32, kphi32, kv32, cos32, pt, npt, aco2, cos75,
                         wn_h, wv_h, wp_h, wn_c, wv_c, phi_c, x, is_lit, qh);
        if (NP < 0) {                                  // PhiCO2 = sum qC * KPhi (code 14)
            acc[0] = acc[0] + phi_c;
            continue;
        }
        // record_factors with the recomputed xlit and lit
        const double af_face_value = af_face[f];       // positive_activity: A of the facet
        double af = af_face_value;                     // the record's activity factor
        // positive_activity: H2O flux floor F = max(q_H2O * A, F_min), F_min = qmax_H * (r_floor * A_orig)
        const double fmin = pt[npt + p] * aco2[nfaces + f];
        const bool at_floor = qh * af_face_value < fmin;
        if (at_floor) {
            const double kn = (double)kn32[r], kphi = (double)kphi32[r], kv = (double)kv32[r];
            wn_h = fmin * kn;                          // the record's H2O factors with F_min
            wv_h = fmin * pt[4 * npt + p] * kv;        // (F_min * VB) * KV
            wp_h = fmin * kphi;
            af = 1.0;                                  // x * 1.0 == x: the sums below use F_min
        }
        const double tl = P.Tmax + P.KT * (P.Tmax - P.T75) * x;
        const bool use_lit = is_lit && (tl > P.Tmin);
        const double te = use_lit ? tl : P.Tmin;
        const double u = 1.0 / sqrt(te);
        const double du = -0.5 * u / te;
        const double dl = use_lit ? du : 0.0;
        const double dn = use_lit ? 0.0 : du;
        const double dlx = dl * x;
        if (NP >= 5) {
            // positive_activity (CF only): CO2 flux floor F_C = max(kappa * qc, F_min,C),
            // F_min,C = qmax_C * (r * ACO2); at_floor_c is the constant false when CF is false
            bool at_floor_c = false;
            double fn_c = 0.0, fv_c = 0.0, fp_c = 0.0;
            if (CF) {
                const double fminc = pt[3 * npt + p] * (ratio_c * aco2[f]);
                at_floor_c = kappa_c * co2_flux_m(r, p, f, cos32, pt, npt, aco2) < fminc;
                if (at_floor_c) {
                    const double kn = (double)kn32[r], kphi = (double)kphi32[r], kv = (double)kv32[r];
                    fn_c = fminc * kn;                     // the record's CO2 factors with F_min,C
                    fv_c = fminc * pt[4 * npt + p] * kv;   // (F_min,C * VB) * KV
                    fp_c = fminc * kphi;
                }
            }
            acc[0] = acc[0] + wn_h * af * u;
            acc[1] = acc[1] + wp_h * af;
            acc[2] = acc[2] + wv_h * af * u;
            if (!at_floor_c) {                             // CF: active CO2 records only (x kappa on the host)
                acc[3] = acc[3] + wn_c * u;
                acc[4] = acc[4] + wv_c * u;
            }
            if (CF) {
                if (at_floor_c) {
                    acc[bf + 1] = acc[bf + 1] + fn_c * u;  // code 16
                    acc[bf + 2] = acc[bf + 2] + fv_c * u;  // code 17
                    acc[bf + 3] = acc[bf + 3] + fp_c;      // code 18
                } else {
                    acc[bf] = acc[bf] + phi_c;             // code 15: PhiCO2 of the active records
                }
            }
            if (NP >= 14) {
                const double w0 = !at_floor_c ? P.cNbmH * wn_h * af + P.kcNbmC * wn_c
                                              : P.cNbmH * wn_h * af;   // CF: gN without its CO2 term (19..21)
                const double w1 = P.bmH * wv_h * af;
                const double w2 = P.bmC * wv_c;
                acc[5] = acc[5] + w0 * dl;  acc[6] = acc[6] + w0 * dlx;  acc[7] = acc[7] + w0 * dn;
                acc[8] = acc[8] + w1 * dl;  acc[9] = acc[9] + w1 * dlx;  acc[10] = acc[10] + w1 * dn;
                if (!at_floor_c) {                         // CF: active CO2 records only (x kappa on the host)
                    acc[11] = acc[11] + w2 * dl; acc[12] = acc[12] + w2 * dlx; acc[13] = acc[13] + w2 * dn;
                }
                if (CF && at_floor_c) {
                    const double w3 = P.bmC * fv_c;        // bm_C * ((F_min,C * VB) * KV)
                    acc[bf + 4] = acc[bf + 4] + fn_c * dl;  acc[bf + 5] = acc[bf + 5] + fn_c * dlx;
                    acc[bf + 6] = acc[bf + 6] + fn_c * dn;                                    // codes 19..21
                    acc[bf + 7] = acc[bf + 7] + w3 * dl;    acc[bf + 8] = acc[bf + 8] + w3 * dlx;
                    acc[bf + 9] = acc[bf + 9] + w3 * dn;                                      // codes 22..24
                }
            }
        }
        if (G > 0) {
            // positive_activity: active 0 = clipped, 1 = active (clip), 2 = exp activity
            // (dA/dc = A * basis); a record at the H2O floor does not depend on the activity
            const unsigned char act = active[f];
            const bool on = act != 0 && !at_floor;
            const double wn_u = wn_h * u;
            const double rg = wp_h + P.bmH * wv_h * u;
            const int base = NP > 0 ? NP : 0;
#pragma unroll
            for (int g = 0; g < G; ++g) {
                const double bv = colbasis[(long long)(col0 + g) * nfaces + f];
                const double mb = on ? (act == 2 ? af_face_value * bv : bv) : 0.0 * bv;
                acc[base + g] = acc[base + g] + wn_u * mb;
                acc[base + G + g] = acc[base + G + g] + rg * mb;
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

#define COPS_FUSED_KERNEL_M(NAME, NP, G)                                                       \
extern "C" __global__ void __launch_bounds__(BLOCK) NAME(                                       \
    const long long* __restrict__ blocks, const int* __restrict__ face,                       \
    const float* __restrict__ kn32, const float* __restrict__ kphi32,                         \
    const float* __restrict__ kv32, const float* __restrict__ cos32,                          \
    const double* __restrict__ pt, long long npt,   \
    const double* __restrict__ aco2, double cos75,                                            \
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,             \
    const double* __restrict__ colbasis, long long nfaces, int col0,                          \
    const double* __restrict__ pp, double* __restrict__ partial)                              \
{                                                                                              \
    fused_body_m<NP, G>(blocks, face, kn32, kphi32, kv32, cos32, pt, npt, aco2, cos75,          \
                        af_face, active, colbasis, nfaces, col0, pp, partial);                 \
}

COPS_FUSED_KERNEL_M(eval_fused_m_phic, -1, 0)
COPS_FUSED_KERNEL_M(eval_fused_m_p5, 5, 0)
COPS_FUSED_KERNEL_M(eval_fused_m_p14, 14, 0)
COPS_FUSED_KERNEL_M(eval_fused_m_c1, 0, 1)
COPS_FUSED_KERNEL_M(eval_fused_m_c2, 0, 2)
COPS_FUSED_KERNEL_M(eval_fused_m_c3, 0, 3)
COPS_FUSED_KERNEL_M(eval_fused_m_c4, 0, 4)
COPS_FUSED_KERNEL_M(eval_fused_m_c5, 0, 5)
COPS_FUSED_KERNEL_M(eval_fused_m_c6, 0, 6)
COPS_FUSED_KERNEL_M(eval_fused_m_c7, 0, 7)
COPS_FUSED_KERNEL_M(eval_fused_m_c8, 0, 8)
COPS_FUSED_KERNEL_M(eval_fused_m_c9, 0, 9)
COPS_FUSED_KERNEL_M(eval_fused_m_c10, 0, 10)
COPS_FUSED_KERNEL_M(eval_fused_m_c11, 0, 11)
COPS_FUSED_KERNEL_M(eval_fused_m_c12, 0, 12)
COPS_FUSED_KERNEL_M(eval_fused_m_c13, 0, 13)
COPS_FUSED_KERNEL_M(eval_fused_m_c14, 0, 14)
COPS_FUSED_KERNEL_M(eval_fused_m_c15, 0, 15)
COPS_FUSED_KERNEL_M(eval_fused_m_c16, 0, 16)

// positive_activity (2026-11-08): the CO2-floor kernels (header item 5), launched only by
// positive_activity.py's evaluation with COPS_POS_FLOOR_CO2=1 (outputs: the codes of eval_fused_m_p5
// / eval_fused_m_p14, then 15..18 / 15..24). The activity/template columns need no CO2 variant.
#define COPS_FUSED_KERNEL_MF(NAME, NP, G)                                                      \
extern "C" __global__ void __launch_bounds__(BLOCK) NAME(                                       \
    const long long* __restrict__ blocks, const int* __restrict__ face,                       \
    const float* __restrict__ kn32, const float* __restrict__ kphi32,                         \
    const float* __restrict__ kv32, const float* __restrict__ cos32,                          \
    const double* __restrict__ pt, long long npt,                                             \
    const double* __restrict__ aco2, double cos75,                                            \
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,             \
    const double* __restrict__ colbasis, long long nfaces, int col0,                          \
    const double* __restrict__ pp, double* __restrict__ partial)                              \
{                                                                                              \
    fused_body_m<NP, G, true>(blocks, face, kn32, kphi32, kv32, cos32, pt, npt, aco2, cos75,    \
                              af_face, active, colbasis, nfaces, col0, pp, partial);           \
}

COPS_FUSED_KERNEL_MF(eval_fused_m_p5f, 5, 0)
COPS_FUSED_KERNEL_MF(eval_fused_m_p14f, 14, 0)
