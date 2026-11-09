// ==== version-controlled copy: header added 2026-11-08, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/model_ext/kernels/finetune.cu (on Voyager: /data001/highresolu_vt/cops_output/model_ext/kernels/finetune.cu)
// SHA256 of the run copy = SHA256 of this file without these header lines: b4aa9009400d01b0bc0939c4df0dbd293e64e2c2939a04584be76c917ac0bb15
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: CUDA source. ft_jtu_m (J^T u scatter over the records, H2O part, active records only) and ft_seen_m (seen
//   facets per point) for positive_finetune.py; derived from kernels/pressure.cu sha 714b9d90 (the driver refuses
//   another). user decisions 2026-11-07/08.
// ==== end of header ====
// =============================================================================================
// model_ext/kernels/finetune.cu -- per-element fine tuning of the positive-activity COPS fits
// (model_ext/positive_finetune.py, user decisions 2026-11-07). NEW kernels; nothing in the git
// tree and no other model_ext kernel is changed. Compiled by positive_finetune.py with cupy
// RawModule and the reference NVRTC options (--std=c++14 --fmad=false) in every GPU worker.
//
// Derived from model_ext/kernels/pressure.cu (positive_activity.py copy of the evaluator kernels)
//     sha256 714b9d908fee29c8135931ecd8a5209fba5ef1dd6905bec5b219cf5a0a4280f3
// positive_finetune.py refuses to run when the pressure.cu in use has another SHA-256: the
// transpose below must be the transpose of THAT forward/Jacobian code (static_factors_m and the
// H2O floor / exp-activity logic of fused_body_m are copied verbatim, expression by expression).
//
// Parameterization (per facet e): A_e = exp(s_e), s = s_fit + held + delta; d A_e / d delta_e = A_e.
// The forward pressures and J v use the UNCHANGED evaluator kernels (eval_fused_m_p5 and
// eval_fused_m_c<k> with colbasis column = v and active == 2: column factor mb = A_f * v_f, records
// at the H2O floor contribute 0). This file adds what the evaluator cannot do:
//
//   ft_jtu_m   J^T u for the H2O activity block. For every record r (local point p, device facet f)
//              of ONE point (the host launches the points one after another on one stream, so the
//              per-facet sums are accumulated in point order: deterministic):
//                  g[f] += (cN * wn_h*u + cR * (wp_h + bm_H*wv_h*u)) * mbf,
//              mbf = A_f (active 2, exp) | 1 (active 1, clip) ; 0 for active 0 or a record at the floor,
//              cN = u_NG(p) * cNG_H * bm_H, cR = u_RG(p) * cRG_H(p) (host), i.e. exactly the terms whose
//              point sums eval_fused_m_c1 forms for J v: JNG(p) = cNG_H bm_H sum_r wn_h u mb,
//              JRG(p) = cRG_H(p) sum_r (wp_h + bm_H wv_h u) mb, mb = A_f v_f. CO2 does not depend on the
//              H2O activity. atomicAdd is used for safety (a point holds each facet once; then no two
//              threads of a launch touch the same facet and the result does not depend on scheduling).
//              2026-11-08 (speed-up): use_atomic = 0 replaces atomicAdd by g[f] = g[f] + term, the SAME
//              IEEE addition in the same order (per facet: point order); used only when ft_dup_m found no
//              facet twice within any point (then no two threads of a launch touch one facet).
//   ft_dup_m   once per bind: facets listed twice within one point (atomicExch of the point index into
//              last[f], initialised to -1; the launches of the points are sequential).
//   ft_map     dst[dst_idx[i]] = src[src_idx[i]]: device facets <-> parameters (2026-11-08 speed-up).
//   ft_seen_m  facets seen from the selected points: the reference predicate of fit_cops
//              Model.seen_faces / masks.cu seen_positive_f32 (Kn > 0 | KPhi > 0 | KV > 0), restricted
//              to points with mask[p] != 0 (benign same-value stores); structural != 0 (2026-11-08): the
//              structural-contribution predicate of the old late-window fine tuning
//              (../LiouvilleCOPS/FINE_TUNING_STEP.md section 2: Kn != 0 | KPhi != 0 | KV != 0).
// =============================================================================================

#define BLOCK 256

struct Params {
    double Tmax, T75, Tmin, KT;      // temperature law
    double cNbmH;                    // cNG_H * bm_H
    double kcNbmC;                   // (kappa * cNG_C) * bm_C
    double bmH, bmC;                 // bm_H, bm_C
};

__device__ __forceinline__ Params load_params(const double* __restrict__ pp)
{
    Params P;
    P.Tmax = pp[0]; P.T75 = pp[1]; P.Tmin = pp[2]; P.KT = pp[3];
    P.cNbmH = pp[4]; P.kcNbmC = pp[5]; P.bmH = pp[6]; P.bmC = pp[7];
    return P;
}

// verbatim copy of static_factors_m (model_ext/kernels/pressure.cu)
__device__ __forceinline__ void static_factors_m(
    long long r, long long p, int f,
    const float* __restrict__ kn32, const float* __restrict__ kphi32,
    const float* __restrict__ kv32, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt,
    const double* __restrict__ aco2, double cos75,
    double& wn_h, double& wv_h, double& wp_h, double& wn_c, double& wv_c, double& phi_c,
    double& x, bool& is_lit, double& qh_out)
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
    qh_out = qh;
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

// J^T u, H2O activity block, records of ONE local point p (blocks = the point's block-table rows).
extern "C" __global__ void __launch_bounds__(BLOCK) ft_jtu_m(
    const long long* __restrict__ blocks, const int* __restrict__ face,
    const float* __restrict__ kn32, const float* __restrict__ kphi32,
    const float* __restrict__ kv32, const float* __restrict__ cos32,
    const double* __restrict__ pt, long long npt,
    const double* __restrict__ aco2, double cos75,
    const double* __restrict__ af_face, const unsigned char* __restrict__ active,
    long long nfaces, const double* __restrict__ pp,
    long long p_expected, double cN, double cR,
    double* __restrict__ g, unsigned long long* __restrict__ bad, int use_atomic)
{
    const Params P = load_params(pp);
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    if (p != p_expected) {                          // block table / launch mismatch: reported, no sums
        if (threadIdx.x == 0) atomicAdd(bad, 1ULL);
        return;
    }
    for (long long r = begin + threadIdx.x; r < end; r += BLOCK) {
        const int f = face[r];
        double wn_h, wv_h, wp_h, wn_c, wv_c, phi_c, x, qh;
        bool is_lit;
        static_factors_m(r, p, f, kn32, kphi32, kv32, cos32, pt, npt, aco2, cos75,
                         wn_h, wv_h, wp_h, wn_c, wv_c, phi_c, x, is_lit, qh);
        // fused_body_m: H2O flux floor and the activity column factor
        const double af_face_value = af_face[f];
        const double fmin = pt[npt + p] * aco2[nfaces + f];
        const bool at_floor = qh * af_face_value < fmin;
        const unsigned char act = active[f];
        if (act == 0 || at_floor) continue;           // column factor 0 (mb = 0 * basis)
        const double tl = P.Tmax + P.KT * (P.Tmax - P.T75) * x;
        const bool use_lit = is_lit && (tl > P.Tmin);
        const double te = use_lit ? tl : P.Tmin;
        const double u = 1.0 / sqrt(te);
        const double wn_u = wn_h * u;
        const double rg = wp_h + P.bmH * wv_h * u;
        const double mbf = act == 2 ? af_face_value : 1.0;
        const double term = (cN * wn_u + cR * rg) * mbf;
        if (use_atomic) atomicAdd(&g[f], term);
        else g[f] = g[f] + term;
    }
}

// Facets listed twice within one point (records of ONE local point p; launches sequential over p).
extern "C" __global__ void ft_dup_m(
    const long long* __restrict__ blocks, const int* __restrict__ face, long long p_expected,
    int* __restrict__ last, unsigned long long* __restrict__ dup)
{
    const long long b = blockIdx.x;
    const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
    if (p != p_expected) {
        if (threadIdx.x == 0) atomicAdd(dup, 1ULL);
        return;
    }
    for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
        const int old = atomicExch(&last[face[r]], (int)p);
        if (old == (int)p) atomicAdd(dup, 1ULL);
    }
}

// dst[dst_idx[i]] = src[src_idx[i]] for i < n: the map between device facets and parameters (both index arrays
// unique). Used for the direction (parameter space -> colbasis) and for J^T u (device facets -> parameter space).
extern "C" __global__ void ft_map(
    long long n, const long long* __restrict__ dst_idx, double* __restrict__ dst,
    const long long* __restrict__ src_idx, const double* __restrict__ src)
{
    for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += (long long)gridDim.x * blockDim.x)
        dst[dst_idx[i]] = src[src_idx[i]];
}

// Facets seen from the points with mask[p] != 0 (all blocks of the device in one launch).
extern "C" __global__ void ft_seen_m(
    const long long* __restrict__ blocks, long long nblocks, const int* __restrict__ face,
    const float* __restrict__ kn32, const float* __restrict__ kphi32, const float* __restrict__ kv32,
    const unsigned char* __restrict__ mask, unsigned char* __restrict__ seen, int structural)
{
    for (long long b = blockIdx.x; b < nblocks; b += gridDim.x) {
        const long long p = blocks[3 * b], begin = blocks[3 * b + 1], end = blocks[3 * b + 2];
        if (!mask[p]) continue;
        for (long long r = begin + threadIdx.x; r < end; r += blockDim.x) {
            const bool hit = structural ? (kn32[r] != 0.0f || kphi32[r] != 0.0f || kv32[r] != 0.0f)
                                        : (kn32[r] > 0.0f || kphi32[r] > 0.0f || kv32[r] > 0.0f);
            if (hit) seen[face[r]] = 1;
        }
    }
}
