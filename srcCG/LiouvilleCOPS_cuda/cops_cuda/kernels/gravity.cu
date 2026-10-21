// Exact homogeneous-polyhedron gravity (Tsoulis 2012) in FP64 on CUDA.
//
// SPDX-License-Identifier: GPL-3.0-only
//
// This file is a translation to CUDA C++ of parts of polyhedral-gravity-model
// 3.3 (Copyright (C) Jonas Schuhmacher and contributors,
// https://github.com/esa/polyhedral-gravity-model), which is distributed under
// the GNU General Public License version 3.  This file is distributed under the
// same license; the license text is in
// third_party/src_pypi/polyhedral_gravity-3.3/LICENSE.
//
// Derived from (sdist polyhedral_gravity-3.3.tar.gz, sha256 1476ab9d...7715):
//   src/polyhedralGravity/model/GravityEvaluable.cpp   evaluate, evaluateFace
//   src/polyhedralGravity/model/GravityModelDetail.cpp all per-face steps
//   src/polyhedralGravity/model/Polyhedron.cpp/.h      transformIterator, scaling
//   src/polyhedralGravity/util/UtilityContainer.h      cross, dot, norm, det, sgn
//   src/polyhedralGravity/util/UtilityConstants.h      PI, PI2, PI_2, G
// Modifications (2026-10-19, LiouvilleCOPS_cuda C03):
//   * scalar CUDA device code instead of std::array/thrust/xsimd containers;
//   * only the potential and the acceleration are accumulated (the
//     gradiometric tensor, its sum1T and singularity vector are not needed by
//     the campaign and are omitted; they do not influence the other outputs);
//   * per-face precomputed terms (segment vectors G, plane normal N, segment
//     normals n) are supplied by the host (cops_cuda/gravity.py, a translation
//     of GravityEvaluable::prepare) together with the face's three vertices;
//   * the final scaling by G*density is applied by the host;
//   * the HEAL/VERIFY orientation test is not translated (GRAVITY_REFERENCE.md
//     section 5.1); the caller admits a body by an independent topology audit.
//
// Please cite Schuhmacher et al. (2024), JOSS, doi:10.21105/joss.06384;
// Tsoulis (2012), Geophysics 77(2) F1-F11, doi:10.1190/geo2010-0334.1; and
// Tsoulis & Gavriilidou (2021), Geophys. Prospect. 69, doi:10.1111/1365-2478.13134.
//
// ---------------------------------------------------------------------------
// Numerical contract (GRAVITY_REFERENCE.md sections 3 and 4)
//
// Every arithmetic expression below is written in the evaluation order of the
// C++ source, and the module is compiled with --fmad=false, so each product
// and sum is rounded exactly as in the reference (which targets baseline
// x86-64 without FMA).  sqrt and division are correctly rounded IEEE
// operations on both sides.  The only intended differences are the
// transcendental functions: CUDA libdevice log/atan/acos versus glibc log/acos
// and xsimd atan (each within a few ulp).
//
// "Translated coordinates": the package moves the computation point P to the
// origin by subtracting P from every vertex (Polyhedron::transformIterator)
// and evaluates every face in those coordinates; the same subtraction is done
// here, per face and per target, on the original vertex values.  G, N and n
// are computed once from the ORIGINAL (untranslated) vertices, as in prepare().
//
// Summation order (gravity_fold): one thread per target folds the faces in
// index order, result = result + face_result starting from 0, exactly the
// sequential thrust::host transform_reduce of evaluate<false>.  Faces are
// staged through shared memory in tiles but consumed in order, and a launch
// may cover only a face range [face_begin, face_end) because the running sums
// are carried through global memory between launches, so the fold is the same
// for any tiling or launch split.  gravity_partial + gravity_reduce instead
// split the faces into chunks with a fixed two-level order; that changes the
// rounding (bounded by the C00 budget) but uses many more threads when there
// are few targets.
//
// Outputs per target (k = row of acc[NSUMS][nt], unscaled):
//   0 potential sum, 1..3 acceleration sum, 4 term magnitude M.
// Row 4 is not part of the package.  M is the sum, over faces, of the
// magnitudes of every term that carries a transcendental result into S:
//   sum_q |sigma_pq h_pq LN_pq| + h_p sum_q (|atan(t0_q)| + |atan(t1_q)|) + |sing_A|
// (|N S| = |S| because |N| = 1).  libdevice and glibc/xsimd results differ by a
// few ulp of each such term, so |g_gpu - g_ref| <= c eps M for a small c,
// independent of cancellation.  Cancellation is severe far from the body
// (each AN is the sum of two atans close to +-pi/2) and at the centre of a
// symmetric body; there |g| << M and a budget relative to |g| would be
// meaningless.  The C03 budget is stated relative to M (VALIDATION.md).
// Row 4 is accumulated separately and does not change the rounding of 0..3.
// ---------------------------------------------------------------------------

#define FACE_STRIDE 30      // v0[3] v1[3] v2[3] G0[3] G1[3] G2[3] N[3] n0[3] n1[3] n2[3]
#define NSUMS 5             // potential, ax, ay, az, M
#define EPSILON_ZERO_OFFSET 1e-14
#ifndef GRAVITY_BLOCK
#define GRAVITY_BLOCK 128
#endif
#ifndef GRAVITY_TILE
#define GRAVITY_TILE 32
#endif

// util::PI, PI2, PI_2 as written in UtilityConstants.h (identical doubles).
#define PG_PI   3.1415926535897932384626433832795028841971693993751058209749445923
#define PG_PI2  6.2831853071795864769252867665590057683943387987502116419498891846
#define PG_PI_2 1.5707963267948966192313216916397514420985846996875529104874722961

struct v3 { double x, y, z; };

__device__ __forceinline__ v3 make(double x, double y, double z)
{
    v3 r; r.x = x; r.y = y; r.z = z; return r;
}

__device__ __forceinline__ v3 sub(const v3 a, const v3 b)
{
    return make(a.x - b.x, a.y - b.y, a.z - b.z);
}

// util::cross: result[0] = l1*r2 - l2*r1, result[1] = l2*r0 - l0*r2, result[2] = l0*r1 - l1*r0
__device__ __forceinline__ v3 cross(const v3 l, const v3 r)
{
    return make(l.y * r.z - l.z * r.y, l.z * r.x - l.x * r.z, l.x * r.y - l.y * r.x);
}

// util::dot: l0*r0 + l1*r1 + l2*r2, left to right
__device__ __forceinline__ double dot(const v3 l, const v3 r)
{
    return l.x * r.x + l.y * r.y + l.z * r.z;
}

// util::euclideanNorm: sqrt(inner_product(x, x, 0.0)).  0.0 + x*x == x*x
// exactly (x*x is never -0), so the sum is (x*x + y*y) + z*z.
__device__ __forceinline__ double norm(const v3 a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// util::sgn with the absolute cut-off EPSILON_ZERO_OFFSET.
__device__ __forceinline__ int sgn(const double value)
{
    return value < -EPSILON_ZERO_OFFSET ? -1 : (value > EPSILON_ZERO_OFFSET ? 1 : 0);
}

// util::det (Sarrus) of the matrix whose ROWS are a, b, c:
// m00 m11 m22 + m01 m12 m20 + m02 m10 m21 - m02 m11 m20 - m00 m12 m21 - m01 m10 m22
__device__ __forceinline__ double det3(const v3 a, const v3 b, const v3 c)
{
    return a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y -
           a.z * b.y * c.x - a.x * b.z * c.y - a.y * b.x * c.z;
}

// projectPointOrthogonallyOntoSegment (eqs. 24-26): Cramer's rule on the
// transposed row matrix, exactly as written in the package.  Note that the
// segment vector here is recomputed from the TRANSLATED vertices.
__device__ __forceinline__ v3 project_onto_segment(const v3 v1, const v3 v2, const v3 pp)
{
    const v3 row1 = sub(v2, v1);
    const v3 row2 = cross(sub(v1, pp), row1);
    const v3 row3 = cross(row2, row1);
    const v3 d = make(dot(row1, pp), dot(row2, pp), dot(row3, v1));
    // columnMatrix = transpose({row1, row2, row3})
    const v3 c0 = make(row1.x, row2.x, row3.x);
    const v3 c1 = make(row1.y, row2.y, row3.y);
    const v3 c2 = make(row1.z, row2.z, row3.z);
    const double determinant = det3(c0, c1, c2);
    return make(det3(d, c1, c2) / determinant, det3(c0, d, c2) / determinant,
                det3(c0, c1, d) / determinant);
}

// GravityEvaluable::evaluateFace for one face and one computation point P,
// returning the face's (unscaled) potential and acceleration terms.
__device__ void face_contribution(const double* __restrict__ f,
                                  const double px, const double py, const double pz,
                                  double& potential, double& ax, double& ay, double& az,
                                  double& magnitude)
{
    // Translated vertices v_q - P (transformIterator) and precomputed terms.
    const v3 v[3] = {make(f[0] - px, f[1] - py, f[2] - pz),
                     make(f[3] - px, f[4] - py, f[5] - pz),
                     make(f[6] - px, f[7] - py, f[8] - pz)};
    const v3 G[3] = {make(f[9], f[10], f[11]), make(f[12], f[13], f[14]),
                     make(f[15], f[16], f[17])};
    const v3 N = make(f[18], f[19], f[20]);
    const v3 n[3] = {make(f[21], f[22], f[23]), make(f[24], f[25], f[26]),
                     make(f[27], f[28], f[29])};

    // 1-04 sigma_p = sgn(N . v0)  (computeUnitNormalOfPlaneDirection)
    const double sigma_p = (double)sgn(dot(N, v[0]));

    // 1-05 Hessian plane: (a,b,c) = (v0 - v1) x (v0 - v2), d = sum((0 - v0) * (a,b,c))
    const v3 abc = cross(sub(v[0], v[1]), sub(v[0], v[2]));
    const double hd = (0.0 - v[0].x) * abc.x + (0.0 - v[0].y) * abc.y + (0.0 - v[0].z) * abc.z;

    // 1-06 h_p = |d / sqrt(a^2 + b^2 + c^2)|
    const double h_p = fabs(hd / sqrt(abc.x * abc.x + abc.y * abc.y + abc.z * abc.z));

    // 1-07 P' = N h_p with the package's literal sign rules from the axis
    // intercepts d/a, d/b, d/c (0 when the denominator is exactly 0).
    const double Nc[3] = {N.x, N.y, N.z};
    const double coef[3] = {abc.x, abc.y, abc.z};
    double pc[3] = {N.x * h_p, N.y * h_p, N.z * h_p};
    for (int i = 0; i < 3; ++i) {
        const double intersection = coef[i] == 0.0 ? 0.0 : hd / coef[i];
        if (intersection < 0) {
            pc[i] = fabs(pc[i]);
        } else if (Nc[i] > 0) {
            pc[i] = -1.0 * pc[i];
        }
    }
    const v3 pp = make(pc[0], pc[1], pc[2]);

    // 1-08 sigma_pq = sgn(n_q . (P' - v_q)) * -1.0.  An int 0 times -1.0 is
    // -0.0, as in the package; every later test compares with == 0.0 or uses
    // fabs, so the sign of zero only matters where the package's does.
    double sigma[3];
    for (int q = 0; q < 3; ++q) {
        sigma[q] = (double)sgn(dot(n[q], sub(pp, v[q]))) * -1.0;
    }

    // 1-09 P''_q (P' itself when sigma_pq == 0) and 1-10 h_pq = |P''_q - P'|
    v3 ppp[3];
    double hpq[3];
    for (int q = 0; q < 3; ++q) {
        ppp[q] = sigma[q] == 0.0 ? pp : project_onto_segment(v[q], v[(q + 1) % 3], pp);
        hpq[q] = norm(sub(ppp[q], pp));
    }

    // 1-11 l1, l2, s1, s2 with the sign options of Tsoulis & Gavriilidou (2021)
    double l1[3], l2[3], s1[3], s2[3];
    for (int q = 0; q < 3; ++q) {
        const int q1 = (q + 1) % 3;
        double dl1 = norm(v[q]), dl2 = norm(v[q1]);
        double ds1 = norm(sub(ppp[q], v[q])), ds2 = norm(sub(ppp[q], v[q1]));
        if (fabs(ds1 - dl1) < EPSILON_ZERO_OFFSET && fabs(ds2 - dl2) < EPSILON_ZERO_OFFSET) {
            // Option 4: P on the segment's line
            if (ds2 < ds1) {                      // case 2
                ds1 *= -1.0; ds2 *= -1.0; dl1 *= -1.0; dl2 *= -1.0;
            } else if (fabs(ds2 - ds1) < EPSILON_ZERO_OFFSET) {   // case 1
                ds1 *= -1.0; dl1 *= -1.0;
            }                                     // case 3: unchanged
        } else {
            const double g = norm(G[q]);          // precomputed segment vector
            if (ds1 < g && ds2 < g) {             // option 1: P'' inside the segment
                ds1 *= -1.0;
            } else if (ds2 < ds1) {               // option 2: right of the segment
                ds1 *= -1.0; ds2 *= -1.0;
            }                                     // option 3: unchanged
        }
        l1[q] = dl1; l2[q] = dl2; s1[q] = ds1; s2[q] = ds2;
    }

    // 1-12 |P' - v_k|
    const double r[3] = {norm(sub(pp, v[0])), norm(sub(pp, v[1])), norm(sub(pp, v[2]))};

    // 1-13 LN_pq (eq. 14) and AN_pq (eq. 15)
    double ln[3], an[3], an_mag[3];
    for (int q = 0; q < 3; ++q) {
        const double r1n = r[(q + 1) % 3];
        const double r2n = r[q];
        if ((sigma[q] == 0.0 && (r1n < EPSILON_ZERO_OFFSET || r2n < EPSILON_ZERO_OFFSET)) ||
            (fabs(s1[q] + s2[q]) < EPSILON_ZERO_OFFSET && fabs(l1[q] + l2[q]) < EPSILON_ZERO_OFFSET)) {
            ln[q] = 0.0;
        } else {
            ln[q] = log((s2[q] + l2[q]) / (s1[q] + l1[q]));
        }
        if (h_p < EPSILON_ZERO_OFFSET || hpq[q] < EPSILON_ZERO_OFFSET) {
            an[q] = 0.0;
            an_mag[q] = 0.0;
        } else {
            // xsimd batch (s2, s1)*h_p / ((l2, -l1)*h_pq), atan, reduce_add.
            const double t0 = (s2[q] * h_p) / (l2[q] * hpq[q]);
            const double t1 = (s1[q] * h_p) / ((-l1[q]) * hpq[q]);
            const double a0 = atan(t0), a1 = atan(t1);
            an[q] = a0 + a1;
            an_mag[q] = fabs(a0) + fabs(a1);
        }
    }

    // 1-14 singularity term sing_A (computeSingularityTerms); first match wins.
    double sing_a;
    if (sigma[0] == 1.0 && sigma[1] == 1.0 && sigma[2] == 1.0) {
        // case 1: P' inside the face
        sing_a = -1.0 * PG_PI2 * h_p;
    } else {
        bool on_edge = false;
        for (int q = 0; q < 3 && !on_edge; ++q) {
            if (fabs(sigma[q]) > EPSILON_ZERO_OFFSET) continue;
            const double gn = norm(G[q]);
            on_edge = r[(q + 1) % 3] < gn && r[q] < gn;
        }
        if (on_edge) {
            // case 2: P' on a segment, not on a vertex
            sing_a = -1.0 * PG_PI * h_p;
        } else {
            // case 3: P' at a vertex; the package selects the first such q and
            // then tests r1 == 0.0 EXACTLY (not the eps test) to pick g1, g2.
            int j = -1;
            double r1sel = 0.0;
            for (int q = 0; q < 3; ++q) {
                if (fabs(sigma[q]) > EPSILON_ZERO_OFFSET) continue;
                const double r1n = r[(q + 1) % 3];
                const double r2n = r[q];
                if (r1n < EPSILON_ZERO_OFFSET || r2n < EPSILON_ZERO_OFFSET) {
                    j = q; r1sel = r1n; break;
                }
            }
            if (j >= 0) {
                const v3 g1 = r1sel == 0.0 ? G[j] : G[(j + 2) % 3];
                const v3 g2 = r1sel == 0.0 ? G[(j + 1) % 3] : G[j];
                const double gdot = dot(make(g1.x * -1.0, g1.y * -1.0, g1.z * -1.0), g2);
                const double theta = gdot == 0.0 ? PG_PI_2 : acos(gdot / (norm(g1) * norm(g2)));
                sing_a = -1.0 * theta * h_p;
            } else {
                sing_a = 0.0;                     // case 4
            }
        }
    }

    // Steps 2, 4, 5: sums in q order (std::accumulate from 0.0)
    double sum1 = 0.0, sum2 = 0.0;
    for (int q = 0; q < 3; ++q) sum1 = sum1 + sigma[q] * hpq[q] * ln[q];
    for (int q = 0; q < 3; ++q) sum2 = sum2 + sigma[q] * an[q];
    const double plane_sum = sum1 + h_p * sum2 + sing_a;

    // Step 7: potential sigma_p h_p S, acceleration N S (before scaling)
    potential = sigma_p * h_p * plane_sum;
    ax = N.x * plane_sum;
    ay = N.y * plane_sum;
    az = N.z * plane_sum;

    // Term magnitude (diagnostic, see the header); |sigma_pq| <= 1.
    double term = 0.0, an_term = 0.0;
    for (int q = 0; q < 3; ++q) term = term + fabs(sigma[q] * hpq[q] * ln[q]);
    for (int q = 0; q < 3; ++q) an_term = an_term + fabs(sigma[q]) * an_mag[q];
    magnitude = term + h_p * an_term + fabs(sing_a);
}

// Cooperative load of faces [base, base+count) into shared memory.
__device__ __forceinline__ void load_tile(double* tile, const double* __restrict__ faces,
                                          long long base, int count)
{
    const double* src = faces + base * FACE_STRIDE;
    for (int k = threadIdx.x; k < count * FACE_STRIDE; k += blockDim.x) tile[k] = src[k];
}

// One thread per target; continues the running sums acc[NSUMS][nt] over faces
// [face_begin, face_end) in index order.  acc must be zero before the first
// launch (the package's result{} initial value).
extern "C" __global__ void __launch_bounds__(GRAVITY_BLOCK)
gravity_fold(const double* __restrict__ faces, long long face_begin, long long face_end,
             const double* __restrict__ targets, long long ntargets, double* __restrict__ acc)
{
    __shared__ double tile[GRAVITY_TILE * FACE_STRIDE];
    const long long t = (long long)blockIdx.x * GRAVITY_BLOCK + threadIdx.x;
    const bool active = t < ntargets;
    double px = 0.0, py = 0.0, pz = 0.0, sp = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, sm = 0.0;
    if (active) {
        px = targets[t]; py = targets[ntargets + t]; pz = targets[2 * ntargets + t];
        sp = acc[t]; sx = acc[ntargets + t]; sy = acc[2 * ntargets + t]; sz = acc[3 * ntargets + t];
        sm = acc[4 * ntargets + t];
    }
    for (long long base = face_begin; base < face_end; base += GRAVITY_TILE) {
        const int count = (int)(face_end - base < GRAVITY_TILE ? face_end - base : GRAVITY_TILE);
        __syncthreads();                         // previous tile fully consumed
        load_tile(tile, faces, base, count);
        __syncthreads();
        if (active) {
            for (int i = 0; i < count; ++i) {
                double p, a, b, c, m;
                face_contribution(tile + i * FACE_STRIDE, px, py, pz, p, a, b, c, m);
                // util::operator+ on the result tuple: element-wise, in face order
                sp = sp + p; sx = sx + a; sy = sy + b; sz = sz + c; sm = sm + m;
            }
        }
    }
    if (active) {
        acc[t] = sp; acc[ntargets + t] = sx; acc[2 * ntargets + t] = sy; acc[3 * ntargets + t] = sz;
        acc[4 * ntargets + t] = sm;
    }
}

// Face-split variant, level 1: block (x, y) folds face chunk
// c = chunk_offset + y for its targets from 0 and writes partial[c][NSUMS][nt].
// chunk_offset lets the host spread the chunks over several short launches.
extern "C" __global__ void __launch_bounds__(GRAVITY_BLOCK)
gravity_partial(const double* __restrict__ faces, long long nfaces, long long chunk_faces,
                long long chunk_offset, const double* __restrict__ targets, long long ntargets,
                double* __restrict__ partial)
{
    __shared__ double tile[GRAVITY_TILE * FACE_STRIDE];
    const long long t = (long long)blockIdx.x * GRAVITY_BLOCK + threadIdx.x;
    const long long chunk = chunk_offset + blockIdx.y;
    const long long face_begin = chunk * chunk_faces;
    const long long face_end = face_begin + chunk_faces < nfaces ? face_begin + chunk_faces : nfaces;
    const bool active = t < ntargets;
    double px = 0.0, py = 0.0, pz = 0.0, sp = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, sm = 0.0;
    if (active) {
        px = targets[t]; py = targets[ntargets + t]; pz = targets[2 * ntargets + t];
    }
    for (long long base = face_begin; base < face_end; base += GRAVITY_TILE) {
        const int count = (int)(face_end - base < GRAVITY_TILE ? face_end - base : GRAVITY_TILE);
        __syncthreads();
        load_tile(tile, faces, base, count);
        __syncthreads();
        if (active) {
            for (int i = 0; i < count; ++i) {
                double p, a, b, c, m;
                face_contribution(tile + i * FACE_STRIDE, px, py, pz, p, a, b, c, m);
                sp = sp + p; sx = sx + a; sy = sy + b; sz = sz + c; sm = sm + m;
            }
        }
    }
    if (active) {
        double* out = partial + chunk * NSUMS * ntargets;
        out[t] = sp; out[ntargets + t] = sx; out[2 * ntargets + t] = sy; out[3 * ntargets + t] = sz;
        out[4 * ntargets + t] = sm;
    }
}

// Face-split variant, level 2: out[k][t] = (((0 + p_0) + p_1) + ...) in chunk
// order -- a fixed order, so the result is independent of scheduling.
extern "C" __global__
void gravity_reduce(const double* __restrict__ partial, long long nchunks, long long ntargets,
                    double* __restrict__ out)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;   // i = k*nt + t
    if (i >= NSUMS * ntargets) return;
    double s = 0.0;
    for (long long c = 0; c < nchunks; ++c) s = s + partial[c * NSUMS * ntargets + i];
    out[i] = s;
}

// ---------------------------------------------------------------------------
// MIXED PRECISION (experiment, user request 2026-10-21; not the C03/C04 reference)
//
// The exact plane sum S_p of face_contribution is the face integral
//     S_p = integral over face p of dS / |r - P|
// (the line-integral form of Tsoulis 2012), and the face's terms are
//     potential = sigma_p h_p S_p = (N . (v - P)) S_p,   acceleration = N S_p.
// For a face far from P the closed form cancels badly (each AN is the sum of two atans
// near +-pi/2), so FP32 cannot evaluate it. The integral itself is smooth there, and a
// 3-point (degree-2) triangle quadrature in FP32 approximates it with an error of order
// (size/distance)^3:
//     S_p ~ (A/3) sum_k 1/|q_k - P|,  q_k at barycentric (2/3, 1/6, 1/6) and permutations.
// gravity_mixed_fold uses, per face and target:
//   * near (|c - P|^2 < rnear2 = (near_factor * circumradius)^2): the exact FP64
//     face_contribution above, unchanged;
//   * far: the FP32 quadrature, the face's (potential, acceleration, magnitude) in FP32,
//     summed in FP32 over one tile of MIXED_TILE faces, and the tile sums added to the
//     FP64 running sums.
// The magnitude row M of far faces is S_p (|N S_p| = S_p), so M stays an upper bound
// of the summed terms. near_count[t] counts the exact faces of target t.
// faces32 row (FACE32_STRIDE floats): c[3] q0[3] q1[3] q2[3] N[3] wA rnear2 pad pad pad.

#define FACE32_STRIDE 20
#ifndef MIXED_TILE
#define MIXED_TILE 64
#endif

extern "C" __global__ void __launch_bounds__(GRAVITY_BLOCK)
gravity_mixed_fold(const double* __restrict__ faces, const float* __restrict__ faces32,
                   long long face_begin, long long face_end,
                   const double* __restrict__ targets, long long ntargets,
                   double* __restrict__ acc, unsigned int* __restrict__ near_count)
{
    __shared__ float tile[MIXED_TILE * FACE32_STRIDE];
    const long long t = (long long)blockIdx.x * GRAVITY_BLOCK + threadIdx.x;
    const bool active = t < ntargets;
    double px = 0.0, py = 0.0, pz = 0.0, sp = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, sm = 0.0;
    unsigned int nnear = 0;
    if (active) {
        px = targets[t]; py = targets[ntargets + t]; pz = targets[2 * ntargets + t];
        sp = acc[t]; sx = acc[ntargets + t]; sy = acc[2 * ntargets + t]; sz = acc[3 * ntargets + t];
        sm = acc[4 * ntargets + t];
    }
    const float fx = (float)px, fy = (float)py, fz = (float)pz;
    for (long long base = face_begin; base < face_end; base += MIXED_TILE) {
        const int count = (int)(face_end - base < MIXED_TILE ? face_end - base : MIXED_TILE);
        __syncthreads();
        const float* src = faces32 + base * FACE32_STRIDE;
        for (int k = threadIdx.x; k < count * FACE32_STRIDE; k += blockDim.x) tile[k] = src[k];
        __syncthreads();
        if (active) {
            float tp = 0.0f, tx = 0.0f, ty = 0.0f, tz = 0.0f, tm = 0.0f;
            for (int i = 0; i < count; ++i) {
                const float* f = tile + i * FACE32_STRIDE;
                const float dx = f[0] - fx, dy = f[1] - fy, dz = f[2] - fz;
                if (dx * dx + dy * dy + dz * dz < f[16]) {
                    double p, a, b, c, m;
                    face_contribution(faces + (base + i) * FACE_STRIDE, px, py, pz, p, a, b, c, m);
                    sp = sp + p; sx = sx + a; sy = sy + b; sz = sz + c; sm = sm + m;
                    ++nnear;
                } else {
                    float s = 0.0f;
                    for (int k = 0; k < 3; ++k) {
                        const float qx = f[3 + 3 * k] - fx, qy = f[4 + 3 * k] - fy, qz = f[5 + 3 * k] - fz;
                        s = s + rsqrtf(qx * qx + qy * qy + qz * qz);
                    }
                    s = s * f[15];
                    const float h = f[12] * dx + f[13] * dy + f[14] * dz;     // N . (c - P)
                    tp = tp + h * s; tx = tx + f[12] * s; ty = ty + f[13] * s; tz = tz + f[14] * s;
                    tm = tm + s;
                }
            }
            sp = sp + (double)tp; sx = sx + (double)tx; sy = sy + (double)ty; sz = sz + (double)tz;
            sm = sm + (double)tm;
        }
    }
    if (active) {
        acc[t] = sp; acc[ntargets + t] = sx; acc[2 * ntargets + t] = sy; acc[3 * ntargets + t] = sz;
        acc[4 * ntargets + t] = sm;
        near_count[t] += nnear;
    }
}
