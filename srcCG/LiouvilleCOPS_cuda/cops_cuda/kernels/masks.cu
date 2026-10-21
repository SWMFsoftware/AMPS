// Geometry-only visibility masks (ROADMAP C05, section 3.2).
//
// Per record r of a selected observation point (same point for all four
// predicates, exactly as build_physical_activity_templates.late_gauge_template):
//   point_ng = Kn > 0
//   point_rg = (KPhi != 0) | (KV != 0)
// and the facet face[r] is marked in
//   bit 0  ng_seen         point_ng
//   bit 1  rg_seen         point_rg
//   bit 2  late_ng_not_rg  point_ng & !point_rg
//   bit 3  late_rg_not_ng  point_rg & !point_ng
// The union over points is formed with atomicOr on 32-bit words, one word per
// facet holding the four bits.  Integer OR is associative and idempotent, so
// the result is independent of thread scheduling and duplicates; there are no
// unsynchronized Boolean stores (ROADMAP section 3.2).
//
// The predicates read the ORIGINAL FP64 kernel values; no rounding or
// thresholding precedes them, so a tiny positive value can never be lost.

extern "C" __global__ void visibility_union(
    const int* __restrict__ face, const double* __restrict__ kn,
    const double* __restrict__ kphi, const double* __restrict__ kv,
    long long nrecords, unsigned int* __restrict__ bits)
{
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x; r < nrecords; r += stride) {
        const bool ng = kn[r] > 0.0;
        const bool rg = (kphi[r] != 0.0) || (kv[r] != 0.0);
        const unsigned int value = (ng ? 1u : 0u) | (rg ? 2u : 0u) |
                                   ((ng && !rg) ? 4u : 0u) | ((rg && !ng) ? 8u : 0u);
        if (value) atomicOr(&bits[face[r]], value);
    }
}

// Facets that contribute to at least one loaded point, exactly as
// fit_cops.Model.seen_faces:  m[face[(Kn > 0) | (KPhi > 0) | (KV > 0)]] = True.
// Note the predicate is "> 0" for all three kernels (not the "!= 0" of the RG
// visibility mask above).  Every store writes the same value 1, so concurrent
// stores to one facet are benign and the result is schedule-independent.
extern "C" __global__ void seen_positive(
    const int* __restrict__ face, const double* __restrict__ kn,
    const double* __restrict__ kphi, const double* __restrict__ kv,
    long long nrecords, unsigned char* __restrict__ seen)
{
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x; r < nrecords; r += stride) {
        if (kn[r] > 0.0 || kphi[r] > 0.0 || kv[r] > 0.0) seen[face[r]] = 1;
    }
}

// FP32 record storage (mixed precision): the same predicates on float fields. Zero-ness is
// exactly that of the FP64 file values (the loader refuses any record where FP32 changes it).
extern "C" __global__ void visibility_union_f32(
    const int* __restrict__ face, const float* __restrict__ kn,
    const float* __restrict__ kphi, const float* __restrict__ kv,
    long long nrecords, unsigned int* __restrict__ bits)
{
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x; r < nrecords; r += stride) {
        const bool ng = kn[r] > 0.0f;
        const bool rg = (kphi[r] != 0.0f) || (kv[r] != 0.0f);
        const unsigned int value = (ng ? 1u : 0u) | (rg ? 2u : 0u) |
                                   ((ng && !rg) ? 4u : 0u) | ((rg && !ng) ? 8u : 0u);
        if (value) atomicOr(&bits[face[r]], value);
    }
}

extern "C" __global__ void seen_positive_f32(
    const int* __restrict__ face, const float* __restrict__ kn,
    const float* __restrict__ kphi, const float* __restrict__ kv,
    long long nrecords, unsigned char* __restrict__ seen)
{
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x; r < nrecords; r += stride) {
        if (kn[r] > 0.0f || kphi[r] > 0.0f || kv[r] > 0.0f) seen[face[r]] = 1;
    }
}
