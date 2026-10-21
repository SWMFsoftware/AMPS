// Device-side audit of a resident transport-kernel store (ROADMAP C02 gate).
//
// The host validates every record while streaming it from disk; that check is
// authoritative.  This audit repeats the record invariants on the bytes that
// actually reside in device memory after the upload, so an indexing error in
// the staging/transfer code (a record copied to the wrong offset, a field
// written into another field, an uninitialised tail) is detected instead of
// being read later by the pressure kernels.
//
// Determinism: the outputs are integer COUNTS combined with atomicAdd on
// unsigned 64-bit counters.  Integer addition is associative, so the result is
// independent of thread scheduling (unlike floating-point atomics).  Each
// thread accumulates privately and issues at most one atomic per counter.
//
// Non-finite test: a double is NaN or +-inf exactly when all 11 exponent bits
// are set.  Testing the bits avoids any dependence on math-library macros and
// cannot be affected by compiler floating-point options.
//
// Counter layout (must match device_store.AUDIT_COUNTERS):
//   0  records whose face index is outside [0, nfaces)
//   1  records with a non-finite Kn, KPhi, KV or cos value
//   2  records whose shadow flag is not 0 or 1
//   3  points whose CSR offsets decrease (offsets[p+1] < offsets[p])
//   4  CSR endpoints wrong (offsets[0] != 0 or offsets[npt] != nrecords)

__device__ __forceinline__ bool nonfinite(double value)
{
    const unsigned long long bits = (unsigned long long)__double_as_longlong(value);
    return (bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL;
}

extern "C" __global__
void audit_records(const int* __restrict__ face,
                   const double* __restrict__ kn,
                   const double* __restrict__ kphi,
                   const double* __restrict__ kv,
                   const double* __restrict__ cosine,
                   const unsigned char* __restrict__ shadow,
                   long long nrecords, int nfaces,
                   unsigned long long* counters)
{
    unsigned long long bad_face = 0, bad_value = 0, bad_shadow = 0;
    // 64-bit grid-stride loop: the campaign store has 2.46e8 records and
    // larger caches may exceed 2^31, so the index must not be a 32-bit int.
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x;
         r < nrecords; r += stride) {
        const int f = face[r];
        bad_face += (f < 0 || f >= nfaces);
        bad_value += (nonfinite(kn[r]) || nonfinite(kphi[r]) ||
                      nonfinite(kv[r]) || nonfinite(cosine[r]));
        bad_shadow += (shadow[r] > 1);
    }
    if (bad_face) atomicAdd(&counters[0], bad_face);
    if (bad_value) atomicAdd(&counters[1], bad_value);
    if (bad_shadow) atomicAdd(&counters[2], bad_shadow);
}

extern "C" __global__
void audit_offsets(const long long* __restrict__ offsets, long long npoints,
                   long long nrecords, unsigned long long* counters)
{
    unsigned long long decreasing = 0;
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long p = (long long)blockIdx.x * blockDim.x + threadIdx.x;
         p < npoints; p += stride) {
        decreasing += (offsets[p + 1] < offsets[p]);
    }
    if (decreasing) atomicAdd(&counters[3], decreasing);
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        if (offsets[0] != 0 || offsets[npoints] != nrecords) atomicAdd(&counters[4], 1ULL);
    }
}

// FP32 record storage (mixed precision, precision.py "fp32"): the same audit on float fields;
// the shadow flag is folded into cos there (no shadow array; it was checked on the host).
__device__ __forceinline__ bool nonfinite_f(float value)
{
    const unsigned int bits = (unsigned int)__float_as_int(value);
    return (bits & 0x7f800000U) == 0x7f800000U;
}

extern "C" __global__
void audit_records_f32(const int* __restrict__ face,
                       const float* __restrict__ kn,
                       const float* __restrict__ kphi,
                       const float* __restrict__ kv,
                       const float* __restrict__ cosine,
                       long long nrecords, int nfaces,
                       unsigned long long* counters)
{
    unsigned long long bad_face = 0, bad_value = 0, bad_shadow = 0;
    const long long stride = (long long)gridDim.x * blockDim.x;
    for (long long r = (long long)blockIdx.x * blockDim.x + threadIdx.x;
         r < nrecords; r += stride) {
        const int f = face[r];
        bad_face += (f < 0 || f >= nfaces);
        bad_value += (nonfinite_f(kn[r]) || nonfinite_f(kphi[r]) ||
                      nonfinite_f(kv[r]) || nonfinite_f(cosine[r]));
    }
    if (bad_face) atomicAdd(&counters[0], bad_face);
    if (bad_value) atomicAdd(&counters[1], bad_value);
    if (bad_shadow) atomicAdd(&counters[2], bad_shadow);
}
