// GPU implementation of the two expensive loops of the LiouvilleCOPS forward model (ROADMAP C11).
//
// liouville_cops_cuda = the UNCHANGED ../LiouvilleCOPS/src/*.cpp + this file + access_shim.cpp.
// In the translation units of Liouville.cpp and Source.cpp the CPU definitions of
//   cLiouville::EvaluateKernels  (renamed EvaluateKernelsCPU by -D in native/Makefile)
//   cSourceModel::Evaluate       (renamed EvaluateCPU)
// are set aside, and the definitions below take their place.  Everything else -- SPICE
// geometry, input parsing, the NG shield-table Monte Carlo, the source model, ApplySource,
// every writer, the output audit and the run manifest -- is the original C++ code.
//
// Translation rules
// -----------------
// * Every expression of Liouville.cpp (EvaluateKernels), BVH.cpp (FirstHit, AnyHit,
//   RayTriangle, RayBox), Mesh.cpp (PointTriangleDistance), Gauges.cpp (NGWeight,
//   RGAcceptance, ShieldFraction) and Source.cpp (Evaluate) is written in the same
//   operand order.  FP64 throughout; compiled with --fmad=false (IEEE operations as written).
//   The reference binary is compiled by GCC with -O3 -march=native, which contracts some
//   multiply-adds into FMA (objdump: 40 FMAs in EvaluateKernels, 96 in the BVH routines),
//   so results are close but NOT bitwise; see the C11 validation policy (ROADMAP).
// * The BVH is the very cBVH that main() built (nodes and FaceOrder read through
//   access_shim), traversed with the same stack discipline (push left, push right, pop
//   right first), so the same leaves are tested in the same order.
// * Samples of a facet are enumerated exactly as in EvaluateKernels (i, j, first then second
//   sub-triangle).  To bound kernel run time (display watchdog) a facet's samples are
//   split into consecutive chunks; each chunk is summed in sample order and the chunk sums
//   of a facet are added in chunk order (deterministic).  A facet with <= CHUNK samples
//   (every facet beyond a few tens of metres) is summed exactly in the CPU order.
// * Supported: SurfaceBC flat|rough; NGAcceptance srccg_shield|hemisphere|full;
//   RGAcceptance hemisphere; FlightTimeRotation none (version-1 kernels).  ik_fov is
//   refused with a clear error (it needs the FOV polygons; not used by the campaign).
//
// Device: environment variable COPS_CUDA_DEVICE (default 0).  Input files are unchanged
// (an unknown input key would be an error in the reference parser).

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include <chrono>

#include <cuda_runtime.h>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_reduce.cuh>

#include "Liouville.h"
#include "Source.h"
#include "access_shim.h"

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "Error: CUDA %s at %s:%d: %s\n", #call, __FILE__, __LINE__, cudaGetErrorString(e_)); \
    exit(EXIT_FAILURE); } } while (0)

namespace {

const int CHUNK = 32;              // samples per work item
const int STACK = 64;              // BVH traversal stack (tree depth ~ 20 for 788k faces)
const long long ITEMS_PER_LAUNCH = 1LL << 18;
const double PI = 3.14159265358979323846;   // Const::Pi

// ------------------------------------------------------------------ device vector algebra (Vec3.h)
struct V3 { double x, y, z; };
__host__ __device__ __forceinline__ V3 v3(double a, double b, double c) { V3 r; r.x = a; r.y = b; r.z = c; return r; }
__device__ __forceinline__ V3 add(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
__device__ __forceinline__ V3 sub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
__device__ __forceinline__ V3 neg(V3 a) { return v3(-a.x, -a.y, -a.z); }
__device__ __forceinline__ V3 mul(V3 a, double s) { return v3(s * a.x, s * a.y, s * a.z); }   // Vec3*double: a*x
__device__ __forceinline__ V3 dvd(V3 a, double s) { return v3(a.x / s, a.y / s, a.z / s); }
__device__ __forceinline__ double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__device__ __forceinline__ V3 cross(V3 a, V3 b) {
  return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
__device__ __forceinline__ double len(V3 a) { return sqrt(dot(a, a)); }
__device__ __forceinline__ double comp(V3 a, int i) { return (i == 0) ? a.x : ((i == 1) ? a.y : a.z); }

// ------------------------------------------------------------------ device mesh and BVH
struct DMesh {
  const double *ax, *ay, *az, *bx, *by, *bz, *cx, *cy, *cz;   // vertices of each face
  const double *nx, *ny, *nz, *gx, *gy, *gz;                  // unit normal, centroid
  const double *area, *maxEdge;
  const double *bmin, *bmax;                                  // BVH node boxes (3 per node)
  const int *left, *right, *first, *count, *faceOrder;
};

__device__ __forceinline__ V3 A(const DMesh& m, int f) { return v3(m.ax[f], m.ay[f], m.az[f]); }
__device__ __forceinline__ V3 B(const DMesh& m, int f) { return v3(m.bx[f], m.by[f], m.bz[f]); }
__device__ __forceinline__ V3 C(const DMesh& m, int f) { return v3(m.cx[f], m.cy[f], m.cz[f]); }
__device__ __forceinline__ V3 N(const DMesh& m, int f) { return v3(m.nx[f], m.ny[f], m.nz[f]); }

// BVH.cpp RayTriangle (Moller-Trumbore)
__device__ bool rayTriangle(const DMesh& m, int f, V3 x0, V3 d, double tMin, double tMax, double& t) {
  V3 a = A(m, f), e1 = sub(B(m, f), a), e2 = sub(C(m, f), a);
  V3 p = cross(d, e2);
  double det = dot(e1, p);
  if (det == 0.0) return false;
  double inv = 1.0 / det;
  V3 s = sub(x0, a);
  double u = dot(s, p) * inv;
  if ((u < 0.0) || (u > 1.0)) return false;
  V3 q = cross(s, e1);
  double v = dot(d, q) * inv;
  if ((v < 0.0) || (u + v > 1.0)) return false;
  t = dot(e2, q) * inv;
  return (t > tMin) && (t < tMax);
}

// BVH.cpp RayBox
__device__ __forceinline__ bool rayBox(const DMesh& m, int node, V3 x0, V3 invd, double tMin, double tMax) {
  for (int i = 0; i < 3; i++) {
    double t0 = (m.bmin[3 * node + i] - comp(x0, i)) * comp(invd, i);
    double t1 = (m.bmax[3 * node + i] - comp(x0, i)) * comp(invd, i);
    if (comp(invd, i) < 0.0) { double tmp = t0; t0 = t1; t1 = tmp; }
    tMin = fmax(tMin, t0);
    tMax = fmin(tMax, t1);
    if (tMax < tMin) return false;
  }
  return true;
}

// BVH.cpp FirstHit; *overflow is set if the stack would overflow (never for tree depth < STACK-1)
__device__ int firstHit(const DMesh& m, V3 x0, V3 d, double tMin, double tMax, double& tHit, int ignore, int* overflow) {
  V3 invd = v3(1.0 / d.x, 1.0 / d.y, 1.0 / d.z);
  int stack[STACK], nstack = 0, res = -1;
  double best = tMax, t;
  stack[nstack++] = 0;
  while (nstack > 0) {
    int n = stack[--nstack];
    if (!rayBox(m, n, x0, invd, tMin, best)) continue;
    if (m.left[n] < 0) {
      for (int i = m.first[n]; i < m.first[n] + m.count[n]; i++) {
        int iface = m.faceOrder[i];
        if ((iface != ignore) && rayTriangle(m, iface, x0, d, tMin, best, t)) best = t, res = iface;
      }
    } else {
      if (nstack + 2 > STACK) { atomicExch(overflow, 1); return -2; }
      stack[nstack++] = m.left[n];
      stack[nstack++] = m.right[n];
    }
  }
  if (res != -1) tHit = best;
  return res;
}

// BVH.cpp AnyHit
__device__ bool anyHit(const DMesh& m, V3 x0, V3 d, double tMin, double tMax, int ignore, int* overflow) {
  V3 invd = v3(1.0 / d.x, 1.0 / d.y, 1.0 / d.z);
  int stack[STACK], nstack = 0;
  double t;
  stack[nstack++] = 0;
  while (nstack > 0) {
    int n = stack[--nstack];
    if (!rayBox(m, n, x0, invd, tMin, tMax)) continue;
    if (m.left[n] < 0) {
      for (int i = m.first[n]; i < m.first[n] + m.count[n]; i++) {
        int iface = m.faceOrder[i];
        if ((iface != ignore) && rayTriangle(m, iface, x0, d, tMin, tMax, t)) return true;
      }
    } else {
      if (nstack + 2 > STACK) { atomicExch(overflow, 1); return true; }
      stack[nstack++] = m.left[n];
      stack[nstack++] = m.right[n];
    }
  }
  return false;
}

// Mesh.cpp PointTriangleDistance
__device__ double pointTriangleDistance(V3 p, V3 a, V3 b, V3 c) {
  V3 ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
  double d1 = dot(ab, ap), d2 = dot(ac, ap);
  if ((d1 <= 0.0) && (d2 <= 0.0)) return len(sub(p, a));
  V3 bp = sub(p, b);
  double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if ((d3 >= 0.0) && (d4 <= d3)) return len(sub(p, b));
  double vc = d1 * d4 - d3 * d2;
  if ((vc <= 0.0) && (d1 >= 0.0) && (d3 <= 0.0)) {
    double v = d1 / (d1 - d3);
    return len(sub(p, add(a, mul(ab, v))));
  }
  V3 cp = sub(p, c);
  double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if ((d6 >= 0.0) && (d5 <= d6)) return len(sub(p, c));
  double vb = d5 * d2 - d1 * d6;
  if ((vb <= 0.0) && (d2 >= 0.0) && (d6 <= 0.0)) {
    double w = d2 / (d2 - d6);
    return len(sub(p, add(a, mul(ac, w))));
  }
  double va = d3 * d6 - d5 * d4;
  if ((va <= 0.0) && ((d4 - d3) >= 0.0) && ((d5 - d6) >= 0.0)) {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return len(sub(p, add(b, mul(sub(c, b), w))));
  }
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom, w = vc * denom;
  return len(sub(p, add(add(a, mul(ab, v)), mul(ac, w))));
}

// ------------------------------------------------------------------ gauges (Gauges.cpp)
enum { NG_SHIELD = 0, NG_HEMISPHERE = 1, NG_FULL = 2 };

struct DGauge {
  V3 ngB, rgB;            // NG and RG boresights (body frame)
  int ngModel, headOn, flat;
  const double* shield;   // cGaugeModel::ShieldTable
  int nShield;
};

__device__ __forceinline__ double shieldFraction(const DGauge& g, double mu) {
  double x = fmax(0.0, fmin(1.0, mu)) * (g.nShield - 1);
  int i = min((int)x, g.nShield - 2);
  double w = x - i;
  return (1.0 - w) * g.shield[i] + w * g.shield[i + 1];
}

__device__ __forceinline__ double ngWeight(const DGauge& g, V3 l) {
  double mu = dot(l, g.ngB);
  double acceptance;
  if (g.ngModel == NG_FULL) acceptance = 1.0;
  else if (g.ngModel == NG_HEMISPHERE) acceptance = (mu < 0.0) ? 1.0 : 0.0;
  else {
    if (mu < 0.0) acceptance = 1.0;
    else acceptance = (dot(l, g.rgB) < 0.0) ? shieldFraction(g, mu) : 0.0;
  }
  if (g.headOn && (mu < 0.0)) acceptance *= 1.0 - mu;
  return acceptance;
}

__device__ __forceinline__ double rgAcceptance(const DGauge& g, V3 l) {   // hemisphere
  if (dot(l, g.rgB) >= 0.0) return 0.0;
  return 1.0;
}

// ------------------------------------------------------------------ kernels
// Pass 1 per face: distance (all faces), front test, subdivision level, sample count.
__global__ void k_pass1(DMesh m, int nFaces, V3 xg, double accuracy, int maxLevel,
                        double* dist, int* level, long long* nsamples, long long* itemCount) {
  int f = blockIdx.x * blockDim.x + threadIdx.x;
  if (f >= nFaces) return;
  V3 a = A(m, f), b = B(m, f), c = C(m, f);
  double dd = pointTriangleDistance(xg, a, b, c);
  dist[f] = dd;
  if (dot(N(m, f), sub(xg, a)) <= 0.0) { level[f] = -1; nsamples[f] = 0; itemCount[f] = 0; return; }
  int lev = 0;
  double h = m.maxEdge[f];
  while ((lev < maxLevel) && (h > accuracy * dd)) h *= 0.5, lev++;
  long long n = 1LL << lev;
  level[f] = lev;
  nsamples[f] = n * n;
  itemCount[f] = (n * n + CHUNK - 1) / CHUNK;
}

// Work items of each face (consecutive CHUNK-sample ranges), at the face's scanned offset.
__global__ void k_fill_items(int nFaces, const long long* faceItem0, const long long* nsamples,
                             int* itemFace, long long* itemS0, long long* itemS1) {
  int f = blockIdx.x * blockDim.x + threadIdx.x;
  if (f >= nFaces) return;
  long long base = faceItem0[f], ns = nsamples[f];
  for (long long k = 0; base + k < faceItem0[f + 1]; k++) {
    itemFace[base + k] = f;
    itemS0[base + k] = k * CHUNK;
    itemS1[base + k] = min(ns, (k + 1) * CHUNK);
  }
}

// Pass 2 per work item (face, [s0, s1)): partial Kn, KPhi, KV, KnBelow in sample order.
__global__ void k_pass2(DMesh m, DGauge g, V3 xg, long long nItems, const int* itemFace,
                        const long long* itemS0, const long long* itemS1, const int* level,
                        double* pKn, double* pKphi, double* pKv, double* pKnb, int* overflow) {
  long long it = (long long)blockIdx.x * blockDim.x + threadIdx.x;
  if (it >= nItems) return;
  const int iface = itemFace[it];
  const int n = 1 << level[iface];
  V3 a = A(m, iface), b = B(m, iface), c = C(m, iface), nrm = N(m, iface);
  double subArea = m.area[iface] / ((double)n * n);
  V3 ea = dvd(sub(b, a), (double)n), eb = dvd(sub(c, a), (double)n);
  double kn = 0.0, kphi = 0.0, kv = 0.0, knb = 0.0;
  // Locate (i, j, second) of sample s0: row i starts at i*(2n-i); within the row the
  // sequence is (j,first),(j,second),... and the last j has only its first sample.
  long long s = itemS0[it], s1 = itemS1[it];
  int i = 0;
  while ((long long)(i + 1) * (2 * n - (i + 1)) <= s) i++;
  long long o = s - (long long)i * (2 * n - i);
  int j = (int)(o / 2), second = (int)(o % 2);
  for (; s < s1; s++) {
    V3 p = second ? add(add(a, mul(ea, i + 2.0 / 3.0)), mul(eb, j + 2.0 / 3.0))
                  : add(add(a, mul(ea, i + 1.0 / 3.0)), mul(eb, j + 1.0 / 3.0));
    // ---- the Sample lambda of cLiouville::EvaluateKernels
    V3 d = sub(p, xg);
    double t;
    if (firstHit(m, xg, d, 0.0, 1.0 + 1.0E-9, t, -1, overflow) == iface) {
      double r2 = dot(d, d), r = sqrt(r2);
      V3 l = dvd(neg(d), r);
      double cosTheta = dot(l, nrm);
      if (cosTheta > 0.0) {
        double w = subArea / r2 * (g.flat ? cosTheta / PI : 1.0 / (2.0 * PI));
        double wng = w * ngWeight(g, l);
        kn += wng;
        if (dot(l, g.rgB) > 0.0) knb += wng;     // Zsc = RG boresight
        double acc = rgAcceptance(g, l);
        if (acc > 0.0) {
          kphi += w * acc * (-dot(l, g.rgB));
          kv += w * acc;
        }
      }
    }
    // ---- next sample in the CPU loop order
    if (!second && (i + j < n - 1)) second = 1;
    else { second = 0; j++; if (i + j >= n) { i++; j = 0; } }
  }
  pKn[it] = kn; pKphi[it] = kphi; pKv[it] = kv; pKnb[it] = knb;
}

// Pass 3 per face: add its chunk partials in order (faces without samples get zeros).
__global__ void k_pass3(int nFaces, const long long* faceItem0, const double* pKn, const double* pKphi,
                        const double* pKv, const double* pKnb, double* Kn, double* Kphi, double* Kv, double* Knb,
                        int* nonzero) {
  int f = blockIdx.x * blockDim.x + threadIdx.x;
  if (f >= nFaces) return;
  double kn = 0.0, kphi = 0.0, kv = 0.0, knb = 0.0;
  for (long long k = faceItem0[f]; k < faceItem0[f + 1]; k++) {
    kn = kn + pKn[k]; kphi = kphi + pKphi[k]; kv = kv + pKv[k]; knb = knb + pKnb[k];
  }
  Kn[f] = kn; Kphi[f] = kphi; Kv[f] = kv; Knb[f] = knb;
  nonzero[f] = ((kn != 0.0) || (kphi != 0.0) || (kv != 0.0)) ? 1 : 0;   // EvaluateKernels' record test
}

// Records of the nonzero faces, in face order (scanned positions).
__global__ void k_compact(int nFaces, const int* nonzero, const int* pos, const double* Kn, const double* Kphi,
                          const double* Kv, const double* Knb, int* cFace, double* cKn, double* cKphi,
                          double* cKv, double* cKnb) {
  int f = blockIdx.x * blockDim.x + threadIdx.x;
  if (f >= nFaces || !nonzero[f]) return;
  int k = pos[f];
  cFace[k] = f; cKn[k] = Kn[f]; cKphi[k] = Kphi[f]; cKv[k] = Kv[f]; cKnb[k] = Knb[f];
}

// Source.cpp Evaluate: cast shadow of the listed (sunlit) faces.
__global__ void k_shadow(DMesh m, int nFaces, const unsigned char* lit, V3 sun, unsigned char* shadow, int* overflow) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nFaces) return;
  if (!lit[i]) { shadow[i] = 0; return; }
  V3 x0 = add(v3(m.gx[i], m.gy[i], m.gz[i]), mul(N(m, i), 1.0E-6 * m.maxEdge[i]));
  shadow[i] = anyHit(m, x0, sun, 0.0, 1.0E12, i, overflow) ? 1 : 0;
}

// Self-check kernels (run once): FirstHit and NGWeight against the CPU classes.
__global__ void k_test_hits(DMesh m, int nRays, const double* o, const double* d, int* face, double* t, int* overflow) {
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= nRays) return;
  double th = -1.0;
  face[k] = firstHit(m, v3(o[3 * k], o[3 * k + 1], o[3 * k + 2]), v3(d[3 * k], d[3 * k + 1], d[3 * k + 2]),
                     0.0, 1.0E30, th, -1, overflow);
  t[k] = th;
}

__global__ void k_test_ng(DGauge g, int n, const double* l, double* out) {
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;
  out[k] = ngWeight(g, v3(l[3 * k], l[3 * k + 1], l[3 * k + 2]));
}

// ------------------------------------------------------------------ host context
template <class T> T* upload(const std::vector<T>& v) {
  T* p = nullptr;
  CK(cudaMalloc(&p, std::max<size_t>(1, v.size()) * sizeof(T)));
  if (!v.empty()) CK(cudaMemcpy(p, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice));
  return p;
}

template <class T> struct DBuf {
  T* p = nullptr; size_t n = 0;
  void ensure(size_t want) {
    if (want <= n) return;
    if (p) CK(cudaFree(p));
    n = std::max<size_t>(want, 1);
    CK(cudaMalloc(&p, n * sizeof(T)));
  }
};

struct Context {
  const cMesh* mesh = nullptr; const cBVH* bvh = nullptr;
  const void* nodesData = nullptr; const void* facesData = nullptr; size_t nFaces = 0;
  DMesh dm{};
  std::vector<void*> owned;
  int* overflow = nullptr;
  bool announced = false;
  // per-point buffers
  DBuf<double> dist, Kn, Kphi, Kv, Knb, pKn, pKphi, pKv, pKnb, shield;
  DBuf<int> level, itemFace, nonzero, pos, cFace; DBuf<long long> nsamp, itemCount, itemS0, itemS1, faceItem0;
  DBuf<double> cKn, cKphi, cKv, cKnb, dmin;
  DBuf<unsigned char> lit, shadow, scratch;
  const cGaugeModel* gauges = nullptr;
  bool gaugesChecked = false;
} ctx;

int Device() {
  const char* s = getenv("COPS_CUDA_DEVICE");
  return s ? atoi(s) : 0;
}

void CheckOverflow() {
  int ov = 0;
  CK(cudaMemcpy(&ov, ctx.overflow, sizeof(int), cudaMemcpyDeviceToHost));
  if (ov) { fprintf(stderr, "Error: GPU BVH traversal stack overflow (tree deeper than %d)\n", STACK - 2); exit(EXIT_FAILURE); }
}

void Bind(const cMesh& mesh, const cBVH& bvh) {
  if (ctx.mesh == &mesh && ctx.bvh == &bvh && ctx.nodesData == mesh.Nodes.data() &&
      ctx.facesData == mesh.Faces.data() && ctx.nFaces == mesh.Faces.size()) return;
  CK(cudaSetDevice(Device()));
  for (void* p : ctx.owned) CK(cudaFree(p));
  ctx.owned.clear();
  const size_t nf = mesh.Faces.size();
  std::vector<double> ax(nf), ay(nf), az(nf), bx(nf), by(nf), bz(nf), cx(nf), cy(nf), cz(nf),
      nx(nf), ny(nf), nz(nf), gx(nf), gy(nf), gz(nf), area(nf), maxEdge(nf);
  for (size_t i = 0; i < nf; i++) {
    const cTriangle& f = mesh.Faces[i];
    const Vec3 &a = mesh.Nodes[f.node[0]], &b = mesh.Nodes[f.node[1]], &c = mesh.Nodes[f.node[2]];
    ax[i] = a.x; ay[i] = a.y; az[i] = a.z; bx[i] = b.x; by[i] = b.y; bz[i] = b.z;
    cx[i] = c.x; cy[i] = c.y; cz[i] = c.z;
    nx[i] = f.Normal.x; ny[i] = f.Normal.y; nz[i] = f.Normal.z;
    gx[i] = f.Centroid.x; gy[i] = f.Centroid.y; gz[i] = f.Centroid.z;
    area[i] = f.Area; maxEdge[i] = f.MaxEdge;
  }
  std::vector<ShimNode> nodes;
  std::vector<int> order;
  ShimBVH(bvh, nodes, order);
  if (order.size() != nf || nodes.empty()) { fprintf(stderr, "Error: BVH does not match the mesh\n"); exit(EXIT_FAILURE); }
  std::vector<double> bmin(3 * nodes.size()), bmax(3 * nodes.size());
  std::vector<int> left(nodes.size()), right(nodes.size()), first(nodes.size()), count(nodes.size());
  for (size_t k = 0; k < nodes.size(); k++) {
    for (int i = 0; i < 3; i++) bmin[3 * k + i] = nodes[k].bmin[i], bmax[3 * k + i] = nodes[k].bmax[i];
    left[k] = nodes[k].left; right[k] = nodes[k].right; first[k] = nodes[k].first; count[k] = nodes[k].count;
  }
  auto keep = [&](auto* p) { ctx.owned.push_back((void*)p); return p; };
  DMesh& m = ctx.dm;
  m.ax = keep(upload(ax)); m.ay = keep(upload(ay)); m.az = keep(upload(az));
  m.bx = keep(upload(bx)); m.by = keep(upload(by)); m.bz = keep(upload(bz));
  m.cx = keep(upload(cx)); m.cy = keep(upload(cy)); m.cz = keep(upload(cz));
  m.nx = keep(upload(nx)); m.ny = keep(upload(ny)); m.nz = keep(upload(nz));
  m.gx = keep(upload(gx)); m.gy = keep(upload(gy)); m.gz = keep(upload(gz));
  m.area = keep(upload(area)); m.maxEdge = keep(upload(maxEdge));
  m.bmin = keep(upload(bmin)); m.bmax = keep(upload(bmax));
  m.left = keep(upload(left)); m.right = keep(upload(right)); m.first = keep(upload(first));
  m.count = keep(upload(count)); m.faceOrder = keep(upload(order));
  if (!ctx.overflow) CK(cudaMalloc(&ctx.overflow, sizeof(int)));
  CK(cudaMemset(ctx.overflow, 0, sizeof(int)));
  ctx.mesh = &mesh; ctx.bvh = &bvh; ctx.nodesData = mesh.Nodes.data(); ctx.facesData = mesh.Faces.data();
  ctx.nFaces = nf;
  if (!ctx.announced) {
    cudaDeviceProp prop;
    CK(cudaGetDeviceProperties(&prop, Device()));
    printf("liouville_cops_cuda: GPU kernels on device %d (%s); BVH %zu nodes; FP64, --fmad=false\n",
           Device(), prop.name, nodes.size());
    ctx.announced = true;
  }

  // Self-check: device FirstHit against cBVH::FirstHit (the public CPU method) on random rays.
  const int nRays = 4096;
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  Vec3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
  for (const Vec3& v : mesh.Nodes) lo = Min(lo, v), hi = Max(hi, v);
  Vec3 mid = (lo + hi) * 0.5, half = (hi - lo) * 0.75;
  std::vector<double> o(3 * nRays), d(3 * nRays);
  for (int k = 0; k < nRays; k++) {
    Vec3 x(mid.x + half.x * U(rng), mid.y + half.y * U(rng), mid.z + half.z * U(rng));
    Vec3 dir = Normalized(Vec3(U(rng), U(rng), U(rng)));
    o[3 * k] = x.x; o[3 * k + 1] = x.y; o[3 * k + 2] = x.z; d[3 * k] = dir.x; d[3 * k + 1] = dir.y; d[3 * k + 2] = dir.z;
  }
  double* dO = upload(o); double* dD = upload(d);
  int* dF = nullptr; double* dT = nullptr;
  CK(cudaMalloc(&dF, nRays * sizeof(int))); CK(cudaMalloc(&dT, nRays * sizeof(double)));
  k_test_hits<<<(nRays + 127) / 128, 128>>>(ctx.dm, nRays, dO, dD, dF, dT, ctx.overflow);
  CK(cudaGetLastError()); CK(cudaDeviceSynchronize()); CheckOverflow();
  std::vector<int> gf(nRays); std::vector<double> gt(nRays);
  CK(cudaMemcpy(gf.data(), dF, nRays * sizeof(int), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(gt.data(), dT, nRays * sizeof(double), cudaMemcpyDeviceToHost));
  CK(cudaFree(dO)); CK(cudaFree(dD)); CK(cudaFree(dF)); CK(cudaFree(dT));
  int mismatch = 0, hits = 0;
  double maxRel = 0.0;
  for (int k = 0; k < nRays; k++) {
    double t = -1.0;
    int f = bvh.FirstHit(Vec3(o[3 * k], o[3 * k + 1], o[3 * k + 2]), Vec3(d[3 * k], d[3 * k + 1], d[3 * k + 2]), 0.0, 1.0E30, t);
    if (f != gf[k]) mismatch++;
    else if (f >= 0) { hits++; maxRel = std::max(maxRel, std::fabs(gt[k] - t) / std::max(std::fabs(t), 1e-300)); }
  }
  printf("liouville_cops_cuda: self-check FirstHit GPU vs CPU: %d rays, %d hits, %d face mismatches, max rel |dt| %.2e\n",
         nRays, hits, mismatch, maxRel);
  if (mismatch > nRays / 1000 || maxRel > 1e-9) {
    fprintf(stderr, "Error: GPU BVH traversal disagrees with cBVH::FirstHit (layout or translation error)\n");
    exit(EXIT_FAILURE);
  }
}

DGauge BindGauges(const cGaugeModel& gauges, const cObservationGeometry& g, bool flat) {
  if (gauges.RGAcceptanceModel != "hemisphere" ||
      !(gauges.NGAcceptanceModel == "srccg_shield" || gauges.NGAcceptanceModel == "hemisphere" ||
        gauges.NGAcceptanceModel == "full")) {
    fprintf(stderr, "Error: liouville_cops_cuda supports NGAcceptance=srccg_shield|hemisphere|full and "
                    "RGAcceptance=hemisphere (got %s/%s); use liouville_cops\n",
            gauges.NGAcceptanceModel.c_str(), gauges.RGAcceptanceModel.c_str());
    exit(EXIT_FAILURE);
  }
  if (ctx.gauges != &gauges) {
    const std::vector<double>& table = ShimShieldTable(gauges);
    ctx.shield.ensure(table.size());
    if (!table.empty()) CK(cudaMemcpy(ctx.shield.p, table.data(), table.size() * sizeof(double), cudaMemcpyHostToDevice));
    ctx.gauges = &gauges;
    ctx.gaugesChecked = false;
  }
  DGauge d;
  d.ngB = v3(g.NG.Boresight.x, g.NG.Boresight.y, g.NG.Boresight.z);
  d.rgB = v3(g.RG.Boresight.x, g.RG.Boresight.y, g.RG.Boresight.z);
  d.ngModel = (gauges.NGAcceptanceModel == "srccg_shield") ? NG_SHIELD :
              ((gauges.NGAcceptanceModel == "hemisphere") ? NG_HEMISPHERE : NG_FULL);
  d.headOn = gauges.NGHeadOnWeight ? 1 : 0;
  d.flat = flat ? 1 : 0;
  d.shield = ctx.shield.p;
  d.nShield = (int)ShimShieldTable(gauges).size();
  if (d.ngModel == NG_SHIELD && d.nShield < 2) { fprintf(stderr, "Error: NG shield table missing\n"); exit(EXIT_FAILURE); }
  if (!ctx.gaugesChecked) {
    // Self-check: device NGWeight against cGaugeModel::NGWeight for random directions.
    const int n = 20000;
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    std::vector<double> l(3 * n), out(n);
    for (int k = 0; k < n; k++) {
      Vec3 v = Normalized(Vec3(U(rng), U(rng), U(rng)));
      l[3 * k] = v.x; l[3 * k + 1] = v.y; l[3 * k + 2] = v.z;
    }
    double* dl = upload(l); double* dout = nullptr;
    CK(cudaMalloc(&dout, n * sizeof(double)));
    k_test_ng<<<(n + 127) / 128, 128>>>(d, n, dl, dout);
    CK(cudaGetLastError());
    CK(cudaMemcpy(out.data(), dout, n * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaFree(dl)); CK(cudaFree(dout));
    double maxAbs = 0.0;
    for (int k = 0; k < n; k++)
      maxAbs = std::max(maxAbs, std::fabs(out[k] - gauges.NGWeight(Vec3(l[3 * k], l[3 * k + 1], l[3 * k + 2]), g)));
    printf("liouville_cops_cuda: self-check NGWeight GPU vs CPU: %d directions, max |diff| %.2e\n", n, maxAbs);
    if (maxAbs > 1e-12) { fprintf(stderr, "Error: GPU NG weight disagrees with cGaugeModel::NGWeight\n"); exit(EXIT_FAILURE); }
    ctx.gaugesChecked = true;
  }
  return d;
}

unsigned Blocks(long long n, int threads) { return (unsigned)((n + threads - 1) / threads); }

// Optional stage timing (COPS_CUDA_TIMING=1): cumulative seconds, printed every 100 points.
struct Timing {
  bool on = getenv("COPS_CUDA_TIMING") != nullptr;
  double kernel_pass1 = 0, host_items = 0, kernel_pass2 = 0, pass3_copy = 0, host_tail = 0,
         shadow = 0, source_host = 0, total_evaluate = 0;
  long long calls = 0, items = 0;
  std::chrono::steady_clock::time_point t;
  void start() { if (on) { cudaDeviceSynchronize(); t = std::chrono::steady_clock::now(); } }
  void lap(double& acc) {
    if (!on) return;
    cudaDeviceSynchronize();
    auto now = std::chrono::steady_clock::now();
    acc += std::chrono::duration<double>(now - t).count();
    t = now;
  }
  void report() {
    if (!on || calls % 100) return;
    printf("liouville_cops_cuda timing after %lld points: pass1 %.2f s, items %.2f s, pass2 %.2f s, "
           "pass3+copy %.2f s, tail %.2f s, shadow %.2f s, source host %.2f s; work items %lld\n",
           calls, kernel_pass1, host_items, kernel_pass2, pass3_copy, host_tail, shadow, source_host, items);
  }
} timing;

}  // namespace

// =====================================================================================================
namespace {
// CUB device scan / reduction with a reusable scratch buffer.
template <class In, class Out> void ExclusiveSum(In in, Out out, int n) {
  size_t bytes = 0;
  CK(cub::DeviceScan::ExclusiveSum(nullptr, bytes, in, out, n));
  ctx.scratch.ensure(bytes);
  CK(cub::DeviceScan::ExclusiveSum(ctx.scratch.p, bytes, in, out, n));
}
double DeviceMin(const double* in, int n) {
  size_t bytes = 0;
  ctx.dmin.ensure(1);
  CK(cub::DeviceReduce::Min(nullptr, bytes, in, ctx.dmin.p, n));
  ctx.scratch.ensure(bytes);
  CK(cub::DeviceReduce::Min(ctx.scratch.p, bytes, in, ctx.dmin.p, n));
  double v;
  CK(cudaMemcpy(&v, ctx.dmin.p, sizeof(double), cudaMemcpyDeviceToHost));
  return v;
}
}  // namespace

void cLiouville::EvaluateKernels(const cMesh& mesh, const cBVH& bvh, const cGaugeModel& gauges,
                                 const cObservationGeometry& g, cLiouvilleResult& res) const {
  const int nFaces = (int)mesh.Faces.size();
  const Vec3 xg = g.xSC;
  Bind(mesh, bvh);
  DGauge dg = BindGauges(gauges, g, SurfaceBC == "flat");

  res.GaugeInsideNucleus = PointInsideSurface(mesh, bvh, xg);       // CPU (one ray)
  res.VdotB_RG = SpacecraftVelocity ? Dot(g.vSC, g.RG.Boresight) : 0.0;

  timing.start();
  ctx.dist.ensure(nFaces); ctx.level.ensure(nFaces); ctx.nsamp.ensure(nFaces);
  ctx.itemCount.ensure(nFaces + 1); ctx.faceItem0.ensure(nFaces + 1);
  V3 x = {xg.x, xg.y, xg.z};
  CK(cudaMemset(ctx.itemCount.p + nFaces, 0, sizeof(long long)));
  k_pass1<<<Blocks(nFaces, 256), 256>>>(ctx.dm, nFaces, x, SubdivisionAccuracy, MaxSubdivisionLevel,
                                        ctx.dist.p, ctx.level.p, ctx.nsamp.p, ctx.itemCount.p);
  CK(cudaGetLastError());
  // res.Altitude = min over all faces of the point-triangle distance (min is order independent: exact).
  const double altitude = DeviceMin(ctx.dist.p, nFaces);
  timing.lap(timing.kernel_pass1);

  // Work items: consecutive CHUNK-sample ranges per face, faces in index order (device scan).
  ExclusiveSum(ctx.itemCount.p, ctx.faceItem0.p, nFaces + 1);
  long long nItems = 0;
  CK(cudaMemcpy(&nItems, ctx.faceItem0.p + nFaces, sizeof(long long), cudaMemcpyDeviceToHost));
  ctx.itemFace.ensure(nItems); ctx.itemS0.ensure(nItems); ctx.itemS1.ensure(nItems);
  ctx.pKn.ensure(nItems); ctx.pKphi.ensure(nItems); ctx.pKv.ensure(nItems); ctx.pKnb.ensure(nItems);
  k_fill_items<<<Blocks(nFaces, 256), 256>>>(nFaces, ctx.faceItem0.p, ctx.nsamp.p, ctx.itemFace.p,
                                             ctx.itemS0.p, ctx.itemS1.p);
  CK(cudaGetLastError());
  timing.lap(timing.host_items);
  timing.items += nItems;
  for (long long start = 0; start < nItems; start += ITEMS_PER_LAUNCH) {
    long long count = std::min(ITEMS_PER_LAUNCH, nItems - start);
    k_pass2<<<Blocks(count, 128), 128>>>(ctx.dm, dg, x, count, ctx.itemFace.p + start, ctx.itemS0.p + start,
                                         ctx.itemS1.p + start, ctx.level.p, ctx.pKn.p + start,
                                         ctx.pKphi.p + start, ctx.pKv.p + start, ctx.pKnb.p + start, ctx.overflow);
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());
  }
  CheckOverflow();
  timing.lap(timing.kernel_pass2);

  // Per-face sums, then the records of the nonzero faces in face order (device compaction).
  ctx.Kn.ensure(nFaces); ctx.Kphi.ensure(nFaces); ctx.Kv.ensure(nFaces); ctx.Knb.ensure(nFaces);
  ctx.nonzero.ensure(nFaces + 1); ctx.pos.ensure(nFaces + 1);
  CK(cudaMemset(ctx.nonzero.p + nFaces, 0, sizeof(int)));
  k_pass3<<<Blocks(nFaces, 256), 256>>>(nFaces, ctx.faceItem0.p, ctx.pKn.p, ctx.pKphi.p, ctx.pKv.p, ctx.pKnb.p,
                                        ctx.Kn.p, ctx.Kphi.p, ctx.Kv.p, ctx.Knb.p, ctx.nonzero.p);
  CK(cudaGetLastError());
  ExclusiveSum(ctx.nonzero.p, ctx.pos.p, nFaces + 1);
  int nRec = 0;
  CK(cudaMemcpy(&nRec, ctx.pos.p + nFaces, sizeof(int), cudaMemcpyDeviceToHost));
  ctx.cFace.ensure(nRec); ctx.cKn.ensure(nRec); ctx.cKphi.ensure(nRec); ctx.cKv.ensure(nRec); ctx.cKnb.ensure(nRec);
  k_compact<<<Blocks(nFaces, 256), 256>>>(nFaces, ctx.nonzero.p, ctx.pos.p, ctx.Kn.p, ctx.Kphi.p, ctx.Kv.p, ctx.Knb.p,
                                          ctx.cFace.p, ctx.cKn.p, ctx.cKphi.p, ctx.cKv.p, ctx.cKnb.p);
  CK(cudaGetLastError());
  std::vector<int> face(nRec);
  std::vector<double> Kn(nRec), KPhi(nRec), KV(nRec), KnBelow(nRec);
  if (nRec) {
    CK(cudaMemcpy(face.data(), ctx.cFace.p, nRec * sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(Kn.data(), ctx.cKn.p, nRec * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(KPhi.data(), ctx.cKphi.p, nRec * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(KV.data(), ctx.cKv.p, nRec * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(KnBelow.data(), ctx.cKnb.p, nRec * sizeof(double), cudaMemcpyDeviceToHost));
  }
  timing.lap(timing.pass3_copy);

  // The tail of cLiouville::EvaluateKernels: altitude, records of the nonzero faces in face order.
  res.Altitude = altitude;
  res.Kernels.resize(nRec);
  for (int k = 0; k < nRec; k++) res.Kernels[k] = {face[k], Kn[k], KPhi[k], KV[k], KnBelow[k]};
  res.nVisibleElements = (int)res.Kernels.size();
  timing.lap(timing.host_tail);
  timing.calls++;
  timing.report();
}

// =====================================================================================================
void cSourceModel::Evaluate(const cMesh& mesh, const cBVH& bvh, const Vec3& SunDirection, double HeliocentricDistance,
                            std::vector<cSurfaceState>& state) const {
  timing.start();
  state.resize(mesh.Faces.size());
  const int nFaces = (int)mesh.Faces.size();
  std::vector<unsigned char> lit(nFaces), shadow(nFaces, 0);
  bool anyLit = false;
  #pragma omp parallel for schedule(static) reduction(||:anyLit)
  for (int i = 0; i < nFaces; i++) {
    cSurfaceState& st = state[i];
    st.CosSZA = Dot(mesh.Faces[i].Normal, SunDirection);      // as Source.cpp (host)
    st.Shadow = false;
    lit[i] = ((st.CosSZA > 0.0) && CastShadows) ? 1 : 0;
    anyLit = anyLit || lit[i];
  }
  if (anyLit) {
    Bind(mesh, bvh);
    ctx.lit.ensure(nFaces); ctx.shadow.ensure(nFaces);
    CK(cudaMemcpy(ctx.lit.p, lit.data(), nFaces, cudaMemcpyHostToDevice));
    V3 sun = {SunDirection.x, SunDirection.y, SunDirection.z};
    k_shadow<<<Blocks(nFaces, 128), 128>>>(ctx.dm, nFaces, ctx.lit.p, sun, ctx.shadow.p, ctx.overflow);
    CK(cudaGetLastError());
    CK(cudaMemcpy(shadow.data(), ctx.shadow.p, nFaces, cudaMemcpyDeviceToHost));
    CheckOverflow();
  }
  timing.lap(timing.shadow);
  #pragma omp parallel for schedule(dynamic,256)
  for (int i = 0; i < nFaces; i++) {
    cSurfaceState& st = state[i];
    st.Shadow = shadow[i] != 0;
    EvaluateElement(i, st.CosSZA, (st.CosSZA <= 0.0) || st.Shadow, HeliocentricDistance, st.Flux, st.Temperature);
  }
  timing.lap(timing.source_host);
}
