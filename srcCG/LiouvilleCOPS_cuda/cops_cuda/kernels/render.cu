// Z-buffer rasterization of a triangle mesh for figures (cops_cuda/figures.py).
//
// This kernel draws pictures; it computes no model quantity.  The caller projects every vertex
// to image coordinates (u, v) in pixels and a depth w (smaller = nearer the viewer, w >= 0):
//   * orthographic views: u, v along two screen axes, w along the viewing direction;
//   * longitude/latitude maps: u = longitude, v = latitude, w = R_max - |r| (the outermost surface
//     wins where overhangs put several surface layers at one longitude/latitude).
// One thread per triangle loops over the pixel centres (i + 0.5, j + 0.5) of its bounding box and
// keeps those inside the triangle (edge functions, both windings accepted).  The nearest triangle
// of each pixel wins through atomicMin on a 64-bit key:
//   key = (IEEE bits of (float)w) << 32 | face index.
// For w >= 0 the float bit pattern increases with w, so the minimum key is the nearest face; ties
// in w go to the lower face index.  Integer atomics are order independent, so the image does not
// depend on thread scheduling (kernels/README.md, determinism).

extern "C" __global__ void raster_faces(
    const double* __restrict__ u, const double* __restrict__ v, const double* __restrict__ w,
    const long long* __restrict__ faces, const unsigned char* __restrict__ draw,
    long long nfaces, int width, int height, int max_box,
    unsigned long long* __restrict__ zbuf)
{
    for (long long f = blockIdx.x * (long long)blockDim.x + threadIdx.x; f < nfaces;
         f += (long long)gridDim.x * blockDim.x) {
        if (!draw[f]) continue;
        const long long a = faces[3 * f], b = faces[3 * f + 1], c = faces[3 * f + 2];
        const double ua = u[a], va = v[a], ub = u[b], vb = v[b], uc = u[c], vc = v[c];
        const double area = (ub - ua) * (vc - va) - (uc - ua) * (vb - va);
        if (area == 0.0) continue;                                  // degenerate in projection
        int i0 = (int)floor(fmin(ua, fmin(ub, uc)) - 0.5), i1 = (int)ceil(fmax(ua, fmax(ub, uc)) - 0.5);
        int j0 = (int)floor(fmin(va, fmin(vb, vc)) - 0.5), j1 = (int)ceil(fmax(va, fmax(vb, vc)) - 0.5);
        if (i1 < 0 || j1 < 0 || i0 >= width || j0 >= height) continue;
        i0 = max(i0, 0); j0 = max(j0, 0); i1 = min(i1, width - 1); j1 = min(j1, height - 1);
        // A face covering a huge box is a projection artefact (e.g. a longitude seam); skip it.
        if (i1 - i0 > max_box || j1 - j0 > max_box) continue;
        const double inv = 1.0 / area;
        for (int j = j0; j <= j1; ++j) {
            const double py = j + 0.5;
            for (int i = i0; i <= i1; ++i) {
                const double px = i + 0.5;
                // barycentric weights of a, b, c (sign-normalized by the signed area)
                const double la = ((ub - px) * (vc - py) - (uc - px) * (vb - py)) * inv;
                const double lb = ((uc - px) * (va - py) - (ua - px) * (vc - py)) * inv;
                const double lc = 1.0 - la - lb;
                if (la < 0.0 || lb < 0.0 || lc < 0.0) continue;
                const double depth = la * w[a] + lb * w[b] + lc * w[c];
                const unsigned int bits = __float_as_uint((float)fmax(depth, 0.0));
                const unsigned long long key = ((unsigned long long)bits << 32) | (unsigned long long)(f & 0xffffffffLL);
                atomicMin(&zbuf[(long long)j * width + i], key);
            }
        }
    }
}

// First intersection of rays with the mesh (figures: the sub-spacecraft intercept point).
// Ray r starts at orig[r] and runs along dir[r]; every (ray, face) pair is tested with the
// Moller-Trumbore algorithm in FP64.  Grid: x over faces (grid-stride), y over rays.  The nearest
// hit (smallest t > 0) wins through atomicMin on (IEEE bits of (float)t) << 32 | face; the caller
// recomputes the exact intersection with the winning face in FP64.  Rays that hit nothing keep
// the initial key 0xffffffffffffffff.
extern "C" __global__ void ray_first_hit(
    const double* __restrict__ vert, const long long* __restrict__ faces, long long nfaces,
    const double* __restrict__ orig, const double* __restrict__ dir,
    unsigned long long* __restrict__ best)
{
    const int r = blockIdx.y;
    const double ox = orig[3 * r], oy = orig[3 * r + 1], oz = orig[3 * r + 2];
    const double dx = dir[3 * r], dy = dir[3 * r + 1], dz = dir[3 * r + 2];
    for (long long f = blockIdx.x * (long long)blockDim.x + threadIdx.x; f < nfaces;
         f += (long long)gridDim.x * blockDim.x) {
        const double* A = vert + 3 * faces[3 * f];
        const double* B = vert + 3 * faces[3 * f + 1];
        const double* C = vert + 3 * faces[3 * f + 2];
        const double e1x = B[0] - A[0], e1y = B[1] - A[1], e1z = B[2] - A[2];
        const double e2x = C[0] - A[0], e2y = C[1] - A[1], e2z = C[2] - A[2];
        const double hx = dy * e2z - dz * e2y, hy = dz * e2x - dx * e2z, hz = dx * e2y - dy * e2x;
        const double det = e1x * hx + e1y * hy + e1z * hz;
        if (fabs(det) < 1e-30) continue;
        const double inv = 1.0 / det;
        const double sx = ox - A[0], sy = oy - A[1], sz = oz - A[2];
        const double u = (sx * hx + sy * hy + sz * hz) * inv;
        if (u < 0.0 || u > 1.0) continue;
        const double qx = sy * e1z - sz * e1y, qy = sz * e1x - sx * e1z, qz = sx * e1y - sy * e1x;
        const double v = (dx * qx + dy * qy + dz * qz) * inv;
        if (v < 0.0 || u + v > 1.0) continue;
        const double t = (e2x * qx + e2y * qy + e2z * qz) * inv;
        if (t <= 0.0) continue;
        const unsigned long long key =
            ((unsigned long long)__float_as_uint((float)t) << 32) | (unsigned long long)(f & 0xffffffffLL);
        atomicMin(&best[r], key);
    }
}
