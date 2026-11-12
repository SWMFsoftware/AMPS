// ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_maps_shadow.cpp (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_maps_shadow.cpp)
// SHA256 of the run copy = SHA256 of this file without these header lines: 492781396726663324bb66c3fc7befc7b846a45c7118d12344f128bef2df6022
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: C++ source. Cast shadow for every facet (ray from centroid + 1e-6 MaxEdge normal toward the Sun, only for cos SZA > 0) as the base
//   builder's k_shadow; compiled with -ffp-contract=off (the builder uses --fmad=false). Validated: 0 shadow and 0
//   night mismatches and bitwise cos on 55,100,036 kernel records at 17 points.
// USED BY: nonneg_maps.py geometry (jobs 124318-124325).
// ==== end of header ====
// nonneg_maps_shadow.cpp -- cast shadow of every sunlit facet (part of nonneg_maps_geom; see nonneg_maps_geom.cpp).
// Compiled with -ffp-contract=off (as ../LiouvilleCOPS/src/BVH.cpp in jobs/nonneg_maps.sbatch): the base builder
// evaluated exactly these operations on the GPU with --fmad=false (gpu_liouville.cu k_shadow):
//   x0 = Centroid + Normal * (1e-6 MaxEdge);  shadow = AnyHit(x0, SunDirection, 0, 1e12, ignore = the facet itself)
// for the facets with cos SZA > 0 when CastShadows = on (Source.cpp cSourceModel::Evaluate); 0 elsewhere.
#include <vector>

#include "Mesh.h"
#include "BVH.h"

void NonnegMapsShadow(const cMesh& mesh,const cBVH& bvh,const Vec3& SunDirection,const std::vector<double>& CosSZA,
                      bool CastShadows,std::vector<unsigned char>& shadow) {
  const long n=(long)mesh.Faces.size();
  shadow.assign(n,0);
  #pragma omp parallel for schedule(dynamic,1024)
  for (long i=0;i<n;i++) {
    if (!((CosSZA[i]>0.0)&&CastShadows)) continue;
    const cTriangle& f=mesh.Faces[i];
    Vec3 x0=f.Centroid+f.Normal*(1.0E-6*f.MaxEdge);
    shadow[i]=bvh.AnyHit(x0,SunDirection,0.0,1.0E12,(int)i) ? 1 : 0;
  }
}
