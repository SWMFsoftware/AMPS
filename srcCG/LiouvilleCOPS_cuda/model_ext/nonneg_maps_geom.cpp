// ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_maps_geom.cpp (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_maps_geom.cpp)
// SHA256 of the run copy = SHA256 of this file without these header lines: 16d9971e392660425f6e7e6df9adfca1bf74b5e86aa5ace00e72d40fb6af2297
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: C++ source. All-facet solar geometry of a LiouvilleCOPS base at given UTC times for nonneg_maps.py: Sun = spkpos_c(SUN, et,
//   67P/C-G_CK, LT+S) as Geometry.cpp, cos SZA = normal . Sun direction, r_h; compiled by nonneg_maps.py geometry
//   against the UNCHANGED reference sources of ../LiouvilleCOPS (no edit of the tree) and CSPICE.
// USED BY: nonneg_maps.py geometry (jobs 124318-124325).
// ==== end of header ====
// nonneg_maps_geom.cpp -- all-facet solar geometry (cos SZA, cast shadow) of a LiouvilleCOPS fitting base at given
// UTC times, for the activity / source-rate maps of model_ext/nonneg_maps.py (user request 2026-11-11: "save the
// activity and source rate maps for visualization"; maps at 00:00 .. 10:00, 10:33:10, 10:39:26 UTC for ALL facets).
//
// Local, untracked model extension (no file of the git tree is edited). It is compiled together with the UNCHANGED
// reference sources ../LiouvilleCOPS/src/*.cpp (except main.cpp) by jobs/nonneg_maps.sbatch.
//
// The geometry is the base builder's (native/liouville_cops_cuda, which wrote runs/.../base):
//  * SPICE: the base's own input file (boundary.in); cSpiceGeometry::Init loads the same kernels. The Sun vector is the
//    code of Geometry.cpp cSpiceGeometry::Evaluate: et = utc2et_c(UTC); spkpos_c("SUN", et, BodyFrame, "LT+S",
//    BodyName) [km] x 1e3 -> xSun [m]; rh = Length(xSun); SunDirection = xSun / rh. (The spacecraft part of Evaluate is
//    not needed and not called: the map times need no spacecraft or CK coverage of the spacecraft.)
//  * mesh and BVH: cMesh::ReadNastran / ReadDSK + FinalizeGeometry, cBVH::Build (as main.cpp).
//  * cos SZA = Dot(face normal, SunDirection): the host part of gpu_liouville.cu cSourceModel::Evaluate (= Source.cpp).
//  * cast shadow (CastShadows=on): for cos > 0, x0 = Centroid + Normal * (1e-6 MaxEdge), shadow = AnyHit(x0,
//    SunDirection, 0, 1e12, self) (k_shadow of gpu_liouville.cu = Source.cpp Evaluate). The builder ran this on the GPU
//    with --fmad=false (IEEE operations as written); here nonneg_maps_shadow.cpp and BVH.cpp are compiled with
//    -ffp-contract=off for the same arithmetic.
//  * flux (validation only): cSourceModel::EvaluateElement (Source.cpp) with the base's own activity (srcCG map of the
//    input) and Q_s = sum_e F_s,e area_e in the order of Liouville.cpp ApplySource = the Q_H2O, Q_CO2 columns of the
//    base's COPS_model.dat.
//
// usage: nonneg_maps_geom <base input file> <times file> <output dir (must not exist)> <threads>
//   times file: one "label UTC" pair per line ('#' comments)
// products in <output dir>:
//   mesh_nodes.f8       nNodes x 3 float64 [m] (cMesh::Nodes)
//   mesh_faces.i4       nFaces x 3 int32, 0-based node indices in the order of FinalizeGeometry (outward normals)
//   mesh_facegeom.f8    nFaces x 8 float64: centroid[3] [m], normal[3], area [m^2], MaxEdge [m]
//   <label>.cos.f8      nFaces float64 cos SZA
//   <label>.shadow.u1   nFaces uint8: 1 = cast shadow (set only where cos > 0 and CastShadows = on)
//   times.txt           label UTC et sun_x sun_y sun_z rh_m n_cos_pos n_shadow Q_H2O_base Q_CO2_base (%.17g)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <omp.h>

#include "SpiceUsr.h"

#include "Input.h"
#include "Mesh.h"
#include "BVH.h"
#include "Geometry.h"
#include "Source.h"
#include "Constants.h"

// nonneg_maps_shadow.cpp (compiled with -ffp-contract=off)
void NonnegMapsShadow(const cMesh& mesh,const cBVH& bvh,const Vec3& SunDirection,const std::vector<double>& CosSZA,
                      bool CastShadows,std::vector<unsigned char>& shadow);

static void Check(const char* where) {
  if (failed_c()) {
    SpiceChar msg[1841];
    getmsg_c("LONG",sizeof(msg),msg);
    fprintf(stderr,"SPICE error in %s: %s\n",where,msg);
    exit(EXIT_FAILURE);
  }
}

template <class T> static void WriteArray(const std::string& fname,const T* data,size_t n) {
  FILE* f=fopen(fname.c_str(),"wb");
  if (f==NULL) {fprintf(stderr,"Error: cannot write '%s'\n",fname.c_str()); exit(EXIT_FAILURE);}
  if (fwrite(data,sizeof(T),n,f)!=n) {fprintf(stderr,"Error: short write '%s'\n",fname.c_str()); exit(EXIT_FAILURE);}
  if (fclose(f)!=0) {fprintf(stderr,"Error: cannot close '%s'\n",fname.c_str()); exit(EXIT_FAILURE);}
}

int main(int argc,char** argv) {
  if (argc!=5) {
    fprintf(stderr,"usage: %s <base input file> <times file> <output dir (new)> <threads>\n",argv[0]);
    return EXIT_FAILURE;
  }
  const std::string InputFile=argv[1],TimesFile=argv[2],OutDir=argv[3];
  const int threads=atoi(argv[4]);
  if (threads<1) {fprintf(stderr,"Error: threads must be >= 1\n"); return EXIT_FAILURE;}
  omp_set_num_threads(threads);
  if (std::filesystem::exists(OutDir)) {
    fprintf(stderr,"Error: the output directory exists (never overwritten): %s\n",OutDir.c_str());
    return EXIT_FAILURE;
  }

  // times
  std::vector<std::string> labels,utcs;
  {
    std::ifstream fin(TimesFile);
    std::string line;
    if (!fin) {fprintf(stderr,"Error: cannot read '%s'\n",TimesFile.c_str()); return EXIT_FAILURE;}
    while (std::getline(fin,line)) {
      size_t c=line.find('#');
      if (c!=std::string::npos) line=line.substr(0,c);
      std::istringstream ss(line);
      std::string l,u,extra;
      if (!(ss >> l)) continue;
      if (!(ss >> u)||(ss >> extra)) {fprintf(stderr,"Error: times file line must be 'label UTC': %s\n",line.c_str()); return EXIT_FAILURE;}
      labels.push_back(l),utcs.push_back(u);
    }
  }
  if (labels.empty()) {fprintf(stderr,"Error: no times in '%s'\n",TimesFile.c_str()); return EXIT_FAILURE;}

  double t0=omp_get_wtime();
  cInput input;
  input.Read(InputFile.c_str());

  cSpiceGeometry spice;
  spice.Init(input);

  cMesh mesh;
  std::string ShapeFormat=input.GetOption("ShapeModelFormat",{"nastran","dsk"},"nastran");
  if (ShapeFormat=="nastran") mesh.ReadNastran(input.GetString("ShapeModel"));
  else mesh.ReadDSK(input.GetString("ShapeModel"),spice.BodyFrame);
  mesh.FinalizeGeometry(true,input.GetBool("RequireClosedSurface",true));

  cBVH bvh;
  bvh.Build(&mesh);

  cSourceModel source;
  source.Init(input,mesh,spice.BodyFrame);
  printf("nonneg_maps_geom: SZAMode=%s CastShadows=%s HeliocentricModel=%s ActivityModel=%s; setup %.1f s; %d threads\n",
         source.SZAMode.c_str(),source.CastShadows ? "on" : "off",source.HeliocentricModel.c_str(),
         source.ActivityModel.c_str(),omp_get_wtime()-t0,threads);
  fflush(stdout);

  std::filesystem::create_directories(OutDir);
  const size_t nF=mesh.Faces.size(),nN=mesh.Nodes.size();
  {
    std::vector<double> xyz(3*nN);
    for (size_t i=0;i<nN;i++) xyz[3*i]=mesh.Nodes[i].x,xyz[3*i+1]=mesh.Nodes[i].y,xyz[3*i+2]=mesh.Nodes[i].z;
    WriteArray(OutDir+"/mesh_nodes.f8",xyz.data(),xyz.size());
    std::vector<int> fc(3*nF);
    std::vector<double> fg(8*nF);
    for (size_t i=0;i<nF;i++) {
      const cTriangle& f=mesh.Faces[i];
      for (int k=0;k<3;k++) fc[3*i+k]=f.node[k];
      double v[8]={f.Centroid.x,f.Centroid.y,f.Centroid.z,f.Normal.x,f.Normal.y,f.Normal.z,f.Area,f.MaxEdge};
      for (int k=0;k<8;k++) fg[8*i+k]=v[k];
    }
    WriteArray(OutDir+"/mesh_faces.i4",fc.data(),fc.size());
    WriteArray(OutDir+"/mesh_facegeom.f8",fg.data(),fg.size());
  }

  FILE* ftab=fopen((OutDir+"/times.txt").c_str(),"w");
  if (ftab==NULL) {fprintf(stderr,"Error: cannot write times.txt\n"); return EXIT_FAILURE;}
  fprintf(ftab,"# nonneg_maps_geom: input %s; shape %s (%zu faces, %zu nodes); SZAMode=%s CastShadows=%s HeliocentricModel=%s\n",
          InputFile.c_str(),input.GetString("ShapeModel").c_str(),nF,nN,source.SZAMode.c_str(),
          source.CastShadows ? "on" : "off",source.HeliocentricModel.c_str());
  fprintf(ftab,"# sun: spkpos_c(SUN, et, %s, LT+S, %s) (Geometry.cpp); Q_*_base = sum_e F_e area_e with the input's activity (Liouville.cpp ApplySource order)\n",
          spice.BodyFrame.c_str(),spice.BodyName.c_str());
  fprintf(ftab,"#label UTC et sun_x sun_y sun_z rh_m n_cos_pos n_shadow Q_H2O_base Q_CO2_base\n");

  std::vector<double> CosSZA(nF);
  std::vector<unsigned char> shadow(nF);
  std::vector<double> FluxH2O(nF),FluxCO2(nF);

  for (size_t k=0;k<labels.size();k++) {
    double t1=omp_get_wtime();
    SpiceDouble et,pos[3],lt;
    utc2et_c(utcs[k].c_str(),&et);
    Check("utc2et_c");
    // Geometry.cpp cSpiceGeometry::Evaluate: Sun as seen from the nucleus (apparent direction: light time and aberration)
    spkpos_c("SUN",et,spice.BodyFrame.c_str(),"LT+S",spice.BodyName.c_str(),pos,&lt);
    Check("spkpos_c(SUN)");
    Vec3 xSun(1.0E3*pos[0],1.0E3*pos[1],1.0E3*pos[2]);
    const double HeliocentricDistance=Length(xSun);
    const Vec3 SunDirection=xSun/HeliocentricDistance;

    // gpu_liouville.cu cSourceModel::Evaluate, host part (as Source.cpp)
    #pragma omp parallel for schedule(static)
    for (long i=0;i<(long)nF;i++) CosSZA[i]=Dot(mesh.Faces[i].Normal,SunDirection);
    NonnegMapsShadow(mesh,bvh,SunDirection,CosSZA,source.CastShadows,shadow);
    #pragma omp parallel for schedule(dynamic,256)
    for (long i=0;i<(long)nF;i++) {
      double Flux[Const::nSpecies],T;
      source.EvaluateElement((int)i,CosSZA[i],(CosSZA[i]<=0.0)||(shadow[i]!=0),HeliocentricDistance,Flux,T);
      FluxH2O[i]=Flux[Const::H2O],FluxCO2[i]=Flux[Const::CO2];
    }
    // Liouville.cpp ApplySource: Q+=state[i].Flux[s]*mesh.Faces[i].Area in face order
    double Q[2]={0.0,0.0};
    for (size_t i=0;i<nF;i++) Q[0]+=FluxH2O[i]*mesh.Faces[i].Area;
    for (size_t i=0;i<nF;i++) Q[1]+=FluxCO2[i]*mesh.Faces[i].Area;
    long nPos=0,nShadow=0;
    for (size_t i=0;i<nF;i++) nPos+=(CosSZA[i]>0.0),nShadow+=(shadow[i]!=0);

    WriteArray(OutDir+"/"+labels[k]+".cos.f8",CosSZA.data(),nF);
    WriteArray(OutDir+"/"+labels[k]+".shadow.u1",shadow.data(),nF);
    fprintf(ftab,"%s %s %.17g %.17g %.17g %.17g %.17g %ld %ld %.17g %.17g\n",labels[k].c_str(),utcs[k].c_str(),et,
            SunDirection.x,SunDirection.y,SunDirection.z,HeliocentricDistance,nPos,nShadow,Q[0],Q[1]);
    fflush(ftab);
    printf("nonneg_maps_geom: %s %s et=%.6f rh=%.9f AU cos>0 %ld shadow %ld Q_H2O=%.6e Q_CO2=%.6e [%.1f s]\n",
           labels[k].c_str(),utcs[k].c_str(),et,HeliocentricDistance/Const::AU,nPos,nShadow,Q[0],Q[1],omp_get_wtime()-t1);
    fflush(stdout);
  }
  if (fclose(ftab)!=0) {fprintf(stderr,"Error: cannot close times.txt\n"); return EXIT_FAILURE;}
  printf("nonneg_maps_geom: done, %zu times, %.1f s\n",labels.size(),omp_get_wtime()-t0);
  return EXIT_SUCCESS;
}
