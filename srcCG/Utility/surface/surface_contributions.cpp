//Surface elements seen by the COPS nude gauge (NG) only or by the ram gauge (RG) only in a time/altitude interval
//
//usage: surface_contributions <input file> <output directory> <epsilon> <temperature [K]> <stride>
//           <utc start|-> <utc end|-> <altitude min [m]|-> <altitude max [m]|-> <write facets 0|1>
//
//Built and run by surface_contributions.py (same directory), which writes the input file (LiouvilleCOPS keys:
//SPICE kernels, shape model, observation file, gauge and transport settings) and compiles this file together with
//the unchanged LiouvilleCOPS sources (../../LiouvilleCOPS/src, all except main.cpp).
//
//Contribution of a surface element f at the sample time t (user decision 2026-10-27): the gauge pressure response to a
//unit uniform H2O flux (F=1 m^-2 s^-1) emitted with the temperature T (command line), computed from the geometric
//kernels of cLiouville::EvaluateKernels and the gauge conversions of cGaugeModel exactly as cLiouville::ApplySource:
//   H_NG(t,f)=NGPressure(H2O, F beta sqrt(pi)/2 Kn_f)
//   H_RG(t,f)=RGPressure(H2O, F (KPhi_f + (V.b_RG) beta sqrt(pi)/2 KV_f), T_RG(t)),   beta=sqrt(m_H2O/(2 k T))
//Classification at one sample, with epsilon (default 1e-3; the rule of the fine tuning, FINE_TUNING_STEP.md §3):
//   seen by NG: H_NG>0;  seen by RG: |H_RG|>0
//   NG only   : H_NG>0  and |H_RG| <= epsilon H_NG
//   RG only   : |H_RG|>0 and  H_NG  <= epsilon |H_RG|
//Over the N selected samples:
//   at least once        : NG only at >=1 sample
//   all the time, strict : NG only at every one of the N samples (hence seen by NG at every sample)
//   all the time, loose  : NG only at every sample at which either gauge sees the element (seen at >=1 sample)
//Samples: the observation times of the observation file inside the UTC and/or altitude interval (every stride-th);
//samples with the gauge inside the nucleus are skipped and reported.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <omp.h>

#include "SpiceUsr.h"

#include "Input.h"
#include "Mesh.h"
#include "BVH.h"
#include "Observations.h"
#include "Geometry.h"
#include "Gauges.h"
#include "Liouville.h"

int main(int argc,char **argv) {
  if (argc!=11) {
    fprintf(stderr,"usage: %s <input> <outdir> <epsilon> <T> <stride> <utc0|-> <utc1|-> <altmin|-> <altmax|-> <facets 0|1>\n",argv[0]);
    return EXIT_FAILURE;
  }

  std::string OutDir=argv[2];
  double eps=atof(argv[3]),T=atof(argv[4]);
  int stride=atoi(argv[5]);
  std::string utc0=argv[6],utc1=argv[7],amin=argv[8],amax=argv[9];
  bool WriteFacets=atoi(argv[10])!=0;

  if (!(eps>0.0)||!(eps<1.0)||!(T>0.0)||(stride<1)) {
    fprintf(stderr,"Error: require 0<epsilon<1, T>0, stride>=1\n");
    return EXIT_FAILURE;
  }

  cInput input;
  input.Read(argv[1]);
  if (input.Has("Threads")) omp_set_num_threads(input.GetInt("Threads",1));
  printf("OpenMP threads: %i\n",omp_get_max_threads());

  cSpiceGeometry spice;
  spice.Init(input);

  cMesh mesh;
  std::string ShapeFormat=input.GetOption("ShapeModelFormat",{"nastran","dsk"},"nastran");
  if (ShapeFormat=="nastran") mesh.ReadNastran(input.GetString("ShapeModel"));
  else mesh.ReadDSK(input.GetString("ShapeModel"),spice.BodyFrame);
  mesh.FinalizeGeometry(true,input.GetBool("RequireClosedSurface",true));

  cBVH bvh;
  bvh.Build(&mesh);

  std::vector<cObservation> obs;
  ReadObservations(input.GetString("ObservationFile"),obs);

  cGaugeModel gauges;
  gauges.Init(input);
  gauges.PrintSummary();

  cLiouville lv;
  lv.Init(input);
  if (lv.FlightTimeRotation!="none") {
    fprintf(stderr,"Error: FlightTimeRotation must be none (static geometric kernels)\n");
    return EXIT_FAILURE;
  }

  //interval of the samples
  double et0=-1.0E300,et1=1.0E300,a0=-1.0E300,a1=1.0E300;
  if (utc0!="-") {utc2et_c(utc0.c_str(),&et0);}
  if (utc1!="-") {utc2et_c(utc1.c_str(),&et1);}
  if (failed_c()) {
    fprintf(stderr,"Error: cannot convert the UTC interval '%s' .. '%s'\n",utc0.c_str(),utc1.c_str());
    return EXIT_FAILURE;
  }
  if (amin!="-") a0=atof(amin.c_str());
  if (amax!="-") a1=atof(amax.c_str());
  if ((et1<et0)||(a1<a0)) {
    fprintf(stderr,"Error: empty interval\n");
    return EXIT_FAILURE;
  }

  //selection: observation times in the UTC interval, then the altitude interval (closest distance to the surface)
  std::vector<int> candidates;
  int nInside=0;
  for (int ip=0;ip<(int)obs.size();ip++) {
    SpiceDouble et;
    utc2et_c(obs[ip].UTC.c_str(),&et);
    if ((et<et0)||(et>et1)) continue;

    if ((amin!="-")||(amax!="-")) {
      cObservationGeometry g;
      spice.Evaluate(obs[ip].UTC,g);
      double alt=bvh.NearestDistance(g.xSC);
      if ((alt<a0)||(alt>a1)) continue;
    }
    candidates.push_back(ip);
  }

  std::vector<int> selected;
  for (size_t i=0;i<candidates.size();i+=stride) selected.push_back(candidates[i]);

  if (selected.empty()) {
    fprintf(stderr,"Error: no observation time in the interval\n");
    return EXIT_FAILURE;
  }
  printf("samples: %zu observation times in the interval, %zu selected (stride %i)\n",candidates.size(),selected.size(),stride);

  const int nFaces=(int)mesh.Faces.size();
  const int s=Const::H2O;
  const double beta=sqrt(Const::SpeciesMass[s]/(2.0*Const::kBoltzmann*T));
  std::vector<uint32_t> nSeenNG(nFaces,0),nSeenRG(nFaces,0),nSeenAny(nFaces,0),nNGonly(nFaces,0),nRGonly(nFaces,0);

  FILE *fs=fopen((OutDir+"/samples.dat").c_str(),"w");
  if (fs==NULL) {
    fprintf(stderr,"Error: cannot write '%s/samples.dat'\n",OutDir.c_str());
    return EXIT_FAILURE;
  }
  fprintf(fs,"# per sample: facets seen by NG / RG / either, NG only, RG only, both above epsilon (mixed); areas [m^2]\n");
  fprintf(fs,"# epsilon=%.6g T=%.6g K (unit uniform H2O flux); RG temperature from the gauge model\n",eps,T);
  fprintf(fs,"#iPoint UTC et altitude_m T_RG VdotB_RG nSeenNG nSeenRG nSeenAny nNGonly nRGonly nMixed areaNGonly areaRGonly areaMixed\n");

  int nUsed=0;
  double t0=omp_get_wtime();

  for (int ip : selected) {
    cObservationGeometry g;
    cLiouvilleResult res={};
    double t1=omp_get_wtime();

    spice.Evaluate(obs[ip].UTC,g);
    lv.EvaluateKernels(mesh,bvh,gauges,g,res);

    if (res.GaugeInsideNucleus) {
      nInside++;
      printf("point %4i %s: gauge inside the nucleus, skipped\n",ip,obs[ip].UTC.c_str());
      continue;
    }

    double TRG=gauges.TemperatureRGAt(g.et);
    double cNG=gauges.NGPressure(s,1.0)*beta*Const::SqrtPi/2.0;   //Pa per unit Kn (linear in n)
    double cRG=gauges.RGPressure(s,1.0,TRG);                         //Pa per unit PhiIn (linear in PhiIn)
    long nNG=0,nRG=0,nAny=0,nN=0,nR=0,nM=0;
    double aN=0.0,aR=0.0,aM=0.0;

    for (auto& k : res.Kernels) {
      double hNG=fabs(cNG*k.Kn);
      double hRG=fabs(cRG*(k.KPhi+res.VdotB_RG*beta*Const::SqrtPi/2.0*k.KV));
      bool seenNG=hNG>0.0,seenRG=hRG>0.0;
      if (!seenNG&&!seenRG) continue;

      bool ngOnly=seenNG&&(hRG<=eps*hNG),rgOnly=seenRG&&(hNG<=eps*hRG);
      double area=mesh.Faces[k.Face].Area;

      nSeenAny[k.Face]++,nAny++;
      if (seenNG) nSeenNG[k.Face]++,nNG++;
      if (seenRG) nSeenRG[k.Face]++,nRG++;
      if (ngOnly) nNGonly[k.Face]++,nN++,aN+=area;
      else if (rgOnly) nRGonly[k.Face]++,nR++,aR+=area;
      else nM++,aM+=area;
    }

    nUsed++;
    fprintf(fs,"%5i %s %.6f %.6e %.4f %.6e %ld %ld %ld %ld %ld %ld %.6e %.6e %.6e\n",ip,obs[ip].UTC.c_str(),g.et,res.Altitude,TRG,res.VdotB_RG,
        nNG,nRG,nAny,nN,nR,nM,aN,aR,aM);
    fflush(fs);
    printf("point %4i %s alt=%9.2f m seen NG %ld RG %ld; NG only %ld, RG only %ld, mixed %ld [%.1f s]\n",ip,obs[ip].UTC.c_str(),res.Altitude,
        nNG,nRG,nN,nR,nM,omp_get_wtime()-t1);
    fflush(stdout);
  }
  fclose(fs);

  if (nUsed==0) {
    fprintf(stderr,"Error: every selected sample has the gauge inside the nucleus\n");
    return EXIT_FAILURE;
  }

  //classification over the interval
  enum {SeenNG,SeenRG,SeenBoth,SeenAny,NGonce,RGonce,NGstrict,RGstrict,NGloose,RGloose,nClass};
  const char* name[nClass]={"seen_by_NG_at_least_once","seen_by_RG_at_least_once","seen_by_both_at_least_once","seen_by_either_at_least_once",
      "NG_only_at_least_once","RG_only_at_least_once","NG_only_all_the_time_strict","RG_only_all_the_time_strict",
      "NG_only_all_the_time_loose","RG_only_all_the_time_loose"};
  long count[nClass]={0};
  double area[nClass]={0.0};
  std::vector<uint8_t> flags(WriteFacets ? nFaces : 0,0);
  const uint32_t N=(uint32_t)nUsed;

  for (int f=0;f<nFaces;f++) {
    bool c[nClass];
    c[SeenNG]=nSeenNG[f]>0,c[SeenRG]=nSeenRG[f]>0,c[SeenBoth]=c[SeenNG]&&c[SeenRG],c[SeenAny]=nSeenAny[f]>0;
    c[NGonce]=nNGonly[f]>0,c[RGonce]=nRGonly[f]>0;
    c[NGstrict]=nNGonly[f]==N,c[RGstrict]=nRGonly[f]==N;
    c[NGloose]=c[SeenAny]&&(nNGonly[f]==nSeenAny[f]),c[RGloose]=c[SeenAny]&&(nRGonly[f]==nSeenAny[f]);

    for (int k=0;k<nClass;k++) if (c[k]) count[k]++,area[k]+=mesh.Faces[f].Area;
    if (WriteFacets) flags[f]=(c[NGonce]?1:0)|(c[RGonce]?2:0)|(c[NGstrict]?4:0)|(c[RGstrict]?8:0)|(c[NGloose]?16:0)|(c[RGloose]?32:0)|(c[SeenNG]?64:0)|(c[SeenRG]?128:0);
  }

  FILE *fj=fopen((OutDir+"/counts.json").c_str(),"w");
  fprintf(fj,"{\n \"epsilon\": %.17g,\n \"temperature_K\": %.17g,\n \"n_facets\": %d,\n \"total_area_m2\": %.10e,\n",eps,T,nFaces,mesh.TotalArea);
  fprintf(fj," \"n_samples\": %d,\n \"n_samples_in_interval\": %zu,\n \"n_samples_gauge_inside_nucleus\": %d,\n \"stride\": %d,\n",nUsed,candidates.size(),nInside,stride);
  fprintf(fj," \"first_sample_point\": %d,\n \"last_sample_point\": %d,\n \"first_sample_utc\": \"%s\",\n \"last_sample_utc\": \"%s\",\n",
      selected.front(),selected.back(),obs[selected.front()].UTC.c_str(),obs[selected.back()].UTC.c_str());
  fprintf(fj," \"counts\": {");
  for (int k=0;k<nClass;k++) fprintf(fj,"%s\n  \"%s\": %ld",k ? "," : "",name[k],count[k]);
  fprintf(fj,"\n },\n \"areas_m2\": {");
  for (int k=0;k<nClass;k++) fprintf(fj,"%s\n  \"%s\": %.10e",k ? "," : "",name[k],area[k]);
  fprintf(fj,"\n }\n}\n");
  fclose(fj);

  if (WriteFacets) {
    //facets.bin: int32 nFaces, int32 nSamples, then uint32 nSeenNG, nSeenRG, nSeenAny, nNGonly, nRGonly [nFaces each], uint8 flags [nFaces]
    FILE *fb=fopen((OutDir+"/facets.bin").c_str(),"wb");
    int hdr[2]={nFaces,nUsed};
    fwrite(hdr,sizeof(int),2,fb);
    for (auto* v : {&nSeenNG,&nSeenRG,&nSeenAny,&nNGonly,&nRGonly}) fwrite(v->data(),sizeof(uint32_t),nFaces,fb);
    fwrite(flags.data(),1,nFaces,fb);
    fclose(fb);
  }

  printf("Done: %d samples, %.1f s\n",nUsed,omp_get_wtime()-t0);
  return EXIT_SUCCESS;
}
