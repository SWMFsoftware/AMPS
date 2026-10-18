//LiouvilleCOPS: collisionless model of the ROSINA/COPS nude and ram gauge pressures
//during the Rosetta end-of-mission descent (2016-09-30)
//
//usage: liouville_cops <input file>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <sys/stat.h>
#include <omp.h>

#include "SpiceUsr.h"

#include "Input.h"
#include "Mesh.h"
#include "BVH.h"
#include "Observations.h"
#include "Geometry.h"
#include "Source.h"
#include "Gauges.h"
#include "Liouville.h"
#include "Paths.h"

static void MakeDirectory(const std::string& dir) {
  std::error_code error;
  std::filesystem::create_directories(dir,error);
  if (error) {
    fprintf(stderr,"Error: cannot create the directory '%s': %s\n",dir.c_str(),error.message().c_str());
    exit(EXIT_FAILURE);
  }
}

//binary kernel file for the later fitting of the source rate and temperature (format described in README.md)
static void WriteKernels(const std::string& fname,int iPoint,const cObservationGeometry& g,const cLiouvilleResult& res,const std::vector<cSurfaceState>& state) {
  FILE *f=fopen(fname.c_str(),"wb");

  if (f==NULL) {
    fprintf(stderr,"Error: cannot write '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  int version=1,nK=(int)res.Kernels.size();
  double hdr[10]={g.et,res.VdotB_RG,g.xSC.x,g.xSC.y,g.xSC.z,g.SunDirection.x,g.SunDirection.y,g.SunDirection.z,g.HeliocentricDistance,res.Altitude};

  fwrite(&version,sizeof(int),1,f);
  fwrite(&iPoint,sizeof(int),1,f);
  fwrite(hdr,sizeof(double),10,f);
  fwrite(&nK,sizeof(int),1,f);

  for (auto& k : res.Kernels) {
    double v[4]={k.Kn,k.KPhi,k.KV,state[k.Face].CosSZA};
    int shadow=state[k.Face].Shadow ? 1 : 0;

    fwrite(&k.Face,sizeof(int),1,f);
    fwrite(v,sizeof(double),4,f);
    fwrite(&shadow,sizeof(int),1,f);
  }

  fclose(f);
}

//kernel file version 2 (FlightTimeRotation=exact): per-species contributions per unit emitted flux
static void WriteKernelsExact(const std::string& fname,int iPoint,const cObservationGeometry& g,const cLiouvilleResult& res,const std::vector<cSurfaceState>& state) {
  FILE *f=fopen(fname.c_str(),"wb");

  if (f==NULL) {
    fprintf(stderr,"Error: cannot write '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  int version=2,nK=(int)res.KernelsExact.size(),nS=Const::nSpecies;
  double hdr[13]={g.et,res.VdotB_RG,g.xSC.x,g.xSC.y,g.xSC.z,g.SunDirection.x,g.SunDirection.y,g.SunDirection.z,g.HeliocentricDistance,res.Altitude,g.Omega.x,g.Omega.y,g.Omega.z};

  fwrite(&version,sizeof(int),1,f);
  fwrite(&iPoint,sizeof(int),1,f);
  fwrite(hdr,sizeof(double),13,f);
  fwrite(&nS,sizeof(int),1,f);
  fwrite(&nK,sizeof(int),1,f);

  for (auto& k : res.KernelsExact) {
    double v[2*Const::nSpecies+2];
    int shadow=state[k.Face].Shadow ? 1 : 0;

    for (int s=0;s<Const::nSpecies;s++) v[s]=k.Gn[s],v[Const::nSpecies+s]=k.GPhi[s];
    v[2*Const::nSpecies]=state[k.Face].CosSZA;
    v[2*Const::nSpecies+1]=state[k.Face].Temperature;

    fwrite(&k.Face,sizeof(int),1,f);
    fwrite(v,sizeof(double),2*Const::nSpecies+2,f);
    fwrite(&shadow,sizeof(int),1,f);
  }

  fclose(f);
}

//=====================================================================================================
//verification against analytic solutions for a uniformly emitting sphere (convex body, constant F and T)
//   flat : n=sqrt(pi) beta F (1-sqrt(1-a^2/r^2))
//   rough: n=sqrt(pi) beta F a/(4r) ln((r+a)/(r-a))
//   both : radial flux Phi=Q/(4 pi r^2)
static int SelfTest(const cInput& input) {
  double a=input.GetDouble("SelfTestRadius",1000.0);
  int nRef=input.GetInt("SelfTestRefinement",6);
  double F=1.0E18,T=150.0;
  int s=Const::H2O;
  double beta=sqrt(Const::SpeciesMass[s]/(2.0*Const::kBoltzmann*T));
  int nFail=0;
  const double tol=input.GetDouble("SelfTestTolerance",0.01);

  cMesh mesh;
  mesh.MakeIcosphere(a,nRef);
  mesh.FinalizeGeometry(true);

  cBVH bvh;
  bvh.Build(&mesh);

  //1. BVH vs brute force ray tracing
  {
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> U(-1.0,1.0);
    int nMismatch=0,nTest=input.GetInt("SelfTestRays",2000);

    for (int i=0;i<nTest;i++) {
      Vec3 x0(3*a*U(rng),3*a*U(rng),3*a*U(rng)),d(U(rng),U(rng),U(rng));
      double t1=-1,t2=-1;
      int i1=bvh.FirstHit(x0,d,0.0,1.0E30,t1),i2=bvh.FirstHitBruteForce(x0,d,0.0,1.0E30,t2);

      if ((i1!=i2)||((i1>=0)&&(fabs(t1-t2)>1.0E-9*fabs(t2)))) nMismatch++;
    }

    printf("SelfTest BVH: %i random rays, %i mismatches with brute force -> %s\n",nTest,nMismatch,(nMismatch==0) ? "PASS" : "FAIL");
    if (nMismatch!=0) nFail++;
  }

  //1b. BVH nearest-surface distance vs. a brute-force scan
  {
    std::mt19937_64 rng(2);
    std::uniform_real_distribution<double> U(-1.0,1.0);
    int nMismatch=0,nTest=200;
    double maxError=0.0;

    for (int i=0;i<nTest;i++) {
      Vec3 x(3*a*U(rng),3*a*U(rng),3*a*U(rng));
      double brute=1.0E300;

      for (const auto& f : mesh.Faces) brute=std::min(brute,PointTriangleDistance(
          x,mesh.Nodes[f.node[0]],mesh.Nodes[f.node[1]],mesh.Nodes[f.node[2]]));

      double tree=bvh.NearestDistance(x);
      double error=fabs(tree-brute);
      maxError=std::max(maxError,error);
      if (error>1.0E-10*std::max(a,brute)) nMismatch++;
    }

    printf("SelfTest BVH distance: %i random points, %i mismatches, max |tree-brute|=%.3e m -> %s\n",
        nTest,nMismatch,maxError,(nMismatch==0) ? "PASS" : "FAIL");
    if (nMismatch!=0) nFail++;
  }

  //2. sphere density and flux
  std::vector<cSurfaceState> state(mesh.Faces.size());

  for (auto& st : state) {
    st.CosSZA=1.0,st.Shadow=false,st.Temperature=T;
    for (int k=0;k<Const::nSpecies;k++) st.Flux[k]=(k==s) ? F : 0.0;
  }

  double Q=F*mesh.TotalArea;
  //the effective radius of the polyhedron: same surface area as the sphere
  double aEff=sqrt(mesh.TotalArea/(4.0*Const::Pi));

  cGaugeModel gauges;
  gauges.Init(input);
  gauges.NGAcceptanceModel="full",gauges.RGAcceptanceModel="hemisphere",gauges.NGHeadOnWeight=false;

  const double rr[]={1.02,1.1,1.5,3.0,10.0};

  for (std::string bc : {"flat","rough"}) {
    cLiouville lv;
    lv.Init(input);
    lv.SurfaceBC=bc,lv.SpacecraftVelocity=false;

    for (double ratio : rr) {
      cObservationGeometry g;
      Vec3 dir=Normalized(Vec3(0.3,-0.5,0.81));

      g.et=0.0,g.xSC=dir*(ratio*a),g.vSC=Vec3();
      g.NG.Boresight=dir,g.RG.Boresight=-dir;   //RG looks at the nucleus centre
      g.NG.FOVShape=g.RG.FOVShape="RECTANGLE";

      cLiouvilleResult res;
      lv.EvaluateKernels(mesh,bvh,gauges,g,res);
      lv.ApplySource(mesh,gauges,state,gauges.TemperatureRG,res);

      double r=ratio*a,x=aEff/r;
      double nExact=(bc=="flat") ? Const::SqrtPi*beta*F*(1.0-sqrt(1.0-x*x)) : Const::SqrtPi*beta*F*aEff/(4.0*r)*log((r+aEff)/(r-aEff));
      double phiExact=Q/(4.0*Const::Pi*r*r);
      double en=res.n[s]/nExact-1.0,ef=res.PhiIn[s]/phiExact-1.0;
      bool ok=(fabs(en)<tol)&&(fabs(ef)<tol);

      printf("SelfTest sphere %-5s r/a=%5.2f: n=%.6e (exact %.6e, err %+.2e)  Phi=%.6e (exact %.6e, err %+.2e) -> %s\n",
          bc.c_str(),ratio,res.n[s],nExact,en,res.PhiIn[s],phiExact,ef,ok ? "PASS" : "FAIL");

      if (!ok) nFail++;
    }
  }

  //2b. FlightTimeRotation=exact
  //  0: omega=0 (srcCG source)           -> identical to the static solver
  //  1: uniform source, gauge on the rotation axis -> the static solution: the rotation maps the sphere onto itself
  //     and the first-order effect of the surface velocity (azimuthal) cancels by symmetry
  //  2: uniform source, gauge off the axis -> physical effect of the surface motion (reported)
  //  3: srcCG source, gauge off the axis   -> physical effect of the rotation during the flight (reported)
  {
    for (int test=0;test<4;test++) for (std::string bc : {"flat","rough"}) for (double ratio : {1.1,10.0}) {
      cSourceModel src;
      src.Init(input,mesh);
      if ((test==1)||(test==2)) src.UniformTestFlux=F;

      double omega=(test==0) ? 0.0 : 1.447E-4;
      Vec3 axis=Normalized(Vec3(0.2,0.1,1.0));
      Vec3 dir=(test==1) ? axis : Normalized(Vec3(0.3,-0.5,0.81));

      cLiouville lv;
      lv.Init(input);
      lv.SurfaceBC=bc,lv.SpacecraftVelocity=false,lv.FlightTimeRotation="exact";

      cObservationGeometry g;
      g.et=0.0,g.xSC=dir*(ratio*a),g.vSC=Vec3();
      g.NG.Boresight=dir,g.RG.Boresight=-dir;
      g.NG.FOVShape=g.RG.FOVShape="RECTANGLE";
      g.SunDirection=Normalized(Vec3(1.0,0.3,0.2)),g.HeliocentricDistance=3.8*Const::AU;
      g.Omega=axis*omega;

      std::vector<cSurfaceState> st;
      src.Evaluate(mesh,bvh,g.SunDirection,g.HeliocentricDistance,st);

      cLiouvilleResult rs,re;
      cLiouville ls=lv;
      ls.FlightTimeRotation="none";
      ls.EvaluateKernels(mesh,bvh,gauges,g,rs);
      ls.ApplySource(mesh,gauges,st,gauges.TemperatureRG,rs);
      lv.EvaluateExact(mesh,bvh,gauges,src,st,g,gauges.TemperatureRG,re);

      double en=re.n[s]/rs.n[s]-1.0,ef=re.PhiIn[s]/rs.PhiIn[s]-1.0;
      const char* name[4]={"omega=0, srcCG source       ","uniform source, on the axis ","uniform source, off the axis","srcCG source, off the axis  "};
      double tol=(test==0) ? 1.0E-6 : 2.0E-5;
      bool reported=(test>=2),ok=reported || ((fabs(en)<tol)&&(fabs(ef)<tol));

      printf("SelfTest exact rotation (%s) %-5s r/a=%5.1f: n/n_static-1=%+.2e Phi/Phi_static-1=%+.2e -> %s\n",name[test],bc.c_str(),ratio,en,ef,
          reported ? "(physical effect, reported)" : (ok ? "PASS" : "FAIL"));
      if (!ok) nFail++;
    }
  }

  //3. spacecraft velocity correction: for a gauge that sees the whole nucleus in front of the RG and with the NG
  //acceptance "full", the RG flux change must equal (V.b) n exactly
  {
    cLiouville lv;
    lv.Init(input);
    lv.SurfaceBC="flat",lv.SpacecraftVelocity=true;

    cObservationGeometry g;
    Vec3 dir=Normalized(Vec3(0.3,-0.5,0.81));
    cLiouvilleResult r0,r1;

    g.et=0.0,g.xSC=dir*(3.0*a),g.vSC=Vec3();
    g.NG.Boresight=dir,g.RG.Boresight=-dir;
    g.NG.FOVShape=g.RG.FOVShape="RECTANGLE";
    lv.EvaluateKernels(mesh,bvh,gauges,g,r0);
    lv.ApplySource(mesh,gauges,state,gauges.TemperatureRG,r0);

    g.vSC=-dir*2.0;   //2 m/s toward the nucleus (along the RG boresight)
    lv.EvaluateKernels(mesh,bvh,gauges,g,r1);
    lv.ApplySource(mesh,gauges,state,gauges.TemperatureRG,r1);

    double expected=2.0*r0.n[s],dphi=r1.PhiIn[s]-r0.PhiIn[s],err=dphi/expected-1.0;
    bool ok=fabs(err)<1.0E-10;

    printf("SelfTest spacecraft velocity: dPhi=%.6e, (V.b) n=%.6e, err %+.2e -> %s\n",dphi,expected,err,ok ? "PASS" : "FAIL");
    if (!ok) nFail++;
  }

  //4. cast shadows: a convex body has no shadowed sunlit elements
  {
    cInput in2=input;
    cSourceModel src;
    std::vector<cSurfaceState> st2;

    src.Init(in2,mesh);
    src.Evaluate(mesh,bvh,Normalized(Vec3(1.0,0.2,-0.3)),3.8*Const::AU,st2);

    int nShadow=0,nLit=0;
    for (auto& x : st2) {
      if (x.CosSZA>0.0) nLit++;
      if (x.Shadow) nShadow++;
    }

    printf("SelfTest cast shadows (convex sphere): %i sunlit elements, %i shadowed -> %s\n",nLit,nShadow,(nShadow==0) ? "PASS" : "FAIL");
    if (nShadow!=0) nFail++;
  }

  printf("SelfTest: %s (%i failures)\n",(nFail==0) ? "ALL PASSED" : "FAILED",nFail);
  return (nFail==0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

//=====================================================================================================
int main(int argc,char **argv) {
  if (argc!=2) {
    fprintf(stderr,"usage: %s <input file>\n",argv[0]);
    return EXIT_FAILURE;
  }

  cInput input;
  std::string InputFile=argv[1];
  if (!std::filesystem::exists(InputFile)) InputFile=CopsPaths::ProductInput(InputFile);
  input.Read(InputFile.c_str());

  std::string mode=input.GetOption("Mode",{"observations","selftest"},"observations");

  if (input.Has("Threads")) omp_set_num_threads(input.GetInt("Threads",1));
  printf("OpenMP threads: %i\n",omp_get_max_threads());

  if (mode=="selftest") return SelfTest(input);

  //Resolve and validate the destination before loading expensive mission data.
  std::string OutDir=CopsPaths::OutputPath(input.GetString("OutputDirectory","output"));
  printf("Output directory: %s\n",OutDir.c_str());

  //SPICE first: the DSK shape models need the frame kernel to resolve the frame of the plates
  cSpiceGeometry spice;
  spice.Init(input);

  //nucleus shape
  cMesh mesh;
  std::string ShapeFormat=input.GetOption("ShapeModelFormat",{"nastran","dsk"},"nastran");

  if (ShapeFormat=="nastran") mesh.ReadNastran(input.GetString("ShapeModel"));
  else mesh.ReadDSK(input.GetString("ShapeModel"),spice.BodyFrame);

  mesh.FinalizeGeometry(true,input.GetBool("RequireClosedSurface",true));

  cBVH bvh;
  bvh.Build(&mesh);

  //observations
  std::vector<cObservation> obs;
  ReadObservations(input.GetString("ObservationFile"),obs);

  std::string RefUnit=input.GetOption("ReferenceDataUnit",{"Pa","mbar"},"Pa");
  double RefToPa=(RefUnit=="mbar") ? 100.0 : 1.0;

  cSourceModel source;
  source.Init(input,mesh);

  cGaugeModel gauges;
  gauges.Init(input);
  gauges.PrintSummary();

  cLiouville lv;
  lv.Init(input);

  int iFirst=input.GetInt("FirstPoint",0);
  int iLast=input.GetInt("LastPoint",(int)obs.size()-1);
  int iStep=input.GetInt("PointStep",1);
  bool KernelOutput=input.GetBool("WriteKernels",true);
  bool GeometryOnly=input.GetBool("GeometryOnly",false);
  iLast=std::min(iLast,(int)obs.size()-1);

  if (GeometryOnly&&KernelOutput) {
    fprintf(stderr,"Error: GeometryOnly=on requires WriteKernels=off (transport kernels are not evaluated)\n");
    return EXIT_FAILURE;
  }

  //intervals of the reference data excluded from the comparison/fit: "Exclude<NG|RG> = <UTC begin> <UTC end>"
  std::vector<std::pair<double,double> > ExcludeNG,ExcludeRG;

  auto ReadIntervals=[&](const char* key,std::vector<std::pair<double,double> >& list) {
    for (auto& v : input.GetAll(key)) {
      char t0[64],t1[64];
      SpiceDouble e0,e1;

      if (sscanf(v.c_str(),"%63s %63s",t0,t1)!=2) {
        fprintf(stderr,"Error: %s must be '<UTC begin> <UTC end>': %s\n",key,v.c_str());
        exit(EXIT_FAILURE);
      }

      utc2et_c(t0,&e0),utc2et_c(t1,&e1);
      if (failed_c()||(e1<=e0)) {
        fprintf(stderr,"Error: invalid interval in %s: %s\n",key,v.c_str());
        exit(EXIT_FAILURE);
      }

      list.push_back({e0,e1});
      printf("%s: %s .. %s\n",key,t0,t1);
    }
  };

  ReadIntervals("ExcludeNG",ExcludeNG);
  ReadIntervals("ExcludeRG",ExcludeRG);

  //time interval of the data used for the fit: "FitInterval = <UTC begin> <UTC end>" (default: all observations)
  std::vector<std::pair<double,double> > FitInterval;
  ReadIntervals("FitInterval",FitInterval);

  if (FitInterval.size()>1) {
    fprintf(stderr,"Error: FitInterval can be given only once\n");
    exit(EXIT_FAILURE);
  }

  auto InFitInterval=[&](double et) {
    return FitInterval.empty() || ((et>=FitInterval[0].first)&&(et<=FitInterval[0].second));
  };

  auto Excluded=[](const std::vector<std::pair<double,double> >& list,double et) {
    for (auto& i : list) if ((et>=i.first)&&(et<=i.second)) return true;
    return false;
  };

  for (auto& k : input.UnusedKeys()) {
    fprintf(stderr,"Error: unknown key '%s' in the input file\n",k.c_str());
    return EXIT_FAILURE;
  }

  MakeDirectory(OutDir);
  if (KernelOutput) {
    MakeDirectory(OutDir+"/kernels");
    source.WriteElements(OutDir+"/elements.bin",mesh);
  }

  std::string fname=OutDir+"/COPS_model.dat";
  FILE *fout=fopen(fname.c_str(),"w");

  fprintf(fout,"# LiouvilleCOPS %s; input file: %s\n",GeometryOnly ? "geometry-only observation table" : "model of the COPS gauge pressures",InputFile.c_str());
  fprintf(fout,"# GeometryOnly=%s; geometry-only tables have zero transport/source/pressure columns and are not fitting bases\n",GeometryOnly ? "on" : "off");
  fprintf(fout,"# ShapeModel=%s (%s) FlightTimeRotation=%s FitInterval=%s\n",input.GetString("ShapeModel").c_str(),ShapeFormat.c_str(),lv.FlightTimeRotation.c_str(),
      input.Has("FitInterval") ? input.GetString("FitInterval").c_str() : "all");
  fprintf(fout,"# AnalyticalTmax=%.10g AnalyticalT75=%.10g AnalyticalTmin=%.10g HeliocentricModel=%s ActivityFile=%s\n",source.AnalyticalTmax,source.AnalyticalT75,source.AnalyticalTmin,
      source.HeliocentricModel.c_str(),(source.ActivityModel=="file") ? input.GetString("ActivityFile").c_str() : "none");
  fprintf(fout,"# SurfaceBC=%s TemperatureModel=%s ActivityModel=%s SZAMode=%s CastShadows=%s SpacecraftVelocity=%s\n",
      lv.SurfaceBC.c_str(),source.TemperatureModel.c_str(),source.ActivityModel.c_str(),source.SZAMode.c_str(),source.CastShadows ? "on" : "off",lv.SpacecraftVelocity ? "on" : "off");
  fprintf(fout,"# NGAcceptance=%s NGHeadOnWeight=%s RGAcceptance=%s RGPressureModel=%s TemperatureNG=%g TemperatureRGMode=%s SpeciesSensitivity=%s (H2O %g, CO2 %g) ReferenceDataUnit=%s\n",
      gauges.NGAcceptanceModel.c_str(),gauges.NGHeadOnWeight ? "on" : "off",gauges.RGAcceptanceModel.c_str(),gauges.RGPressureModel.c_str(),
      gauges.TemperatureNG,(gauges.TemperatureRGMode=="constant") ? (std::to_string(gauges.TemperatureRG)+"K").c_str() : "hk_table",gauges.SpeciesSensitivity ? "on" : "off",gauges.Sensitivity[0],gauges.Sensitivity[1],RefUnit.c_str());
  fprintf(fout,"# units: SI (m, m/s, m^-3, m^-2 s^-1, s^-1, Pa, K); reference values converted to Pa; a reference value <=0 is a data gap\n# useNG/useRG=1: the point is used in the fit (reference value >0, gauge outside the nucleus, not in Exclude<NG|RG>, inside FitInterval)\n");
  fprintf(fout,"#iPoint UTC et x y z altitude inside rh_AU VdotB_RG nVisible Q_H2O Q_CO2 n_H2O n_CO2 PhiRG_H2O PhiRG_CO2 pNG_H2O pNG_CO2 pNG pRG_H2O pRG_CO2 pRG refNG refRG NGBelowDeckFraction T_RG useNG useRG\n");

  double t0=omp_get_wtime();

  for (int ip=iFirst;ip<=iLast;ip+=iStep) {
    cObservationGeometry g;
    std::vector<cSurfaceState> state;
    cLiouvilleResult res={};
    double t1=omp_get_wtime();

    spice.Evaluate(obs[ip].UTC,g);
    if (GeometryOnly) {
      res.Altitude=bvh.NearestDistance(g.xSC);
      res.GaugeInsideNucleus=PointInsideSurface(mesh,bvh,g.xSC);
      res.VdotB_RG=lv.SpacecraftVelocity ? Dot(g.vSC,g.RG.Boresight) : 0.0;
      res.TemperatureRG=gauges.TemperatureRGAt(g.et);
    }
    else {
      source.Evaluate(mesh,bvh,g.SunDirection,g.HeliocentricDistance,state);
      if (lv.FlightTimeRotation=="exact") lv.EvaluateExact(mesh,bvh,gauges,source,state,g,gauges.TemperatureRGAt(g.et),res);
      else {
        lv.EvaluateKernels(mesh,bvh,gauges,g,res);
        lv.ApplySource(mesh,gauges,state,gauges.TemperatureRGAt(g.et),res);
      }
    }

    double pNG=res.pNG[0]+res.pNG[1],pRG=res.pRG[0]+res.pRG[1];
    int useNG=((obs[ip].NG>0.0)&&(!res.GaugeInsideNucleus)&&(!Excluded(ExcludeNG,g.et))&&InFitInterval(g.et)) ? 1 : 0;
    int useRG=((obs[ip].RG>0.0)&&(!res.GaugeInsideNucleus)&&(!Excluded(ExcludeRG,g.et))&&InFitInterval(g.et)) ? 1 : 0;

    fprintf(fout,"%5i %s %.6f %.6e %.6e %.6e %.6e %i %.8f %.6e %7i %.6e %.6e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.6e %.6e %.4e %.6f %i %i\n",
        ip,obs[ip].UTC.c_str(),g.et,g.xSC.x,g.xSC.y,g.xSC.z,res.Altitude,res.GaugeInsideNucleus ? 1 : 0,g.HeliocentricDistance/Const::AU,res.VdotB_RG,res.nVisibleElements,
        res.TotalProduction[0],res.TotalProduction[1],res.n[0],res.n[1],res.PhiIn[0],res.PhiIn[1],
        res.pNG[0],res.pNG[1],pNG,res.pRG[0],res.pRG[1],pRG,obs[ip].NG*RefToPa,obs[ip].RG*RefToPa,res.NGBelowDeckFraction,res.TemperatureRG,useNG,useRG);
    fflush(fout);

    if (KernelOutput) {
      char kname[64];

      sprintf(kname,"/kernels/point=%04i.bin",ip);
      if (lv.FlightTimeRotation=="exact") WriteKernelsExact(OutDir+kname,ip,g,res,state);
      else WriteKernels(OutDir+kname,ip,g,res,state);
    }

    printf("point %4i %s alt=%9.2f m visible=%6i pNG=%.3e pRG=%.3e (ref %.3e %.3e) %s%s[%.1f s]\n",ip,obs[ip].UTC.c_str(),res.Altitude,res.nVisibleElements,
        pNG,pRG,obs[ip].NG*RefToPa,obs[ip].RG*RefToPa,GeometryOnly ? "GEOMETRY ONLY " : "",res.GaugeInsideNucleus ? "INSIDE NUCLEUS " : "",omp_get_wtime()-t1);
    fflush(stdout);
  }

  fclose(fout);
  printf("Done: %s (%.1f s)\n",fname.c_str(),omp_get_wtime()-t0);
  return EXIT_SUCCESS;
}
