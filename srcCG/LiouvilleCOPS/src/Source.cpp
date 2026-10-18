#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "Source.h"
#include "SrcCGTables.h"
#include "Paths.h"
#include "Provenance.h"

//real spherical harmonics as defined in srcCG (Comet.cpp, sphericalHarmonic()); copied verbatim
static double SphericalHarmonic(int i,double colatitude,double longitude) {
  switch (i) {
  case 0:  return 1.0;
  case 1:  return sin(colatitude)*sin(longitude);
  case 2:  return cos(colatitude);
  case 3:  return sin(colatitude)*cos(longitude);
  case 4:  return pow(sin(colatitude),2.0)*sin(2*longitude);
  case 5:  return cos(colatitude)*sin(colatitude)*sin(longitude);
  case 6:  return 3*pow(cos(colatitude),2.0)-1.0;
  case 7:  return cos(colatitude)*sin(colatitude)*cos(longitude);
  case 8:  return pow(sin(colatitude),2.0)*cos(2*longitude);
  case 9:  return pow(sin(colatitude),3.0)*sin(3*longitude);
  case 10: return pow(sin(colatitude),2.0)*cos(colatitude)*sin(2*longitude);
  case 11: return 5*(pow(cos(colatitude),2.0)-1)*sin(colatitude)*sin(longitude);
  case 12: return 5*pow(cos(colatitude),3.0)-3*cos(colatitude);
  case 13: return 5*(pow(cos(colatitude),2.0)-1)*sin(colatitude)*cos(longitude);
  case 14: return pow(sin(colatitude),2.0)*cos(colatitude)*cos(2*longitude);
  case 15: return pow(sin(colatitude),3.0)*cos(3*longitude);
  case 16: return pow(sin(colatitude),4.0)*sin(4*longitude);
  case 17: return pow(sin(colatitude),3.0)*cos(colatitude)*sin(3*longitude);
  case 18: return pow(sin(colatitude),2.0)*(7*pow(cos(colatitude),2.0)-1)*sin(2*longitude);
  case 19: return sin(colatitude)*(7*pow(cos(colatitude),3.0)-3*cos(colatitude))*sin(longitude);
  case 20: return 35*pow(cos(colatitude),4.0)-30*pow(cos(colatitude),2.0)+3;
  case 21: return sin(colatitude)*(7*pow(cos(colatitude),3.0)-3*cos(colatitude))*cos(longitude);
  case 22: return pow(sin(colatitude),2.0)*(7*pow(cos(colatitude),2.0)-1)*cos(2*longitude);
  case 23: return pow(sin(colatitude),3.0)*cos(colatitude)*cos(3*longitude);
  case 24: return pow(sin(colatitude),4.0)*cos(4*longitude);
  }

  fprintf(stderr,"Error: order of the spherical harmonic is too high\n");
  exit(EXIT_FAILURE);
}

void cSourceModel::Init(const cInput& input,const cMesh& mesh,const std::string& bodyFrame) {
  ActivityModel=input.GetOption("ActivityModel",{"shap5","shap5_1","file"},"shap5");
  TemperatureModel=input.GetOption("TemperatureModel",{"analytical","constant","bjorn"},"analytical");
  SZAMode=input.GetOption("SZAMode",{"continuous","srccg_table"},"continuous");
  CastShadows=input.GetBool("CastShadows",true);

  ConstantTemperature=input.GetDouble("ConstantTemperature",180.0);
  TemperatureCorrectionFactor=input.GetDouble("SurfaceTemperatureCorrectionFactor",1.0);
  BjornDistanceIndex=input.GetInt("BjornTemperatureDistanceIndex",4);

  if ((BjornDistanceIndex<0)||(BjornDistanceIndex>5)) {
    fprintf(stderr,"Error: BjornTemperatureDistanceIndex must be in [0,5]\n");
    exit(EXIT_FAILURE);
  }

  //activity coefficients: srcCG tables or a file with 3 rows (H2O, CO2, CO) of 25 coefficients (e.g. written by tools/fit_cops.py)
  if (ActivityModel=="file") {
    std::string fname=CopsPaths::ProductInput(input.GetString("ActivityFile"));
    std::ifstream fin(fname);
    std::string line,all;

    if (!fin) {
      fprintf(stderr,"Error: cannot open the activity file '%s'\n",fname.c_str());
      exit(EXIT_FAILURE);
    }

    while (std::getline(fin,line)) {
      size_t c=line.find('#');
      if (c!=std::string::npos) line=line.substr(0,c);
      for (auto& ch : line) if (ch==',') ch=' ';
      all+=line+" ";
    }

    std::istringstream ss(all);
    for (int s=0;s<3;s++) for (int j=0;j<25;j++) if (!(ss >> Coefficient[s][j])) {
      fprintf(stderr,"Error: the activity file '%s' must contain 3x25 coefficients\n",fname.c_str());
      exit(EXIT_FAILURE);
    }

    //heliocentric dependence of the H2O production that goes with the coefficients
    HeliocentricModel=input.GetOption("ActivityFileHeliocentricModel",{"shap5","shap5_1"},"shap5");
    try {LegacyActivityIdentity="sha256:"+ProvenanceFileSHA256(fname);}
    catch (const std::exception& exception) {
      fprintf(stderr,"Error: cannot identify the legacy activity file: %s\n",exception.what());
      exit(EXIT_FAILURE);
    }
    printf("Activity coefficients from '%s' (heliocentric model %s)\n",fname.c_str(),HeliocentricModel.c_str());
  }
  else {
    const double (*Activity)[25]=(ActivityModel=="shap5") ? srccg::Activity_SHAP5 : srccg::Activity_SHAP5_1;

    for (int s=0;s<3;s++) for (int j=0;j<25;j++) Coefficient[s][j]=Activity[s][j];
    HeliocentricModel=ActivityModel;
    LegacyActivityIdentity="builtin:"+ActivityModel;
  }

  //The enhancement is a separate H2O residual, never an extension of the
  //legacy three-species coefficient file.  Its default is strictly disabled.
  const std::string enhancementModel=input.GetOption("ActivityEnhancementModel",
      {"none",ActivityEnhancement::ModelId},"none");
  if (enhancementModel=="none") {
    //Degree 4 is the explicit zero-length control.  Consume it when supplied
    //so a runner can state the control setting without triggering UnusedKeys.
    const double requestedDegree=input.GetDouble("ActivityEnhancementDegree",4.0);
    if (!std::isfinite(requestedDegree)||requestedDegree!=4.0) {
      fprintf(stderr,"Error: ActivityEnhancementDegree must be 4 when ActivityEnhancementModel=none\n");
      exit(EXIT_FAILURE);
    }
    if (input.Has("ActivityEnhancementFile")||input.Has("ActivityEnhancementScale")) {
      //Read supplied values solely to mark them as recognized before issuing
      //the more useful incompatibility message.
      if (input.Has("ActivityEnhancementFile")) input.GetString("ActivityEnhancementFile");
      if (input.Has("ActivityEnhancementScale")) input.GetDouble("ActivityEnhancementScale");
      fprintf(stderr,"Error: ActivityEnhancementFile/Scale require ActivityEnhancementModel=%s\n",
          ActivityEnhancement::ModelId);
      exit(EXIT_FAILURE);
    }
    Enhancement=ActivityEnhancement::Configuration();
  }
  else {
    const std::string filename=CopsPaths::ProductInput(input.GetString("ActivityEnhancementFile"));
    try {Enhancement=ActivityEnhancement::ReadFile(filename);}
    catch (const std::exception& exception) {
      fprintf(stderr,"Error: %s\n",exception.what());
      exit(EXIT_FAILURE);
    }
    try {ActivityEnhancement::ValidateCompatibility(Enhancement,bodyFrame,LegacyActivityIdentity);}
    catch (const std::exception& exception) {
      fprintf(stderr,"Error: %s\n",exception.what());
      exit(EXIT_FAILURE);
    }
    const double requestedDegree=input.GetDouble("ActivityEnhancementDegree",
        static_cast<double>(Enhancement.MaximumDegree));
    if (!std::isfinite(requestedDegree)||std::floor(requestedDegree)!=requestedDegree||
        requestedDegree<4.0||requestedDegree>8.0||
        static_cast<int>(requestedDegree)!=Enhancement.MaximumDegree) {
      fprintf(stderr,"Error: ActivityEnhancementDegree must be the integer %d supplied by the enhancement file\n",
          Enhancement.MaximumDegree);
      exit(EXIT_FAILURE);
    }
    const double requestedScale=input.GetDouble("ActivityEnhancementScale",Enhancement.Scale);
    if (!std::isfinite(requestedScale)||requestedScale<=0.0||requestedScale!=Enhancement.Scale) {
      fprintf(stderr,"Error: ActivityEnhancementScale must exactly match the file value %.17g\n",Enhancement.Scale);
      exit(EXIT_FAILURE);
    }
    printf("Activity enhancement: model=%s species=H2O basis=%s maximum_degree=%d modes=%zu scale=%.17g\n",
        Enhancement.Model.c_str(),Enhancement.BasisId.c_str(),Enhancement.MaximumDegree,
        Enhancement.Coefficients.size(),Enhancement.Scale);
  }

  AnalyticalTmax=input.GetDouble("AnalyticalTmax",182.09922);
  AnalyticalT75=input.GetDouble("AnalyticalT75",136.11149);
  AnalyticalTmin=input.GetDouble("AnalyticalTmin",133.0);

  //activity factor: evaluated at the element centroid (longitude/colatitude in the body-fixed frame), negative values set to zero (as in srcCG)
  Basis.resize(25*mesh.Faces.size());

  for (size_t i=0;i<mesh.Faces.size();i++) {
    const Vec3& c=mesh.Faces[i].Centroid;
    double colatitude=acos(c.z/Length(c)),longitude,rxy=sqrt(c.x*c.x+c.y*c.y);

    if (rxy==0.0) longitude=0.0;
    else if (c.y>0.0) longitude=acos(c.x/rxy);
    else longitude=-acos(c.x/rxy);

    for (int j=0;j<25;j++) Basis[25*i+j]=SphericalHarmonic(j,colatitude,longitude);
  }

  if (Enhancement.Enabled()) {
    LegacyLinearH2O.resize(mesh.Faces.size());
    EnhancementH2O.resize(mesh.Faces.size());
    TotalLinearH2O.resize(mesh.Faces.size());
  }

  for (int s=0;s<Const::nSpecies;s++) {
    ActivityFactor[s].resize(mesh.Faces.size());

    for (size_t i=0;i<mesh.Faces.size();i++) {
      double factor=0.0;
      for (int j=0;j<25;j++) factor+=Coefficient[s][j]*Basis[25*i+j];
      if ((s==Const::H2O)&&Enhancement.Enabled()) {
        LegacyLinearH2O[i]=factor;
        ActivityFactor[s][i]=ActivityEnhancement::Apply(Enhancement,
            mesh.Faces[i].Centroid,factor,&EnhancementH2O[i],&TotalLinearH2O[i]);
      }
      else ActivityFactor[s][i]=std::max(factor,0.0);
    }
  }

  printf("Source model: activity=%s, temperature=%s, SZA=%s, cast shadows=%s\n",ActivityModel.c_str(),TemperatureModel.c_str(),SZAMode.c_str(),CastShadows ? "on" : "off");
  if (Enhancement.Enabled()) {
    const auto residualRange=std::minmax_element(EnhancementH2O.begin(),EnhancementH2O.end());
    const std::size_t clipped=std::count_if(TotalLinearH2O.begin(),TotalLinearH2O.end(),
        [](double value) {return value<0.0;});
    printf("Activity enhancement field: residual=[%.6e, %.6e], clipped facets=%zu/%zu, rule=%s\n",
        *residualRange.first,*residualRange.second,clipped,mesh.Faces.size(),Enhancement.Clipping.c_str());
  }
}

//srcCG: Comet::GetTotalProductionRateBjornNASTRAN()
void cSourceModel::GetQ(int spec,double rh_AU,double& Qmin,double& Qmax) const {
  if (spec==Const::H2O) {
    double p=(HeliocentricModel=="shap5") ? 4.2143229 : 4.8;

    Qmin=10.0*0.02/pow(rh_AU,p)*600;
    Qmax=1.0/pow(rh_AU,p)*600;
  }
  else {
    Qmin=0.1/pow(rh_AU,2.0)*600;
    Qmax=1.0/pow(rh_AU,2.0)*600;
  }
}

//srcCG: Exosphere::GetSurfaceTemperature() (Comet.cpp)
double cSourceModel::Temperature(double CosSZA) const {
  double res;

  if (TemperatureModel=="constant") {
    res=ConstantTemperature;
  }
  else if (TemperatureModel=="analytical") {
    const double Tmax=AnalyticalTmax,T75=AnalyticalT75,Tmin=AnalyticalTmin;
    double a=(Tmax-T75)/(1-1/cos(75*Const::Pi/180.0));
    double b=Tmax-a;
    double angle=acos(std::max(-1.0,std::min(1.0,CosSZA)));

    res=(angle<=75*Const::Pi/180.0) ? std::max(Tmin,a/cos(angle)+b) : Tmin;
  }
  else {
    //"bjorn": table SurfaceTemp[angle][1+distance index] (srcCG, CG-BC.h); distance index fixed by the input (srcCG: Comet::ndist=4)
    const double DistanceFromTheSun[6]={1.3,2.0,2.7,3.0,3.25,3.5};
    const double minTemp[6]={172.0,163.0,150.0,145.0,139.0,133.0};
    int nd=BjornDistanceIndex,angle;

    if (CosSZA<0.0) res=minTemp[nd];
    else if (nd<5) {
      angle=(int)(acos(CosSZA)*180.0/Const::Pi);
      if (angle>89) angle=89;
      res=std::max(srccg::SurfaceTemp[angle][nd+1],minTemp[nd]);
    }
    else {
      angle=(int)(acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[nd],2.0)*CosSZA)*180/Const::Pi);
      if (angle>89) angle=89;
      res=std::max(srccg::SurfaceTemp[angle][4+1],minTemp[nd]);
    }
  }

  return TemperatureCorrectionFactor*res;
}

void cSourceModel::EvaluateElement(int iface,double CosSZA,bool night,double HeliocentricDistance,double* Flux,double& T) const {
  double Qmin,Qmax;

  if (UniformTestFlux>0.0) {
    T=ConstantTemperature;
    for (int s=0;s<Const::nSpecies;s++) Flux[s]=UniformTestFlux;
    return;
  }

  //night side and shadowed elements: night-side production and night temperature (cos SZA=-1, as in srcCG)
  T=Temperature(night ? -1.0 : CosSZA);

  for (int s=0;s<Const::nSpecies;s++) {
    double flux;

    GetQ(s,HeliocentricDistance/Const::AU,Qmin,Qmax);

    if (night) flux=Qmin;
    else if (SZAMode=="continuous") flux=Qmin+(Qmax-Qmin)*CosSZA;
    else {
      //srcCG tabulates the flux at integer degrees and truncates the angle: fluxBjorn[(int)angle]
      int angle=(int)(acos(std::min(1.0,CosSZA))*180.0/Const::Pi);
      if (angle>89) angle=89;
      flux=Qmin+(Qmax-Qmin)*cos(angle*Const::Pi/180.0);
    }

    Flux[s]=flux*ActivityFactor[s][iface];
  }
}

void cSourceModel::Evaluate(const cMesh& mesh,const cBVH& bvh,const Vec3& SunDirection,double HeliocentricDistance,std::vector<cSurfaceState>& state) const {
  state.resize(mesh.Faces.size());

  #pragma omp parallel for schedule(dynamic,256)
  for (size_t i=0;i<mesh.Faces.size();i++) {
    const cTriangle& f=mesh.Faces[i];
    cSurfaceState& st=state[i];

    st.CosSZA=Dot(f.Normal,SunDirection);
    st.Shadow=false;

    if ((st.CosSZA>0.0)&&CastShadows) {
      //ray from the element toward the Sun (the Sun is at infinity at the nucleus scale)
      Vec3 x0=f.Centroid+f.Normal*(1.0E-6*f.MaxEdge);

      st.Shadow=bvh.AnyHit(x0,SunDirection,0.0,1.0E12,(int)i);
    }

    EvaluateElement((int)i,st.CosSZA,(st.CosSZA<=0.0)||st.Shadow,HeliocentricDistance,st.Flux,st.Temperature);
  }
}

//binary element file: int32 version(=1), int32 nFaces, int32 nBasis(=25), then for each element
//double centroid[3], normal[3], area, Y[25], ActivityFactor_H2O, ActivityFactor_CO2
//The layout and Y[25] remain legacy-compatible.  If Stage A3 is enabled, the
//H2O factor is the final clipped total; its components are kept in the separate
//ActivityEnhancementFacets.dat diagnostic rather than expanding this cache.
void cSourceModel::WriteElements(const std::string& fname,const cMesh& mesh) const {
  FILE *f=fopen(fname.c_str(),"wb");

  if (f==NULL) {
    fprintf(stderr,"Error: cannot write '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  int version=1,n=(int)mesh.Faces.size(),nb=25;

  fwrite(&version,sizeof(int),1,f);
  fwrite(&n,sizeof(int),1,f);
  fwrite(&nb,sizeof(int),1,f);

  for (int i=0;i<n;i++) {
    const cTriangle& t=mesh.Faces[i];
    double v[7+25+2]={t.Centroid.x,t.Centroid.y,t.Centroid.z,t.Normal.x,t.Normal.y,t.Normal.z,t.Area};

    for (int j=0;j<25;j++) v[7+j]=Basis[25*i+j];
    v[32]=ActivityFactor[Const::H2O][i],v[33]=ActivityFactor[Const::CO2][i];
    fwrite(v,sizeof(double),34,f);
  }

  fclose(f);
}

void cSourceModel::WriteActivityEnhancementDiagnostics(const std::string& fname,
                                                       const cMesh& mesh) const {
  if (!Enhancement.Enabled()) return;
  if (mesh.Faces.size()!=LegacyLinearH2O.size()||
      mesh.Faces.size()!=EnhancementH2O.size()||
      mesh.Faces.size()!=TotalLinearH2O.size()) {
    fprintf(stderr,"Error: internal activity-enhancement diagnostic size mismatch\n");
    exit(EXIT_FAILURE);
  }
  FILE* output=fopen(fname.c_str(),"w");
  if (output==nullptr) {
    fprintf(stderr,"Error: cannot write '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }
  fprintf(output,"# LiouvilleCOPS activity-enhancement facet diagnostic v1\n");
  fprintf(output,"# Model=%s BasisID=%s Species=H2O BodyFrame=%s MinimumDegree=5 MaximumDegree=%d Scale=%.17g\n",
      Enhancement.Model.c_str(),Enhancement.BasisId.c_str(),Enhancement.BodyFrame.c_str(),
      Enhancement.MaximumDegree,Enhancement.Scale);
  fprintf(output,"# Ordering=%s Clipping=%s LegacyActivityIdentity=%s\n",
      Enhancement.CoefficientOrdering.c_str(),Enhancement.Clipping.c_str(),
      Enhancement.LegacyActivityIdentity.c_str());
  fprintf(output,"# Units: centroid/area are m/m^2; all activity columns are dimensionless\n");
  fprintf(output,"# face x y z area legacy_linear_h2o enhancement_h2o total_linear_h2o activity_h2o clipped\n");
  for (std::size_t i=0;i<mesh.Faces.size();i++) {
    const cTriangle& face=mesh.Faces[i];
    fprintf(output,"%zu %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %d\n",
        i,face.Centroid.x,face.Centroid.y,face.Centroid.z,face.Area,
        LegacyLinearH2O[i],EnhancementH2O[i],TotalLinearH2O[i],
        ActivityFactor[Const::H2O][i],TotalLinearH2O[i]<0.0 ? 1 : 0);
  }
  if (ferror(output)||fclose(output)!=0) {
    fprintf(stderr,"Error: failed while writing '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }
}
