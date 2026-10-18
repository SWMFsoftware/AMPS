#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>
#include <fstream>
#include <sstream>

#include "SpiceUsr.h"

#include "Gauges.h"
#include "Paths.h"

void cGaugeModel::Init(const cInput& input) {
  NGAcceptanceModel=input.GetOption("NGAcceptance",{"srccg_shield","hemisphere","ik_fov","full"},"srccg_shield");
  RGAcceptanceModel=input.GetOption("RGAcceptance",{"hemisphere","ik_fov"},"hemisphere");
  RGPressureModel=input.GetOption("RGPressureModel",{"equilibrium_sphere","calibration_srccg","dsmc_srccg"},"calibration_srccg");

  NGHeadOnWeight=input.GetBool("NGHeadOnWeight",true);
  SpeciesSensitivity=input.GetBool("SpeciesSensitivity",true);
  TemperatureNG=input.GetDouble("TemperatureNG",293.0);
  TemperatureRG=input.GetDouble("TemperatureRG",293.0);
  TemperatureRGMode=input.GetOption("TemperatureRGMode",{"constant","hk_table"},"constant");
  Sensitivity[Const::H2O]=input.GetDouble("SensitivityH2O",0.893);
  Sensitivity[Const::CO2]=input.GetDouble("SensitivityCO2",0.704);

  ShieldRmax=input.GetDouble("NGShieldRmax",10.3);
  ShieldHmax=input.GetDouble("NGShieldHmax",33.0);
  ShieldRtube=input.GetDouble("NGShieldRtube",20.5);

  if (TemperatureRGMode=="hk_table") LoadTemperatureRGTable(CopsPaths::ProductInput(input.GetString("TemperatureRGTable")),input.GetInt("TemperatureRGTableColumn",6));

  if (NGAcceptanceModel=="srccg_shield") BuildShieldTable(input.GetInt("NGShieldTableBins",2001),input.GetInt("NGShieldTableSamples",20000));
}

void cGaugeModel::PrintSummary() const {
  printf("Gauges: NG acceptance=%s, head-on weight=%s, T_NG=%g K; RG acceptance=%s, RG pressure=%s, T_RG=%s; species sensitivity=%s (H2O %g, CO2 %g)\n",
      NGAcceptanceModel.c_str(),NGHeadOnWeight ? "on" : "off",TemperatureNG,RGAcceptanceModel.c_str(),RGPressureModel.c_str(),
      (TemperatureRGMode=="constant") ? (std::to_string(TemperatureRG)+" K").c_str() : "housekeeping table",
      SpeciesSensitivity ? "on" : "off",Sensitivity[0],Sensitivity[1]);
}

//Monte Carlo evaluation identical to srcCG: random point in the ionization cylinder (uniform in the disk area and height);
//follow the direction toward the source (-l) back to the base plane of the cylinder; the point is not shielded when
//the crossing is outside the tube radius
void cGaugeModel::BuildShieldTable(int nBins,int nSamples) {
  std::mt19937_64 rng(20160930);
  std::uniform_real_distribution<double> U(0.0,1.0);

  ShieldTable.resize(nBins);

  for (int ib=0;ib<nBins;ib++) {
    double mu=(double)ib/(nBins-1);     //mu=l.b>=0: molecule arrives from behind the gauge
    double sinth=sqrt(std::max(0.0,1.0-mu*mu));
    //direction toward the source in the gauge frame (e2 = boresight); azimuth irrelevant (axial symmetry)
    double sx=sinth,sz=-mu;
    int cnt=0;

    for (int i=0;i<nSamples;i++) {
      double r=sqrt(U(rng))*ShieldRmax,phi=2.0*Const::Pi*U(rng);
      double px=r*cos(phi),py=r*sin(phi),pz=U(rng)*ShieldHmax;

      if (sz>=0.0) {cnt++; continue;}   //grazing direction (mu=0): never crosses the base plane

      double t=-pz/sz;
      double cx=px+sx*t,cy=py;

      if (cx*cx+cy*cy>ShieldRtube*ShieldRtube) cnt++;
    }

    ShieldTable[ib]=(double)cnt/nSamples;
  }
}

double cGaugeModel::ShieldFraction(double mu) const {
  double x=std::max(0.0,std::min(1.0,mu))*(ShieldTable.size()-1);
  int i=std::min((int)x,(int)ShieldTable.size()-2);
  double w=x-i;

  return (1.0-w)*ShieldTable[i]+w*ShieldTable[i+1];
}

double cGaugeModel::NGWeight(const Vec3& l,const cObservationGeometry& g) const {
  double mu=Dot(l,g.NG.Boresight);
  double acceptance;

  if (NGAcceptanceModel=="full") acceptance=1.0;
  else if (NGAcceptanceModel=="hemisphere") acceptance=(mu<0.0) ? 1.0 : 0.0;
  else if (NGAcceptanceModel=="ik_fov") acceptance=InsideFOV(g.NG,-l) ? 1.0 : 0.0;
  else {
    //srccg_shield: full acceptance from the front; from behind only through the shield opening and
    //(as in srcCG) only for molecules that arrive from the front of the RG
    if (mu<0.0) acceptance=1.0;
    else acceptance=(Dot(l,g.RG.Boresight)<0.0) ? ShieldFraction(mu) : 0.0;
  }

  //optional head-on weight of srcCG ("GasFluxQuality"): 1-l.b for molecules arriving from the front
  if (NGHeadOnWeight&&(mu<0.0)) acceptance*=1.0-mu;

  return acceptance;
}

double cGaugeModel::RGAcceptance(const Vec3& l,const cObservationGeometry& g) const {
  if (Dot(l,g.RG.Boresight)>=0.0) return 0.0;  //the molecule must move into the aperture
  if (RGAcceptanceModel=="ik_fov") return InsideFOV(g.RG,-l) ? 1.0 : 0.0;
  return 1.0;
}

double cGaugeModel::NGPressure(int spec,double n) const {
  double p=n*Const::kBoltzmann*TemperatureNG;

  return SpeciesSensitivity ? p/Sensitivity[spec] : p;
}

//table of the RG temperature: column 'column' (1-based) of the housekeeping table written by tools/cops_hk_table.py;
//columns 1-3: UTC of the row (middle of the interval), start and stop of the measurement interval; '#' starts a comment
void cGaugeModel::LoadTemperatureRGTable(const std::string& fname,int column) {
  std::ifstream fin(fname);
  std::string line;

  if (!fin) {
    fprintf(stderr,"Error: cannot open the RG temperature table '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  while (std::getline(fin,line)) {
    if (line.empty()||(line[0]=='#')) continue;

    std::istringstream ss(line);
    std::vector<std::string> f;
    std::string w;

    while (ss >> w) f.push_back(w);

    if ((int)f.size()<column) {
      fprintf(stderr,"Error: the RG temperature table '%s' has fewer than %i columns\n",fname.c_str(),column);
      exit(EXIT_FAILURE);
    }

    SpiceDouble et,et0,et1;
    utc2et_c(f[0].c_str(),&et);
    utc2et_c(f[1].c_str(),&et0);
    utc2et_c(f[2].c_str(),&et1);

    if (TableEt.empty()) TableStartEt=et0;
    TableStopEt=et1;

    if (failed_c()) {
      fprintf(stderr,"Error: cannot convert the time '%s' in '%s'\n",f[0].c_str(),fname.c_str());
      exit(EXIT_FAILURE);
    }

    TableEt.push_back(et);
    TableT.push_back(atof(f[column-1].c_str()));
  }

  for (size_t i=1;i<TableEt.size();i++) if (TableEt[i]<=TableEt[i-1]) {
    fprintf(stderr,"Error: the RG temperature table '%s' is not strictly ordered in time (row %zu)\n",fname.c_str(),i+1);
    exit(EXIT_FAILURE);
  }

  printf("RG temperature table '%s' (column %i): %zu rows\n",fname.c_str(),column,TableEt.size());
}

double cGaugeModel::TemperatureRGAt(double et) const {
  if (TemperatureRGMode=="constant") return TemperatureRG;

  if ((et<TableStartEt)||(et>TableStopEt)) {
    fprintf(stderr,"Error: time outside of the measurement intervals of the RG temperature table (et=%f)\n",et);
    exit(EXIT_FAILURE);
  }

  //within the first/last measurement interval but outside the range of the row times: the measured value of that interval
  if (et<=TableEt.front()) return TableT.front();
  if (et>=TableEt.back()) return TableT.back();

  size_t i=std::upper_bound(TableEt.begin(),TableEt.end(),et)-TableEt.begin();
  if (i==TableEt.size()) i--;
  if (i==0) i=1;

  double w=(et-TableEt[i-1])/(TableEt[i]-TableEt[i-1]);
  return (1.0-w)*TableT[i-1]+w*TableT[i];
}

//PhiIn: molecular flux into the RG aperture [m^-2 s^-1]
double cGaugeModel::RGPressure(int spec,double PhiIn,double TemperatureRG) const {
  double m=Const::SpeciesMass[spec],kT=Const::kBoltzmann*TemperatureRG,p;

  if (RGPressureModel=="equilibrium_sphere") {
    //flux balance of an equilibrium sphere at T_RG: PhiIn=n_s*v_s/4, v_s=sqrt(8kT/(pi m)) -> p=n_s kT=PhiIn*sqrt(2 pi m kT)
    p=PhiIn*sqrt(2.0*Const::Pi*m*kT);
  }
  else if (RGPressureModel=="calibration_srccg") {
    //srcCG PressureCalculationRG_ModeCalibrationPressure: p=PhiIn*sqrt(pi m kT/2)
    p=PhiIn*sqrt(Const::Pi*m*kT/2.0);
  }
  else {
    //srcCG DSMC output (RosinaMeasurements.cpp): p=kT*pi*beta/2*PhiIn, beta=sqrt(m/(2kT))
    p=kT*Const::Pi*sqrt(m/(2.0*kT))/2.0*PhiIn;
  }

  return SpeciesSensitivity ? p/Sensitivity[spec] : p;
}
