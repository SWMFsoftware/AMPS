//surface source model: production rate and surface temperature (models of srcCG)
#ifndef SOURCE_H
#define SOURCE_H

#include <string>
#include <vector>

#include "Constants.h"
#include "Input.h"
#include "Mesh.h"
#include "BVH.h"

//state of a surface element for a given position of the Sun
struct cSurfaceState {
  double CosSZA;                     //cosine of the solar zenith angle (local normal vs. Sun direction)
  bool Shadow;                       //the element faces the Sun but is shadowed by another part of the nucleus
  double Temperature;                //surface (gas emission) temperature [K]
  double Flux[Const::nSpecies];      //emitted molecular flux [m^-2 s^-1]
};

class cSourceModel {
public:
  void Init(const cInput& input,const cMesh& mesh);

  //evaluate the flux and temperature of all surface elements
  void Evaluate(const cMesh& mesh,const cBVH& bvh,const Vec3& SunDirection,double HeliocentricDistance,std::vector<cSurfaceState>& state) const;

  //flux and temperature of the element 'iface' for a given cos(SZA); 'night' marks night side or shadow
  void EvaluateElement(int iface,double CosSZA,bool night,double HeliocentricDistance,double* Flux,double& Temperature) const;

  //element file for the fit (tools/fit_cops.py): centroid, normal, area, spherical-harmonic basis, activity factors
  void WriteElements(const std::string& fname,const cMesh& mesh) const;

  std::string TemperatureModel,ActivityModel,SZAMode,HeliocentricModel;
  double AnalyticalTmax,AnalyticalT75,AnalyticalTmin;   //TemperatureModel=analytical
  double Coefficient[3][25];                           //activity coefficients in use (rows: H2O, CO2, CO)

  //verification only: when >0, every element emits this flux of every species at ConstantTemperature
  double UniformTestFlux=0.0;
  bool CastShadows;

private:
  std::vector<double> ActivityFactor[Const::nSpecies];   //spherical-harmonic activity factor of each element (>=0)
  std::vector<double> Basis;                             //Y_j of each element (nFaces x 25)

  //temperature models
  double ConstantTemperature,TemperatureCorrectionFactor;
  int BjornDistanceIndex;
  double Temperature(double CosSZA) const;

  //production rate of srcCG (Comet::GetTotalProductionRateBjornNASTRAN): Qmin, Qmax at the heliocentric distance
  void GetQ(int spec,double rh_AU,double& Qmin,double& Qmax) const;
};

#endif
