//COPS gauge response: angular acceptance and conversion of density/flux into the gauge pressure
//
//Notation: l = unit direction of the molecule velocity (from the surface element toward the gauge),
//          b = unit boresight of the gauge (body-fixed frame). A molecule arrives from the front of the gauge when l.b<0.
#ifndef GAUGES_H
#define GAUGES_H

#include <string>
#include <vector>

#include "Input.h"
#include "Geometry.h"
#include "Constants.h"

class cGaugeModel {
public:
  void Init(const cInput& input);

  //weight with which a molecule moving in direction l contributes to the NG density (acceptance x optional head-on weight)
  double NGWeight(const Vec3& l,const cObservationGeometry& g) const;

  //acceptance of the RG aperture for a molecule moving in direction l (the flux projection factor -l.b is applied separately)
  double RGAcceptance(const Vec3& l,const cObservationGeometry& g) const;

  //gauge pressures [Pa] from the density/flux of a species
  double NGPressure(int spec,double n) const;
  double RGPressure(int spec,double PhiIn,double TemperatureRG) const;

  //RG temperature at the ephemeris time et: constant (TemperatureRG) or interpolated in the housekeeping table
  double TemperatureRGAt(double et) const;
  void LoadTemperatureRGTable(const std::string& fname,int column);   //requires the SPICE leapseconds kernel

  std::string NGAcceptanceModel,RGAcceptanceModel,RGPressureModel,TemperatureRGMode;
  bool NGHeadOnWeight,SpeciesSensitivity;
  double TemperatureNG,TemperatureRG;
  double Sensitivity[Const::nSpecies];   //p = n k T / Sensitivity (srcCG convention: BetaFactor)

  void PrintSummary() const;

private:
  //srcCG nude gauge shielding model (RosinaMeasurements_Liouville.cpp): the fraction of the ionization volume
  //(cylinder of radius Rmax and height Hmax along the boresight) that sees a molecule arriving from behind
  //the gauge through the opening of the tube of radius Rtube; tabulated vs. mu=l.b in [0,1]
  std::vector<double> TableEt,TableT;   //RG temperature table (time of a row = middle of its measurement interval)
  double TableStartEt,TableStopEt;      //start of the first and end of the last measurement interval

  double ShieldRmax,ShieldHmax,ShieldRtube;
  std::vector<double> ShieldTable;
  void BuildShieldTable(int nBins,int nSamples);
  double ShieldFraction(double mu) const;
};

#endif
