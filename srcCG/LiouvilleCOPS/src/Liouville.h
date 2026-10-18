//Collisionless (Liouville theorem) evaluation of the gas density and flux at the COPS gauges
//
//Each surface element emits molecules with the speed distribution of a flux-weighted Maxwellian at the local
//surface temperature T,  h(v)=2 beta^4 v^3 exp(-beta^2 v^2),  beta=sqrt(m/(2kT)),  int h dv=1,
//with the angular distribution of the emitted flux (per unit element area and solid angle)
//   flat  : dF/dOmega = F cos(theta)/pi     (Lambertian; isotropic half-Maxwellian velocity distribution)
//   rough : dF/dOmega = F/(2 pi)            (uniform over the hemisphere defined by the external normal)
//where F is the flux per unit (planar) element area and theta the emission angle from the external normal.
//Without collisions and with straight trajectories, the distribution function is conserved along the
//line of sight, which gives for a visible element of area dS at distance r:
//   number flux at the gauge   dPhi = (dF/dOmega) dS/r^2
//   number density             dn   = (dF/dOmega) dS/r^2 * int h(v)/v dv = (dF/dOmega) dS/r^2 * sqrt(pi) beta/2
//The geometric part does not depend on F and T, so it is stored per element ("kernels"):
//   Kn   = sum (dS/r^2) A(theta) W_NG(l)
//   KPhi = sum (dS/r^2) A(theta) Acc_RG(l) (-l.b_RG)
//   KV   = sum (dS/r^2) A(theta) Acc_RG(l)
//with A=cos(theta)/pi (flat) or 1/(2 pi) (rough); the sums run over the visible parts of the element. Then
//   n_s    = sum_f F_sf beta_sf sqrt(pi)/2 Kn_f
//   Phi_s  = sum_f F_sf [ KPhi_f + (V.b_RG) beta_sf sqrt(pi)/2 KV_f ]
//The second term in Phi_s is the exact correction for the spacecraft velocity V relative to the nucleus
//(flux into the aperture of a molecule with velocity v l in the spacecraft frame: -(v l - V).b).
#ifndef LIOUVILLE_H
#define LIOUVILLE_H

#include <vector>
#include <string>

#include "Mesh.h"
#include "BVH.h"
#include "Gauges.h"
#include "Source.h"
#include "Geometry.h"

struct cElementKernel {
  int Face;
  double Kn,KPhi,KV;
  double KnBelowDeck;                    //part of Kn from directions below the +Z deck plane
};

//FlightTimeRotation=exact: per-species contributions of an element per unit emitted flux (flux at the reception time)
struct cElementKernelExact {
  int Face;
  double Gn[Const::nSpecies];      //NG-weighted density per unit flux [m^-3 / (m^-2 s^-1)]
  double GPhi[Const::nSpecies];    //RG aperture flux per unit flux (including the spacecraft velocity term)
};

struct cLiouvilleResult {
  double Altitude;                       //closest distance from the gauge to the surface [m]
  bool GaugeInsideNucleus;
  int nVisibleElements;
  double VdotB_RG;                       //spacecraft velocity component along the RG boresight [m/s]
  double NGBelowDeckFraction;            //fraction of the NG density carried by molecules arriving from below the +Z deck plane
                                         //of the spacecraft (l.Z_sc>0), where the bus, solar arrays and HGA are located
  double n[Const::nSpecies];             //NG-weighted number density [m^-3]
  double PhiIn[Const::nSpecies];         //flux into the RG aperture [m^-2 s^-1]
  double pNG[Const::nSpecies],pRG[Const::nSpecies];
  double TemperatureRG;                  //RG temperature used in the flux->pressure conversion [K]
  double TotalProduction[Const::nSpecies];  //total production rate of the nucleus [s^-1]
  std::vector<cElementKernel> Kernels;
  std::vector<cElementKernelExact> KernelsExact;
};

class cLiouville {
public:
  void Init(const cInput& input);

  //geometric kernels for the gauge location/orientation in 'g'
  void EvaluateKernels(const cMesh& mesh,const cBVH& bvh,const cGaugeModel& gauges,const cObservationGeometry& g,cLiouvilleResult& res) const;

  //apply the source model (flux, temperature) to the kernels
  void ApplySource(const cMesh& mesh,const cGaugeModel& gauges,const std::vector<cSurfaceState>& state,double TemperatureRG,cLiouvilleResult& res) const;

  //FlightTimeRotation=exact: rotation of the nucleus during the flight of the molecules (see Liouville.cpp);
  //evaluates the densities, fluxes and pressures directly (the geometry depends on the molecular speed)
  void EvaluateExact(const cMesh& mesh,const cBVH& bvh,const cGaugeModel& gauges,const cSourceModel& source,
      const std::vector<cSurfaceState>& state,const cObservationGeometry& g,double TemperatureRG,cLiouvilleResult& res) const;

  std::string SurfaceBC,FlightTimeRotation;
  int nSpeedNodes;
  double PathTolerance;
  std::vector<double> SpeedNode,SpeedWeight;   //Gauss-Legendre rule for x=beta*v on [0,SpeedRangeMax]
  bool SpacecraftVelocity;
  double SubdivisionAccuracy;
  int MaxSubdivisionLevel;
};

//true when the point is inside the closed surface (ray parity test)
bool PointInsideSurface(const cMesh& mesh,const cBVH& bvh,const Vec3& x);

#endif
