//$Id$

#ifndef _SEP3D_H_
#define _SEP3D_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

// SPICE is NOT used for this application: include empty substitutes
#include "SpiceEmptyDefinitions.h"


namespace SEP3D {
  // it is based on a generic application named Exosphere
  using namespace Exosphere;

  //--------------------------------------------------------------------------
  //particle tracking condition
  namespace ParticleTracker {
    inline bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
      //only those particles are traced, which are close to the magnetic island
      if (x[0]<9.0E8||x[2]>0.8E7||x[2]<-0.8E7) return false;

      return PIC::ParticleTracker::TrajectoryTrackingCondition_default(x,v,spec,ParticleData);
    }
  }

  //--------------------------------------------------------------------------
  // specific initialization procedures
  void Init_BeforeParser();

  //--------------------------------------------------------------------------
  // methods for the reaction processing specific for the application
  namespace Physics {
    //------------------------------------------------------------------------
    // particle mover
    int Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, bool FirstBoundaryFlag);
    int Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    //------------------------------------------------------------------------
    // self-explanatory
    void inline TotalParticleAcceleration(double *accl, int spec, long int ptr,
					  double *x, double *v, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
      //......................................................................
      // local values of position, velocity and acceleration
      double x_LOCAL[3], v_LOCAL[3], accl_LOCAL_CLASSIC[3]={0.0};
      double accl_LOCAL_RELATIVISTIC[3]={0.0};
      memcpy(x_LOCAL,x,3*sizeof(double));
      memcpy(v_LOCAL,v,3*sizeof(double));

      //......................................................................
      // calculate Lorentz force if needed
//#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
      // find electro-magnetic field
      double E[3],B[3];
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      // if no coupler is used -> use characteristic values
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere::swB_Typical,3*sizeof(double));
#else 
      // if coupler is used -> get values from it
      PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
      PIC::CPLR::GetBackgroundElectricField(E);
      PIC::CPLR::GetBackgroundMagneticField(B);
#endif//_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 

      //......................................................................
      // calculate acceleraton due to Lorentz force
      double ElCharge=PIC::MolecularData::GetElectricCharge(spec);
      double mass=PIC::MolecularData::GetMass(spec);
      accl_LOCAL_CLASSIC[0]+=ElCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
      accl_LOCAL_CLASSIC[1]+=ElCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
      accl_LOCAL_CLASSIC[2]+=ElCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;
//#endif//_FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

      //......................................................................
      // calculate gravitational force if needed
//#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
      // find the distance to the Sun's center, here it's {0, 0, 0} 
      double r, r2;
      r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
      r=sqrt(r2);
      // calculate acceleration due to gravity
      for (int idim=0;idim<DIM;idim++) {
	accl_LOCAL_CLASSIC[idim]-=GravityConstant*_SUN__MASS_/r2*x_LOCAL[idim]/r;
      }
//#endif//_FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_

      //......................................................................
      // convert classic acceleration to relativistic
      double speed2=0.0, VelDotAccl=0.0;
      double c2=SpeedOfLight*SpeedOfLight;
      for(int idim=0; idim<DIM;idim++){
	speed2    += v_LOCAL[idim]*v_LOCAL[idim];
	VelDotAccl+= v_LOCAL[idim]*accl_LOCAL_CLASSIC[idim]; 
      }
      if(speed2 > c2) exit(__LINE__,__FILE__,"Error: superluminous particle");
      double RelGamma = pow(1.0 - speed2/c2, -0.5);
      for(int idim=0; idim<DIM;idim++)
	accl_LOCAL_RELATIVISTIC[idim]= 
	  (accl_LOCAL_CLASSIC[idim]-VelDotAccl*v_LOCAL[idim]/c2)/RelGamma;

      //copy the local values of the acceleration to the global ones
      memcpy(accl,accl_LOCAL_RELATIVISTIC,3*sizeof(double));
    }
  }
}

#endif//_SEP3D_H_
