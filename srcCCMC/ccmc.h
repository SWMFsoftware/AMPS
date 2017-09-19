
#ifndef CCMC_H
#define CCMC_H

//$Id$

#include <cmath>

#include "ccmc.dfn"

#include "pic.h"
#include "Exosphere.h"
#include "constants.h"

namespace CCMC {
  using namespace Exosphere;

  //the size of the computational domain
  namespace Domain {
    static const double xmin[]={0.0,0.0,0.0};
    static const double xmax[]={0.0,0.0,0.0};
  }

  //the condition of the particle trajectory tracking
  namespace ParticleTracker {
    inline bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {

/*      //only those solar wind ions are traced, which trajectories are close to the surface of Mercury
      if (spec==_H_PLUS_SPEC_) {
        if (x[1]*x[1]+x[2]*x[2]>pow(2.0*_RADIUS_(_TARGET_),2)) return false;
      }*/

      return PIC::ParticleTracker::TrajectoryTrackingCondition_default(x,v,spec,ParticleData);
    }
  }

  //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

/*
    memcpy(accl,accl_LOCAL,3*sizeof(double));
    return;
*/


    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

    //the Lorentz force
    double elCharge;

    if ((elCharge=PIC::MolecularData::GetElectricCharge(spec))>0.0) {
      int i,j,k;
      double E[3],B[3];

      if (PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false)==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

      #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
      #endif


#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere_swB_Typical,3*sizeof(double));
#else 
      startNode=PIC::Mesh::Search::FindBlock(x_LOCAL);
      PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
      PIC::CPLR::GetBackgroundElectricField(E);
      PIC::CPLR::GetBackgroundMagneticField(B);
#endif

      elCharge/=PIC::MolecularData::GetMass(spec);

      accl_LOCAL[0]+=elCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1]);
      accl_LOCAL[1]+=elCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0]);
      accl_LOCAL[2]+=elCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0]);

    }


    if (std::isnan(accl_LOCAL[0])||std::isnan(accl_LOCAL[1])||std::isnan(accl_LOCAL[2])) exit(__LINE__,__FILE__,"Error in calculation of the acceleration");


/*    //account for the planetary rotation around the Sun
    double aCen[3],aCorr[3],t3,t7,t12;

    t3 = RotationVector_SO_FROZEN[0] * x_LOCAL[1] - RotationVector_SO_FROZEN[1] * x_LOCAL[0];
    t7 = RotationVector_SO_FROZEN[2] * x_LOCAL[0] - RotationVector_SO_FROZEN[0] * x_LOCAL[2];
    t12 = RotationVector_SO_FROZEN[1] * x_LOCAL[2] - RotationVector_SO_FROZEN[2] * x_LOCAL[1];

    aCen[0] = -RotationVector_SO_FROZEN[1] * t3 + RotationVector_SO_FROZEN[2] * t7;
    aCen[1] = -RotationVector_SO_FROZEN[2] * t12 + RotationVector_SO_FROZEN[0] * t3;
    aCen[2] = -RotationVector_SO_FROZEN[0] * t7 + RotationVector_SO_FROZEN[1] * t12;


    aCorr[0] = -2.0*(RotationVector_SO_FROZEN[1] * v_LOCAL[2] - RotationVector_SO_FROZEN[2] * v_LOCAL[1]);
    aCorr[1] = -2.0*(RotationVector_SO_FROZEN[2] * v_LOCAL[0] - RotationVector_SO_FROZEN[0] * v_LOCAL[2]);
    aCorr[2] = -2.0*(RotationVector_SO_FROZEN[0] * v_LOCAL[1] - RotationVector_SO_FROZEN[1] * v_LOCAL[0]);

    accl_LOCAL[0]+=aCen[0]+aCorr[0];
    accl_LOCAL[1]+=aCen[1]+aCorr[1];
    accl_LOCAL[2]+=aCen[2]+aCorr[2];*/

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }

}


#endif
