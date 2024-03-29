//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

#ifndef MERCURY_H_
#define MERCURY_H_

#include <vector>

#include "Exosphere.h"

namespace Moon {
  using namespace Exosphere;

  extern bool UseKaguya;

  //electron impact ionozation probability 
  double ElectronImpactIonizationRate(PIC::ParticleBuffer::byte* ParticleData,int& ResultSpeciesIndex,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //iteration counter
  extern int nIterationCounter;

  //init the model
  void Init_AfterParser();

  //output the column integrals in the anti-solar direction
  namespace AntiSolarDirectionColumnMap {
    const int nAzimuthPoints=150;
    const double maxZenithAngle=Pi/4.0;
    const double dZenithAngleMin=0.001*maxZenithAngle,dZenithAngleMax=0.02*maxZenithAngle;


    void Print(int DataOutputFileNumber);
  }



  //sampling procedure uniqueck for the model of the lunar exosphere
  namespace Sampling {
    using namespace Exosphere::Sampling;

    //sample the velocity distribution in the anti-sunward direction
    namespace VelocityDistribution {
      static const int nVelocitySamplePoints=500;
      static const double maxVelocityLimit=20.0E3;
      static const double VelocityBinWidth=2.0*maxVelocityLimit/nVelocitySamplePoints;

      static const int nAzimuthPoints=150;
      static const double maxZenithAngle=Pi/4.0;
      static const double dZenithAngleMin=0.001*maxZenithAngle,dZenithAngleMax=0.02*maxZenithAngle;

      extern int nZenithPoints;

      class cVelocitySampleBuffer {
      public:
        double VelocityLineOfSight[PIC::nTotalSpecies][nVelocitySamplePoints];
        double VelocityRadialHeliocentric[PIC::nTotalSpecies][nVelocitySamplePoints];
        double Speed[PIC::nTotalSpecies][nVelocitySamplePoints];
        double meanVelocityLineOfSight[PIC::nTotalSpecies],meanVelocityRadialHeliocentric[PIC::nTotalSpecies],ColumnDensityIntegral[PIC::nTotalSpecies],meanSpeed[PIC::nTotalSpecies];
        double Brightness[PIC::nTotalSpecies];


        SpiceDouble lGSE[6];

        cVelocitySampleBuffer() {
          for (int s=0;s<PIC::nTotalSpecies;s++) {
            meanVelocityLineOfSight[s]=0.0,meanVelocityRadialHeliocentric[s]=0.0,ColumnDensityIntegral[s]=0.0,meanSpeed[s]=0.0,Brightness[s]=0.0;

            for (int n=0;n<nVelocitySamplePoints;n++) VelocityLineOfSight[s][n]=0.0,VelocityRadialHeliocentric[s][n]=0.0,Speed[s][n]=0.0;
          }

          for (int i=0;i<6;i++) lGSE[i]=0.0;
        }
      };

      extern cVelocitySampleBuffer *SampleBuffer;
      extern int nTotalSampleDirections;

      void Init();
      void Sampling();
      void OutputSampledData(int DataOutputFileNumber);
    }


    //calcualte the integrals alonf the direction of observation of TVIS instrument on Kaguya
    namespace Kaguya {

      namespace TVIS {
        struct cTvisOrientation {
          double Epoch;
          char UTC[100];
          double EquatorialCrossingTime;
          double fovRightAscension__J2000;
          double fovDeclination__J2000;
          double fov[3];
          double scPosition[3];
        };

        //scan continuos observations
        struct cTvisOrientationListElement {
          cTvisOrientation *TvisOrientation;
          int nTotalTvisOrientationElements;
        };

        extern vector<cTvisOrientationListElement> TvisOrientationVector;

        //scan individual points
        extern cTvisOrientation individualPointsTvisOrientation[];
        extern int nTotalIndividualPointsTvisOrientationElements;

        void OutputModelData(int DataOutputFileNumber);
        void Init();
      }

      inline void Init () {
        if (!UseKaguya) return;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
        SpiceDouble et,lt,State[6],xform[6][6];
        SpiceDouble lJ2000[6]={0.0,0.0,0.0,0.0,0.0,0.0};
#endif

        SpiceDouble lLSO[6]={1.0,0.0,0.0,0.0,0.0,0.0};

        if (PIC::ThisThread==0) cout << "$PREFIX: Moon::Sampling::Kaguya::TVIS - init the model" << endl;

        TVIS::Init();

        for (unsigned int DataSetCnt=0;DataSetCnt<1+TVIS::TvisOrientationVector.size();DataSetCnt++) {
          unsigned int DataSetLength;
          TVIS::cTvisOrientation *DataSet;

          if (DataSetCnt<TVIS::TvisOrientationVector.size()) DataSetLength=TVIS::TvisOrientationVector[DataSetCnt].nTotalTvisOrientationElements,DataSet=TVIS::TvisOrientationVector[DataSetCnt].TvisOrientation;
          else DataSetLength=TVIS::nTotalIndividualPointsTvisOrientationElements,DataSet=TVIS::individualPointsTvisOrientation;


          for (unsigned int i=0;i<DataSetLength;i++,DataSet++) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
            utc2et_c(DataSet->UTC,&et);
            spkezr_c("Kaguya",et,"LSO","none","Moon",State,&lt);

            for (int idim=0;idim<3;idim++) DataSet->scPosition[idim]=1.0E3*State[idim];

            DataSet->Epoch=et;

            //calcualte the pointing direction of the instrument
            lJ2000[0]=cos(DataSet->fovRightAscension__J2000)*cos(DataSet->fovDeclination__J2000);
            lJ2000[1]=sin(DataSet->fovRightAscension__J2000)*cos(DataSet->fovDeclination__J2000);
            lJ2000[2]=sin(DataSet->fovDeclination__J2000);

            //get pointing direction in the LSO frame
            sxform_c ("J2000","LSO",et,xform);
            mxvg_c(xform,lJ2000,6,6,lLSO);
#endif

            memcpy(DataSet->fov,lLSO,3*sizeof(double));
          }
        }
      }

    }

    //calculate the column integrals along the limb direction
    namespace SubsolarLimbColumnIntegrals {
      //sample values of the column integrals at the subsolar limb as a function of the phase angle
      const int nPhaseAngleIntervals=50;
      const double dPhaseAngle=Pi/nPhaseAngleIntervals;

      //sample the altitude variatin of the column integral in the solar direction from the limb
      const double maxAltitude=50.0; //altititude in radii of the object
      const double dAltmin=0.001,dAltmax=0.001*maxAltitude; //the minimum and maximum resolution along the 'altitude' line
      extern double rr;
      extern long int nSampleAltitudeDistributionPoints;

      //offsets in the sample for individual parameters of separate species
      extern int _NA_EMISSION_5891_58A_SAMPLE_OFFSET_,_NA_EMISSION_5897_56A_SAMPLE_OFFSET_,_NA_COLUMN_DENSITY_OFFSET_;

/*
      struct cSampleAltitudeDistrubutionBufferElement {
        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
      };

      extern cSampleAltitudeDistrubutionBufferElement *SampleAltitudeDistrubutionBuffer;
*/


      extern SpiceDouble etSampleBegin;
      extern int SamplingPhase;
      extern int firstPhaseRadialVelocityDirection;

      extern int nOutputFile;

      class cSampleBufferElement {
      public:
//        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
//        cSampleAltitudeDistrubutionBufferElement *AltitudeDistrubutionBuffer;

        double *SampleColumnIntegrals;
        double **AltitudeDistributionColumnIntegrals;

        int nSamples;
        double JulianDate;

        cSampleBufferElement() {
          SampleColumnIntegrals=NULL,AltitudeDistributionColumnIntegrals=NULL,nSamples=0,JulianDate=0;
        }
      };

      //sample buffer for the current sampling cycle
      extern cSampleBufferElement SampleBuffer_AntiSunwardMotion[nPhaseAngleIntervals];
      extern cSampleBufferElement SampleBuffer_SunwardMotion[nPhaseAngleIntervals];

      //sample buffer for the lifetime of the simulation
      extern cSampleBufferElement SampleBuffer_AntiSunwardMotion__TotalModelRun[nPhaseAngleIntervals];
      extern cSampleBufferElement SampleBuffer_SunwardMotion__TotalModelRun[nPhaseAngleIntervals];

      void EmptyFunction();

      void init();
      void CollectSample(int DataOutputFileNumber); //the function will ba called by that part of the core that prints output files
      void PrintDataFile();
    }

  }


  //check if the point is in the Earth shadow: the function returns 'true' if it is in the shadow, and 'false' if the point is outside of the Earth shadow
  bool inline EarthShadowCheck(double *x_LOCAL_SO_OBJECT) {
/*    double lPerp[3],lSun2=0.0,c=0.0,xSun_LOCAL[3],xEarth_LOCAL[3];
    int idim;

    memcpy(xSun_LOCAL,xSun_SO,3*sizeof(double));
    memcpy(xEarth_LOCAL,xEarth_SO,3*sizeof(double));

    for (idim=0;idim<3;idim++) {
      xSun_LOCAL[idim]-=x_LOCAL_SO_OBJECT[idim];
      xEarth_LOCAL[idim]-=x_LOCAL_SO_OBJECT[idim];

      lSun2+=pow(xSun_LOCAL[idim],2);
      c+=xSun_LOCAL[idim]*xEarth_LOCAL[idim];
    }

    for (idim=0;idim<3;idim++) {
      lPerp[idim]=xEarth_LOCAL[idim]-c*xSun_LOCAL[idim]/lSun2;
    }

    return (lPerp[0]*lPerp[0]+lPerp[1]*lPerp[1]+lPerp[2]*lPerp[2]<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) ? true : false;*/

    int idim;
    double e0[3],c,lParallel=0.0,lPerp2=0.0;

    for (c=0.0,idim=0;idim<3;idim++) {
      e0[idim]=xSun_SO[idim]-xEarth_SO[idim];
      c+=e0[idim]*e0[idim];
    }

    for (c=sqrt(c),idim=0;idim<3;idim++) {
      e0[idim]/=c;
      lParallel+=(x_LOCAL_SO_OBJECT[idim]-xEarth_SO[idim])*e0[idim];
    }

    if (lParallel>0.0) return false; //the point is between the Earth and the Sun

    for (idim=0;idim<3;idim++) {
      lPerp2+=pow(x_LOCAL_SO_OBJECT[idim]-xEarth_SO[idim]-lParallel*e0[idim],2);
    }

    return (lPerp2<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) ? true : false;
  }


    //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};



/*    //Test: no acceleration:
    memcpy(accl,accl_LOCAL,3*sizeof(double));
    return;*/


    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));


    //get the radiation pressure acceleration
    if (spec==_NA_SPEC_) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
      if ( ((x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x_LOCAL[0]>0.0)) && (EarthShadowCheck(x_LOCAL)==false) ) { //calculate the radiation pressure force
        double rHeliocentric,vHeliocentric,radiationPressureAcceleration;

        //calcualte velocity of the particle in the "Frozen SO Frame";
        double v_LOCAL_SO_FROZEN[3];

        v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+v_LOCAL[0]+
            Exosphere::RotationVector_SO_FROZEN[1]*x_LOCAL[2]-Exosphere::RotationVector_SO_FROZEN[2]*x_LOCAL[1];

        v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+v_LOCAL[1]-
            Exosphere::RotationVector_SO_FROZEN[0]*x_LOCAL[2]+Exosphere::RotationVector_SO_FROZEN[2]*x_LOCAL[0];

        v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+v_LOCAL[2]+
            Exosphere::RotationVector_SO_FROZEN[0]*x_LOCAL[1]-Exosphere::RotationVector_SO_FROZEN[1]*x_LOCAL[0];


        rHeliocentric=sqrt(pow(x_LOCAL[0]-Exosphere::xObjectRadial,2)+(x_LOCAL[1]*x_LOCAL[1])+(x_LOCAL[2]*x_LOCAL[2]));
        vHeliocentric=(
            (v_LOCAL_SO_FROZEN[0]*(x_LOCAL[0]-Exosphere::xObjectRadial))+
            (v_LOCAL_SO_FROZEN[1]*x_LOCAL[1])+(v_LOCAL_SO_FROZEN[2]*x_LOCAL[2]))/rHeliocentric;

        radiationPressureAcceleration=SodiumRadiationPressureAcceleration__Combi_1997_icarus(vHeliocentric,rHeliocentric);

        accl_LOCAL[0]+=radiationPressureAcceleration*(x_LOCAL[0]-Exosphere::xObjectRadial)/rHeliocentric;
        accl_LOCAL[1]+=radiationPressureAcceleration*x_LOCAL[1]/rHeliocentric;
        accl_LOCAL[2]+=radiationPressureAcceleration*x_LOCAL[2]/rHeliocentric;
      }
#endif
    }

    if (PIC::MolecularData::ElectricChargeTable[spec]!=0.0) { //the Lorentz force
      double E[3],B[3];

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  #endif


      if (_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__OFF_) {
	 for (int idim=0;idim<3;idim++) B[idim]=Exosphere::swB_Typical[idim],E[idim]=Exosphere::swE_Typical[idim]; 
      }
      else {
        PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
        PIC::CPLR::GetBackgroundFieldsVector(E,B);
      }


      accl_LOCAL[0]+=PIC::MolecularData::ElectricChargeTable[spec]*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/PIC::MolecularData::MolMass[spec];
      accl_LOCAL[1]+=PIC::MolecularData::ElectricChargeTable[spec]*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/PIC::MolecularData::MolMass[spec];
      accl_LOCAL[2]+=PIC::MolecularData::ElectricChargeTable[spec]*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/PIC::MolecularData::MolMass[spec];

    }



    //the gravity force
    double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
    double r=sqrt(r2);
    int idim;

    for (idim=0;idim<DIM;idim++) {
      accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
    }


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    //correct the gravity acceleration: accout for solar gravity of the particle location
    //correction of the solar gravity

    double rSun2Moon,rSun2Particle;

    rSun2Moon=xObjectRadial;
    rSun2Particle=sqrt(pow(x_LOCAL[0]-xObjectRadial,2)+pow(x_LOCAL[1],2)+pow(x_LOCAL[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_SUN_)*((x_LOCAL[0]-xObjectRadial)/pow(rSun2Particle,3)+xObjectRadial/pow(rSun2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[1]/pow(rSun2Particle,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[2]/pow(rSun2Particle,3));

    //correction of the Earth gravity
    double rEarth2Moon,rEarth2Particle;

    rEarth2Moon=sqrt(pow(xEarth_SO[0],2)+pow(xEarth_SO[1],2)+pow(xEarth_SO[2],2));
    rEarth2Particle=sqrt(pow(x_LOCAL[0]-xEarth_SO[0],2)+pow(x_LOCAL[1]-xEarth_SO[1],2)+pow(x_LOCAL[2]-xEarth_SO[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[0]-xEarth_SO[0])/pow(rEarth2Particle,3)+xEarth_SO[0]/pow(rEarth2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[1]-xEarth_SO[1])/pow(rEarth2Particle,3)+xEarth_SO[1]/pow(rEarth2Moon,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[2]-xEarth_SO[2])/pow(rEarth2Particle,3)+xEarth_SO[2]/pow(rEarth2Moon,3));

    //account for the planetary rotation around the Sun
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
    accl_LOCAL[2]+=aCen[2]+aCorr[2];
#endif

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }



  inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    static const double LifeTime=3600.0*5.8/pow(0.4,2);


    //only sodium can be ionized
    if (spec!=_NA_SPEC_) {
      PhotolyticReactionAllowedFlag=false;
      return -1.0;
    }

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    double res,r2=x[1]*x[1]+x[2]*x[2];

    //check if the particle is outside of the Earth and lunar shadows
    if ( ((r2>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x[0]<0.0)) && (Moon::EarthShadowCheck(x)==false) ) {
      res=LifeTime,PhotolyticReactionAllowedFlag=true;
    }
    else {
      res=-1.0,PhotolyticReactionAllowedFlag=false;
    }

    //check if the particle intersect the surface of the Earth
    if (pow(x[0]-Moon::xEarth_SO[0],2)+pow(x[1]-Moon::xEarth_SO[1],2)+pow(x[2]-Moon::xEarth_SO[2],2)<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) {
      res=1.0E-10*LifeTime,PhotolyticReactionAllowedFlag=true;
    }
#else
    double res=LifeTime;
    PhotolyticReactionAllowedFlag=true;
#endif

    return res;
  }

  inline int ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {

    if ((spec==_NA_SPEC_)&&(_NA_PLUS_SPEC_>=0)) {
      PIC::ParticleBuffer::SetI(_NA_PLUS_SPEC_,ParticleData);
      return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;
    }

    return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
  }

}

#endif /* MERCURY_H_ */
