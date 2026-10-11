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

  /**
   * Coordinate-frame and force conventions used by the lunar application.
   *
   * Keeping these names in one production namespace prevents a copied frame
   * literal from silently reintroducing the Mercury-era ``MSGR_HCI`` frame.
   * The LSO frame is defined by ``Kernels/OTHER/Moon.LSO.tf``: +X points from
   * the Moon to the Sun, -Y follows the Moon's heliocentric orbital velocity,
   * and +Z completes a right-handed triad.  ``MOON_ME_DE421`` is the mean-
   * Earth/polar-axis frame associated with the DE421 lunar ephemeris and is
   * also the coordinate convention declared by the LDEM_4 LOLA label.
   */
  namespace Frames {
    static const char Inertial[]="J2000";
    static const char SolarOrbital[]="LSO";
    static const char BodyFixed[]="MOON_ME_DE421";
    static const char ForceAberrationCorrection[]="NONE";
  }

  namespace OrbitalDynamics {
    /** Add lunar point-mass gravity, -GM r/|r|^3, in SI units. */
    inline void AddLunarPointMassAcceleration(
        double *accelerationMPerS2,const double *particlePositionM) {
      double radius2=0.0;
      for (int i=0;i<3;i++) {
        radius2+=particlePositionM[i]*particlePositionM[i];
      }
      const double radius3=radius2*sqrt(radius2);
      const double lunarGravitationalParameter=
          GravityConstant*_MASS_(_MOON_);
      for (int i=0;i<3;i++) {
        accelerationMPerS2[i]-=lunarGravitationalParameter*
            particlePositionM[i]/radius3;
      }
    }

    /**
     * Add the differential acceleration of a point-mass attractor.
     *
     * All vectors are expressed in the same Moon-centred frame and use metres;
     * ``attractorPositionM`` points from the Moon to the attracting body.
     * The returned contribution is
     *
     *   GM [ (R-r)/|R-r|^3 - R/|R|^3 ],
     *
     * which subtracts the acceleration of the lunar origin.  This is the
     * production kernel used for both the Sun and Earth, and is intentionally
     * exposed so U04 can compare it with an independent analytical evaluator.
     * The caller must not supply a particle at the attractor or an attractor at
     * the origin; those singular physical states retain the existing failure
     * behavior (non-finite floating-point output).
     */
    inline void AddDifferentialPointMassAcceleration(
        double *accelerationMPerS2,const double *particlePositionM,
        const double *attractorPositionM,double attractorMassKg) {
      double moonToAttractor2=0.0,particleToAttractor2=0.0;

      for (int i=0;i<3;i++) {
        moonToAttractor2+=attractorPositionM[i]*attractorPositionM[i];
        const double separation=attractorPositionM[i]-particlePositionM[i];
        particleToAttractor2+=separation*separation;
      }

      const double moonToAttractor3=
          moonToAttractor2*sqrt(moonToAttractor2);
      const double particleToAttractor3=
          particleToAttractor2*sqrt(particleToAttractor2);
      const double gravitationalParameter=GravityConstant*attractorMassKg;

      for (int i=0;i<3;i++) {
        accelerationMPerS2[i]+=gravitationalParameter*(
            (attractorPositionM[i]-particlePositionM[i])/
                particleToAttractor3-
            attractorPositionM[i]/moonToAttractor3);
      }
    }

    /**
     * Add centrifugal and Coriolis acceleration in the rotating LSO frame.
     *
     * Inputs use metres, metres per second, and radians per second, all with
     * LSO components.  Euler acceleration is absent because the legacy Moon
     * mover treats the angular velocity as frozen during one particle step;
     * I10 must quantify the resulting time-discretization error.
     */
    inline void AddRotatingFrameAcceleration(
        double *accelerationMPerS2,const double *particlePositionM,
        const double *particleVelocityMPerS,
        const double *angularVelocityRadPerS) {
      double omegaCrossPosition[3],omegaCrossOmegaCrossPosition[3];
      double omegaCrossVelocity[3];

      omegaCrossPosition[0]=angularVelocityRadPerS[1]*particlePositionM[2]-
          angularVelocityRadPerS[2]*particlePositionM[1];
      omegaCrossPosition[1]=angularVelocityRadPerS[2]*particlePositionM[0]-
          angularVelocityRadPerS[0]*particlePositionM[2];
      omegaCrossPosition[2]=angularVelocityRadPerS[0]*particlePositionM[1]-
          angularVelocityRadPerS[1]*particlePositionM[0];

      omegaCrossOmegaCrossPosition[0]=
          angularVelocityRadPerS[1]*omegaCrossPosition[2]-
          angularVelocityRadPerS[2]*omegaCrossPosition[1];
      omegaCrossOmegaCrossPosition[1]=
          angularVelocityRadPerS[2]*omegaCrossPosition[0]-
          angularVelocityRadPerS[0]*omegaCrossPosition[2];
      omegaCrossOmegaCrossPosition[2]=
          angularVelocityRadPerS[0]*omegaCrossPosition[1]-
          angularVelocityRadPerS[1]*omegaCrossPosition[0];

      omegaCrossVelocity[0]=angularVelocityRadPerS[1]*particleVelocityMPerS[2]-
          angularVelocityRadPerS[2]*particleVelocityMPerS[1];
      omegaCrossVelocity[1]=angularVelocityRadPerS[2]*particleVelocityMPerS[0]-
          angularVelocityRadPerS[0]*particleVelocityMPerS[2];
      omegaCrossVelocity[2]=angularVelocityRadPerS[0]*particleVelocityMPerS[1]-
          angularVelocityRadPerS[1]*particleVelocityMPerS[0];

      for (int i=0;i<3;i++) {
        accelerationMPerS2[i]-=omegaCrossOmegaCrossPosition[i]+2.0*
            omegaCrossVelocity[i];
      }
    }

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    /**
     * Return the instantaneous angular velocity of LSO relative to J2000.
     *
     * CSPICE ``xf2rav`` applied to the LSO-to-J2000 state transform returns
     * the angular velocity of J2000 relative to LSO, resolved in LSO.  The
     * required rotating-frame vector has the opposite sign, hence the explicit
     * negation below.  Using the state-transform derivative avoids the
     * ``acos``/``sin(angle)`` singularity of the former finite-step extraction.
     * Ephemeris time is TDB seconds past J2000; output is rad s^-1 in LSO.
     * CSPICE's configured error action controls failure for missing frames or
     * kernels, so an unavailable authoritative frame cannot silently fall back.
     */
    inline void GetSolarOrbitalAngularVelocityLSO(
        SpiceDouble ephemerisTime,double *angularVelocityRadPerS) {
      SpiceDouble stateTransform[6][6],rotation[3][3],inverseAngularVelocity[3];

      sxform_c(Frames::SolarOrbital,Frames::Inertial,ephemerisTime,
          stateTransform);
      xf2rav_c(stateTransform,rotation,inverseAngularVelocity);
      for (int i=0;i<3;i++) {
        angularVelocityRadPerS[i]=-inverseAngularVelocity[i];
      }
    }
#endif
  }

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



    // Lunar point gravity is evaluated by the same callable production kernel
    // used by U03; _TARGET_ is fixed to _MOON_ by the srcMoon configuration.
    OrbitalDynamics::AddLunarPointMassAcceleration(accl_LOCAL,x_LOCAL);


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    // LSO +X points from Moon to Sun, so the solar position is exactly the
    // positive radial distance.  Earth already has LSO components in metres.
    const double sunPositionLSOM[3]={xObjectRadial,0.0,0.0};
    OrbitalDynamics::AddDifferentialPointMassAcceleration(accl_LOCAL,x_LOCAL,
        sunPositionLSOM,_MASS_(_SUN_));
    OrbitalDynamics::AddDifferentialPointMassAcceleration(accl_LOCAL,x_LOCAL,
        xEarth_SO,_MASS_(_EARTH_));

    // Apply the two fictitious terms for a frame whose angular velocity is
    // frozen for this particle step.  The production time-step hook populates
    // the vector directly from the derivative of the SPICE state transform.
    OrbitalDynamics::AddRotatingFrameAcceleration(accl_LOCAL,x_LOCAL,v_LOCAL,
        RotationVector_SO_FROZEN);
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

    // LSO +X points from the Moon toward the Sun.  The cylindrical lunar
    // shadow therefore occupies x<0 with transverse radius below R_Moon;
    // points on the sunward side (x>0) or outside that cylinder are lit.  This
    // sign must match the radiation-pressure gate above.  The previous x<0
    // test incorrectly enabled photoionization in the near-lunar nightside.
    if ( ((r2>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x[0]>0.0)) && (Moon::EarthShadowCheck(x)==false) ) {
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
