//========================================================================
//$Id: Mars.h,v 1.11 2016/09/11 20:39:17 yunilee Exp $
//========================================================================
#ifndef _MARS_EXOSPEHRE_
#define _MARS_EXOSPEHRE_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>
#include <dirent.h>

#include "pic.h"
#include "Exosphere.h"

#include "specfunc.h"
#include "ifileopr.h"
#include "MTGCM.h"
#include "m-gitm.h"
#include "Mars.dfn"

#include "quadrature.h"
#include "SphericalVolumeMesh.h"

//Background atmosphere from J. Fox
#include "MarsBackgroundAtmosphereFox.h"

//the forward scattering cross section
#include "Kharchenko-2000-fig-3.h"
#include "NewOCO2-allE-total-smoothed.h"
#include "NewOCO2-table_E.h"

//List of types of O2+ DR coefficients
#define Hodges 0
#define Mehr 1
#define Peverall 2

#define DRO2p Mehr

//List of types of CO+ DR coefficients
#define Rosen 3
#define Cloutier 4
#define Mitchell 5

#define DRCOp Rosen

using namespace std;

const double massO=26.56E-27; //kg 
const double massC=19.93E-27; //kg
const double massO2=2*26.56E-27; //kg
const double massCO=massC+massO; //kg
const double k=1.3806503E-23; //m2 kg s-2 k-1
const double Semimajoraxis=1.523679; //AU, 2227,939,100km

//the offset in the cells' data to the local maximum of the cells' background densit
extern int maxLocalBackdroundDensityOffset;

//BEGIN: the global parameters of the forward scattering cross section
#include "DifferentialCrossSection.h"
#include "DifferentialCrossSection_O-CO2-Edependent.h"
extern cDiffCrossSection ForwardCollisionCrossSection;
extern cDiffCrossSection_OCO2 ForwardCollisionCrossSection_OCO2;
// extern cDiffCrossSection_CCO2 ForwardCollisionCrossSection_CCO2;

extern int maxLocalBackdroundDensityOffset;

extern double TotalIntegratedForwardScatteringCrossSection;
extern int CumulativeDistributionMaskList;
extern int *CumulativeDistributionMask;
extern int nForwardScatteringCrossSectionLines;

#ifndef _cForwardScatteringCrossSection_
#define _cForwardScatteringCrossSection_
struct cForwardScatteringCrossSection {
  double Angle,DifferentialCrossSection,CumulativeDistributionFunction,deltaCumulativeDistributionFunction;
};
#endif

extern cForwardScatteringCrossSection *ForwardScatteringCrossSectionData;
//END: the global parameters of the forward scattering cross section
 
namespace newMars {
inline std::string GetMGITMInputFilePath() {
    return std::string(PIC::UserModelInputDataPath)+"/MTGCM_equinox_SL/MGITM_PERMAX2.dat";
}

    extern cDataSetMGITM Te,Ti,Tn,CO2,O,N2,CO,O2p,E,Un,Vn,Wn;

//offsets of the model sampled data
extern int maxLocalCellOxigenProductionRateOffset,minLocalCellOxigenProductionRateOffset,minLocalBackdroundDensityOffset;
extern int sampledLocalInjectionRateOffset;
extern double *SampledEscapeRate;
extern double *IonizationRate;

//int maxLocalBackdroundDensityOffset;

//read the forward collision cross section data
//extern cDiffCrossSection ForwardCollisionCrossSection;
 
  //spherical mesh for sampling macroscopic parameters
    const int nSampledFieldVariablesSphericalVolumeMesh=2;
    const int LocalDensityOffsetSamplingSphericalVolumeMesh=0;
    const int LocalSourceRateOffsetSamplingSphericalVolumeMesh=1;
    extern cSphericalVolumeMesh SphericalSamplingMesh;
    int PrintVariableListSphericalVolumeMesh(FILE*);
    double GetSampledVariableSphericalSamplingMesh(int,double*,int);
    
    
    //sample particle data on the spherical mesh
    void inline SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec) {
        int el;
        double x[3];
        
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
        if (!isfinite(LocalParticleWeight)) {
            exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
        }
#endif
        
        
        PIC::ParticleBuffer::GetX(x,(PIC::ParticleBuffer::byte*)ParticleData);
        el=SphericalSamplingMesh.GetSamplingElementNumber(x);
        if (el!=-1) SphericalSamplingMesh.SamplingBuffer[LocalDensityOffsetSamplingSphericalVolumeMesh][el]+=LocalParticleWeight;
    }
    
    
//request sampling buffer
int RequestSamplingData(int);
int RequestStaticCellData(int);

inline void ReadMTGCM() {
    const std::string fname=GetMGITMInputFilePath();
    
		Te.PlanetRadius=_RADIUS_(_TARGET_);
		Te.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
		Te.ReadDataFile(_Te_MGITM_,fname.c_str());
		
		Tn.PlanetRadius=_RADIUS_(_TARGET_);
		Tn.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
		Tn.ReadDataFile(_Tn_MGITM_,fname.c_str());
				
		Ti.PlanetRadius=_RADIUS_(_TARGET_);
		Ti.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
		Ti.ReadDataFile(_Ti_MGITM_,fname.c_str());
		
		O2p.PlanetRadius=_RADIUS_(_TARGET_);
		O2p.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
		O2p.ReadDataFile(_nO2P_MGITM_,fname.c_str());
				
		E.PlanetRadius=_RADIUS_(_TARGET_);
		E.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
		E.ReadDataFile(_Ne_MGITM_,fname.c_str());
		
		O.PlanetRadius=_RADIUS_(_TARGET_);
		O.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
		O.ReadDataFile(_nO_MGITM_,fname.c_str());
			
		CO2.PlanetRadius=_RADIUS_(_TARGET_);
		CO2.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
		CO2.ReadDataFile(_nCO2_MGITM_,fname.c_str());

                CO.PlanetRadius=_RADIUS_(_TARGET_);
                CO.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
                CO.ReadDataFile(_nCO_MGITM_,fname.c_str());

                N2.PlanetRadius=_RADIUS_(_TARGET_);
                N2.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_;
                N2.ReadDataFile(_nN2_MGITM_,fname.c_str());

                Un.PlanetRadius=_RADIUS_(_TARGET_);
                Un.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
                Un.ReadDataFile(_UN_MGITM_,fname.c_str());

                Vn.PlanetRadius=_RADIUS_(_TARGET_);
                Vn.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
                Vn.ReadDataFile(_VN_MGITM_,fname.c_str());

                Wn.PlanetRadius=_RADIUS_(_TARGET_);
                Wn.OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_;
                Wn.ReadDataFile(_WN_MGITM_,fname.c_str());	
}
	
  //the total acceleration of the particles in the exosphere
  void TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);

	//print on the screen sampled total escape rate
	void SampleModelData();

  void OutputSampledModelData(int DataOutputFileNumber);

	//process particles that leaves the boundary of the computational domain: calcualte the escape rate
	int ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode); //(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);

	//init the model
  inline void Init_AfterParser() {
    int i,iAngle;
    double da;


    //read the forward collision cross section data
//    ForwardCollisionCrossSection.Set_ScatteringAngleMinimumIntegratedValue_Degrees(1.5);

//    ForwardCollisionCrossSection.AddProfile(0.03*eV2J,_KHANCHERKO_2000__0_03EV_DATA_LENGTH_,_KHANCHERKO_2000__0_03EV_DATA_);
//    ForwardCollisionCrossSection.AddProfile(0.3*eV2J,_KHANCHERKO_2000__0_3EV_DATA_LENGTH_,_KHANCHERKO_2000__0_3EV_DATA_);
//    ForwardCollisionCrossSection.AddProfile(3*eV2J,_KHANCHERKO_2000__3EV_DATA_LENGTH_,_KHANCHERKO_2000__3EV_DATA_);

    ForwardCollisionCrossSection.AddProfile(3*eV2J,_FOX_LENGTH_,_FOX_DATA_);
//    ForwardCollisionCrossSection.AddProfile(3.1*eV2J,_KHANCHERKO_LENGTH_,_KHANCHERKO_DATA_);
    ForwardCollisionCrossSection.InitInternalTable();
    ForwardCollisionCrossSection.PrintCrossSectionPlot();

//New O-CO2 cross section for all collision energy
    ForwardCollisionCrossSection_OCO2.AddProfile(_NEWOCO2_TABLE_E_LENGTH_,_NEWOCO2_TABLE_E_DATA_,_NEWOCO2_allE_LENGTH_,_NEWOCO2_allE_DATA_);
    ForwardCollisionCrossSection_OCO2.InitInternalTable();
    ForwardCollisionCrossSection_OCO2.PrintCrossSectionPlot();

//New C-CO2 cross section
    // ForwardCollisionCrossSection_CCO2.AddProfile(2.0*eV2J,_NEWCCO2_2EV_LENGTH_,_NEWCCO2_2EV_DATA_);
    // ForwardCollisionCrossSection_CCO2.AddProfile(_NEWCCO2_TABLE_E_LENGTH_,_NEWCCO2_TABLE_E_DATA_,_NEWCCO2_allE_LENGTH_,_NEWCCO2_allE_DATA_);
    // ForwardCollisionCrossSection_CCO2.InitInternalTable();
    // ForwardCollisionCrossSection_CCO2.PrintCrossSectionPlot();

    ReadMTGCM();
      //init the spherical sampling mesh
      SphericalSamplingMesh.init(100.0E3+_RADIUS_(_TARGET_),1000.0e3+_RADIUS_(_TARGET_),1000,50,25,nSampledFieldVariablesSphericalVolumeMesh,PrintVariableListSphericalVolumeMesh,GetSampledVariableSphericalSamplingMesh);
      SphericalSamplingMesh.plot("SamplingSphericalMesh.dat");


    //request sampling buffer for the Mars model
    PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);
    PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestStaticCellData);

    //init the buffer for sampling the escape rate
    SampledEscapeRate=new double[PIC::nTotalSpecies];
    for (i=0;i<PIC::nTotalSpecies;i++) SampledEscapeRate[i]=0.0;

    //init the buffer for sampling the ionization rate
    IonizationRate=new double[PIC::nTotalSpecies];
    for (i=0;i<PIC::nTotalSpecies;i++) IonizationRate[i]=0.0;
    
    //set up the functions that calculated and output the escape rate
    PIC::Mover::ProcessOutsideDomainParticles=ProcessOutsideDomainParticles;

    //set up the model sampling procedure
    PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SampleModelData,OutputSampledModelData);

    //init the procedure for calculation of the forward scatternig cross section
    /*TotalIntegratedForwardScatteringCrossSection=0.0;

    for (iAngle=0;iAngle<nForwardScatteringCrossSectionLines;iAngle++) {
      ForwardScatteringCrossSectionData[iAngle].Angle*=Pi/180.0;
      ForwardScatteringCrossSectionData[iAngle].DifferentialCrossSection*=1.0E-16*1.0E-4;
    }

    for (i=0;i<CumulativeDistributionMaskList;i++) CumulativeDistributionMask[i]=-1;

    for (iAngle=0;iAngle<nForwardScatteringCrossSectionLines-1;iAngle++) {
      da=ForwardScatteringCrossSectionData[iAngle+1].Angle-ForwardScatteringCrossSectionData[iAngle].Angle;

      ForwardScatteringCrossSectionData[iAngle].deltaCumulativeDistributionFunction=2.0*Pi* 0.5*da*(
          ForwardScatteringCrossSectionData[iAngle].DifferentialCrossSection*sin(ForwardScatteringCrossSectionData[iAngle].Angle)+
          ForwardScatteringCrossSectionData[iAngle+1].DifferentialCrossSection*sin(ForwardScatteringCrossSectionData[iAngle+1].Angle));

      ForwardScatteringCrossSectionData[iAngle].CumulativeDistributionFunction=TotalIntegratedForwardScatteringCrossSection;
      TotalIntegratedForwardScatteringCrossSection+=ForwardScatteringCrossSectionData[iAngle].deltaCumulativeDistributionFunction;
    }

    //normalize the cumulative distribution function; init the cumulative distribution mask
    for (iAngle=0;iAngle<nForwardScatteringCrossSectionLines;iAngle++) {
      ForwardScatteringCrossSectionData[iAngle].CumulativeDistributionFunction/=TotalIntegratedForwardScatteringCrossSection;
      ForwardScatteringCrossSectionData[iAngle].deltaCumulativeDistributionFunction/=TotalIntegratedForwardScatteringCrossSection;

      i=(int)(CumulativeDistributionMaskList*ForwardScatteringCrossSectionData[iAngle].CumulativeDistributionFunction);
      if (CumulativeDistributionMask[i]==-1) CumulativeDistributionMask[i]=iAngle;
    }

    //fill all possible hole in the cumulative distribution mask
    for (i=1;i<CumulativeDistributionMaskList;i++) if (CumulativeDistributionMask[i]==-1) {
      for (int j=i-1;j>=0;j--) if (CumulativeDistributionMask[j]!=-1) {
        CumulativeDistributionMask[i]=CumulativeDistributionMask[j];
        break;
      }
    }
*/
  }
	
inline double ProductionRateCaluclation_HotO(double *x) {
    double production_rate=0.0;
    if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<pow(O.PlanetRadius+O.minAltitude,2)) return 0.0;
	
    #if DRO2p == Hodges //Hodges [2000]
            production_rate=1.6E-7*pow(300/Te.Interpolate(x),0.55)*O2p.Interpolate(x)*E.Interpolate(x);
    #endif
    #if DRO2p == Mehr //Mehr and Biondi [1969]
            if (Te.Interpolate(x)<=1200.0) {production_rate=1.95E-7*pow(300/Te.Interpolate(x),0.7)*(O2p.Interpolate(x)/1.0e6)*(E.Interpolate(x)/1.0e6);}
            else if (Te.Interpolate(x)>1200.0) {production_rate=7.39E-8*pow(1200/Te.Interpolate(x),0.56)*(O2p.Interpolate(x)/1.0e6)*(E.Interpolate(x)/1.0e6);}
    #endif
    #if DRO2p == Peverall //Peverall [2001]
            production_rate=2.4E-7*pow(300/Te.Interpolate(x),0.7)*O2p.Interpolate(x)*E.Interpolate(x);
    #endif
        
    
    return production_rate*1.0e6;

}

//Below is OLD hot C function which can be deleted//
//Newer version for hot C model is available//
//
//inline double ProductionRateCaluclation_HotC(double *x) {
//    
//    #if _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Dissociative_Recombination_COp_
//    double production_rate=0.0;
//    double t=0.0;
//    //	bool validAlt = false;
//    double z0 = 200E3;
//    //	while (validAlt == false) {
//    if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<pow(O.PlanetRadius+O.minAltitude,2)) return 0.0;
//    
//    double Altitude2,rpeak;
//    double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
//    double r=sqrt(r2);
//    double AtmPressure,AtmPressurePeak,AzimuthalAng,PolarAng=0.0;
//    
//    double Altitude = r-_RADIUS_(_TARGET_);
//    
//    
//    
//    // Hot Carbon
//    double COpScaleHeight=53.8E3;//50.2E3; //Solar Min
//    //double COpScaleHeight=72.6E3;//102.1E3; //Solar Max
//    double peak=210.0E3;//solar min
//    //double peak=240.0E3;//solar max
//    
//    double vectop2[3]={0.0,0.0,0.0},vectop[3]={0.0,0.0,0.0},vectopP[3]={0.0,0.0,0.0};
//    //	double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
//    //	double r=sqrt(r2);
//    int idim;
//    
//    
//    
//    AtmPressure=Kbol*(O.Interpolate(x)+CO2.Interpolate(x))*Tn.Interpolate(x);
//    
//    for (idim=0;idim<DIM;idim++) {vectopP[idim]= (3388.25E3+210.0E3)*(x[idim]/r);}
//    AtmPressurePeak=Kbol*(OEQU.Interpolate(vectopP)+CO2EQU.Interpolate(vectopP))*TnEQU.Interpolate(vectopP);
//    // AtmPressurePeak=1.7E-13;
//    
//    double Pz=0.0;
//    Altitude2=0.0;
//    if (AtmPressure<=AtmPressurePeak+(0.1*AtmPressurePeak) && AtmPressure>=AtmPressurePeak-(0.1*AtmPressurePeak))
//    {Altitude2=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-_RADIUS_(_TARGET_);}
//    else if (Altitude <200.0E3)
//    {Altitude2=peak;}
//    else if (Altitude >= 200.0E3 && Altitude2==0.0) {
//        while (AtmPressure!=AtmPressurePeak)/*(AtmPressure>AtmPressurePeak+(0.1*AtmPressurePeak) || AtmPressure<AtmPressurePeak-(0.1*AtmPressurePeak))*/ {
//            
//            
//            if (AtmPressure<=AtmPressurePeak+(0.1*AtmPressurePeak) && AtmPressure>=AtmPressurePeak-(0.1*AtmPressurePeak))
//            {for (idim=0;idim<DIM;idim++) {vectopP[idim]=(3388.25E3+(r-(Pz+3388.25E3)))*(x[idim]/r);}
//                Altitude2=sqrt(vectopP[0]*vectopP[0]+vectopP[1]*vectopP[1]+vectopP[2]*vectopP[2])-_RADIUS_(_TARGET_);
//                break;}
//            
//            else if (((r-3388.25E3)-Pz)<200.0E3) {
//                Altitude2=peak;
//                break;}
//            
//            for (idim=0;idim<DIM;idim++) {vectopP[idim]=(3388.25E3+(r-(Pz+3388.25E3)))*(x[idim]/r);}
//            AtmPressure=Kbol*(O.Interpolate(vectopP)+CO2.Interpolate(vectopP))*Tn.Interpolate(vectopP);
//            
//            
//            Pz+=1E3;
//        }
//    }
//    
//    
//    
//    AtmPressure=Kbol*(O.Interpolate(x)+CO2.Interpolate(x))*Tn.Interpolate(x);
//    
//    bool validt=false;
//    while (validt==false) {
//        
//        for (idim=0;idim<DIM;idim++) {
//            vectop[idim] = (3388.25E3+z0)*(x[idim]/r);
//            vectop2[idim]= (3388.25E3+(z0-10E3))*(x[idim]/r);
//        }
//        
//        double localH=10E3/(log(COp.Interpolate(vectop)/COp.Interpolate(vectop2)));
//        bool isnan(localH);
//        if (localH<COpScaleHeight/3) {
//            localH=COpScaleHeight/2;
//        }
//        else if (isnan(localH)==true) {
//            localH=COpScaleHeight;
//        }
//        else if (localH>COpScaleHeight) {
//            localH=COpScaleHeight;
//        }
//        
//        if (Altitude>=135.0E3 && Altitude<z0){t=COp.Interpolate(x);}
//        else if (Altitude>=z0 && Altitude<Altitude2){t=COp.Interpolate(vectop)*exp((Altitude-z0)/localH);}
//        // else if (Altitude>=z0 && Altitude<peak){t=COp.Interpolate(vectop)*exp((Altitude-z0)/COpScaleHeight);}
//        else if (Altitude>=Altitude2){t=(COp.Interpolate(vectop)*exp((Altitude2-z0)/localH))*exp(-(Altitude-Altitude2)/COpScaleHeight);}
//        
//        //	if (nz>100E3){exit(__LINE__,__FILE__,"Error");}
//        bool isnan(t);
//        bool isinf(t);
//        
//        double nz=1E3,tt=0.0;
//        while (isnan(t)==true) {
//            if (nz>=Altitude) {t=1.0;}
//            for (idim=0;idim<DIM;idim++) {vectop2[idim]= (3388.25E3+(Altitude-nz))*(x[idim]/r);}
//            t=COp.Interpolate(vectop2);
//            nz+=1E3;
//        }
//        
//        if (Altitude>=1500E3) {
//            if (t>1E10) {z0-=1E3;}
//            else {validt=true;}
//        }
//        
//        else {
//            if (t>1E10 || t<1E-20) {z0-=1E3;}
//            else {validt=true;}
//        }
//        if (z0<135.0E3) {
//            t=1.0;
//            validt=true;}
//    }
//    
//    if (Altitude>200.0E3) {
//    #if DRCOp == Rosen //Rosen [1998]
//            production_rate=2.75E-7*pow(300/(4200-3750*exp((180-(Altitude/1000))/89.6)),0.55)*t*E.Interpolate(x);
//    #endif
//    #if DRCOp == Cloutier //Cloutier and Daniell [1979]
//            production_rate=6.47E-7*pow(300/(4200-3750*exp((180-(Altitude/1000))/89.6)),0.53)*t*E.Interpolate(x);
//    #endif
//    #if DRCOp == Mitchell //Mitchell and Hus [1985]
//            production_rate=2.0E-7*pow(300/(4200-3750*exp((180-(Altitude/1000))/89.6)),0.48)*t*E.Interpolate(x);
//    #endif
//        }
//        else {
//    #if DRCOp == Rosen //Rosen [1998]
//            production_rate=2.75E-7*pow(300/Te.Interpolate(x),0.55)*t*E.Interpolate(x);
//    #endif
//    #if DRCOp == Cloutier //Cloutier and Daniell [1979]
//            production_rate=6.47E-7*pow(300/Te.Interpolate(x),0.53)*t*E.Interpolate(x);
//    #endif
//    #if DRCOp == Mitchell //Mitchell and Hus [1985]
//            production_rate=2.0E-7*pow(300/Te.Interpolate(x),0.48)*t*E.Interpolate(x);
//    #endif
//    }
//    
//    return production_rate*1.0e6;
//    
//    #elif _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Photodissociation_CO_
//
//    double production_rate=0.0;
//    if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<pow(O.PlanetRadius+O.minAltitude,2)) return 0.0;
//    
//        
//    //Hot Carbon Photodissociation of CO
//            
//        //Photodissociation of CO
//        double frequency=(4.4E-7)/pow(Semimajoraxis,2); //s^-1, Solar low, Fox[89] //Huebner et al.[92]
//        //double frequency=(1.21E-6)/pow(Semimajoraxis,2); //s^-1, Solar High, Fox[89] // Huebner et al.[92]
//            
//        //determine the positon of subsolar point (tile angle = 25.19 deg)
//        double Lat,Lon,xsun[3],xunit[3],x_proj;
//        Lat=0.0*Pi/180; //Aphelion=25.19, Perihelion=-25.19, Equinox=0.0
//        Lon=180.0*Pi/180;
//    
//        //unit vector for anti-subsolar point
//        xunit[0]=-cos(Lat)*cos(Lon);
//        xunit[1]=-cos(Lat)*sin(Lon);
//        xunit[2]=-sin(Lat);
//    
//        //x_proj_to_anti-subsolar=unit_vec(unit_vec.x)
//        xsun[0]=xunit[0]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
//        xsun[1]=xunit[1]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
//        xsun[2]=xunit[2]*(xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2]);
//    
//        x_proj=sqrt(pow((x[0]-xsun[0]),2)+pow((x[1]-xsun[1]),2)+pow((x[2]-xsun[2]),2));
//    
//        if (((xunit[0]*x[0]+xunit[1]*x[1]+xunit[2]*x[2])>0.0)&&(x_proj<_RADIUS_(_TARGET_))) frequency=0.0; //No photodissociation on the nightside
//            
//        production_rate=CO.Interpolate(x)*frequency;
//            
//        
//    return production_rate;
//    
//#else
//    exit(__LINE__,__FILE__,"Error: HotC source is not defined");
//#endif
//
//}
    

	void ProductionRateCaluclation(bool *InjectionFlag,double *Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
	void SpeciesProductionRateCaluclation(int spec,bool &InjectionFlag,double &Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //print the model parameters
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

  //Ionization of hot corona
    double PhotoIonizationHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData);
    double ChargeExchangeHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData);
    double ElectronImpactHotOFrequency(PIC::ParticleBuffer::byte *modelParticleData);
    double IonizationHotO(PIC::ParticleBuffer::byte *modelParticleData);
    double TotalLifeTime(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    void PhotochemicalModelProcessor(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
    
    
    
	///////////////////////////////////

	inline void ProductionRateSliceMap() {
		
		
		double x[3],production=0.0,O2pdensity=0.0,Edensity=0.0,Etemp=0.0,Tneutral=0.0;
		
		const int nPoints=300;
		const double R=190E3+Te.PlanetRadius; //determine the altitude of interest
		const double dLon=2.0*Pi/(nPoints-1),dLat=Pi/(nPoints-1);
		
		FILE *fout=fopen("Production_rate.dat","w");
		fprintf(fout,"VARIABLES = \"Lon\", \"Lat\", \"Production Rate\",\"O2+ Density\",\"electron Density\",\"e temp\",\"neutral temp\"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",nPoints,nPoints);
		
		int i,j;
		double Lon,Lat;
		
		for (i=0;i<nPoints;i++) {
			Lat=-Pi/2.0+i*dLat;
			
			for (j=0;j<nPoints;j++) {
				Lon=j*dLon;
				
				x[0]=R*cos(Lat)*cos(Lon);
				x[1]=R*cos(Lat)*sin(Lon);
				x[2]=R*sin(Lat);
								
				production=ProductionRateCaluclation_HotO(x);
				O2pdensity=O2p.Interpolate(x);
				Tneutral=Tn.Interpolate(x);
				Edensity=E.Interpolate(x);
				Etemp=Te.Interpolate(x);
				
				if (production>1E10) {
					std::cout << __LINE__ << __FILE__ << std::endl;
				}
								
				fprintf(fout,"%e  %e  %e  %e  %e  %e  %e\n",Lon,Lat,production,O2pdensity,Edensity,Etemp,Tneutral);
								
			}
		}
		fclose(fout);
	}
	/*	
	double maxvalue(double array[]) {
		int length=array.length()l; //size of array
		int max=array[0]; //max=1st element
		
		for (int i=1;i<length;i++) {
			if (array[i]>max) max=array[i];}
		return max;
	}*/
	
	/*void CartesianToSpherical(double *Spherical, double *k) {
		double X=k[0],Y=k[1],Z=k[2];
		double R=sqrt(X*X+Y*Y+Z*Z);
		double LAT=asin(Z/R)/Pi*180.0;
		double LONG=0.0;
		if (Y>0) {
			if (X>0) {LONG=atan(Y/X)/Pi*180.0;}
			else if (X<0) {LONG=180.0+atan(Y/X)/Pi*180.0;}
			else if (X==0) {LONG=90.0;}
		}
		else if (Y<0) {
			if (X>0) {LONG=-atan(-Y/X)/Pi*180.0;}
			else if (X<0) {LONG=-(180.0+atan(-Y/X)/Pi*180.0);}
			else if (X==0) {LONG=-90.0;}
		}
		else if (Y==0) {
			if (X>0) {LONG=0.0;}
			if (X<0) {LONG=180.0;}
		}
		Spherical[0]=R;
		Spherical[1]=LONG;
		Spherical[2]=LAT;
	}*/
	
	inline void BGMeanFlowVelocity(double *velocity, double *x) {
		for (int idim=0;idim<3;idim++) velocity[idim]=0.0;
        //  #if WindsAndRotation==on
		double RotationSpeed=240.0;
        
        //	CartesianToSpherical(sph,x);
        //	double R=sph[0],LON=sph[1],LAT=sph[2];
		double R=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
		double LON=0.0, LAT=0.0;
		LAT=asin(x[2]/R)/Pi*180.0;
		if (x[1]>1E-10) {
			if (x[0]>1E-10) {LON=atan(x[1]/x[0])/Pi*180.0;}
			else if (x[0]<1E-10) {LON=180.0+atan(x[1]/x[0])/Pi*180.0;}
			else {LON=90.0;}
		}
		else if (x[1]<1E-10) {
			if (x[0]>1E-10) {LON=-atan(-x[1]/x[0])/Pi*180.0;}
			else if (x[0]<1E-10) {LON=-(180.0+atan(-x[1]/x[0])/Pi*180.0);}
			else {LON=-90.0;}
		}
		else {
			if (x[0]>1E-10) {LON=0.0;}
			if (x[0]<1E-10) {LON=180.0;}
		}
        
		double Z=Un.Interpolate(x)+RotationSpeed*cos(LAT/180.0*Pi);
		double M=Vn.Interpolate(x);
		double V=Wn.Interpolate(x);
		velocity[0]=V*cos(LAT/180.0*Pi)*cos(LON/180.0*Pi)-M*sin(LAT/180.0*Pi)*cos(LON/180.0*Pi)-Z*sin(LON/180.0*Pi);
		velocity[1]=V*cos(LAT/180.0*Pi)*sin(LON/180.0*Pi)-M*sin(LAT/180.0*Pi)*sin(LON/180.0*Pi)+Z*cos(LON/180.0*Pi);
		velocity[2]=V*sin(LAT/180.0*Pi)+M*cos(LAT/180.0*Pi);
        //	#endif

	}
	
namespace HotOxygen {

		//KE = (eV from each channel) * (eV to J) 
		
		extern float KineticEnergy1;
		extern float KineticEnergy2;
		extern float KineticEnergy3;
		extern float KineticEnergy4;
                extern float KineticEnergy;
  inline void init() {
       KineticEnergy1=0.83*1.60217653E-19;
       KineticEnergy2=3.05*1.60217653E-19;
       KineticEnergy3=5.02*1.60217653E-19;
       KineticEnergy4=6.98*1.60217653E-19;
       KineticEnergy=KineticEnergy4; 
   } 
        //const float KineticEnergy=KineticEnergy4;


		 long int HotOProduction(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);


}
namespace HotCarbon {
    
#if _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Dissociative_Recombination_COp_

    //Hot Carbon DR
    
    const float KineticEnergy1=2.90*1.60217653E-19;
    const float KineticEnergy2=1.64*1.60217653E-19;
    const float KineticEnergy3=0.94*1.60217653E-19;
#elif _HotCarbon_Source_Reaction_==_HotCarbon_Source_Reaction__Photodissociation_CO_

    //Hot Carbon Photodissociation
    //        const float KineticEnergy=2.56*1.60217652E-19; //solar low
    		const float KineticEnergy=2.58*1.60217652E-19; //solar high
    
#else
    exit(__LINE__,__FILE__,"Error: HotC source is not defined");
#endif
    
    long int HotCProduction(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

}
    
    double LocalTimeStep(int spec,bool& TimeStepLimitationImposed, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    //wrapper function
    long int HotAtomProduction_wrapper(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    
    
    
}



#endif




