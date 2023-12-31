/*
 * pic_kameleon.cpp
 *
 *  Created on: Jan 1, 2015
 *      Author: vtenishe
 */
//$Id$
//AMPS' coupler functions that call CCMC's data reader Kameleon for interpolatino of the pre-calculated plasma model results stodered in the .cdf format

#include "pic.h"


//init the reader

void PIC::CPLR::DATAFILE::KAMELEON::Init() {
  PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::ElectricField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.allocate=true;
}

double PIC::CPLR::DATAFILE::KAMELEON::PlasmaSpeciesAtomicMass=1.0*_AMU_;
double PIC::CPLR::DATAFILE::KAMELEON::cdfDataFile2m_ConversionFactor=_RADIUS_(_EARTH_);


#ifdef _PIC_COMPILE__KAMELEON_ //complile the KAMELEON part only when _PIC_COMPILE__KAMELEON_ is defined; The macro is defined in the Makefile only when 'ifneq ($(KAMELEON),nokameleon)' 
#undef Pi //Kameleon uses its own cpnstant 'Pi'; The definition of the "global" 'Pi' is remove to eliminated the naming conflict with Kameleon library
#undef DIM

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <ccmc/LFM.h>
#include <ccmc/LFMInterpolator.h>
#include <ccmc/Interpolator.h>
#include <ccmc/Point3f.h>
#include <ccmc/Kameleon.h>
#include <ccmc/Tracer.h>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <vector>

void PIC::CPLR::DATAFILE::KAMELEON::GetDomainLimits(double *xmin,double *xmax,const char *fname) {
  ccmc::Kameleon kameleon;
  char DataFileFullName[_MAX_STRING_LENGTH_PIC_];

  sprintf(DataFileFullName,"%s/%s",PIC::CPLR::DATAFILE::path,fname);
  long status = kameleon.open(DataFileFullName);

  if(status == ccmc::FileReader::OK){
    //proceed and get domain limits
    xmin[0]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("x","actual_min").getAttributeFloat();
    xmin[1]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("y","actual_min").getAttributeFloat();
    xmin[2]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("z","actual_min").getAttributeFloat();
    
    xmax[0]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("x","actual_max").getAttributeFloat();
    xmax[1]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("y","actual_max").getAttributeFloat();
    xmax[2]=cdfDataFile2m_ConversionFactor*kameleon.getVariableAttribute("z","actual_max").getAttributeFloat();

    //check wether xmin and xmax are consistent
    for (int idim=0;idim<3;idim++) if (xmax[idim]-xmin[idim]<=0.0) exit(__LINE__,__FILE__,"Error: Particle tracking domain limites are incorrect (xmin exeeds xmax)"); 

    // close file
    kameleon.close();
  }
  else
    exit(__LINE__, __FILE__, "ERROR: couldn't open Kameleon file");
}

void PIC::CPLR::DATAFILE::KAMELEON::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  static ccmc::Kameleon kameleon;
  //  create interpolator
  static ccmc::Interpolator *interpolator;//=kameleon.createNewInterpolator();

  double ConversionVelocity = 1E3;
  double ConversionBField   = 1E-9;
  double ConversionEField   = -1E-6;
  double ConversionDensity  = 1E6;
  double ConversionPressure = 1E-9;

  if (startNode==PIC::Mesh::mesh->rootTree) {

    //open the file
    char DataFileFullName[_MAX_STRING_LENGTH_PIC_];

    sprintf(DataFileFullName,"%s/%s",PIC::CPLR::DATAFILE::path,fname);
    long status = kameleon.open(DataFileFullName);

    if (PIC::ThisThread==0) printf("$PREFIX: Loading file: %s\n",DataFileFullName);

    if (status != ccmc::FileReader::OK) exit(__LINE__, __FILE__, "ERROR: could'nt open Kameleon file");
    //allocate the interpolator and read the data file

    // load variables
    kameleon.loadVariable("p");
    kameleon.loadVariable("rho");
    kameleon.loadVariable("bx");
    kameleon.loadVariable("by");
    kameleon.loadVariable("bz");
    kameleon.loadVariable("ux");
    kameleon.loadVariable("uy");
    kameleon.loadVariable("uz");

    //  create interpolator
    interpolator=kameleon.createNewInterpolator();
  }

  //symbolic names of the interpolated variables
  static const std::string bVariables[] = {"bx", "by", "bz"};
  static const std::string eVariables[] = {"ex", "ey", "ez"};
  static const std::string vVariables[] = {"ux", "uy", "uz"};

  //perform the interpolation loop
  int i,j,k,nd,idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;



  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double x[3],T,n,p;

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
      x[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
      x[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);

      //transfor the coordinated into the LFM frame
      for (idim=0;idim<3;idim++) x[idim]/=_RADIUS_(_EARTH_);

      //locate the cell
      nd=_getCenterNodeLocalNumber(i,j,k);
      if (startNode->block==NULL) continue;
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
      offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;

      //save the interpolated values
      for (idim=0;idim<3;idim++) {
        *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset))=ConversionBField*interpolator->interpolate(bVariables[idim],x[0],x[1],x[2]);
        *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))=ConversionEField*interpolator->interpolate(eVariables[idim],x[0],x[1],x[2]);
        *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.RelativeOffset))=ConversionVelocity*interpolator->interpolate(vVariables[idim],x[0],x[1],x[2]);
      }

      p=ConversionPressure*interpolator->interpolate("p",x[0],x[1],x[2]);
      n=ConversionDensity *interpolator->interpolate("rho",x[0],x[1],x[2])/PlasmaSpeciesAtomicMass;
      T=(n>0.0) ? p/(n*Kbol) : 0.0;

      *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.RelativeOffset))=p;
      *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.RelativeOffset))=n;
      *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.RelativeOffset))=T;
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(fname,startNode->downNode[nDownNode]);
  }


  if (startNode==PIC::Mesh::mesh->rootTree) {
    //de-allocate the interpolation procedure
    //delete lfm;
    delete interpolator;

    interpolator=NULL;

    //close file
    kameleon.close();
    std::cout << "Kameleon file "<< fname << " was loaded successfully\n";

    //initialize derived data
    if (PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate==true) {
      if (_PIC_COUPLER__INTERPOLATION_MODE_==_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_) { 
        exit(__LINE__,__FILE__,"ERROR: magnetic field gradient can't be computed with 0th order interpolation method");
      } 

      for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        PIC::CPLR::DATAFILE::GenerateMagneticFieldGradient(node);
      }

      //Exchange derived data betwenn the boundary nodes
      PIC::Mesh::mesh->ParallelBlockDataExchange();
    }

  }
}

#endif

