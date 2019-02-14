//$Id$


#include "pic.h"
#include "constants.h"

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
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>

#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "meshAMRgeneric.h"

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"
#include "PeriodicBCTest.dfn"
#include "RandNum.h"
//for lapenta mover

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"

FILE * fMatVec;
FILE * fRhs;

double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {return 0.0;}
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {return 0.0;}
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EPHIOD";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EPHIOD";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Europa";
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;}

//const double DomainLength[3]={1.0E8,1.0E8,1.0E8},DomainCenterOffset[3]={0.0,0.0,0.0};

#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif

double rho0,rho1,ux0,uy0,uz0,ux1,uy1,uz1,p0,p1,ppar0,ppar1;
double uth0,uth1;
double bx0,by0,bz0,ex0,ey0,ez0;

int nCells[3] ={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
 
int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

void GetGlobalCellIndex(int * index ,double * x, double * dx, double * xmin){
  //global cell index starts from 1 for true cells
  index[0] = round((x[0]-xmin[0])/dx[0]+0.5);
  index[1] = round((x[1]-xmin[1])/dx[1]+0.5);
  index[2] = round((x[2]-xmin[2])/dx[2]+0.5);

}

void GetGlobalCornerIndex(int * index ,double * x, double * dx, double * xmin){
  
  index[0] = round((x[0]-xmin[0])/dx[0]);
  index[1] = round((x[1]-xmin[1])/dx[1]);
  index[2] = round((x[2]-xmin[2])/dx[2]);

}


bool isBoundaryCell(double *x, double *dx, double * xmin, double * xmax, int minIndex, int maxIndex){

  if ( maxIndex < minIndex)
    exit(__LINE__,__FILE__,"Error: minIndex is greater than maxIndex");
  
  int indexBoundary[3]; //index count from the boundary

  for (int idim=0; idim<3; idim++) {
    indexBoundary[idim]=
      (fabs(x[idim]-xmin[idim])<fabs(x[idim]-xmax[idim])?
       round((x[idim]-xmin[idim])/dx[idim]-0.5):round((xmax[idim]-x[idim])/dx[idim]-0.5));
    //minus value means outside the domain
    //positive value and zero mean inside
    for (int idx=minIndex;idx<=maxIndex; idx++){
      if (indexBoundary[idim]==idx) return true;
    }
  }

  return false;
}


// bool isBoundaryCorner(double *x, double *dx, double * xmin, double * xmax, int minIndex, int maxIndex){

//   if ( maxIndex < minIndex)
//     exit(__LINE__,__FILE__,"Error: minIndex is greater than maxIndex");
  
//   int indexBoundary[3]; //index count from the boundary

//   for (int idim=0; idim<3; idim++) {
//     indexBoundary[idim]=
//       (fabs(x[idim]-xmin[idim])<fabs(x[idim]-xmax[idim])?
//        round((x[idim]-xmin[idim])/dx[idim]):round((xmax[idim]-x[idim])/dx[idim]));
//     //minus value means outside the domain
//     //positive value inside
//     for (int idx=minIndex;idx<=maxIndex; idx++){
//       if (indexBoundary[idim]==idx) return true;
//     }
//   }

//   return false;
// }



double CellInterpolatedVar(std::string var,PIC::InterpolationRoutines::CellCentered::cStencil * centerStencilPtr){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  double value =0;
  if (var.substr(0, 2) == "Bx"){    

    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double * B_ptr = (double *)(centerStencilPtr->cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset);
      value+=centerStencilPtr->Weight[iStencil]*B_ptr[0];      
    }

  }else if (var.substr(0, 2) == "By"){

    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double * B_ptr = (double *)(centerStencilPtr->cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset);
      value+=centerStencilPtr->Weight[iStencil]*B_ptr[1];      
    } 

  }else if (var.substr(0, 2) == "Bz"){

    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double * B_ptr = (double *)(centerStencilPtr->cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset+CurrentBOffset);
      value+=centerStencilPtr->Weight[iStencil]*B_ptr[2];      
    } 

  }

  return value;

}

double CellInterpolatedVar(std::string var,PIC::InterpolationRoutines::CellCentered::cStencil * centerStencilPtr,int iSp){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  double value =0;
  if (var.substr(0, 3) == "Rho"){    
    //number density
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      value += centerStencilPtr->Weight[iStencil]*
        centerStencilPtr->cell[iStencil]->GetDatumAverage(PIC::Mesh::DatumNumberDensity, iSp);
      double xTemp[3];
      centerStencilPtr->cell[iStencil]->GetX(xTemp);
    }
  }else if (var.substr(0, 2) == "Ux"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3];
      centerStencilPtr->cell[iStencil]->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      value += centerStencilPtr->Weight[iStencil]*vTemp[0];
    }
  }else if (var.substr(0, 2) == "Uy"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3];
      centerStencilPtr->cell[iStencil]->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      value += centerStencilPtr->Weight[iStencil]*vTemp[1];
    }
  }else if (var.substr(0, 2) == "Uz"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3];
      centerStencilPtr->cell[iStencil]->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      value += centerStencilPtr->Weight[iStencil]*vTemp[2];
    }
  }else if (var.substr(0, 3) == "Pxx"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2Temp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2,v2Temp,iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2Temp[0]-vTemp[0]*vTemp[0]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }else if (var.substr(0, 3) == "Pyy"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2Temp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2,v2Temp,iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2Temp[1]-vTemp[1]*vTemp[1]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }else if (var.substr(0, 3) == "Pzz"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2Temp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2,v2Temp,iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2Temp[2]-vTemp[2]*vTemp[2]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }else if (var.substr(0, 3) == "Pxy"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2TensorTemp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2Tensor, v2TensorTemp, iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2TensorTemp[0]-vTemp[0]*vTemp[1]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }else if (var.substr(0, 3) == "Pyz"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2TensorTemp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2Tensor, v2TensorTemp, iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2TensorTemp[1]-vTemp[1]*vTemp[2]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }else if (var.substr(0, 3) == "Pxz"){
    for (int iStencil=0;iStencil<centerStencilPtr->Length;iStencil++) {
      double vTemp[3],v2TensorTemp[3];
      PIC::Mesh::cDataCenterNode *cell= centerStencilPtr->cell[iStencil]; 
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity,vTemp,iSp);
      cell->GetDatumAverage(PIC::Mesh::DatumParticleVelocity2Tensor, v2TensorTemp, iSp);
      value += centerStencilPtr->Weight[iStencil]*(v2TensorTemp[2]-vTemp[0]*vTemp[2]);
    }
    value *= CellInterpolatedVar("Rho",centerStencilPtr,iSp)
      /fabs(PIC::CPLR::FLUID::FluidInterface.get_qom(iSp));
  }

  return value;
}


double GetCornerVar(std::string var,char * DataPtr,int iSp){

  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  double value =0;
  if (var.substr(0, 3) == "Rho"){    
    value =((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+Rho_];
  }else if (var.substr(0, 2) == "Ux"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUx_];
    value /= ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+Rho_];

  }else if (var.substr(0, 2) == "Mx"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUx_];
  }else if (var.substr(0, 2) == "Uy"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUy_];
    value /= ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+Rho_];

  }else if (var.substr(0, 2) == "My"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUy_];
  }else if (var.substr(0, 2) == "Uz"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUz_];
    value /= ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+Rho_];
  }else if (var.substr(0, 2) == "Mz"){
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUz_];
  }else if (var.substr(0, 3) == "Pxx"){
    double uTemp = GetCornerVar("Ux",DataPtr,iSp);
    double rhoTemp = GetCornerVar("Rho",DataPtr,iSp);
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUxUx_] - rhoTemp*uTemp*uTemp;
  }else if (var.substr(0, 3) == "Pyy"){
    double uTemp = GetCornerVar("Uy",DataPtr,iSp);
    double rhoTemp = GetCornerVar("Rho",DataPtr,iSp);
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUyUy_] - rhoTemp*uTemp*uTemp;
  }else if (var.substr(0, 3) == "Pzz"){
    double uTemp = GetCornerVar("Uz",DataPtr,iSp);
    double rhoTemp = GetCornerVar("Rho",DataPtr,iSp);
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUzUz_] - rhoTemp*uTemp*uTemp;
  }else if (var.substr(0, 3) == "Pxy"){
    double uTemp = GetCornerVar("Uy",DataPtr,iSp);    
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUxUy_] - 
      ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUx_]*uTemp;


  }else if (var.substr(0, 3) == "Pyz"){
    double uTemp = GetCornerVar("Uz",DataPtr,iSp);    
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUyUz_] - 
      ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUy_]*uTemp;

  }else if (var.substr(0, 3) == "Pxz"){
    double uTemp = GetCornerVar("Uz",DataPtr,iSp);    
    value = ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUxUz_] - 
      ((double *)(DataPtr+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset))
      [SpeciesDataIndex[iSp]+RhoUx_]*uTemp;
  }

  return value;
}


void SendDataToFluid(char *NameVar, int *nVarIn, int *nDimIn, int *nPointIn, double *Xyz_I, double *data_I) {
    
  FluidPicInterface * col= &PIC::CPLR::FLUID::FluidInterface; 
    
  int nDim = *nDimIn;
  int nPoint = *nPointIn;
  int nVar = *nVarIn;

  const int ns = col->get_nS();

  double *dataPIC_I;
  // (rho + 3*Moment + 6*p)*ns + 3*E + 3*B;
  const int nVarSpecies = 10;
  int nVarPIC = ns * nVarSpecies + 6;
  const int iRho = 0, iMx = 1, iMy = 2, iMz = 3, iPxx = 4, iPyy = 5, iPzz = 6,
            iPxy = 7, iPxz = 8, iPyz = 9;
  int iBx = ns * nVarSpecies, iBy = iBx + 1, iBz = iBy + 1, iEx = iBz + 1,
      iEy = iEx + 1, iEz = iEy + 1;
  
  dataPIC_I = new double[nVarPIC];


  double xmin, ymin, zmin;
  xmin = col->getphyMin(0);
  ymin = col->getphyMin(1);
  zmin = col->getphyMin(2);

  double mhd_D[3], pic_D[3];
  double xp, yp, zp;


  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    mhd_D[0] = Xyz_I[iPoint * nDim] * col->getSi2NoL();
    mhd_D[1] = Xyz_I[iPoint * nDim + 1] * col->getSi2NoL();
    mhd_D[2] = Xyz_I[iPoint * nDim + 2] * col->getSi2NoL();
    col->mhd_to_Pic_Vec(mhd_D, pic_D);
    
    for (int ii=0;ii<nVarPIC;ii++) dataPIC_I[ii]=0;
    //  if (col->isThisRun(pic_D)) {
    if (true) {
      xp = pic_D[0];
      yp = (nDim > 1) ? pic_D[1] : 0.0;
      zp = (nDim > 2) ? pic_D[2] : 0.0;

      // double weight_I[8];
      int ix, iy, iz;
      // grid->getInterpolateWeight(xp, yp, zp, ix, iy, iz, weight_I);

      
      double xLoc[3]={xp,yp,zp};
      // double xLoc[3]={1,1.5,1};
      PIC::InterpolationRoutines::CornerBased::cStencil CornerStencil(false);
      CornerStencil=*(PIC::InterpolationRoutines::CornerBased::InitStencil(xLoc,NULL));
     
      char * DataPtr_I[8];
      for (int iSt=0; iSt<CornerStencil.Length;iSt++)
        DataPtr_I[iSt] = CornerStencil.cell[iSt]->GetAssociatedDataBufferPointer();


      for (int iSpecies = 0; iSpecies < ns; iSpecies++) {
        const int iStart = iSpecies * nVarSpecies;
        //   const double QoM = col->getQOM(iSpecies);
        // rho
        for (int iStencil=0;iStencil<CornerStencil.Length;iStencil++) {
          double tempWeight= CornerStencil.Weight[iStencil];
          dataPIC_I[iRho+iStart] +=
            tempWeight*GetCornerVar("Rho",DataPtr_I[iStencil],iSpecies);
      
          dataPIC_I[iMx+iStart]  += 
            tempWeight*GetCornerVar("Mx",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iMy+iStart]  += 
            tempWeight*GetCornerVar("My",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iMz+iStart]  +=
            tempWeight*GetCornerVar("Mz",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPxx+iStart]+=
             tempWeight*GetCornerVar("Pxx",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPyy+iStart]+=
             tempWeight*GetCornerVar("Pyy",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPzz+iStart]+=
             tempWeight*GetCornerVar("Pzz",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPxy+iStart]+=
             tempWeight*GetCornerVar("Pxy",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPxz+iStart]+=
             tempWeight*GetCornerVar("Pxz",DataPtr_I[iStencil],iSpecies);
          
          dataPIC_I[iPyz+iStart]+=
             tempWeight*GetCornerVar("Pyz",DataPtr_I[iStencil],iSpecies);
         
        }
      } // iSpecies


      for (int iStencil=0;iStencil<CornerStencil.Length;iStencil++) {
        double * tempB = 
          (double *)(DataPtr_I[iStencil]+
                     PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset 
                     +PIC::FieldSolver::Electromagnetic::ECSIM::OffsetB_corner
                     +PIC::FieldSolver::Electromagnetic::ECSIM::CurrentBOffset);

        double * tempE = 
          (double *)(DataPtr_I[iStencil]+
                     PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset 
                     +PIC::FieldSolver::Electromagnetic::ECSIM::CurrentEOffset);
        
        double tempWeight=CornerStencil.Weight[iStencil];
        // Bx
        dataPIC_I[iBx] += tempB[0]*tempWeight;
                         
        // By
        dataPIC_I[iBy] += tempB[1]*tempWeight;                         

        // Bz
        dataPIC_I[iBz] += tempB[2]*tempWeight;
                         
        // Ex
        dataPIC_I[iEx] += tempE[0]*tempWeight;
                         
        // Ey
        dataPIC_I[iEy] += tempE[1]*tempWeight;

        // Ez
        dataPIC_I[iEz] += tempE[2]*tempWeight;
      }
      // Combine PIC plasma data into MHD fluid data.
      col->CalcFluidState(dataPIC_I, &data_I[iPoint * nVar]);

    } // isThisRun

  } // iPoint

  delete[] dataPIC_I;
  

}

void SetParticleForCell(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,int iBlock,
                        int iCell,int jCell,int kCell,double * dx,double * xminBlock,
                        double * ParticleWeight,double CellVolume,long & nLocalInjectedParticles){

  double cellParNumPerSp = PIC::CPLR::FLUID::npcelx[0]*PIC::CPLR::FLUID::npcely[0]*PIC::CPLR::FLUID::npcelz[0];
  int nSubCells[3]={PIC::CPLR::FLUID::npcelx[0],PIC::CPLR::FLUID::npcely[0],PIC::CPLR::FLUID::npcelz[0]};
  int ind[3]={iCell,jCell,kCell};
  double xCenter[3],xCorner[3];
  for (int idim=0; idim<3; idim++) {
    xCenter[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
    xCorner[idim]=xCenter[idim]-0.5*dx[idim];
  }    
  // double BulkVel[PIC::nTotalSpecies][3];

  /*
  for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++){
    
    BulkVel[iSp][0]= PIC::CPLR::FLUID::FluidInterface.getPICUx(iBlock,
                                                              xCenter[0],xCenter[1],xCenter[2],iSp);
    
    BulkVel[iSp][1]= PIC::CPLR::FLUID::FluidInterface.getPICUy(iBlock,
                                                              xCenter[0],xCenter[1],xCenter[2],iSp);
    
    BulkVel[iSp][2]= PIC::CPLR::FLUID::FluidInterface.getPICUz(iBlock,
                                                              xCenter[0],xCenter[1],xCenter[2],iSp);
    
    */

  /*
    if (iSp==0 && PIC::CPLR::FLUID::iCycle==0){
      ux0=BulkVel[iSp][0];
      uy0=BulkVel[iSp][1];
      uz0=BulkVel[iSp][2];
    }
    
    if (iSp==1 && PIC::CPLR::FLUID::iCycle==0){
      ux1=BulkVel[iSp][0];
      uy1=BulkVel[iSp][1];
      uz1=BulkVel[iSp][2];
    }
    
    
  } 
  */
  nLocalInjectedParticles += cellParNumPerSp*PIC::nTotalSpecies;
  for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++){
    
    RandNum rndNum;
    int ig, jg, kg, nxcg, nycg, nzcg, npcel, nRandom = 7;
    int index_G[3]; 
    //assuming it is uniform
    nxcg=(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0])/dx[0];
    nycg=(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1])/dx[1]; 
    nzcg=(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2])/dx[2];
    
    
    GetGlobalCellIndex(index_G , xCenter, dx, PIC::Mesh::mesh.xGlobalMin);
    npcel = cellParNumPerSp;
    rndNum.set_seed((iSp + 3) * nRandom * npcel *
                            (nxcg * nycg * nzcg * PIC::CPLR::FLUID::iCycle + 
                             nycg * nzcg * index_G[0] + nzcg * index_G[1] + index_G[2]));
    
    //    double NumberDensity=fabs(PIC::CPLR::FLUID::FluidInterface.getPICRhoNum(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp));
    /*

    if (iSp==0 && iCell==1 && jCell==1 && kCell==1) {
      rho0=NumberDensity;
      p0 = PIC::CPLR::FLUID::FluidInterface.getPICP(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
      ppar0 =  PIC::CPLR::FLUID::FluidInterface.getPICPpar(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
      uth0 = PIC::CPLR::FLUID::FluidInterface.getPICUth(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
    } 
    
    if (iSp==1 && iCell==1 && jCell==1 && kCell==1) {
      
      rho1=NumberDensity;
      p1 = PIC::CPLR::FLUID::FluidInterface.getPICP(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
      uth1 = PIC::CPLR::FLUID::FluidInterface.getPICUth(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
      ppar1 =  PIC::CPLR::FLUID::FluidInterface.getPICPpar(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
      
    }
      
    */
    //    double weightCorrection=NumberDensity*CellVolume/ParticleWeight[iSp]/cellParNumPerSp;
    
    int npart = cellParNumPerSp;
    
    for (int ii = 0; ii < PIC::CPLR::FLUID::npcelx[0]; ii++){
      for (int jj = 0; jj < PIC::CPLR::FLUID::npcely[0]; jj++){
        for (int kk = 0; kk < PIC::CPLR::FLUID::npcelz[0]; kk++){
          
          double xPar[3];
          int index_sub[3]={ii,jj,kk};
          for (int idim=0; idim<3; idim++){
            
            xPar[idim] = (index_sub[idim] + rndNum()) * (dx[idim]/nSubCells[idim])
              + xCorner[idim];
          }
  
          double NumberDensity=fabs(PIC::CPLR::FLUID::FluidInterface.getPICRhoNum(iBlock,xPar[0],xPar[1],xPar[2],iSp));
  
          double weightCorrection=NumberDensity*CellVolume/ParticleWeight[iSp]/cellParNumPerSp;
    
          
          double ParVel_D[3];
          double uth[3];
          
          if ((iSp!=0 || PIC::CPLR::FLUID::FluidInterface.get_useElectronFluid()) && PIC::CPLR::FLUID::FluidInterface.getUseAnisoP()){
            
            double rand[4];
            for (int irand=0; irand<4; irand++) rand[irand]=rndNum();
            
            PIC::CPLR::FLUID::FluidInterface.setPICAnisoUth(
                                                         iBlock,xPar[0],xPar[1],xPar[2],
                                                         uth,uth+1,uth+2,
                                                         rand[0], rand[1], rand[2], rand[3],iSp);
            
          }else{
            double rand[4];
            for (int irand=0; irand<4; irand++) rand[irand]=rndNum();
            
            PIC::CPLR::FLUID::FluidInterface.setPICIsoUth(
                                                         iBlock,xPar[0],xPar[1],xPar[2],
                                                         uth,uth+1,uth+2,
                                                         rand[0], rand[1], rand[2], rand[3],iSp);        
          }

          double BulkVel[3];
            
            BulkVel[0]= PIC::CPLR::FLUID::FluidInterface.getPICUx(iBlock,
                                                                      xPar[0],xPar[1],xPar[2],iSp);
    
            BulkVel[1]= PIC::CPLR::FLUID::FluidInterface.getPICUy(iBlock,
                                                                      xPar[0],xPar[1],xPar[2],iSp);
            
            BulkVel[2]= PIC::CPLR::FLUID::FluidInterface.getPICUz(iBlock,
                                                                      xPar[0],xPar[1],xPar[2],iSp);
    



          for (int idim=0;idim<3;idim++)
            ParVel_D[idim] = BulkVel[idim]+uth[idim];
          /*
          static int icnt=0;
          if (icnt==0 && iSp==0) {
            printf("test populate xPar:%e,%e,%e, BulkVel:%e,%e,%e, uth:%e,%e,%e \n",
                   xPar[0],xPar[1],xPar[2],BulkVel[0],BulkVel[1],BulkVel[2],
                   uth[0], uth[1], uth[2]
                   );
            icnt++;
          }
          */
          //electronVelocity[idim]=uth_e* sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd())+electronBulkVelocity[idim];
          //   ionVelocity[idim]=uth_i*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd())+ionBulkVelocity[idim];   
           
          /*
          if (iSp==0){
          for (int idim=0;idim<3;idim++){
            //ParVel_D[idim] = BulkVel[iSp][idim]+uth[idim];
            ParVel_D[0] = ux0+ 
              uth0*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
            ParVel_D[1] = uy0+ 
              uth0*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
            ParVel_D[2] = uz0+
              uth0*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
          }
          }
          
          if (iSp==1){
          for (int idim=0;idim<3;idim++){
            //ParVel_D[idim] = BulkVel[iSp][idim]+uth[idim];
            ParVel_D[0] = ux1+
              uth1*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
            ParVel_D[1] = uy1+
              uth1*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
            ParVel_D[2] = uz1+
              uth1*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd());
         
          }
          }
          */
          /*
          if (PIC::CPLR::FLUID::iCycle==1)
          fprintf(fPartData,"%d %d %d %d %e %e %e %e %e %e %e %e\n",
                  index_G[0],index_G[1],index_G[2],iSp,
                  xPar[0],xPar[1],xPar[2],
                  ParVel_D[0],ParVel_D[1],ParVel_D[2], BulkVel[0], uth[1]);
          */
          PIC::ParticleBuffer::InitiateParticle(xPar, ParVel_D,&weightCorrection,&iSp,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
          
        }
      }
    }
  }
  
  //  fclose(fPartData);

}


long int setFixedParticle_BC(){
  

  printf("setFixedBC, iCycle:%ld\n",PIC::CPLR::FLUID::iCycle);

  // if (PIC::CPLR::FLUID::iCycle==0) return 0;

  static int cnt=0;
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;  
  //set field for one layer ghost corners 

  char fullname[STRING_LENGTH];
  sprintf(fullname,"test_fixedbc_before_iter=%d.dat",cnt);
  //PIC::Mesh::mesh.outputMeshDataTECPLOT(fullname,0);                                                                         
  

  const int nI =_BLOCK_CELLS_X_, nJ = _BLOCK_CELLS_Y_, nK=  _BLOCK_CELLS_Z_; 
  
  long nParticleDeleted=0, nParticleCreated=0;
    
  int ** countedFaceBC;
  countedFaceBC = new int * [PIC::DomainBlockDecomposition::nLocalBlocks];
  countedFaceBC[0] = new int [PIC::DomainBlockDecomposition::nLocalBlocks*6]; 

  for (int iBlk=0; iBlk<PIC::DomainBlockDecomposition::nLocalBlocks;iBlk++){
    countedFaceBC[iBlk] =  countedFaceBC[0]+6*iBlk;
    for (int iFace=0; iFace<6; iFace++) countedFaceBC[iBlk][iFace]=0;
  }  
  //printf("fixed Bc is called!\n");
  printf("PIC::DomainBlockDecomposition::nLocalBlocks:%d\n",PIC::DomainBlockDecomposition::nLocalBlocks);

          
  for (int iFace=0;iFace<6;iFace++){
    
    //int iBlock=0;
    //    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*   node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    for (int iBlk=0; iBlk<PIC::DomainBlockDecomposition::nLocalBlocks;iBlk++){

      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = PIC::DomainBlockDecomposition::BlockTable[iBlk];


    if (!node->block) continue;
    if (node->GetNeibFace(iFace,0,0)!=NULL) continue;
   
    //if (node->GetNeibFace(iFace,0,0)->Thread!=-1) continue;

    int ind_beg[3]={0,0,0}, ind_end[3] = {nCells[0]-1, nCells[1]-1,nCells[2]-1};

    
    int * iFaceMin[6]={&ind_beg[0],&ind_end[0],&ind_beg[0],&ind_beg[0],&ind_beg[0],&ind_beg[0]};
    int * iFaceMax[6]={&ind_beg[0],&ind_end[0],&ind_end[0],&ind_end[0],&ind_end[0],&ind_end[0]};
    int * jFaceMin[6]={&ind_beg[1],&ind_beg[1],&ind_beg[1],&ind_end[1],&ind_beg[1],&ind_beg[1]};
    int * jFaceMax[6]={&ind_end[1],&ind_end[1],&ind_beg[1],&ind_end[1],&ind_end[1],&ind_end[1]};
    int * kFaceMin[6]={&ind_beg[2],&ind_beg[2],&ind_beg[2],&ind_beg[2],&ind_beg[2],&ind_end[2]};
    int * kFaceMax[6]={&ind_end[2],&ind_end[2],&ind_end[2],&ind_end[2],&ind_beg[2],&ind_end[2]};
  
    
  
    double dx[3];
    double *xminBlock= node->xmin, *xmaxBlock= node->xmax;
          
    for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
    PIC::Mesh::cDataBlockAMR *block = node->block;
 
    long int *  FirstCellParticleTable=node->block->FirstCellParticleTable;
    double CellVolume=1;
    for (int idim=0;idim<3;idim++) CellVolume *= dx[idim];
    double ParticleWeight[PIC::nTotalSpecies];
    for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++)
      ParticleWeight[iSp]=node->block->GetLocalParticleWeight(iSp);
        
    //    for (int iFace=0;iFace<6;iFace++){
    //  if (node->neibNodeFace[iFace]!=NULL && node->neibNodeFace[iFace]->Thread!=-1) continue;

    for (int jface=0; jface<iFace; jface+=2) ind_beg[jface/2]+=countedFaceBC[iBlk][jface]; 
    for (int jface=1; jface<iFace; jface+=2) ind_end[(jface-1)/2]-=countedFaceBC[iBlk][jface]; 
    
    for (int i=*iFaceMin[iFace];i<=*iFaceMax[iFace];i++)
      for (int j=*jFaceMin[iFace];j<=*jFaceMax[iFace];j++)
        for (int k=*kFaceMin[iFace];k<=*kFaceMax[iFace];k++){
          
          
            //for (int i=0;i<nCells[0];i++) for (int j=0;j<nCells[1];j++) for (int k=0;k<nCells[2];k++) {
      
            double x[3];
            int ind[3]={i,j,k};
            
            for (int idim=0; idim<3; idim++) {
              x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
            }
 
            
            if (!isBoundaryCell(x, dx, PIC::Mesh::mesh.xGlobalMin, PIC::Mesh::mesh.xGlobalMax, 0, 0)) continue;
           
            //printf("fixRho at x:%e,%e,%e\n", x[0],x[1],x[2]);
              
            if (FirstCellParticleTable!=NULL){
              long int * ptr=FirstCellParticleTable+(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));
              while ((*ptr)!=-1) {
                PIC::ParticleBuffer::DeleteParticle(*ptr,*ptr);
                nParticleDeleted++;
              }
            }
              
            SetParticleForCell(node,iBlk,i,j,k,dx,xminBlock,ParticleWeight,CellVolume, nParticleCreated);
          }//for (int i=0;i<nCells[0];i++) for (int j=0;j<nCells[1];j++) for (int k=0;k<nCells[2];k++)
          
    countedFaceBC[iBlk][iFace]++;
    //iBlock++;
    }

  }


  sprintf(fullname,"test_fixedbc_after_iter=%d.dat",cnt);
  //PIC::Mesh::mesh.outputMeshDataTECPLOT(fullname,0);                                                                cnt++;         
  printf("setfixed bc done\n");

  //fclose(fPartData);

  delete [] countedFaceBC[0];
  delete [] countedFaceBC;
  return  nParticleCreated-nParticleDeleted;
}



void setFixedE_BC_half(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  int iBlock=0;

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {

    if (!node->block) continue;
        
    double dx[3];
    double *xminBlock= node->xmin, *xmaxBlock= node->xmax;
    // for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
       
    for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
    PIC::Mesh::cDataBlockAMR *block = node->block;
    
    for (int i=-1;i<=nCells[0];i++) for (int j=-1;j<=nCells[1];j++) for (int k=-1;k<=nCells[2];k++) {
              
          //if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
          
          //}
          double x[3];
          int ind[3]={i,j,k};
             
          for (int idim=0; idim<3; idim++) {
            x[idim]=xminBlock[idim]+ind[idim]*dx[idim];
          }


          if (!PIC::CPLR::FLUID::isBoundaryCorner(x, dx, PIC::Mesh::mesh.xGlobalMin, PIC::Mesh::mesh.xGlobalMax, 0, 0)) continue;
              
          PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
          if (CornerNode!=NULL){
            char *  offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                
            double Ex,Ey,Ez;
                 
            Ex = PIC::CPLR::FLUID::FluidInterface.getEx(iBlock,x[0],x[1],x[2]);
            Ey = PIC::CPLR::FLUID::FluidInterface.getEy(iBlock,x[0],x[1],x[2]);
            Ez = PIC::CPLR::FLUID::FluidInterface.getEz(iBlock,x[0],x[1],x[2]);
            
            ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=Ex;
            ((double*)(offset+OffsetE_HalfTimeStep))[EyOffsetIndex]=Ey;
            ((double*)(offset+OffsetE_HalfTimeStep))[EzOffsetIndex]=Ez;
           
          }
        }// for (int i=iFaceMin_n[iface];i<=iFaceMax_n[iface];i++)...

    iBlock++;
  }

}



void setFixedE_BC_curr(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  int iBlock=0;

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {

    if (!node->block) continue;
        
    double dx[3];
    double *xminBlock= node->xmin, *xmaxBlock= node->xmax;
    // for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
       
    for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
    PIC::Mesh::cDataBlockAMR *block = node->block;
    
    for (int i=-1;i<=nCells[0];i++) for (int j=-1;j<=nCells[1];j++) for (int k=-1;k<=nCells[2];k++) {
              
          //if (CornerNode==block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k))) {
          
          //}
          double x[3];
          int ind[3]={i,j,k};
             
          for (int idim=0; idim<3; idim++) {
            x[idim]=xminBlock[idim]+ind[idim]*dx[idim];
          }


          if (!PIC::CPLR::FLUID::isBoundaryCorner(x, dx, PIC::Mesh::mesh.xGlobalMin, PIC::Mesh::mesh.xGlobalMax, 0, 0)) continue;
              
          PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
          if (CornerNode!=NULL){
            char *  offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
                
            double Ex,Ey,Ez;
               
            Ex = PIC::CPLR::FLUID::FluidInterface.getEx(iBlock,x[0],x[1],x[2]);
            Ey = PIC::CPLR::FLUID::FluidInterface.getEy(iBlock,x[0],x[1],x[2]);
            Ez = PIC::CPLR::FLUID::FluidInterface.getEz(iBlock,x[0],x[1],x[2]);
       
            if (fabs(x[0]-32)<1e-3 && fabs(x[1]-5)<1e-3 && fabs(x[2]-3)<1e-3)
              printf("fixedBC E i,j,k:%d,%d,%d; x:%e,%e,%e, E:%e,%e,%e\n",i,j,k,x[0],x[1],x[2],Ex,Ey,Ez);

       
            ((double*)(offset+CurrentEOffset))[ExOffsetIndex]=Ex;
            ((double*)(offset+CurrentEOffset))[EyOffsetIndex]=Ey;
            ((double*)(offset+CurrentEOffset))[EzOffsetIndex]=Ez;
	             
          }
        }// for (int i=iFaceMin_n[iface];i<=iFaceMax_n[iface];i++)...

    iBlock++;
  }

}


void setFixedB_center_BC(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  int iBlock=0;

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {

    if (!node->block) continue;
        
    double dx[3];
    double *xminBlock= node->xmin, *xmaxBlock= node->xmax;
    // for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
       
    for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
    PIC::Mesh::cDataBlockAMR *block = node->block;
  
    for (int i=-1;i<=nCells[0];i++) for (int j=-1;j<=nCells[1];j++) for (int k=-1;k<=nCells[2];k++) {
      

          double x[3];
          int ind[3]={i,j,k};
            
          for (int idim=0; idim<3; idim++) {
            x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
          }

          if (!isBoundaryCell(x, dx, PIC::Mesh::mesh.xGlobalMin, PIC::Mesh::mesh.xGlobalMax, 0, 0)) continue;

          PIC::Mesh::cDataCenterNode *CenterNode= node->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
          if (CenterNode!=NULL){
            char *  offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
                
            double Bx,By,Bz;
              
            Bx = PIC::CPLR::FLUID::FluidInterface.getBx(iBlock,x[0],x[1],x[2]);
            By = PIC::CPLR::FLUID::FluidInterface.getBy(iBlock,x[0],x[1],x[2]);
            Bz = PIC::CPLR::FLUID::FluidInterface.getBz(iBlock,x[0],x[1],x[2]);
              
            ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=Bx;
            ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=By;
            ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=Bz;
		
            ((double*)(offset+PrevBOffset))[BxOffsetIndex]=Bx;
            ((double*)(offset+PrevBOffset))[ByOffsetIndex]=By;
            ((double*)(offset+PrevBOffset))[BzOffsetIndex]=Bz;
            
          }
        }
  
    iBlock++;
  }

}

void setFixedB_corner_BC(){
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  int iBlock=0;

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {

    if (!node->block) continue;
        
    double dx[3];
    double *xminBlock= node->xmin, *xmaxBlock= node->xmax;
    // for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
       
    for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
    PIC::Mesh::cDataBlockAMR *block = node->block;
  

    for (int i=-1;i<=nCells[0];i++) for (int j=-1;j<=nCells[1];j++) for (int k=-1;k<=nCells[2];k++) {
          /*   
               for (int i=iFaceMin[iFace];i<=iFaceMax[iFace];i++)
               for (int j=jFaceMin[iFace];j<=jFaceMax[iFace];j++)
               for (int k=kFaceMin[iFace];k<=kFaceMax[iFace];k++){
          */
          double x[3];
          int ind[3]={i,j,k};
            
          for (int idim=0; idim<3; idim++) {
            x[idim]=xminBlock[idim]+(ind[idim])*dx[idim];
          }
              
          if (!PIC::CPLR::FLUID::isBoundaryCorner(x, dx, PIC::Mesh::mesh.xGlobalMin, PIC::Mesh::mesh.xGlobalMax, 0, 0)) continue;

          PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
          if (CornerNode!=NULL){
            char *  offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset+OffsetB_corner;
                
            double Bx,By,Bz;
              
            Bx = PIC::CPLR::FLUID::FluidInterface.getBx(iBlock,x[0],x[1],x[2]);
            By = PIC::CPLR::FLUID::FluidInterface.getBy(iBlock,x[0],x[1],x[2]);
            Bz = PIC::CPLR::FLUID::FluidInterface.getBz(iBlock,x[0],x[1],x[2]);
              
            
            ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=Bx;
            ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=By;
            ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=Bz;
            
            ((double*)(offset+PrevBOffset))[BxOffsetIndex]=Bx;
            ((double*)(offset+PrevBOffset))[ByOffsetIndex]=By;
            ((double*)(offset+PrevBOffset))[BzOffsetIndex]=Bz;
            
          }
        }
    iBlock++;
  }

}


void CleanParticles(){
  
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) if (node->block!=NULL) {
   
     long int *  FirstCellParticleTable=node->block->FirstCellParticleTable;
     if (FirstCellParticleTable==NULL) continue;
     for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
	 for (int i=0;i<_BLOCK_CELLS_X_;i++) {
	   long int * ptr=FirstCellParticleTable+(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));
	   while ((*ptr)!=-1) PIC::ParticleBuffer::DeleteParticle(*ptr,*ptr);

//////
/*
long int next,ptr=FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]; 

while (ptr!=-1) {

  next=PIC::ParticleBuffer::GetNext(ptr);
  PIC::ParticleBuffer::DeleteParticle(ptr);
  ptr=next; 
}

FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]=-1;
*/

/////	   
	 }   
       }
     }
     
  }

}


long int PrepopulateDomain() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
  
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;
 
  int nBlock[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};


  int iBlock=0;
  printf("prepopulation icycle:%ld\n", PIC::CPLR::FLUID::iCycle);

  printf("UseAniso:%s",PIC::CPLR::FLUID::FluidInterface.getUseAnisoP()?"T\n":"F\n" );
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (!node->block) continue;
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
	  //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	  BoundaryBlock=true;
	  break;
	}
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;
 
    double * xminBlock=node->xmin,* xmaxBlock=node->xmax;
    double dx[3];
    double CellVolume=1;
    for (int idim=0;idim<3;idim++) {
      dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nBlock[idim];
      CellVolume *= dx[idim];
    }

    double ParticleWeight[PIC::nTotalSpecies];
    for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++)
      ParticleWeight[iSp]=node->block->GetLocalParticleWeight(iSp);
    PIC::Mesh::cDataBlockAMR *block = node->block;

    for (int iCell=0;iCell<nBlock[0];iCell++) 
      for (int jCell=0;jCell<nBlock[1];jCell++) 
        for (int kCell=0;kCell<nBlock[2];kCell++){
          
          SetParticleForCell(node,iBlock,iCell,jCell,kCell,dx,xminBlock,ParticleWeight,CellVolume,nLocalInjectedParticles);
          /*
          double cellParNumPerSp = npcelx[0]*npcely[0]*npcelz[0];
          int nSubCells[3]={npcelx[0],npcely[0],npcelz[0]};
          int ind[3]={iCell,jCell,kCell};
          double xCenter[3],xCorner[3];
          for (int idim=0; idim<3; idim++) {
            xCenter[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
            xCorner[idim]=xCenter[idim]-0.5*dx[idim];
          }    
          double BulkVel[PIC::nTotalSpecies][3];
          for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++){
         
            BulkVel[iSp][0]= PIC::CPLR::FLUID::FluidInterface.getPICUx(iBlock,
                                                                      xCenter[0],xCenter[1],xCenter[2],iSp);
          
            BulkVel[iSp][1]= PIC::CPLR::FLUID::FluidInterface.getPICUy(iBlock,
                                                                      xCenter[0],xCenter[1],xCenter[2],iSp);
            
            BulkVel[iSp][2]= PIC::CPLR::FLUID::FluidInterface.getPICUz(iBlock,
                                                                      xCenter[0],xCenter[1],xCenter[2],iSp);
            
            if (iSp==0){
              ux0=BulkVel[iSp][0];
              uy0=BulkVel[iSp][1];
              uz0=BulkVel[iSp][2];
            }

            if (iSp==1){
              ux1=BulkVel[iSp][0];
              uy1=BulkVel[iSp][1];
              uz1=BulkVel[iSp][2];
            }


          } 
          
          nLocalInjectedParticles += cellParNumPerSp*PIC::nTotalSpecies;
          for (int iSp=0;iSp<PIC::nTotalSpecies;iSp++){
            
            RandNum rndNum;
            int ig, jg, kg, nxcg, nycg, nzcg, npcel, nRandom = 7;
            int index_G[3]; 
            //assuming it is uniform
            nxcg=(PICxmax[0]-PICxmin[0])/dx[0];
            nycg=(PICxmax[1]-PICxmin[1])/dx[1]; 
            nzcg=(PICxmax[2]-PICxmin[2])/dx[2];

           
            GetGlobalCellIndex(index_G , xCenter, dx, PICxmin);
            npcel = cellParNumPerSp;
            rndNum.set_seed((iSp + 3) * nRandom * npcel *
                            (nxcg * nycg * nzcg * PIC::CPLR::FLUID::iCycle + 
                                nycg * nzcg * index_G[0] + nzcg * index_G[1] + index_G[2]));
            
            double NumberDensity=fabs(PIC::CPLR::FLUID::FluidInterface.getPICRhoNum(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp));
                  
            if (iSp==0 && iCell==1 && jCell==1 && kCell==1) {
              rho0=NumberDensity;
              p0 = PIC::CPLR::FLUID::FluidInterface.getPICP(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
              ppar0 =  PIC::CPLR::FLUID::FluidInterface.getPICPpar(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
              uth0 = PIC::CPLR::FLUID::FluidInterface.getPICUth(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
            } 

            if (iSp==1 && iCell==1 && jCell==1 && kCell==1) {

              rho1=NumberDensity;
              p1 = PIC::CPLR::FLUID::FluidInterface.getPICP(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
              uth1 = PIC::CPLR::FLUID::FluidInterface.getPICUth(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);
              ppar1 =  PIC::CPLR::FLUID::FluidInterface.getPICPpar(iBlock,xCenter[0],xCenter[1],xCenter[2],iSp);

            }
            
            double weightCorrection=NumberDensity*CellVolume/ParticleWeight[iSp]/cellParNumPerSp;

            int npart = cellParNumPerSp;

            for (int ii = 0; ii < npcelx[0]; ii++){
              for (int jj = 0; jj < npcely[0]; jj++){
                for (int kk = 0; kk < npcelz[0]; kk++){
                  
                  double xPar[3];
                  int index_sub[3]={ii,jj,kk};
                  for (int idim=0; idim<3; idim++){
                      
                      xPar[idim] = (index_sub[idim] + rndNum()) * (dx[idim]/nSubCells[idim])
                        + xCorner[idim];
                  }


              double ParVel_D[3];
              double uth[3];
              
              PIC::CPLR::FLUID::FluidInterface.setPICAnisoUth(
                                                             iBlock,xCenter[0],xCenter[1],xCenter[2],
                                                             uth,uth+1,uth+2,
                                                             rndNum(), rndNum(), rndNum(),rndNum(),iSp);
                    

              for (int idim=0;idim<3;idim++)
                ParVel_D[idim] = BulkVel[iSp][idim]+uth[idim];
               
                  
              PIC::ParticleBuffer::InitiateParticle(xPar, ParVel_D,&weightCorrection,&iSp,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);

                }
              }
            }
          }
          */

        }
          
    iBlock++;
  }

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("particles prepopulated!\n");
  //fclose(fPartData);
  return nGlobalInjectedParticles;
}



void SetIC() {
  
    int i,j,k;
    char *offset;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    int iBlock=0;
    double cPi = 3.14159265;
    double waveNumber[3]={0.0,0.0,0.0};
    double lambda=32.0;
   
    waveNumber[0]=2*cPi/lambda;
  
    double x[3];
   
    using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
    int nBreak=0;

    printf("User Set IC called\n");


    int nCells[3] ={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
   
    for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {

      if (!node->block) continue;
      //for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      //node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      
      if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
	bool BoundaryBlock=false;
	
	for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
	    //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	    BoundaryBlock=true;
	    break;
	  }
	
	if (BoundaryBlock==true) continue;
      }
      
      double dx[3];
      double *xminBlock= node->xmin, *xmaxBlock= node->xmax;

      for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nCells[idim];
      //printf("dx:%e,%e,%e\n",dx[0],dx[1],dx[2]);
     
      for (k=-1;k<_BLOCK_CELLS_Z_+1;k++) for (j=-1;j<_BLOCK_CELLS_Y_+1;j++) for (i=-1;i<_BLOCK_CELLS_X_+1;i++) {
	    
            int ind[3]={i,j,k};
            
	    PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
	    if (CornerNode!=NULL){
              
	      offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
	     
              
	      for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+ind[idim]*dx[idim];
	  
              double Ex,Ey,Ez,Bx,By,Bz;
              
              Ex = PIC::CPLR::FLUID::FluidInterface.getEx(iBlock,x[0],x[1],x[2]);
              Ey = PIC::CPLR::FLUID::FluidInterface.getEy(iBlock,x[0],x[1],x[2]);
              Ez = PIC::CPLR::FLUID::FluidInterface.getEz(iBlock,x[0],x[1],x[2]);
          

              Bx = PIC::CPLR::FLUID::FluidInterface.getBx(iBlock,x[0],x[1],x[2]);
              By = PIC::CPLR::FLUID::FluidInterface.getBy(iBlock,x[0],x[1],x[2]);
              Bz = PIC::CPLR::FLUID::FluidInterface.getBz(iBlock,x[0],x[1],x[2]);
          
          

              ex0=Ex, ey0=Ey, ez0=Ez;
              
              //Ex=0, Ey=0, Ez=0;
              ((double*)(offset+CurrentEOffset))[ExOffsetIndex]=Ex;
              ((double*)(offset+CurrentEOffset))[EyOffsetIndex]=Ey;
              ((double*)(offset+CurrentEOffset))[EzOffsetIndex]=Ez;
	  
              ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=Ex;
              ((double*)(offset+OffsetE_HalfTimeStep))[EyOffsetIndex]=Ey;
              ((double*)(offset+OffsetE_HalfTimeStep))[EzOffsetIndex]=Ez;
              
              
              ((double*)(offset+OffsetB_corner+CurrentBOffset))[BxOffsetIndex]=Bx;
              ((double*)(offset+OffsetB_corner+CurrentBOffset))[ByOffsetIndex]=By;
              ((double*)(offset+OffsetB_corner+CurrentBOffset))[BzOffsetIndex]=Bz;
                        
              ((double*)(offset+OffsetB_corner+PrevBOffset))[BxOffsetIndex]=Bx;
              ((double*)(offset+OffsetB_corner+PrevBOffset))[ByOffsetIndex]=By;
              ((double*)(offset+OffsetB_corner+PrevBOffset))[BzOffsetIndex]=Bz;
              


              
	    // ((double*)(offset+CurrentCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;


	    
	    //((double*)(offset+NextCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
	    }//
            else{
              printf("corner node is null!\n");
            }
	  }//for (k=0;k<_BLOCK_CELLS_Z_+1;k++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (i=0;i<_BLOCK_CELLS_X_+1;i++) 
      /*
      for (k=-1;k<_BLOCK_CELLS_Z_+1;k++) for (j=-1;j<_BLOCK_CELLS_Y_+1;j++) for (i=-1;i<_BLOCK_CELLS_X_+1;i++) {
            int ind[3]={i,j,k};
	    PIC::Mesh::cDataCenterNode *CenterNode= node->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
	    if (CenterNode!=NULL){
	      offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

              for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];

              double Bx,By,Bz;
              Bx = PIC::CPLR::FLUID::FluidInterface.getBx(iBlock,x[0],x[1],x[2]);
              By = PIC::CPLR::FLUID::FluidInterface.getBy(iBlock,x[0],x[1],x[2]);
              Bz = PIC::CPLR::FLUID::FluidInterface.getBz(iBlock,x[0],x[1],x[2]);
              
              if (i==-1 && j==-1 && k==-1)
                bx0=Bx,by0=By,bz0=Bz;
              // if (fabs(By-4.02e-2)<1e-4) printf("error2, fluid interface\n");
            

              ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=Bx;
              ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=By;
              ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=Bz;
		
		
              ((double*)(offset+PrevBOffset))[BxOffsetIndex]=Bx;
              ((double*)(offset+PrevBOffset))[ByOffsetIndex]=By;
              ((double*)(offset+PrevBOffset))[BzOffsetIndex]=Bz;
            }// if (CenterNode!=NULL)
	  }//for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) 
      */
      iBlock++;
    }

    PIC::BC::ExternalBoundary::UpdateData();
    InterpolateB_N2C(); 

}


double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  return PIC::CPLR::FLUID::dt;
}


double BulletLocalResolution(double *x) {                                                                                           
  // Assume dx = dy = dz
  double dx = PIC::CPLR::FLUID::FluidInterface.getdx(0);
  // Why use 0.1? How about res = res*(1+1e-6)? --Yuxi
  double res=sqrt(3*dx*dx) + 0.1;
  return res;
}
                       
void amps_init_mesh() {
    
  PIC::InitMPI();
  PIC::Init_BeforeParser();

#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  printf("non-uniform mesh!\n");
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_
  printf("uniform mesh!\n");
#endif


#if _CURRENT_MODE_==_PIC_MODE_ON_
  printf("current on!\n");
#endif

#if _CURRENT_MODE_==_PIC_MODE_OFF_
  printf("current mode off!\n");
#endif


  PIC::Mesh::mesh.PopulateOutsideDomainNodesFlag=true;

  //seed the random number generator
  rnd_seed(100);

  //generate mesh or read from file
  char mesh[_MAX_STRING_LENGTH_PIC_]="none";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  sprintf(mesh,"amr.sig=%s.mesh.bin","test_mesh");


  double xMin[3] = {0, 0, 0};

  double xMax[3] = {PIC::CPLR::FLUID::FluidInterface.getphyMax(0) - 
		    PIC::CPLR::FLUID::FluidInterface.getphyMin(0), 
		    PIC::CPLR::FLUID::FluidInterface.getphyMax(1) - 
		    PIC::CPLR::FLUID::FluidInterface.getphyMin(1), 
		    PIC::CPLR::FLUID::FluidInterface.getphyMax(2) - 
		    PIC::CPLR::FLUID::FluidInterface.getphyMin(2)};


  PIC::Mesh::mesh.AllowBlockAllocation=false;
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
    PIC::BC::ExternalBoundary::Periodic::Init(xMin,xMax,BulletLocalResolution);
  }else{
    PIC::Mesh::mesh.init(xMin,xMax,BulletLocalResolution);
  }
  PIC::Mesh::mesh.memoryAllocationReport();

  //generate mesh or read from file
  bool NewMeshGeneratedFlag=false;

  char fullname[STRING_LENGTH];
  sprintf(fullname,"%s/%s",PIC::UserModelInputDataPath,mesh);

  FILE *fmesh=NULL;

  fmesh=fopen(fullname,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(fullname);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh.ThisThread==0) {
      PIC::Mesh::mesh.buildMesh();
      PIC::Mesh::mesh.saveMeshFile("mesh.msh");
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }
  }


  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();
  PIC::Mesh::mesh.InitCellMeasure();

  //experiment of staircase blocks
  /*  
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)  {
    double xmiddle[3];
    bool isOutside=false;
    for (int idim=0;idim<3;idim++){
      xmiddle[idim]=(node->xmin[idim]+node->xmax[idim])/2;
      if (xmiddle[idim]<PIC::Mesh::mesh.xGlobalMin[idim] || xmiddle[idim]>PIC::Mesh::mesh.xGlobalMax[idim])
        isOutside = true;
    }
    if (isOutside)  {
      printf("deallocate block at xmiddle:%e,%e,%e\n",xmiddle[0],xmiddle[1],xmiddle[2]);
      PIC::Mesh::mesh.DeallocateBlock(node);
      node->Thread = -1;
    }
  }
  */

  //coupling send info from amps to fluid
  PIC::CPLR::FLUID::SendCenterPointData.push_back(SendDataToFluid);
}


void amps_init(){

  //  PIC::BC::UserDefinedParticleInjectionFunction=setFixedBC;
  
  PIC::Init_AfterParser();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  if (PIC::ThisThread==0) printf("test1\n");
  //PIC::Mesh::mesh.outputMeshTECPLOT("mesh_test.dat");
  
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
    PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable();
  }

  if (PIC::ThisThread==0) printf("test2\n");
 
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(1,1);

  PIC::DomainBlockDecomposition::UpdateBlockTable();

  //  PIC::BC::UserDefinedParticleInjectionFunction=setFixedBC;   
  
  PIC::FieldSolver::Electromagnetic::ECSIM::setParticle_BC=
    setFixedParticle_BC;
  PIC::FieldSolver::Electromagnetic::ECSIM::setE_half_BC=
    setFixedE_BC_half;
  PIC::FieldSolver::Electromagnetic::ECSIM::setE_curr_BC=
    setFixedE_BC_curr;
  PIC::FieldSolver::Electromagnetic::ECSIM::setB_center_BC=
    setFixedB_center_BC;
  PIC::FieldSolver::Electromagnetic::ECSIM::setB_corner_BC=
    setFixedB_corner_BC;

  //solve the transport equation
  //set the initial conditions for the transport equation
  //  TransportEquation::SetIC(3);
 

  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange();
    break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
    break;
  }
  //PIC::FieldSolver::Init(); 
  PIC::FieldSolver::Electromagnetic::ECSIM::SetIC=SetIC;
    
  
  //PIC::FieldSolver::Init();
  PIC::FieldSolver::Electromagnetic::ECSIM::Init_IC();

 
     
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange();
    break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
    break;
  }
  //PIC::Mesh::mesh.outputMeshDataTECPLOT("ic.dat",0);
  

  int LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
  int GlobalParticleNumber;
  MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("Before cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
  std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;

  CleanParticles();
  LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
  MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("After cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);

  PrepopulateDomain();

  LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
  MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("After prepopulating, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
  std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;
   
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange();
    break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
    break;
  }
    
  //  PIC::Sampling::Sampling();
  PIC::FieldSolver::Electromagnetic::ECSIM::UpdateJMassMatrix();

  cout<<"Calling set_FluidInterface"<<endl;
  PIC::CPLR::FLUID::set_FluidInterface();
  
}


void amps_time_step(){
  
  static int cnt=0;
  static int niter=0;
  printf("output what read from GM\n");                                                                                                         
  
  char fullname[STRING_LENGTH];
  sprintf(fullname,"test_pic_coupler_data_iter=%d.dat",cnt);
  cnt++;
  //PIC::Mesh::mesh.outputMeshDataTECPLOT(fullname,0);                                                                         
  printf("done: output what read from GM\n");                               
  printf(" Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",
         niter++,
         PIC::RequiredSampleLength,
         PIC::RequiredSampleLength-PIC::CollectingSampleCounter);

 
  PIC::TimeStep();
  
 
  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
    PIC::Mesh::mesh.ParallelBlockDataExchange();
    break;

  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
    break;
  }
  
}


