/*
 * OH_Coupler.cpp
 *
 *  Created on: Mar 20, 2015
 *      Author: vtenishe
 */
//$Id$
//Coupler functions for OH

#include "OH.h"
double OH::Coupling::TimeAfterCoupling[PIC::nTotalSpecies] = {0.0};

void OH::Coupling::Send(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI) {
  int i0=0,i1=0;
  char vname[200];

  int density_source_offset=0;
  int momentum_source_offset=1;
  int energy_source_offset=4;

  int density_sink_offset=0;
  int momentum_sink_offset=1;
  int energy_sink_offset=4;

  int package_length=5;

/*
  int density_sink_offset=5;
  int momentum_sink_offset=6;
  int energy_sink_offset=9;

  int package_length=10;
*/

  //save the data
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;
  PIC::Mesh::cDataCenterNode *cell;
  char *AssociatedDataPointer;
  int i,j,k,LocalCellNumber;

  for (int pt=0;pt<(*nPoint);pt++) {
    double x[3]={0.0,0.0,0.0}; //the location of the point; to be sure that it has all 3 components

    //find the block
    memcpy(x,Xyz_DI+pt*(*nDimIn),(unsigned int)(*nDimIn)*sizeof(double));
    if ((node=PIC::Mesh::mesh.findTreeNode(x,node))==NULL) exit(__LINE__,__FILE__,"Error: can not find the block");
    if (node->Thread!=PIC::ThisThread) exit(__LINE__,__FILE__,"Error: the point data is located on another processor");

    //find the cell
    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");

    cell=node->block->GetCenterNode(LocalCellNumber);
    AssociatedDataPointer=cell->GetAssociatedDataBufferPointer();

    int ifluid=0,fluid_data_offset=0;

    for (ifluid=0;ifluid<PIC::CPLR::SWMF::nCommunicatedIonFluids;ifluid++) {
      if (PIC::LastSampleLength!=0) {
        // copy density source
        Data_VI[fluid_data_offset+density_source_offset+pt*(*nVarIn)]=
          (*(ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceDensityOffset)))/PIC::LastSampleLength;

        // copy momentum source
        for (int idim=0;idim<3;idim++) Data_VI[fluid_data_offset+idim+momentum_source_offset+pt*(*nVarIn)]=
          (*(idim+3*ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceMomentumOffset)))/PIC::LastSampleLength;

        // copy energy source
        Data_VI[fluid_data_offset+energy_source_offset+pt*(*nVarIn)]=
          (*(ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceEnergyOffset)))/PIC::LastSampleLength;


        // copy density sink 
        Data_VI[fluid_data_offset+density_sink_offset+pt*(*nVarIn)]-=
          (*(ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSinkDensityOffset)))/PIC::LastSampleLength;

        // copy momentum sink 
        for (int idim=0;idim<3;idim++) Data_VI[fluid_data_offset+idim+momentum_sink_offset+pt*(*nVarIn)]-=
          (*(idim+3*ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSinkMomentumOffset)))/PIC::LastSampleLength;

        // copy energy sink 
        Data_VI[fluid_data_offset+energy_sink_offset+pt*(*nVarIn)]-=
          (*(ifluid+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSinkEnergyOffset)))/PIC::LastSampleLength;

      }
      else {
        for (int offset=0;offset<package_length;offset++) Data_VI[fluid_data_offset+offset+pt*(*nVarIn)]=0.0; 
      }

      fluid_data_offset+=package_length;
    }
  }

  // reset time after last coupling session
  for(int spec=0; spec < PIC::nTotalSpecies; spec++)
    OH::Coupling::TimeAfterCoupling[spec] = 0.0;
}



