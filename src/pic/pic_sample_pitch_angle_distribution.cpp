//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the function for sampling of the particle's velocity/evergy distribution function

#include "pic.h"

//const int PIC::PitchAngleDistributionSample::_LINEAR_SAMPLING_SCALE_=0,PIC::PitchAngleDistributionSample::_LOGARITHMIC_SAMPLING_SCALE_=1;
//int PIC::PitchAngleDistributionSample::v2SamplingMode=_LINEAR_SAMPLING_SCALE_,PIC::PitchAngleDistributionSample::speedSamplingMode=_LINEAR_SAMPLING_SCALE_;
double PIC::PitchAngleDistributionSample::CosPAMin=-1.0,PIC::PitchAngleDistributionSample::CosPAMax=1.0;
long int PIC::PitchAngleDistributionSample::nSampledFunctionPoints=201;
double** PIC::PitchAngleDistributionSample::SamplingBuffer=NULL;
double PIC::PitchAngleDistributionSample::SamplingLocations[][3]={{0.0,0.0,0.0}};
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** PIC::PitchAngleDistributionSample::SampleNodes=NULL;
double PIC::PitchAngleDistributionSample::dCosPA=0.01;
long int *PIC::PitchAngleDistributionSample::SampleLocalCellNumber=NULL;
int PIC::PitchAngleDistributionSample::nSampleLocations=0;
bool PIC::PitchAngleDistributionSample::SamplingInitializedFlag=false;

//int PIC::PitchAngleDistributionSample::Sample_Velocity_Offset=0,PIC::PitchAngleDistributionSample::Sample_Speed_Offset=0;
int PIC::PitchAngleDistributionSample::Sample_PitchAngle_Offset=0,
PIC::PitchAngleDistributionSample::SampleDataLength=0;

//====================================================
//init the sampling buffers
void PIC::PitchAngleDistributionSample::Init() {//double ProbeLocations[][DIM],int nProbeLocations) {
  int idim,nProbe,i,j,k;

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif


  if (DIM < 3) exit(__LINE__,__FILE__,"Error: DIM < 3, pitch angle calculations are not meaningful ");

  //  nSampleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lenfths of the sampling intervals
  double t0,t1; //tempotary variables to satisfy intel c++ compiler

  dCosPA=1.001*(CosPAMax-CosPAMin)/(nSampledFunctionPoints-1);

  Sample_PitchAngle_Offset=0;
  SampleDataLength=1;

  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nSampleLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nSampleLocations];

  SamplingBuffer=new double* [nSampleLocations];
  SamplingBuffer[0]=new double [nSampleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1)];

  for (nProbe=1;nProbe<nSampleLocations;nProbe++) {
    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nSampleLocations;nProbe++) {
    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void PIC::PitchAngleDistributionSample::flushSamplingBuffers() {
  long int i,TotalDataLength=nSampleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}
//====================================================
//return the offset where the sample data for the particular specie, sampling interval and the sampling point are located
long int PIC::PitchAngleDistributionSample::GetSampleDataOffset(int spec,int SampleVariableOffset) {
  long int offset;

  offset=spec*SampleDataLength*(nSampledFunctionPoints-1);
  offset+=SampleVariableOffset*(nSampledFunctionPoints-1);

  return offset;
}
//====================================================
//Sample the distribution function
void PIC::PitchAngleDistributionSample::SampleDistributionFnction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,idim,offset;
  double LocalParticleWeight, speed, CosPA;

  for (node=SampleNodes[0],nProbe=0;nProbe<nSampleLocations;node=SampleNodes[++nProbe]) if (node->Thread==PIC::ThisThread) {
    double *v;
    double B[3], Bdir[3], Babs;
    PIC::ParticleBuffer::byte *ParticleData;
    int i,j,k;

    PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
    ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    PIC::CPLR::InitInterpolationStencil(SamplingLocations[nProbe],node);
    PIC::CPLR::GetBackgroundMagneticField(B);

    Babs = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5);
    if (Babs==0.0) continue;

    for(idim=0;idim<DIM; ++idim) Bdir[idim] = B[idim] / Babs;

    while (ptr!=-1) {
      ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      spec=PIC::ParticleBuffer::GetI(ParticleData);
      v=PIC::ParticleBuffer::GetV(ParticleData);
      speed = pow(v[0]*v[0]+v[1]*v[1]+v[2]*v[2], 0.5);

      LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
      LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

      for (CosPA=0.0,idim=0;idim<DIM;idim++) CosPA+=Bdir[idim]*v[idim];

      CosPA /= speed;
      i=(int)((CosPA-CosPAMin)/dCosPA);
      offset=GetSampleDataOffset(spec,Sample_PitchAngle_Offset);

      if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight;

      ptr=PIC::ParticleBuffer::GetNext(ParticleData);
    }
  }
}


//====================================================
//print the distribution function into a file
void PIC::PitchAngleDistributionSample::printDistributionFunction(char *fname,int spec) {
  long int idim,nProbe,i,nVariable,thread,offset;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,dInterval=0.0;
  char str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);


  for (nProbe=0;nProbe<nSampleLocations;nProbe++) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      sprintf(str,"%s.nSamplePoint=%ld.dat",fname,nProbe);
      fout=fopen(str,"w");

      fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);

      fprintf(fout,"TITLE=\"Pitch Angle distribution function at x=%e",SamplingLocations[nProbe][0]);
      for (idim=1;idim<DIM;idim++) fprintf(fout,", %e",SamplingLocations[nProbe][idim]);

      fprintf(fout,"\"\nVARIABLES=\"Cos(PitchAngle)\",\"f\"\n");

      //collect the sampled information from other processors
      for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]+=pipe.recv<double>(thread);
      }

      //normalize the distribution functions
      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        norm=0.0;
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) {
          if (nVariable==Sample_PitchAngle_Offset) dInterval=dCosPA;
          else exit(__LINE__,__FILE__,"Error: unknown option");

          norm+=SamplingBuffer[nProbe][i+offset]*dInterval;
        }

        if (fabs(norm)>0.0) for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]/=norm;
      }

      //print the output file
      for (i=0;i<nSampledFunctionPoints-1;i++) {
        double CosPA=0.0;

        CosPA=CosPAMin+i*dCosPA;
        fprintf(fout,"%e ",CosPA);

        offset=GetSampleDataOffset(spec,Sample_PitchAngle_Offset);
        fprintf(fout,"  %e\n",SamplingBuffer[nProbe][i+offset]);
      }

      //close the output file
      fclose(fout);
      fprintf(PIC::DiagnospticMessageStream,"done.\n");
    }
    else {
      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) {
          pipe.send(SamplingBuffer[nProbe][i+offset]);
          SamplingBuffer[nProbe][i+offset]=0.0; //this sampled information is stored by the root processor
        }
      }
    }

  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

