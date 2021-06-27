//functions describing injection at shock

#include "sep.h"

int SEP::ParticleSource::ShockWave::ShockStateFlag_offset=-1;

//condition for presence of a shock in a given cell
bool SEP::ParticleSource::ShockWave::IsShock(PIC::Mesh::cDataCenterNode *CenterNode) {
  double density_current,density_last;
  char *SamplingBuffer; 
  bool flag=false;
  
  const double min_ratio=1.2;

  SamplingBuffer=CenterNode->GetAssociatedDataBufferPointer();

  density_current=*((double*)(SamplingBuffer+PIC::CPLR::SWMF::PlasmaNumberDensityOffset));
  density_last=*((double*)(SamplingBuffer+PIC::CPLR::SWMF::PlasmaNumberDensityOffset_last)); 

  if ((density_last>0.0)&&(density_current>0.0)) {
    flag=(density_current/density_last>min_ratio) ? true : false;
  }

  if (flag==true) {
    *((double*)(SamplingBuffer+ShockStateFlag_offset))=1.0;
  }
  else {
    *((double*)(SamplingBuffer+ShockStateFlag_offset))=0.0;
  }

  return flag;
} 

void SEP::ParticleSource::ShockWave::Output::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double *ptr_flag=(double*)(CenterNode->GetAssociatedDataBufferPointer()+SEP::ParticleSource::ShockWave::ShockStateFlag_offset);

  *ptr_flag=0.0;

  for (int i=0;i<nInterpolationCoeficients;i++) { 
    bool flag=SEP::ParticleSource::ShockWave::IsShock(InterpolationList[i]);

    if (flag==true) {
      *ptr_flag=1.0;
      break;
    }
  }
} 

void SEP::ParticleSource::ShockWave::Output::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double t;

  bool gather_print_data=false;

  if (pipe==NULL) gather_print_data=true;
  else if (pipe->ThisThread==CenterNodeThread) gather_print_data=true;

  if (gather_print_data==true) {
    t=*((double*)(CenterNode->GetAssociatedDataBufferPointer()+SEP::ParticleSource::ShockWave::ShockStateFlag_offset));
  }

  if ((PIC::ThisThread==0)||(pipe==NULL)) {
    if ((CenterNodeThread!=0)&&(pipe!=NULL)) pipe->recv(t,CenterNodeThread);

    fprintf(fout," %e ",t);
  }
  else {
    pipe->send(t);
  }
}


void SEP::ParticleSource::ShockWave::Output::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Is shock\"");
}


