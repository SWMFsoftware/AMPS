
#include "sep.h"

SEP::Sampling::cSamplingBuffer **SEP::Sampling::SamplingBufferTable=NULL;
vector<double> SEP::Sampling::SamplingHeliocentricDistanceList;

void SEP::Sampling::Init() {
  namespace FL=PIC::FieldLine; 

  SEP::OutputAMPS::SamplingParticleData::Init();

  PIC::IndividualModelSampling::SamplingProcedure.push_back(Manager);

  //init the sampling buffer table 
  SamplingBufferTable=new cSamplingBuffer* [FL::nFieldLineMax];

  for (int i=0;i<FL::nFieldLineMax;i++) SamplingBufferTable[i]=NULL;
}  

void SEP::Sampling::InitSingleFieldLineSampling(int iFieldLine) {
  char fname[200];

  if (SamplingBufferTable[iFieldLine]!=NULL) return; 

  if (SamplingHeliocentricDistanceList.size()==0) {
    SamplingBufferTable[iFieldLine]=new cSamplingBuffer [SamplingHeliocentricDistanceTableLength];

    for (int i=0;i<SamplingHeliocentricDistanceTableLength;i++) {
      sprintf(fname,"%s/sample",PIC::OutputDataFileDirectory); 

      SamplingBufferTable[iFieldLine][i].Init(fname,MinSampleEnergy,MaxSampleEnergy,nSampleIntervals,SamplingHeliocentricDistanceTable[i],iFieldLine);
    }
  }
  else {
    int size=SamplingHeliocentricDistanceList.size();
    SamplingBufferTable[iFieldLine]=new cSamplingBuffer [size];

    for (int i=0;i<SamplingHeliocentricDistanceList.size();i++) {
      sprintf(fname,"%s/sample",PIC::OutputDataFileDirectory);

      SamplingBufferTable[iFieldLine][i].Init(fname,MinSampleEnergy,MaxSampleEnergy,nSampleIntervals,SamplingHeliocentricDistanceList[i],iFieldLine);
    }
  }
}


void SEP::Sampling::Manager() {
  namespace FL=PIC::FieldLine;

  static int cnt=0;
  cnt++;

  //return if FL::FieldLinesAll not allocated
  if (FL::FieldLinesAll==NULL) return;

  for (int iFieldLine=0;iFieldLine<FL::nFieldLineMax;iFieldLine++) if (FL::FieldLinesAll[iFieldLine].IsInitialized()==true) {
    if (SamplingBufferTable[iFieldLine]==NULL) {
      InitSingleFieldLineSampling(iFieldLine);
    }  

    //sample the field line data 
    int TableSize=SamplingHeliocentricDistanceList.size();
    if (TableSize==0) TableSize=SamplingHeliocentricDistanceTableLength;

    for (int i=0;i<TableSize;i++) {
      SamplingBufferTable[iFieldLine][i].Sampling();

      //output sampled data
      if (cnt%120==0) {
        SamplingBufferTable[iFieldLine][i].Output();
      }
    }
  }
}



