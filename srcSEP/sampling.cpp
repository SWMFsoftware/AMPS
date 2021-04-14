
#include "sep.h"

int SEP::Sampling::SamplingBufferTableLength=0;
SEP::Sampling::cSamplingBuffer *SEP::Sampling::SamplingBufferTable=NULL;

void SEP::Sampling::Init() {
  PIC::IndividualModelSampling::SamplingProcedure.push_back(Manager);

  //init the samplgin buffer;
  SamplingBufferTableLength=4;

  SamplingBufferTable=new cSamplingBuffer [SamplingBufferTableLength];

  for (int i=0;i<SamplingBufferTableLength;i++) {
    SamplingBufferTable[i].Init("sample",0.1*MeV2J,100.0*MeV2J,5,i*0.25*_AU_,0);
  } 
}  


void SEP::Sampling::Manager() {
  for (int i=0;i<SamplingBufferTableLength;i++) {
    SamplingBufferTable[i].Sampling();
  }

  static int cnt=0;

  cnt++;

  if (cnt%2==0) {
    for (int i=0;i<SamplingBufferTableLength;i++) {
      SamplingBufferTable[i].Output();
    }
  }
}



