
#include "sep.h"

SEP::Sampling::cSamplingBuffer **SEP::Sampling::SamplingBufferTable=NULL;
vector<double> SEP::Sampling::SamplingHeliocentricDistanceList;

array_4d<double>  SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable;
array_3d<double>  SEP::Sampling::PitchAngle::PitchAngleRSamplingTable;

double SEP::Sampling::PitchAngle::emin=0.1*MeV2J;
double SEP::Sampling::PitchAngle::emax=3000.0*MeV2J;
int SEP::Sampling::PitchAngle::nEnergySamplingIntervals=10;
double SEP::Sampling::PitchAngle::dLogE; 

void SEP::Sampling::Init() {
  namespace FL=PIC::FieldLine; 

  SEP::OutputAMPS::SamplingParticleData::Init();

  PIC::IndividualModelSampling::SamplingProcedure.push_back(Manager);

  //init the sampling buffer table 
  SamplingBufferTable=new cSamplingBuffer* [FL::nFieldLineMax];

  for (int i=0;i<FL::nFieldLineMax;i++) SamplingBufferTable[i]=NULL;

  SEP::Sampling::PitchAngle::dLogE=log(SEP::Sampling::PitchAngle::emax/SEP::Sampling::PitchAngle::emin)/SEP::Sampling::PitchAngle::nEnergySamplingIntervals;
  PitchAngle::PitchAngleREnergySamplingTable.init(PitchAngle::nMuIntervals,SEP::Sampling::PitchAngle::nEnergySamplingIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  PitchAngle::PitchAngleREnergySamplingTable=0.0; 

  PitchAngle::PitchAngleRSamplingTable.init(PitchAngle::nMuIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  PitchAngle::PitchAngleRSamplingTable=0.0;

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
  namespace PB = PIC::ParticleBuffer;

  static int cnt=0;
  cnt++;

  //return if FL::FieldLinesAll not allocated
  if (FL::FieldLinesAll==NULL) return;

  for (int iFieldLine=0;iFieldLine<FL::nFieldLineMax;iFieldLine++) if (FL::FieldLinesAll[iFieldLine].IsInitialized()==true) {
    if (SamplingBufferTable[iFieldLine]==NULL) {
      InitSingleFieldLineSampling(iFieldLine);
    }

      //sample pitchaanglesof the partcles for a given field line 
      FL::cFieldLineSegment* Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();
      double speed,e,x[3],v[3],mu,FieldLineCoord;
      long int ptr; 
      int iMu,iR,spec,iE;

      for (;Segment!=NULL;Segment=Segment->GetNext()) {
        ptr=Segment->FirstParticleIndex;

        while (ptr!=-1) {
          PB::byte* ParticleData=PB::GetParticleDataPointer(ptr);
          spec=PB::GetI(ParticleData);

          v[0]=PB::GetVParallel(ParticleData);
          v[1]=PB::GetVNormal(ParticleData);
          v[2]=0.0;

          speed=sqrt(v[0]*v[0]+v[1]*v[1]);
          if (speed>0.99*SpeedOfLight) speed=0.99*SpeedOfLight;
          e=Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec)); 
          iE=log(e/SEP::Sampling::PitchAngle::emin)/SEP::Sampling::PitchAngle::dLogE; 

          if (iE>=SEP::Sampling::PitchAngle::nEnergySamplingIntervals) iE=SEP::Sampling::PitchAngle::nEnergySamplingIntervals-1; 
          if (iE<0)iE=0;

          mu=v[0]/sqrt(v[0]*v[0]+v[1]*v[1]);
          iMu=(int)((mu+1.0)/SEP::Sampling::PitchAngle::dMu);
          if (iMu>=SEP::Sampling::PitchAngle::nMuIntervals) iMu=SEP::Sampling::PitchAngle::nMuIntervals-1; 


          FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
          FL::FieldLinesAll[iFieldLine].GetCartesian(x,FieldLineCoord);

          iR=(int)(Vector3D::Length(x)/SEP::Sampling::PitchAngle::dR);
          if (iR>=SEP::Sampling::PitchAngle::nRadiusIntervals) iR=SEP::Sampling::PitchAngle::nRadiusIntervals-1;

          double ParticleWeight= PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
          ParticleWeight*=PB::GetIndividualStatWeightCorrection(ParticleData);

          SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable(iMu,iE,iR,iFieldLine)+=ParticleWeight;
          SEP::Sampling::PitchAngle::PitchAngleRSamplingTable(iMu,iR,iFieldLine)+=ParticleWeight;

          ptr=PB::GetNext(ptr);
        }
    //  }




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


  if (cnt%120==0) {
    //output sampled pitch angle distribution
    //1. notmalize the distribution
    int iMu,iR,iLine; 
    double summ;

    for (iLine=0;iLine<FL::nFieldLineMax;iLine++) if (FL::FieldLinesAll[iLine].IsInitialized()==true)  for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
      summ=0.0;

      for (iMu=0;iMu<SEP::Sampling::PitchAngle::nMuIntervals;iMu++) {
        summ+=SEP::Sampling::PitchAngle::PitchAngleRSamplingTable(iMu,iR,iLine);
      }

      summ*=SEP::Sampling::PitchAngle::dMu;  
      if (summ==0.0) summ=1.0;

      for (iMu=0;iMu<SEP::Sampling::PitchAngle::nMuIntervals;iMu++) {
        SEP::Sampling::PitchAngle::PitchAngleRSamplingTable(iMu,iR,iLine)/=summ;
      }

      int iE;

      for (iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
        summ=0.0;

        for (iMu=0;iMu<SEP::Sampling::PitchAngle::nMuIntervals;iMu++) {
          summ+=SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable(iMu,iE,iR,iLine);
        }

        summ*=SEP::Sampling::PitchAngle::dMu;
        if (summ==0.0) summ=1.0;

        for (iMu=0;iMu<SEP::Sampling::PitchAngle::nMuIntervals;iMu++) {
          SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable(iMu,iE,iR,iLine)/=summ;
        }
      }
    }

    //output a  file  
    FILE *fout;
    char fname[200];

    sprintf(fname,"%s/PitchAngleRSample.cnt=%i.dat",PIC::OutputDataFileDirectory,cnt);
    fout=fopen(fname,"w");

    fprintf(fout,"VARIABLES=\"R [AU]\", \"Mu\"");

    for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
      if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
        fprintf(fout,", \"f%i\"",iLine); 

        for (int iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) fprintf(fout,", \"f%i(%e-%e)Mev\"",iLine, 
          SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV,
          SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)*J2MeV); 
      }
    }

    fprintf(fout,"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",SEP::Sampling::PitchAngle::nMuIntervals+1,SEP::Sampling::PitchAngle::nRadiusIntervals+1);

    for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals+1;iR++) {
      for (iMu=0;iMu<SEP::Sampling::PitchAngle::nMuIntervals+1;iMu++) {

        fprintf(fout,"%e  %e  ",iR*SEP::Sampling::PitchAngle::dR/_AU_,iMu*SEP::Sampling::PitchAngle::dMu-1.0);

        for (iLine=0;iLine<FL::nFieldLineMax;iLine++) if (FL::FieldLinesAll[iLine].IsInitialized()==true) {
          double f=0.0,base=0.0;
          int di,dj,i,j;
 
          for (di=-1;di<=0;di++) for (dj=-1;dj<=0;dj++) {
            i=iR+di;
            j=iMu+dj;

            if ((i>=0)&&(i<SEP::Sampling::PitchAngle::nRadiusIntervals)&&(j>=0)&&(j<SEP::Sampling::PitchAngle::nMuIntervals)) {
              f+=SEP::Sampling::PitchAngle::PitchAngleRSamplingTable(j,i,iLine);
              base++;
            }  
          }
      
          fprintf(fout,"  %e",f/base);


          for (int iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
            f=0.0,base=0.0; 

            for (di=-1;di<=0;di++) for (dj=-1;dj<=0;dj++) {
              i=iR+di;
              j=iMu+dj;

              if ((i>=0)&&(i<SEP::Sampling::PitchAngle::nRadiusIntervals)&&(j>=0)&&(j<SEP::Sampling::PitchAngle::nMuIntervals)) {
                f+=SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable(j,iE,i,iLine);
                base++;
              }
            }

            fprintf(fout,"  %e",f/base);
          }
       }

       fprintf(fout,"\n");
     }
   }

   fclose(fout);
   SEP::Sampling::PitchAngle::PitchAngleRSamplingTable=0.0;
   SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable=0.0; 

  }


}



