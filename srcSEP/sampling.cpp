
#include "sep.h"

SEP::Sampling::cSamplingBuffer **SEP::Sampling::SamplingBufferTable=NULL;
vector<double> SEP::Sampling::SamplingHeliocentricDistanceList;

array_4d<double>  SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable;
array_3d<double>  SEP::Sampling::PitchAngle::PitchAngleRSamplingTable;
array_5d<double>  SEP::Sampling::PitchAngle::DmumuSamplingTable;

array_3d<double>  SEP::Sampling::RadialDisplacement::DisplacementSamplingTable;
array_4d<double>  SEP::Sampling::RadialDisplacement::DisplacementEnergySamplingTable;

array_3d<double>  SEP::Sampling::Energy::REnergySamplingTable;

array_3d<double>  SEP::Sampling::LarmorRadius::SamplingTable;

double SEP::Sampling::PitchAngle::emin=0.01*MeV2J;
double SEP::Sampling::PitchAngle::emax=3000.0*MeV2J;
int SEP::Sampling::PitchAngle::nEnergySamplingIntervals=70;
double SEP::Sampling::PitchAngle::dLogE; 
double SEP::Sampling::MaxSampleEnergy=3000.0*MeV2J;

//sample particle's mean free path
double SEP::Sampling::MeanFreePath::MaxSampledMeanFreePath=3.0*_AU_;
double SEP::Sampling::MeanFreePath::MinSampledMeanFreePath=1.0E-7*_AU_;

int SEP::Sampling::MeanFreePath::nSampleIntervals=100;
double SEP::Sampling::MeanFreePath::dLogMeanFreePath; 
bool SEP::Sampling::MeanFreePath::active_flag=true;
array_3d<double> SEP::Sampling::MeanFreePath::SamplingTable;


//sample particle's Dxx


//sample particle's D\mu\mu 


void SEP::Sampling::Init() {
  namespace FL=PIC::FieldLine; 

  SEP::OutputAMPS::SamplingParticleData::Init();

  PIC::IndividualModelSampling::SamplingProcedure.push_back(Manager);

  //init the sampling buffer table 
  SamplingBufferTable=new cSamplingBuffer* [FL::nFieldLineMax];

  for (int i=0;i<FL::nFieldLineMax;i++) SamplingBufferTable[i]=NULL;

  if (MeanFreePath::active_flag==true) {
    MeanFreePath::dLogMeanFreePath=log(MeanFreePath::MaxSampledMeanFreePath/MeanFreePath::MinSampledMeanFreePath)/MeanFreePath::nSampleIntervals; 

    MeanFreePath::SamplingTable.init(MeanFreePath::nSampleIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
    MeanFreePath::SamplingTable=0.0;
  }



  SEP::Sampling::PitchAngle::dLogE=log(SEP::Sampling::PitchAngle::emax/SEP::Sampling::PitchAngle::emin)/SEP::Sampling::PitchAngle::nEnergySamplingIntervals;
  PitchAngle::PitchAngleREnergySamplingTable.init(PitchAngle::nMuIntervals,SEP::Sampling::PitchAngle::nEnergySamplingIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  PitchAngle::PitchAngleREnergySamplingTable=0.0; 

  PitchAngle::PitchAngleRSamplingTable.init(PitchAngle::nMuIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  PitchAngle::PitchAngleRSamplingTable=0.0;

  PitchAngle::DmumuSamplingTable.init(4,PitchAngle::nMuIntervals,SEP::Sampling::PitchAngle::nEnergySamplingIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  PitchAngle::DmumuSamplingTable=0.0;

  RadialDisplacement::DisplacementSamplingTable.init(RadialDisplacement::nSampleIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax); 
  RadialDisplacement::DisplacementEnergySamplingTable.init(RadialDisplacement::nSampleIntervals,SEP::Sampling::PitchAngle::nEnergySamplingIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);

  RadialDisplacement::DisplacementSamplingTable=0.0;
  RadialDisplacement::DisplacementEnergySamplingTable=0.0;

   
  Energy::REnergySamplingTable.init(SEP::Sampling::PitchAngle::nEnergySamplingIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  Energy::REnergySamplingTable=0.0;

  SEP::Sampling::LarmorRadius::SamplingTable.init(SEP::Sampling::LarmorRadius::nSampleIntervals,PitchAngle::nRadiusIntervals,FL::nFieldLineMax);
  SEP::Sampling::LarmorRadius::SamplingTable=0.0;
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
      int iL,iMu,iR,spec,iE;
      double rLarmor,MiddleX[3],MiddleB;

      for (;Segment!=NULL;Segment=Segment->GetNext()) {
        ptr=Segment->FirstParticleIndex;

	Segment->GetCartesian(MiddleX,0.5);
	MiddleB=QLT1::B(Vector3D::Length(MiddleX)); 


        while (ptr!=-1) {
          PB::byte* ParticleData=PB::GetParticleDataPointer(ptr);
          spec=PB::GetI(ParticleData);

          v[0]=PB::GetVParallel(ParticleData);
          v[1]=PB::GetVNormal(ParticleData);
          v[2]=0.0;

	  rLarmor= PIC::MolecularData::GetMass(spec)*v[1]/fabs(PIC::MolecularData::GetElectricCharge(spec)*MiddleB); 

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

	  SEP::Sampling::Energy::REnergySamplingTable(iE,iR,iFieldLine)+=ParticleWeight;

	  //sample particle Larmor radius 
	  if (isfinite(rLarmor)==true) { 
	    if (rLarmor<1.0) {
               iL=0;
             }
             else if (rLarmor>=SEP::Sampling::LarmorRadius::rLarmorRadiusMax) {
               iL=SEP::Sampling::LarmorRadius::nSampleIntervals-1;
             }
             else {
               iL=(int)(log(rLarmor)/Sampling::LarmorRadius::dLog);
	     }
            
	    SEP::Sampling::LarmorRadius::SamplingTable(iL,iR,iFieldLine)+=ParticleWeight;  
	  }

	  //sample mean free path 
	  if ((SEP::Offset::MeanFreePath!=-1)&&(SEP::Sampling::MeanFreePath::active_flag==true)) {
            double v=*((double*)(ParticleData+SEP::Offset::MeanFreePath));
	    int i;

	    if (isfinite(v)==true) {
	      if (v<SEP::Sampling::MeanFreePath::MinSampledMeanFreePath) {
	        i=0;
	      }
	      else if (v>=SEP::Sampling::MeanFreePath::MaxSampledMeanFreePath) {
	        i=nSampleIntervals-1;
	      }
	      else {
                i=(int)(log(v/SEP::Sampling::MeanFreePath::MinSampledMeanFreePath)/SEP::Sampling::MeanFreePath::dLogMeanFreePath); 
	      }

	      SEP::Sampling::MeanFreePath::SamplingTable(i,iR,iFieldLine)+=ParticleWeight; 
	    }
          }   

	  //sample displacement of a particle from the magnetic field line 
         if (SEP::Offset::RadialLocation!=-1) {
           double r=*((double*)(ParticleData+SEP::Offset::RadialLocation));
	   int iD;

	   if (isfinite(r)==true) {
	     if (r<1.0) {
               iD=0;
  	     }
	     else if (r>=SEP::Sampling::RadialDisplacement::rDisplacementMax) {
	       iD=SEP::Sampling::RadialDisplacement::nSampleIntervals-1;
             }
             else {
               iD=(int)(log(r)/Sampling::RadialDisplacement::dLogDisplacement); 
	     }	   

	     Sampling::RadialDisplacement::DisplacementSamplingTable(iD,iR,iFieldLine)+=ParticleWeight;
	     Sampling::RadialDisplacement::DisplacementEnergySamplingTable(iD,iE,iR,iFieldLine)+=ParticleWeight;
	   }
	 }


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
      if (cnt%12==0) {
        SamplingBufferTable[iFieldLine][i].Output();
      }
    }
  }


  if (cnt%20==0) {
    //output data in a file 
    SEP::Sampling::Energy::Output(cnt);
    SEP::Sampling::RadialDisplacement::OutputDisplacementSamplingTable(cnt);
    SEP::Sampling::RadialDisplacement::OutputDisplacementEnergySamplingTable(cnt);
    SEP::Sampling::LarmorRadius::Output(cnt);
    SEP::Sampling::MeanFreePath::Output(cnt);

    //output sampled pitch angle distribution
    //1. notmalize the distribution
    int iMu,iR,iLine; 
    double summ;

    SEP::Sampling::PitchAngle::PitchAngleRSamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR); 
    SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR); 
    SEP::Sampling::PitchAngle::DmumuSamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread!=0) goto end;

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

        for (int iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
          fprintf(fout,", \"f%i(%e-%e)Mev\"",iLine, 
          SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV,
          SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)*J2MeV); 

          fprintf(fout,", \"D%i(%e-%e)Mev\"",iLine,
          SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV,
          SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)*J2MeV);

          fprintf(fout,", \"dMudT%i(%e-%e)Mev\"",iLine,
          SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV,
          SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)*J2MeV);

          fprintf(fout,", \"dPdT%i(%e-%e)Mev\"",iLine,
          SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV,
          SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)*J2MeV);
        }
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
            double D=0.0,TotalSampledWeight=0.0;           
            double dPdT=0.0,dMuDt=0.0;

            f=0.0,base=0.0; 

            for (di=-1;di<=0;di++) for (dj=-1;dj<=0;dj++) {
              i=iR+di;
              j=iMu+dj;

              if ((i>=0)&&(i<SEP::Sampling::PitchAngle::nRadiusIntervals)&&(j>=0)&&(j<SEP::Sampling::PitchAngle::nMuIntervals)) {
                f+=SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable(j,iE,i,iLine);
                base++;

                D+=SEP::Sampling::PitchAngle::DmumuSamplingTable(0,j,iE,i,iLine);
                TotalSampledWeight+=SEP::Sampling::PitchAngle::DmumuSamplingTable(1,j,iE,i,iLine);

                dMuDt+=SEP::Sampling::PitchAngle::DmumuSamplingTable(2,j,iE,i,iLine);
                dPdT+=SEP::Sampling::PitchAngle::DmumuSamplingTable(3,j,iE,i,iLine);
              }
            }

            if (TotalSampledWeight==0.0) TotalSampledWeight=1.0;

            fprintf(fout,"  %e  %e  %e  %e ",f/base,D/TotalSampledWeight,dMuDt/TotalSampledWeight,dPdT/TotalSampledWeight);
          }
       }

       fprintf(fout,"\n");
     }
   }

   fclose(fout);

end:
   SEP::Sampling::PitchAngle::PitchAngleRSamplingTable=0.0;
   SEP::Sampling::PitchAngle::PitchAngleREnergySamplingTable=0.0; 
   SEP::Sampling::PitchAngle::DmumuSamplingTable=0.0;


   SEP::Sampling::RadialDisplacement::DisplacementSamplingTable=0.0;
   SEP::Sampling::RadialDisplacement::DisplacementEnergySamplingTable=0.0;

   SEP::Sampling::Energy::REnergySamplingTable=0.0;
   SEP::Sampling::LarmorRadius::SamplingTable=0.0;

   SEP::Sampling::MeanFreePath::SamplingTable=0.0;

  }


}



