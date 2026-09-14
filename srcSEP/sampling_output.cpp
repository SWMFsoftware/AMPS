#include "sep.h"
#include "util/sep_transactional_output.h"

#include <fstream>
#include <sstream>
#include <string>

namespace {

struct StagedSamplingFile {
  FILE* stream = NULL;
  std::string root;
  std::string relative;
  std::string stage;
};

StagedSamplingFile OpenSamplingFile(const std::string& directory,
                                    const std::string& name) {
  StagedSamplingFile file;
  file.root=PIC::OutputDataFileDirectory;
  file.relative=directory+"/"+name;
  file.stage=file.root+"/"+file.relative+".stage";
  const SEP::Transport::Status directoryStatus=
      SEP::Output::EnsureOutputDirectory(file.root,directory);
  if (!directoryStatus.ok())
    exit(__LINE__,__FILE__,directoryStatus.message.c_str());
  file.stream=std::fopen(file.stage.c_str(),"wb");
  if (file.stream==NULL)
    exit(__LINE__,__FILE__,"cannot open staged sampling output");
  return file;
}

void CommitSamplingFile(StagedSamplingFile* file,std::uint64_t records) {
  if (std::fflush(file->stream)!=0 || ::fsync(::fileno(file->stream))!=0 ||
      std::fclose(file->stream)!=0)
    exit(__LINE__,__FILE__,"cannot flush staged sampling output");
  file->stream=NULL;
  std::ifstream input(file->stage.c_str(),std::ios::binary);
  std::ostringstream payload; payload<<input.rdbuf();
  if (!input && !input.eof())
    exit(__LINE__,__FILE__,"cannot read staged sampling output");
  ::unlink(file->stage.c_str());
  SEP::Output::ArtifactMetadata metadata;
  metadata.schemaVersion="srcsep-sampling-v2";
  metadata.recordCount=records;
  metadata.configurationFingerprint=
      SEP::Background::CurrentConfigurationFingerprint();
  const SEP::Output::WriteResult result=SEP::Output::WriteTransactional(
      file->root,file->relative,payload.str(),metadata);
  if (!result.status.ok()) exit(__LINE__,__FILE__,result.status.message.c_str());
}

}  // namespace

void SEP::Sampling::Energy::Output(int cnt) {
  namespace FL=PIC::FieldLine;
  double de,norm;
  int iLine,iE,iR;

  //gather all sampled data
  REnergySamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread!=0) return;

  REnergySamplingTable.find_nan();


  //normalize the energy distribution
  for (iLine=0;iLine<FL::nFieldLineMax;iLine++) if (FL::FieldLinesAll[iLine].IsInitialized()==true)  for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
    norm=0.0;

    for (iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
      de=SEP::Sampling::PitchAngle::emin*exp((iE+1)*SEP::Sampling::PitchAngle::dLogE)-SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE);
      norm+=REnergySamplingTable(iE,iR,iLine)*de;
    }

    if (norm>0.0) {
      for (iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
        REnergySamplingTable(iE,iR,iLine)/=norm;

      }
    }
  }

  REnergySamplingTable.find_nan();

  //output the energy distribution in a file
  const std::string outputName="cnt="+std::to_string(cnt)+".dat";
  StagedSamplingFile output=OpenSamplingFile("EnergyRSample",outputName);
  FILE* fout=output.stream;

  // This legacy table is now explicitly a unit-integral shape dP/dE.  It is
  // never re-scaled to peak one; absolute density/intensity and crossing flux
  // are distinct WP28 products with geometry/time normalizations.
  fprintf(fout,"# schema=srcsep-sampling-v2 product=normalized-shape units=1/MeV\n");
  fprintf(fout,"VARIABLES=\"E [Mev]\"");

   for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
      if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
        for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
          fprintf(fout,", \"F=%i (%.2e-%.2e)AU\"",iLine,
            iR*SEP::Sampling::PitchAngle::dR/_AU_,(iR+1)*SEP::Sampling::PitchAngle::dR/_AU_);
	}
      }
   }

   fprintf(fout,"\n");

   //output the data
   for (iE=0;iE<SEP::Sampling::PitchAngle::nEnergySamplingIntervals;iE++) {
     fprintf(fout," %.2e",SEP::Sampling::PitchAngle::emin*exp(iE*SEP::Sampling::PitchAngle::dLogE)*J2MeV);

     for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
       if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
         for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
           fprintf(fout," %.2e",REnergySamplingTable(iE,iR,iLine));
 	 }
       }
     }

     fprintf(fout,"\n");
   }

   CommitSamplingFile(&output,
       static_cast<std::uint64_t>(SEP::Sampling::PitchAngle::nEnergySamplingIntervals));
}

void SEP::Sampling::LarmorRadius::Output(int cnt) {
  namespace FL=PIC::FieldLine;
  double dL,norm;
  int iLine,iR,iD,iL;

  //gather all sampled data
  SamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread!=0) return;
  SamplingTable.find_nan();

  //normalize the energy distribution
  for (iLine=0;iLine<FL::nFieldLineMax;iLine++) if (FL::FieldLinesAll[iLine].IsInitialized()==true)   for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
    norm=0.0;

    for (iL=0;iL<nSampleIntervals;iL++) {
      dL=exp((iL+1)*dLog)-exp(iL*dLog);

      norm+=SamplingTable(iL,iR,iLine)*dL;
    }


    if (norm>0.0) {
      for (iL=0;iL<nSampleIntervals;iL++) {
        SamplingTable(iL,iR,iLine)/=norm;

      }
    }
  }


  //output the distribution in a file
  const std::string outputName="cnt="+std::to_string(cnt)+".dat";
  StagedSamplingFile output=OpenSamplingFile("LarmorRadius",outputName);
  FILE* fout=output.stream;

  fprintf(fout,"# schema=srcsep-sampling-v2 product=normalized-shape units=1/m\n");
  fprintf(fout,"VARIABLES=\"LarmorRadius [m]\"");

   for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
      if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
        for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
          fprintf(fout,", \"F=%i (%.2e-%.2e)AU\"",iLine,
            iR*SEP::Sampling::PitchAngle::dR/_AU_,(iR+1)*SEP::Sampling::PitchAngle::dR/_AU_);
        }
      }
   }

   fprintf(fout,"\n");

   //output the data
   for (iL=0;iL<nSampleIntervals;iL++) {
     fprintf(fout," %.2e",exp(iL*dLog));

     for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
       if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
         for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
           fprintf(fout," %.2e",SamplingTable(iL,iR,iLine));
         }
       }
     }

     fprintf(fout,"\n");
   }

   CommitSamplingFile(&output,static_cast<std::uint64_t>(nSampleIntervals));
}



void SEP::Sampling::MeanFreePath::Output(int cnt) {
  namespace FL=PIC::FieldLine;
  double d,norm;
  int iLine,iR,i;

  //gather all sampled data
  SamplingTable.reduce(0,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread!=0) return;
  SamplingTable.find_nan();

  //normalize the energy distribution
  for (iLine=0;iLine<FL::nFieldLineMax;iLine++) if (FL::FieldLinesAll[iLine].IsInitialized()==true)   for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
    norm=0.0;

    for (i=0;i<nSampleIntervals;i++) {
      d=MinSampledMeanFreePath*(exp((i+1)*dLogMeanFreePath)-exp(i*dLogMeanFreePath));

      norm+=SamplingTable(i,iR,iLine)*d;
    }


    if (norm>0.0) {
      for (i=0;i<nSampleIntervals;i++) {
        SamplingTable(i,iR,iLine)/=norm;

      }
    }
  }

  //output the distribution in a file
  const std::string outputName="cnt="+std::to_string(cnt)+".dat";
  StagedSamplingFile output=OpenSamplingFile("MeanFreePath",outputName);
  FILE* fout=output.stream;

  fprintf(fout,"# schema=srcsep-sampling-v2 product=normalized-shape units=1/m\n");
  fprintf(fout,"VARIABLES=\"MeanFreePath [AU]\"");

   for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
      if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
        for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
          fprintf(fout,", \"F=%i (%.2e-%.2e)AU\"",iLine,
            iR*SEP::Sampling::PitchAngle::dR/_AU_,(iR+1)*SEP::Sampling::PitchAngle::dR/_AU_);
        }
      }
   }

   fprintf(fout,"\n");

   //output the data
   for (i=0;i<nSampleIntervals;i++) {
     fprintf(fout," %.2e",MinSampledMeanFreePath*exp(i*dLogMeanFreePath)/_AU_);

     for (iLine=0;iLine<FL::nFieldLineMax;iLine++) {
       if (FL::FieldLinesAll[iLine].IsInitialized()==true)  {
         for (iR=0;iR<SEP::Sampling::PitchAngle::nRadiusIntervals;iR++) {
           fprintf(fout," %.2e",SamplingTable(i,iR,iLine));
         }
       }
     }

     fprintf(fout,"\n");
   }

   CommitSamplingFile(&output,static_cast<std::uint64_t>(nSampleIntervals));
}
