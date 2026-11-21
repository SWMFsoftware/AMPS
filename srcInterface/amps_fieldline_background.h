#ifndef AMPS_FIELDLINE_BACKGROUND_IMPLEMENTATION_H
#define AMPS_FIELDLINE_BACKGROUND_IMPLEMENTATION_H

// Included by amps2swmf.cpp only. This owns receive-time immutable derivative
// preparation; particle workers never run collectives or replace its storage.
#include "amps_velocity_gradient.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <map>
#include <vector>

namespace {
std::vector<std::vector<std::array<double,9> > > importedVelocityGradients;
double importedGradientEpochS=std::numeric_limits<double>::quiet_NaN();
std::uint64_t importedVelocityFingerprint=0;
bool importedGradientsReady=false;

// Every PT process participates even when PIC uses MPI_COMM_SELF. Resolve a
// receive/preparation failure collectively before the next data collective or
// publication, so a bad local packet cannot strand peers or partially publish.
void AgreeImportedBackground(const std::string& message) {
  MPI_Comm comm=AMPS2SWMF::FieldLineInputCommunicator;
  if(comm==MPI_COMM_NULL)throw std::logic_error("SWMF field-line input communicator is absent");
  int rank=0;MPI_Comm_rank(comm,&rank);
  int candidate=message.empty()?INT_MAX:rank,first=INT_MAX;
  MPI_Allreduce(&candidate,&first,1,MPI_INT,MPI_MIN,comm);
  if(first==INT_MAX)return;
  std::string reason=rank==first?message:std::string();
  int length=static_cast<int>(reason.size());MPI_Bcast(&length,1,MPI_INT,first,comm);
  reason.resize(length);MPI_Bcast(length?&reason[0]:nullptr,length,MPI_CHAR,first,comm);
  throw std::logic_error(reason);
}

void PrepareImportedVelocityGradients(double epochS,int stride,int lineCount,
    const int* vertexCounts,int variableCount,const double* packet) {
  using Kernel=AMPS::Coupling::VelocityGradient;
  std::string error;std::vector<double> local;
  try {
    if(!std::isfinite(epochS)||stride<=0||lineCount<0||variableCount<13)
      throw std::logic_error("invalid SWMF velocity-gradient packet dimensions/epoch");
    AMPS2SWMF::FieldLineData::SourceLineNamespace(lineCount);
    const double unit=AMPS2SWMF::FieldLineData::CoordinateUnitM();
    for(int line=0;line<lineCount;++line) {
      if(vertexCounts[line]<1||vertexCounts[line]>stride)
        throw std::logic_error("SWMF velocity-gradient vertex count is invalid");
      const auto id=AMPS2SWMF::FieldLineData::SourceGlobalLineId(line,lineCount);
      for(int vertex=0;vertex<vertexCounts[line];++vertex) {
        const std::size_t offset=(static_cast<std::size_t>(line)*stride+vertex)*(variableCount+1);
        // CON_bline declares X/Y/Z in its supplied length unit and U in m/s.
        // Preserve its actual Cartesian frame. Density/temperature do not
        // supply any entry of this full spatial velocity-gradient tensor.
        for(int d=0;d<3;++d)local.push_back(packet[offset+1+d]*unit);
        for(int d=0;d<3;++d)local.push_back(packet[offset+6+d]);
        local.push_back(static_cast<double>(id));local.push_back(packet[offset]);
      }
    }
    if(local.size()>static_cast<std::size_t>(INT_MAX))
      throw std::logic_error("SWMF velocity sample packet exceeds MPI count range");
    for(double value:local)if(!std::isfinite(value))
      throw std::logic_error("SWMF velocity samples contain nonfinite values");
  }catch(const std::exception& exception){error=exception.what();}
  AgreeImportedBackground(error);

  // The spatial fit is an instantaneous field. Reject mixing different host
  // epochs across independent PIC ranks before any global sample assembly.
  double earliest=0.0,latest=0.0;
  MPI_Allreduce(&epochS,&earliest,1,MPI_DOUBLE,MPI_MIN,AMPS2SWMF::FieldLineInputCommunicator);
  MPI_Allreduce(&epochS,&latest,1,MPI_DOUBLE,MPI_MAX,AMPS2SWMF::FieldLineInputCommunicator);
  AgreeImportedBackground(earliest==latest?std::string():"PT ranks supplied different physical background epochs");

  // Gather the actual received spatial samples, not an inferred radial field.
  // Canonical global-line/Lagrangian-node ordering removes MPI decomposition
  // from nearest-neighbor ties and QR summation order for this fixed packet.
  const MPI_Comm comm=AMPS2SWMF::FieldLineInputCommunicator;
  int processes=0;MPI_Comm_size(comm,&processes);
  const int localCount=static_cast<int>(local.size());std::vector<int> counts(processes),offsets(processes);
  MPI_Allgather(&localCount,1,MPI_INT,counts.data(),1,MPI_INT,comm);
  long long total=0;
  for(int rank=0;rank<processes;++rank){offsets[rank]=static_cast<int>(total);total+=counts[rank];
    if(total>INT_MAX){error="global SWMF velocity samples exceed MPI count range";break;}}
  AgreeImportedBackground(error);
  std::vector<double> gathered(static_cast<std::size_t>(total));
  MPI_Allgatherv(local.data(),localCount,MPI_DOUBLE,gathered.data(),counts.data(),offsets.data(),MPI_DOUBLE,comm);
  std::vector<std::array<double,8> > rows(gathered.size()/8);
  for(std::size_t i=0;i<rows.size();++i)std::copy(gathered.begin()+8*i,gathered.begin()+8*i+8,rows[i].begin());
  std::sort(rows.begin(),rows.end(),[](const std::array<double,8>& a,const std::array<double,8>& b){return a[6]<b[6]||(a[6]==b[6]&&a[7]<b[7]);});
  std::uint64_t fingerprint=UINT64_C(1469598103934665603);
  for(const auto& row:rows)for(double value:row){unsigned char bytes[sizeof(double)];std::memcpy(bytes,&value,sizeof(value));
    for(unsigned char byte:bytes){fingerprint^=byte;fingerprint*=UINT64_C(1099511628211);}}
  if(importedGradientsReady&&fingerprint==importedVelocityFingerprint){importedGradientEpochS=epochS;return;}

  std::vector<std::vector<std::array<double,9> > > staged(lineCount);
  try {
    std::vector<Kernel::Sample> samples;std::map<std::pair<double,double>,std::size_t> lookup;
    std::map<double,unsigned> alongLine;
    for(const auto& row:rows) {
      if(!lookup.emplace(std::make_pair(row[6],row[7]),samples.size()).second)
        throw std::logic_error("SWMF velocity packet repeats a global line/Lagrangian node");
      Kernel::Sample sample;sample.line=static_cast<std::uint64_t>(row[6]);
      for(int d=0;d<3;++d){sample.x[d]=row[d];sample.u[d]=row[3+d];}
      samples.push_back(sample);++alongLine[row[6]];
    }
    bool radialExtent=false;for(const auto& item:alongLine)radialExtent=radialExtent||item.second>1;
    // BL_init initially supplies only seed-surface points, before actual line
    // extraction. That bootstrap is not a valid three-dimensional medium.
    // Leave the cache unavailable; an attempted particle query still rejects.
    if(!radialExtent){importedGradientsReady=false;return;}
    if(importedGradientsReady&&epochS==importedGradientEpochS&&AMPS2SWMF::CoupledStep.Sequence()!=0)
      throw std::logic_error("SWMF changed velocity/geometry at the same already-consumed background epoch");
    Kernel kernel(samples);
    for(int line=0;line<lineCount;++line) {
      staged[line].reserve(vertexCounts[line]);
      const double id=static_cast<double>(AMPS2SWMF::FieldLineData::SourceGlobalLineId(line,lineCount));
      for(int vertex=0;vertex<vertexCounts[line];++vertex) {
        const std::size_t offset=(static_cast<std::size_t>(line)*stride+vertex)*(variableCount+1);
        const auto at=lookup.find(std::make_pair(id,packet[offset]));
        if(at==lookup.end())throw std::logic_error("local SWMF velocity node is absent from the collective packet");
        const auto result=kernel.Evaluate(at->second);
        if(!result.valid)throw std::logic_error("cannot reconstruct received Cartesian velocity gradient: "+result.message);
        staged[line].push_back(result.gradient);
      }
    }
  }catch(const std::exception& exception){error=exception.what();}
  AgreeImportedBackground(error);
  importedVelocityGradients.swap(staged);importedGradientEpochS=epochS;
  importedVelocityFingerprint=fingerprint;importedGradientsReady=true;
}
}

bool AMPS2SWMF::FieldLineData::VelocityGradientAt(int line,int vertex,double fraction,
    double epochS,std::array<double,9>* output,std::string* error) {
  // The particle phase only interpolates the prepared immutable endpoint
  // tensors. No fitting, MPI operation, background replacement or epoch
  // extrapolation is allowed from a worker. Failure leaves output untouched.
  if(!output||!importedGradientsReady||epochS!=importedGradientEpochS||line<0||
     static_cast<std::size_t>(line)>=importedVelocityGradients.size()||vertex<0||
     static_cast<std::size_t>(vertex+1)>=importedVelocityGradients[line].size()||
     !std::isfinite(fraction)||fraction<0||fraction>1) {
    if(error)*error="prepared SWMF Cartesian velocity gradient does not match the requested immutable epoch/segment";
    return false;
  }
  std::array<double,9> value;
  for(int i=0;i<9;++i)value[i]=(1-fraction)*importedVelocityGradients[line][vertex][i]+
      fraction*importedVelocityGradients[line][vertex+1][i];
  *output=value;return true;
}
#endif
