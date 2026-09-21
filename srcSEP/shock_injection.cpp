//functions describing injection at shock

#include "sep.h"

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "amps2swmf.h" 
#endif

int SEP::ParticleSource::ShockWave::ShockStateFlag_offset=-1;
double SEP::ParticleSource::ShockWave::MaxLimitCompressionRatio=3.0;

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

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  flag=false;

  if ((PIC::CPLR::SWMF::PlasmaDivUdXOffset>0)&&(AMPS2SWMF::DivUdXShockLocationThrehold>0.0)) {
    flag=(fabs(*((double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::PlasmaDivUdXOffset)))>AMPS2SWMF::DivUdXShockLocationThrehold); 
  }
  #endif


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

void SEP::ParticleSource::PopulateAllFieldLines() {
  namespace FL = PIC::FieldLine;
  int iFieldLine;

  for (iFieldLine=0;iFieldLine<FL::nFieldLine;iFieldLine++) {
    PopulateFieldLine(iFieldLine);
  } 
}

void SEP::ParticleSource::PopulateFieldLine(int iFieldLine) { 
  namespace FL = PIC::FieldLine;
  int iSegment; 
  FL::cFieldLineSegment* Segment;


  for (iSegment=0,Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment(); iSegment<FL::FieldLinesAll[iFieldLine].GetTotalSegmentNumber(); iSegment++,Segment=Segment->GetNext()) {
    double t_sw_end,t_sw_begin,v_sw_end[3],v_sw_begin[3],n_sw_end,n_sw_begin,v[3],NumberDensity,Temperature,Volume;
    auto VertexBegin=Segment->GetBegin();
    auto VertexEnd=Segment->GetEnd();
  
    VertexBegin->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_begin);
    VertexBegin->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_begin);
    VertexBegin->GetPlasmaVelocity(v_sw_begin);

    VertexEnd->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_end);
    VertexEnd->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_end);
    VertexEnd->GetPlasmaVelocity(v_sw_end);

    NumberDensity=0.5*(n_sw_begin+n_sw_end);
    Temperature=0.5*(t_sw_begin+t_sw_end);

    Volume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(Segment,iFieldLine);  

    for (int idim=0;idim<3;idim++) v[idim]=0.5*(v_sw_begin[idim]+v_sw_end[idim]);

    // The coupled prepopulation request applies to the compiled AMPS species
    // table.  Never assume H_PLUS exists or maps to a particular slot: AMPS'
    // build-time SpeciesList is the immutable authority and PopulateSegment
    // records the actual generated index in every created particle.
    for (int species=0;species<PIC::nTotalSpecies;++species) {
      PIC::FieldLine::PopulateSegment(
          species,NumberDensity,Temperature,v,Volume,iSegment,iFieldLine,200);
    }
 }
}
