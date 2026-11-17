#include "standard_batsrus_wave_adapter.h"
#include "transport_common.h"
#include <cmath>
#include <exception>

namespace SEP { namespace Transport { namespace PICAdapter {
Status ReadStandardBatsrusWaveView(const ParticleContext& context,
    Waves::StandardBatsrus::MagneticBranchView* output) {
  namespace FL=PIC::FieldLine;
  if(!output||!FL::FieldLinesAll||context.state.fieldLineId<0||context.state.fieldLineId>=FL::nFieldLine||
      !std::isfinite(context.state.coordinate)||context.state.coordinate<0.0||
      context.state.coordinate>=FL::FieldLinesAll[context.state.fieldLineId].GetTotalSegmentNumber()||
      !context.segment||FL::FieldLinesAll[context.state.fieldLineId].GetSegment(context.state.coordinate)!=context.segment)
    return Status::Error(StatusCode::InvalidArgument,"standard BATSRUS native wave view requires a valid field-line coordinate/segment/output");
  try {
    const auto& snapshot=Background::SnapshotStore::Instance().AcquireForMover();
    const auto source=Coefficient::ValidateSourceAgainstBackground(Coefficient::SourceMode::Swmf,snapshot.provider(),snapshot.ownership());
    if(!source.ok())return source;
    if(context.snapshotGeneration==0||context.snapshotGeneration!=snapshot.field_line_generation())
      return Status::Error(StatusCode::InvalidParticleState,"standard BATSRUS wave context does not match the immutable imported generation");
    const auto& datum=FL::DatumAtVertexPlasmaWaves;
    auto* first=context.segment->GetBegin();auto* last=context.segment->GetEnd();
    if(!first||!last||datum.offset<0||datum.length!=2)
      return Status::Error(StatusCode::InvalidCoefficient,"standard BATSRUS wave view requires allocated Wave1/Wave2 vertex storage");
    const auto fits=[&](FL::cFieldLineVertex* vertex) {
      const int bytes=vertex->AssociatedDataLength();
      return bytes>=0&&datum.offset<=bytes&&2*sizeof(double)<=static_cast<std::size_t>(bytes-datum.offset)&&vertex->GetAssociatedDataBufferPointer();
    };
    if(!fits(first)||!fits(last))return Status::Error(StatusCode::InvalidCoefficient,
        "standard BATSRUS imported wave datum lies outside allocated native vertex storage");
    // PIC exposes a non-const accessor; bind its returned data immediately as
    // const. The imported shared storage remains total energy for other apps.
    const double* begin=first->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);
    const double* end=last->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);
    if(!begin||!end)return Status::Error(StatusCode::InvalidCoefficient,"standard BATSRUS imported wave endpoint is absent");
    Waves::StandardBatsrus::RawTotalAlfvenEndpoint rawBegin,rawEnd;
    for(unsigned branch=0;branch<2;++branch){rawBegin.waveJPerM3[branch]=begin[branch];rawEnd.waveJPerM3[branch]=end[branch];}
    const double fraction=context.state.coordinate-std::floor(context.state.coordinate);
    return Waves::StandardBatsrus::InterpolateMagneticView(rawBegin,rawEnd,fraction,snapshot.field_line_generation(),output);
  }catch(const std::exception& error) {
    return Status::Error(StatusCode::InvalidCoefficient,std::string("standard BATSRUS read-only wave access rejected: ")+error.what());
  }
}
} } }
