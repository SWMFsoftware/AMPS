#include "bg_provider.h"
#include <limits>

namespace SEP3D {
namespace Background {

void CompleteVectorDerivatives(BackgroundSample* s) {
  // Normalize once. Each column j of gradB is the derivative of the full
  // vector with respect to x_j; contracting it with b gives d|B|/dx_j.
  s->absB=s->B.Norm(); s->bHat=s->B/s->absB;
  Core::Vec3 gradAbsB;
  for (int j=0;j<3;++j) {
    double derivative=0;
    const double b[3]={s->bHat.x,s->bHat.y,s->bHat.z};
    for (int i=0;i<3;++i) derivative+=b[i]*s->gradB(i,j);
    if (j==0) gradAbsB.x=derivative;
    if (j==1) gradAbsB.y=derivative;
    if (j==2) gradAbsB.z=derivative;
  }
  const double dln=s->bHat.Dot(gradAbsB)/s->absB;
  // dln=b dot grad(log|B|). A zero derivative has infinite focusing length,
  // represented by +infinity and checked against the tensor by the validator.
  // div(b)=div(B)/|B|-dln includes nonzero div(B) in prescribed regional fields.
  s->focusingLenM=dln==0 ? std::numeric_limits<double>::infinity() : -1/dln;
  s->divBhat=s->gradB.Trace()/s->absB-dln;
  // Project the along-field B derivative perpendicular to b, then divide by
  // |B| to obtain curvature (b dot grad)b [1/m]. Plasma derivatives use the
  // same tensor convention: trace is divU, double contraction is bb:gradU.
  s->curvature=(s->gradB.Apply(s->bHat)-s->bHat*s->bHat.Dot(gradAbsB))/s->absB;
  s->divU=s->gradU.Trace();
  s->fieldAlignedStrain=s->gradU.DoubleContract(s->bHat,s->bHat);
}


Core::Status BackgroundProvider::EvaluateBatchDetailed(
    const double* xM, const double* yM, const double* zM,
    std::size_t count, BackgroundSample* output,
    Core::Status* perSampleStatus) const {
  if ((count != 0) && (xM == nullptr || yM == nullptr || zM == nullptr ||
                       output == nullptr || perSampleStatus == nullptr)) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "batch background arrays must not be null");
  }

  Core::Status aggregate = Core::Status::OK();
  for (std::size_t i = 0; i < count; ++i) {
    const BackgroundSample candidate = Evaluate({xM[i], yM[i], zM[i]});
    Core::Status sampleStatus = candidate.status;
    if (sampleStatus.ok() && !candidate.valid) {
      sampleStatus = Core::Status(
          Core::StatusCode::BackgroundInvalid,
          "provider returned an invalid sample without an error status");
    }
    perSampleStatus[i] = sampleStatus;
    if (sampleStatus.ok() && candidate.valid) {
      output[i] = candidate;
    } else if (aggregate.ok()) {
      aggregate = sampleStatus;
    }
  }
  return aggregate;
}

}  // namespace Background
}  // namespace SEP3D
