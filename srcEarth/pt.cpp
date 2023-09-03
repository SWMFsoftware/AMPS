//functions used fot particle tracking 

#include "pic.h"
#include "Earth.h"

bool Earth::ParticleTracker::TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  static int cnt=0;
  bool res=false;
  PIC::ParticleTracker::cParticleData *ParticleTrajectoryRecord;

  ParticleTrajectoryRecord=(PIC::ParticleTracker::cParticleData*)(PIC::ParticleTracker::ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==true) {
    res=true;
  }
  else {
    if (cnt>40000) res=false;
    else {
      if (Vector3D::Length(x)<5*_EARTH__RADIUS_) { 
        res=true;
        cnt++;
      }
    }
  }

  return res;
}   

