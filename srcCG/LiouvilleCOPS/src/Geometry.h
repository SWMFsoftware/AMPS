//observation geometry derived with SPICE: spacecraft state, gauge pointing and fields of view, Sun direction
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include <vector>

#include "Vec3.h"
#include "Input.h"

//field of view of a gauge as defined in the instrument kernel (getfov_c)
struct cFOV {
  std::string Shape,Frame;
  Vec3 Boresight;               //in the instrument frame
  std::vector<Vec3> Bounds;     //corner vectors (RECTANGLE/POLYGON) in the instrument frame
};

struct cGaugeGeometry {
  Vec3 Boresight;               //unit vector in the body-fixed frame
  std::vector<Vec3> Bounds;     //FOV corner vectors in the body-fixed frame
  std::string FOVShape;
};

struct cObservationGeometry {
  double et;
  Vec3 xSC;                     //gauge location in the body-fixed frame [m]
  Vec3 vSC;                     //spacecraft velocity relative to the body-fixed (rotating) frame [m/s]
  Vec3 SunDirection;            //unit vector from the nucleus centre to the Sun (apparent, LT+S)
  double HeliocentricDistance;  //[m]
  Vec3 Omega;                   //angular velocity of the body-fixed frame w.r.t. J2000, in body-fixed coordinates [rad/s]
  cGaugeGeometry NG,RG;
};

class cSpiceGeometry {
public:
  void Init(const cInput& input);
  void Evaluate(const std::string& utc,cObservationGeometry& g) const;

  std::string BodyFrame,BodyName,SpacecraftName,SpacecraftFrame;
  int NGInstrumentID,RGInstrumentID;
  Vec3 GaugeOffsetSC;           //gauge position relative to the spacecraft reference point in the spacecraft frame [m]

private:
  cFOV FOV_NG,FOV_RG;
  void ReadFOV(int id,cFOV& fov);
};

//true when direction 'd' (unit vector from the gauge outward) is inside the FOV
bool InsideFOV(const cGaugeGeometry& gauge,const Vec3& d);

#endif
