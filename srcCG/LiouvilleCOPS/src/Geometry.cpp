#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

#include "SpiceUsr.h"

#include "Geometry.h"

static void CheckSpice(const char* where) {
  if (failed_c()) {
    SpiceChar msg[1841];

    getmsg_c("LONG",sizeof(msg),msg);
    fprintf(stderr,"SPICE error in %s: %s\n",where,msg);
    exit(EXIT_FAILURE);
  }
}

void cSpiceGeometry::Init(const cInput& input) {
  //SPICE errors are reported by the code (no abort inside the toolkit)
  erract_c("SET",0,(SpiceChar*)"RETURN");
  errprt_c("SET",0,(SpiceChar*)"NONE");

  std::string path=input.GetString("KernelDirectory");

  for (auto& k : input.GetAll("Kernel")) {
    std::string fname=(k.size()>0 && k[0]=='/') ? k : path+"/"+k;
    struct stat buf;

    if (stat(fname.c_str(),&buf)!=0) {
      fprintf(stderr,"Error: SPICE kernel '%s' is not found\n",fname.c_str());
      exit(EXIT_FAILURE);
    }

    furnsh_c(fname.c_str());
    CheckSpice(fname.c_str());
    printf("SPICE: loaded %s\n",fname.c_str());
  }

  BodyFrame=input.GetString("BodyFrame","67P/C-G_CK");
  BodyName=input.GetString("BodyName","CHURYUMOV-GERASIMENKO");
  SpacecraftName=input.GetString("SpacecraftName","ROSETTA");
  SpacecraftFrame=input.GetString("SpacecraftFrame","ROS_SPACECRAFT");
  NGInstrumentID=input.GetInt("NGInstrumentID",-226331);
  RGInstrumentID=input.GetInt("RGInstrumentID",-226332);
  GaugeOffsetSC=input.GetVec3("GaugeOffsetSC",Vec3(0.0,0.0,0.0));

  ReadFOV(NGInstrumentID,FOV_NG);
  ReadFOV(RGInstrumentID,FOV_RG);

  printf("NG (ID %i): FOV %s in frame %s, boresight %g %g %g\n",NGInstrumentID,FOV_NG.Shape.c_str(),FOV_NG.Frame.c_str(),FOV_NG.Boresight.x,FOV_NG.Boresight.y,FOV_NG.Boresight.z);
  printf("RG (ID %i): FOV %s in frame %s, boresight %g %g %g\n",RGInstrumentID,FOV_RG.Shape.c_str(),FOV_RG.Frame.c_str(),FOV_RG.Boresight.x,FOV_RG.Boresight.y,FOV_RG.Boresight.z);
}

void cSpiceGeometry::ReadFOV(int id,cFOV& fov) {
  SpiceChar shape[64],frame[64];
  SpiceDouble bsight[3],bounds[16][3];
  SpiceInt n;

  getfov_c(id,16,sizeof(shape),sizeof(frame),shape,frame,bsight,&n,bounds);
  CheckSpice("getfov_c");

  fov.Shape=shape;
  fov.Frame=frame;
  fov.Boresight=Normalized(Vec3(bsight[0],bsight[1],bsight[2]));
  fov.Bounds.clear();

  for (int i=0;i<n;i++) fov.Bounds.push_back(Normalized(Vec3(bounds[i][0],bounds[i][1],bounds[i][2])));
}

void cSpiceGeometry::Evaluate(const std::string& utc,cObservationGeometry& g) const {
  SpiceDouble et,state[6],lt,pos[3],m[3][3];

  utc2et_c(utc.c_str(),&et);
  CheckSpice("utc2et_c");
  g.et=et;

  //spacecraft state relative to the nucleus in the body-fixed frame (geometric: the gas and the spacecraft are co-located)
  spkezr_c(SpacecraftName.c_str(),et,BodyFrame.c_str(),"NONE",BodyName.c_str(),state,&lt);
  CheckSpice("spkezr_c(spacecraft)");

  Vec3 xRef(1.0E3*state[0],1.0E3*state[1],1.0E3*state[2]);
  g.vSC=Vec3(1.0E3*state[3],1.0E3*state[4],1.0E3*state[5]);

  //gauge position: spacecraft reference point + offset given in the spacecraft frame
  pxform_c(SpacecraftFrame.c_str(),BodyFrame.c_str(),et,m);
  CheckSpice("pxform_c(spacecraft)");

  Vec3 off=GaugeOffsetSC;
  g.xSC=xRef+Vec3(m[0][0]*off.x+m[0][1]*off.y+m[0][2]*off.z,m[1][0]*off.x+m[1][1]*off.y+m[1][2]*off.z,m[2][0]*off.x+m[2][1]*off.y+m[2][2]*off.z);

  //Sun as seen from the nucleus (apparent direction: light time and stellar aberration)
  spkpos_c("SUN",et,BodyFrame.c_str(),"LT+S",BodyName.c_str(),pos,&lt);
  CheckSpice("spkpos_c(SUN)");

  Vec3 xSun(1.0E3*pos[0],1.0E3*pos[1],1.0E3*pos[2]);
  g.HeliocentricDistance=Length(xSun);
  g.SunDirection=xSun/g.HeliocentricDistance;

  //angular velocity of the nucleus: state transformation body-fixed -> J2000, xf=[R 0; dR R];
  //dX/dt=dR x=(dR R^T) X -> omega_J2000 from the skew matrix W=dR R^T, then omega_body=R^T omega_J2000
  SpiceDouble xf[6][6];
  sxform_c(BodyFrame.c_str(),"J2000",et,xf);
  CheckSpice("sxform_c(body)");

  double W[3][3];
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
    W[i][j]=0.0;
    for (int k=0;k<3;k++) W[i][j]+=xf[3+i][k]*xf[j][k];
  }

  Vec3 wJ(0.5*(W[2][1]-W[1][2]),0.5*(W[0][2]-W[2][0]),0.5*(W[1][0]-W[0][1]));
  g.Omega=Vec3(xf[0][0]*wJ.x+xf[1][0]*wJ.y+xf[2][0]*wJ.z,xf[0][1]*wJ.x+xf[1][1]*wJ.y+xf[2][1]*wJ.z,xf[0][2]*wJ.x+xf[1][2]*wJ.y+xf[2][2]*wJ.z);

  //gauge pointing and FOV
  auto Rotate=[&](const cFOV& fov,cGaugeGeometry& gauge) {
    SpiceDouble r[3][3];

    pxform_c(fov.Frame.c_str(),BodyFrame.c_str(),et,r);
    CheckSpice("pxform_c(gauge)");

    auto R=[&](const Vec3& v) {return Vec3(r[0][0]*v.x+r[0][1]*v.y+r[0][2]*v.z,r[1][0]*v.x+r[1][1]*v.y+r[1][2]*v.z,r[2][0]*v.x+r[2][1]*v.y+r[2][2]*v.z);};

    gauge.Boresight=Normalized(R(fov.Boresight));
    gauge.Bounds.clear();
    for (auto& b : fov.Bounds) gauge.Bounds.push_back(Normalized(R(b)));
    gauge.FOVShape=fov.Shape;
  };

  Rotate(FOV_NG,g.NG);
  Rotate(FOV_RG,g.RG);
}

//the FOV pyramid is bounded by the planes through consecutive corner vectors;
//d is inside when it is on the same side of every plane as the boresight
bool InsideFOV(const cGaugeGeometry& gauge,const Vec3& d) {
  if (Dot(d,gauge.Boresight)<=0.0) return false;

  if ((gauge.FOVShape!="RECTANGLE")&&(gauge.FOVShape!="POLYGON")) {
    fprintf(stderr,"Error: FOV shape '%s' is not supported (RECTANGLE and POLYGON only)\n",gauge.FOVShape.c_str());
    exit(EXIT_FAILURE);
  }

  size_t n=gauge.Bounds.size();

  for (size_t i=0;i<n;i++) {
    Vec3 m=Cross(gauge.Bounds[i],gauge.Bounds[(i+1)%n]);

    if (Dot(m,d)*Dot(m,gauge.Boresight)<0.0) return false;
  }

  return true;
}
