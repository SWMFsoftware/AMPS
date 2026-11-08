#include "test_harness.h"

// pic_background_update_mode.h is normally included after global.h.  The
// dependency-light fixture supplies only the annotation macros needed by that
// public enum/storage contract; no device implementation is being emulated.
#define _TARGET_DEVICE_
#define _CUDA_MANAGED_
#define _TARGET_HOST_
#include "pic_background_update_mode.h"
#include "pic_geometry_kernels.h"
#include "pic_integration_kernels.h"
#include "pic_deposition_kernels.h"

#include <array>
#include <cmath>
#include <limits>

// The production header deliberately declares storage in pic_datafile.cpp.
// This dependency-light executable supplies the one definition required to
// test the public policy contract without linking the full PIC/MPI library.
namespace PIC { namespace CPLR { namespace DATAFILE {
BackgroundUpdateMode BackgroundUpdatePolicy=BackgroundUpdateMode::FileSchedule;
}}}

AMPS_TEST_CASE("PIC-AUX-RAY-BOX-001",
               "ray-box clipping returns exact exits, face coordinates, and typed invalid cases") {
  const double lo[3]={-1.0,-2.0,-3.0};
  const double hi[3]={ 3.0, 4.0, 5.0};
  const double origin[3]={0.0,0.0,0.0};
  const double direction[3]={2.0,1.0,-0.5};
  double point[3]={0.0,0.0,0.0}, local[2]={0.0,0.0};
  int face=-1;
  AMPS_EXPECT_TRUE(PIC::GeometryKernels::GetRayBoxExit(
      lo,hi,origin,direction,1.0e-14,point,local,face)==
      PIC::GeometryKernels::RayBoxStatus::Success);
  // x reaches +3 first at t=1.5; y=1.5 and z=-0.75.  On an x-face AMPS
  // stores normalized (y,z) coordinates.
  AMPS_EXPECT_EQ(face,1);
  AMPS_EXPECT_NEAR(point[0],3.0,2.0e-14);
  AMPS_EXPECT_NEAR(point[1],1.5,2.0e-14);
  AMPS_EXPECT_NEAR(point[2],-0.75,2.0e-14);
  AMPS_EXPECT_NEAR(local[0],3.5/6.0,2.0e-14);
  AMPS_EXPECT_NEAR(local[1],2.25/8.0,2.0e-14);

  const double outside[3]={-5.0,0.0,0.0};
  const double along_x[3]={1.0,0.0,0.0};
  const PIC::GeometryKernels::RayBoxInterval interval=
      PIC::GeometryKernels::ClipRayToAxisAlignedBox(
          lo,hi,outside,along_x,0.0,true);
  AMPS_EXPECT_TRUE(interval.status==PIC::GeometryKernels::RayBoxStatus::Success);
  AMPS_EXPECT_NEAR(interval.entry_parameter,4.0,0.0);
  AMPS_EXPECT_NEAR(interval.exit_parameter,8.0,0.0);
  AMPS_EXPECT_EQ(interval.entry_face,0);
  AMPS_EXPECT_EQ(interval.exit_face,1);

  const double zero[3]={0.0,0.0,0.0};
  AMPS_EXPECT_TRUE(PIC::GeometryKernels::ClipRayToAxisAlignedBox(
      lo,hi,origin,zero,0.0,true).status==
      PIC::GeometryKernels::RayBoxStatus::InvalidInput);
  const double miss_origin[3]={-5.0,10.0,0.0};
  AMPS_EXPECT_TRUE(PIC::GeometryKernels::ClipRayToAxisAlignedBox(
      lo,hi,miss_origin,along_x,0.0,true).status==
      PIC::GeometryKernels::RayBoxStatus::NoIntersection);
}

AMPS_TEST_CASE("PIC-AUX-COLUMN-CONSTANT-001",
               "production trapezoid kernel integrates constant vector fields exactly") {
  double integral[3]={0.0,0.0,0.0};
  const double value[3]={2.0,-3.0,0.25};
  const double steps[4]={0.125,0.25,0.5,1.125};
  for (double ds:steps)
    AMPS_EXPECT_TRUE(PIC::IntegrationKernels::TrapezoidAccumulate(
        integral,value,value,3,ds));
  const double length=2.0;
  for (int i=0;i<3;++i)
    AMPS_EXPECT_NEAR(integral[i],value[i]*length,4.0e-15);
}

AMPS_TEST_CASE("PIC-AUX-COLUMN-AFFINE-001",
               "production trapezoid kernel integrates affine path data exactly") {
  // f(s)=a+b*s on s in [0,3].  The exact integral is 3a+9b/2 and every
  // trapezoid partition reproduces it without a production-generated oracle.
  const double a[2]={1.5,-2.0}, b[2]={0.75,4.0};
  const double nodes[5]={0.0,0.2,0.9,1.7,3.0};
  double integral[2]={0.0,0.0};
  for (int n=0;n<4;++n) {
    double left[2],right[2];
    for (int i=0;i<2;++i) {
      left[i]=a[i]+b[i]*nodes[n];
      right[i]=a[i]+b[i]*nodes[n+1];
    }
    AMPS_EXPECT_TRUE(PIC::IntegrationKernels::TrapezoidAccumulate(
        integral,left,right,2,nodes[n+1]-nodes[n]));
  }
  for (int i=0;i<2;++i)
    AMPS_EXPECT_NEAR(integral[i],3.0*a[i]+4.5*b[i],2.0e-14);
  AMPS_EXPECT_TRUE(!PIC::IntegrationKernels::TrapezoidAccumulate(
      integral,a,b,2,-1.0));
}

AMPS_TEST_CASE("PIC-AUX-STOPPING-TABLE-001",
               "stopping-power interpolation accepts knots/interiors and rejects invalid limits") {
  const double energy[4]={1.0,2.0,5.0,11.0};
  const double stopping[4]={4.0,7.0,13.0,25.0};
  double value=-99.0;
  using PIC::IntegrationKernels::InterpolateTable;
  using PIC::IntegrationKernels::TableStatus;
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,1.0,value)==
                   TableStatus::Success);
  AMPS_EXPECT_NEAR(value,4.0,0.0);
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,3.5,value)==
                   TableStatus::Success);
  AMPS_EXPECT_NEAR(value,10.0,2.0e-15);
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,11.0,value)==
                   TableStatus::Success);
  AMPS_EXPECT_NEAR(value,25.0,0.0);
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,0.5,value)==
                   TableStatus::BelowValidityRange);
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,12.0,value)==
                   TableStatus::AboveValidityRange);
  const double bad_energy[3]={1.0,2.0,2.0};
  AMPS_EXPECT_TRUE(InterpolateTable(bad_energy,stopping,3,1.5,value)==
                   TableStatus::InvalidTable);
  AMPS_EXPECT_TRUE(InterpolateTable(energy,stopping,4,
      std::numeric_limits<double>::quiet_NaN(),value)==
      TableStatus::InvalidArgument);
}

AMPS_TEST_CASE("PIC-AUX-BACKGROUND-DISPATCH-001",
               "background update ownership distinguishes file schedules from runtime providers") {
  using namespace PIC::CPLR::DATAFILE;
  BackgroundUpdatePolicy=BackgroundUpdateMode::FileSchedule;
  AMPS_EXPECT_TRUE(UsesFileSchedule());
  BackgroundUpdatePolicy=BackgroundUpdateMode::RuntimeProvider;
  AMPS_EXPECT_TRUE(!UsesFileSchedule());
  // Restore the production default so this case is order independent.
  BackgroundUpdatePolicy=BackgroundUpdateMode::FileSchedule;
}

AMPS_TEST_CASE("PIC-AUX-BACKGROUND-INTERPOLATION-001",
               "trilinear background interpolation reproduces affine interface fixtures") {
  const double xi[3]={0.2,0.65,0.4};
  double weight[8];
  AMPS_EXPECT_TRUE(PIC::DepositionKernels::TrilinearWeights(xi,weight)==
                   PIC::DepositionKernels::Status::Success);
  double value=0.0;
  for (int iz=0;iz<2;++iz) for (int iy=0;iy<2;++iy)
    for (int ix=0;ix<2;++ix) {
      const int corner=ix+2*iy+4*iz;
      const double fixture=2.0+3.0*ix-4.0*iy+5.0*iz;
      value+=weight[corner]*fixture;
    }
  AMPS_EXPECT_NEAR(value,2.0+3.0*xi[0]-4.0*xi[1]+5.0*xi[2],2.0e-15);
}
