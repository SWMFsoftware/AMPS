// Normal production stage, independent actual furnished-SPICE coefficient
// oracle. Empty-population clock controls are not observations/validation.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonTransportClock.h"
#include "MoonSurfaceEnvironment.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
void amps_init();
namespace {int checks=0,failed=0;std::ofstream rows,metrics;
void Need(bool b,const std::string&e){if(!b)throw std::runtime_error(e);}
void Check(bool b,const std::string&name){int a=b?1:0,all=0;MPI_Allreduce(&a,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);++checks;if(!all)++failed;if(PIC::ThisThread==0)rows<<name<<','<<(all?"PASS":"FAIL")<<'\n';}
void Metric(unsigned step,const char* key,unsigned index,double actual,double expected){if(PIC::ThisThread==0)metrics<<step<<','<<key<<','<<index<<','<<actual<<','<<expected<<'\n';}
}
int main(int argc,char**argv){
 if(argc!=5)return 2;const std::string mode=argv[3];if(mode!="frozen"&&mode!="evolving")return 2;
 try{
  Moon::Runtime::Configuration c;std::string error;Need(Moon::Runtime::ParseApplicationInput(argv[1],&c,&error),error);
  c.simulationEnvironmentMode=mode;Need(Moon::Runtime::InstallConfiguration(c,&error),error);
  amps_init();SetExitErrorCode(__LINE__,1);PIC::DomainBlockDecomposition::UpdateBlockTable();
  if(PIC::ThisThread==0){rows.open(std::string(argv[2])+".checks.csv");rows<<"check,status\n";metrics.open(std::string(argv[2])+".metrics.csv");metrics<<std::setprecision(17)<<"step,quantity,index,actual,expected\n";}
  const bool common=Moon::CommonTransportTimeSelected();const double initial=Moon::OrbitalMotion::et,dt=.25;
  for(int s=0;s<PIC::nTotalSpecies;++s)PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=dt;
  for(unsigned i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks;++i){auto*n=PIC::DomainBlockDecomposition::BlockTable[i];if(n&&n->block)for(int s=0;s<PIC::nTotalSpecies;++s)n->block->SetLocalTimeStep(dt,s);}
  for(unsigned step=0;step<3;++step){const std::string prefix="step"+std::to_string(step)+"_";
   const double begin=Moon::OrbitalMotion::et,middle=mode=="frozen"?begin:begin+.5*dt,end=mode=="frozen"?begin:begin+dt;
   const auto before=PIC::ParticleBuffer::GetAllPartNum(),sampling=PIC::CollectingSampleCounter,outputs=PIC::DataOutputFileNumber;
   if(common)Need(Moon::BeginCommonTransportStep(begin,&error),error);
   Moon::PrepareProductionParticleEnvironment();
   double sun[6],earth[6],moon[6],middleMoon[6],light=0,transform[6][6],body[6][6];
   spkezr_c("SUN",middle,Moon::Frames::SolarOrbital,"NONE","MOON",sun,&light);
   spkezr_c("EARTH",middle,Moon::Frames::SolarOrbital,"NONE","MOON",earth,&light);
   spkezr_c("MOON",begin,Moon::Frames::Inertial,"NONE","SUN",moon,&light);
   spkezr_c("MOON",middle,Moon::Frames::Inertial,"NONE","SUN",middleMoon,&light);
   sxform_c(Moon::Frames::Inertial,Moon::Frames::SolarOrbital,middle,transform);
   sxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,end,body);
   Need(!failed_c(),"independent furnished-SPICE stage oracle unavailable");
   bool positions=true,velocities=true,moonState=true,heliocentricVelocity=true,bodyRotation=true,bodyDerivative=true;
   for(unsigned k=0;k<3;++k){positions=positions&&Moon::xSun_SO[k]==1000*sun[k]&&Moon::xEarth_SO[k]==1000*earth[k];
    velocities=velocities&&Moon::vSun_SO[k]==1000*sun[k+3]&&Moon::vEarth_SO[k]==1000*earth[k+3];
    moonState=moonState&&Moon::xObject_HCI[k]==1000*moon[k]&&Moon::vObject_HCI[k]==1000*moon[k+3];
    double expected=0;for(unsigned j=0;j<3;++j)expected+=transform[k][j]*(1000*middleMoon[j+3]);
    heliocentricVelocity=heliocentricVelocity&&Moon::vObject_SO_FROZEN[k]==expected;
    Metric(step,"sun_position_m",k,Moon::xSun_SO[k],1000*sun[k]);Metric(step,"earth_position_m",k,Moon::xEarth_SO[k],1000*earth[k]);
    for(unsigned j=0;j<3;++j){bodyRotation=bodyRotation&&Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[k][j]==body[k][j];
     const double derivative=c.particleTracingMover!="legacy"&&mode=="frozen"?0:body[k+3][j];bodyDerivative=bodyDerivative&&Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[k+3][j]==derivative;}
   }
   Check(positions,prefix+"same_epoch_Sun_Earth_positions_SI");Check(velocities,prefix+"same_epoch_Sun_Earth_velocities_SI");
   Check(moonState,prefix+"begin_heliocentric_state_SI");Check(heliocentricVelocity,prefix+"free_velocity_midpoint_rotation");
   long double spin[3][3]{};for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)for(unsigned k=0;k<3;++k)spin[i][j]+=static_cast<long double>(transform[i+3][k])*transform[j][k];
   const long double omega[3]={.5L*(spin[1][2]-spin[2][1]),.5L*(spin[2][0]-spin[0][2]),.5L*(spin[0][1]-spin[1][0])};
   const long double scale=std::hypot(omega[0],omega[1],omega[2]);bool angular=true;
   for(unsigned k=0;k<3;++k){angular=angular&&std::abs(Moon::RotationVector_SO_FROZEN[k]-omega[k])<=256*std::numeric_limits<double>::epsilon()*scale;Metric(step,"angular_velocity_rad_s",k,Moon::RotationVector_SO_FROZEN[k],static_cast<double>(omega[k]));}
   Check(angular,prefix+"independent_Rdot_Rtranspose_angular_velocity");Check(bodyRotation&&bodyDerivative,prefix+"actual_body_rotation_and_declared_derivative_policy");
   Check(Moon::OrbitalMotion::et==end&&Moon::SurfaceEnvironment::ProductionSourceEpochS()==begin,prefix+"begin_middle_end_epoch_policy");
   Check(PIC::ParticleBuffer::GetAllPartNum()==before&&PIC::CollectingSampleCounter==sampling&&PIC::DataOutputFileNumber==outputs,prefix+"stage_does_not_inject_or_sample");
   if(common)Need(Moon::CommitCommonTransportStep(&error),error);
  }
  if(mode=="frozen"){
    double point[3]={5e6,1e6,2e6},velocity[3]={1000,-500,300};auto* node=PIC::Mesh::mesh->findTreeNode(point);bool valid=node!=nullptr;int localOwner=0;
    if(node&&node->Thread==PIC::ThisThread){++localOwner;valid=node->block!=nullptr;
      if(valid){double acceleration[3];Moon::TotalParticleAcceleration(acceleration,_NA_SPEC_,-1,point,velocity,node);
        valid=std::isfinite(acceleration[0])&&std::isfinite(acceleration[1])&&std::isfinite(acceleration[2])&&std::hypot(acceleration[0],acceleration[1],acceleration[2])>0;}}
    int owners=0;MPI_Allreduce(&localOwner,&owners,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    Check(valid&&owners==1,"actual_normal_Na_force_uses_initialized_coefficients_on_allocated_owner");
  }
  Check((mode=="frozen"&&Moon::OrbitalMotion::et==initial)||(mode=="evolving"&&Moon::OrbitalMotion::et==initial+3*dt),"three_stage_external_epoch_evolution");
  if(PIC::ThisThread==0){std::ofstream kernels(std::string(argv[2])+".spice_files.txt");SpiceInt n=0;ktotal_c("ALL",&n);for(SpiceInt i=0;i<n;++i){SpiceChar file[4096],kind[64],source[4096];SpiceInt handle;SpiceBoolean found;kdata_c(i,"ALL",4096,64,4096,file,kind,source,&handle,&found);Need(found,"actual furnished kernel inventory incomplete");kernels<<file<<'\n';}
   std::cout<<"M3 coefficient stage mode="<<mode<<" checks="<<checks<<" status="<<(failed?"FAIL":"PASS")<<'\n';}
  MPI_Finalize();return failed?1:0;
 }catch(const std::exception&e){std::cerr<<"ERROR coefficient stage: "<<e.what()<<'\n';int init=0;MPI_Initialized(&init);if(init)MPI_Abort(MPI_GLOBAL_COMMUNICATOR,2);return 2;}
}
