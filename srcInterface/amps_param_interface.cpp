//$Id$
//read the SWMF PARAM.in file

#include <iomanip>
#include <sstream>

#include "amps2swmf.h"

using namespace std;

int AMPS2SWMF::PARAMIN::read_paramin(stringstream *param) {//(int argc, char **argv, stringstream *param) {
  string Command;
  bool TestVar=false;

  while (*param) {
    Command="";
    get_next_command(param,&Command);
    if (Command == "") continue;

    if (Command == "#RESTART") {
       PIC::Restart::LoadRestartSWMF=true;

       switch (AMPS2SWMF::ComponentID) {
       case _AMPS_SWMF_PC_:
         sprintf(PIC::Restart::recoverParticleDataRestartFileName,"%s","PC/restartOUT/restart_particle.dat");
         sprintf(PIC::Restart::SamplingData::RestartFileName,"%s","PC/restartOUT/restart_field.dat");
         break;
       case _AMPS_SWMF_PT_:
         sprintf(PIC::Restart::recoverParticleDataRestartFileName,"%s","PT/restartOUT/restart_particle.dat");
         sprintf(PIC::Restart::SamplingData::RestartFileName,"%s","PT/restartOUT/restart_field.dat");
         break;
       default:
         exit(__LINE__,__FILE__,"Error: the option is unlnown");
       }
    }
    else if (Command == "#TEST"){
      read_var(param,"DoTest",   &TestVar);

      if (TestVar==true) {
        PIC::ModelTestRun::mode=true;
      }
    }
    //    if(      Command == "#CASE"){
    //      read_var(param,"Simulation Case",   &Case);
    //      read_var(param,"FieldsInit",        &FieldsInit);
    //      read_var(param,"PartInit",          &PartInit);
    //      read_var(param,"WriteMethod",       &wmethod);
    //      read_var(param,"PoissonCorrection", &PoissonCorrection);
    //      read_var(param,"SimulationName",    &SimName);
    //      read_var(param,"verbose",           &verbose);
    //
    //    }
    //    else if( Command == "#NSYNC"){
    //      int tmp;
    //      read_var(param,"nSync !!!!!!! Not USED !!!!!!!", &tmp);
    //
    //    }
    //    else if( Command == "#UNIFORMSTATE"){
    //      doInitSim = true;
    //      read_var(param,"B0x", &B0x);
    //      read_var(param,"B0y", &B0y);
    //      read_var(param,"B0z", &B0z);
    //      read_var(param,"B1x", &B1x);
    //      read_var(param,"B1y", &B1y);
    //      read_var(param,"B1z", &B1z);
    //
    //    }
    //    else if( Command == "#TIMESTEP"){
    //      doInitSim = true;
    //      read_var(param,"dt", &dt);
    //
    //    }
    else {
      cout<<"Can not find Comand : "<<Command<<endl;
    }

  }


  return 0;
}

