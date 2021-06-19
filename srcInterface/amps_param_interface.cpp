//$Id$
//read the SWMF PARAM.in file

/* Formal of SWMF's PARAM.in COMP PT commands that cabe be parsed here:

#BEGIN_COMP PT ---------------------------------------------------------------

#FIELDLINE
1.5     ROrigin
-175    LonMin
175     LonMax
-85     LatMin
85      LatMax
5       nLon
4       nLat


#TEST 
T


#RESTART


#HELIOSPHERE
2       MinRadius

#END_COMP PT -----------------------------------------------------------------

*/


#include <iomanip>
#include <sstream>

#include "amps2swmf.h"

using namespace std;

int AMPS2SWMF::PARAMIN::read_paramin(list<pair<string,string> >& param_list) {
  string t,Command;
  bool TestVar=false;
  list<pair<string,string> >::iterator it;

  while (param_list.begin()!=param_list.end()) { 
    Command="";

    Command=param_list.front().first;
    cout << param_list.front().second << endl;
    param_list.pop_front();

    if (Command == "#RESTART") {
       PIC::Restart::LoadRestartSWMF=true;

       switch (AMPS2SWMF::ComponentID) {
       case _AMPS_SWMF_PC_:
         sprintf(PIC::Restart::recoverParticleDataRestartFileName,"%s","PC/restartOUT/restart_particle.dat");
         sprintf(PIC::Restart::SamplingData::RestartFileName,"%s","PC/restartOUT/restart_field.dat");
         break;
       case _AMPS_SWMF_PT_:
         sprintf(PIC::Restart::saveParticleDataRestartFileName,"%s","PT/restartOUT/restart_particle.dat");
         sprintf(PIC::Restart::recoverParticleDataRestartFileName,"%s","PT/restartIN/restart_particle.dat");
         sprintf(PIC::Restart::SamplingData::RestartFileName,"%s","PT/restartIN/restart_field.dat");
         break;
       default:
         exit(__LINE__,__FILE__,"Error: the option is unlnown");
       }
    }
    else if (Command == "#HELIOSPHERE") {
      t=param_list.front().first;
      AMPS2SWMF::Heliosphere::rMin=atof(t.c_str())*_RADIUS_(_SUN_);

      cout << param_list.front().second << endl;
      param_list.pop_front();

    }
    else if (Command == "#TEST"){
      t=param_list.front().first;

      cout << param_list.front().second << endl;
      param_list.pop_front();

      TestVar=(t=="T") ? true : false;

      if (TestVar==true) {
        PIC::ModelTestRun::mode=true;
      }
    }
    else if (Command == "#FIELDLINE") { 
      string t;

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::ROrigin=atof(t.c_str());
      cout << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LonMin=atof(t.c_str())*_DEGREE_;
      cout << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LonMax=atof(t.c_str())*_DEGREE_;
      cout << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LatMin=atof(t.c_str())*_DEGREE_;
      cout << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LatMax=atof(t.c_str())*_DEGREE_;
      cout << param_list.front().second << endl;
      param_list.pop_front();


      std::string::size_type sz;


      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::nLon=std::stoi(t,&sz); 
      cout << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::nLat=std::stoi(t,&sz);
      cout << param_list.front().second << endl;
      param_list.pop_front();
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

