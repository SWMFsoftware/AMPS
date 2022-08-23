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

#SAMPLING
F                       Particle data sampling

#SAMPLING_LENGTH
2                       sampling length



#SAMPLE_OUTPUT_CADENCE
20                      the cadence between starting sampled procedure for output AMPS' data file

#TIMESTEP
1.0                     fixedDt



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

  //exit if the command is not recognized
  bool StrictCommandCheck=false;

  while (param_list.begin()!=param_list.end()) { 
    Command="";

    Command=param_list.front().first;
    cout << "PT: "  << param_list.front().second << endl;
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

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

    }

    else if (Command == "#TEST"){
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      TestVar=(t=="T") ? true : false;

      if (TestVar==true) {
        PIC::ModelTestRun::mode=true;
      }
    }

    else if (Command == "#SAMPLING"){
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      TestVar=(t=="T") ? true : false;

      if (TestVar==false) {
        PIC::SamplingMode=_TEMP_DISABLED_SAMPLING_MODE_;
      }
    }

    else if (Command == "#SAMPLING_LENGTH") {
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      PIC::RequiredSampleLength=atoi(t.c_str());
    }

    else if (Command == "#SAMPLING_OUTPUT_CADENCE") {
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      AMPS2SWMF::SamplingOutputCadence=atoi(t.c_str());
    }


    else if (Command == "#TIMESTEP") { 
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      double dt=atof(t.c_str()); 

      if (_SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_) {
        PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;

        for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=dt; 
      }
      else exit(__LINE__,__FILE__,"Error: the option is not defined"); 
    }

    else if (Command == "#STRICT") {  
      StrictCommandCheck=true; 
    }

    else if (Command == "#FIELDLINE") { 
      string t;

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::ROrigin=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LonMin=atof(t.c_str())*_DEGREE_;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LonMax=atof(t.c_str())*_DEGREE_;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LatMin=atof(t.c_str())*_DEGREE_;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::LatMax=atof(t.c_str())*_DEGREE_;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();


      std::string::size_type sz;


      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::nLon=std::stoi(t,&sz); 
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      AMPS2SWMF::FieldLineData::nLat=std::stoi(t,&sz);
      cout << "PT: " << param_list.front().second << endl;
      param_list.pop_front();
    }
    
    //command related to the SEP model
    else if (Command == "#SEP_PITCH_ANGLE_SCATTERING_FUCTION") {
      #ifdef _SEP_MODEL_ON_
      t=param_list.front().first;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
   
      if (t == "none") {
        SEP::Diffusion::GetPitchAngleDiffusionCoefficient=NULL;
      }
      else if (t=="Jokopii") {
        SEP::Diffusion::GetPitchAngleDiffusionCoefficient=SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient;
        SEP::Diffusion::Jokopii1966AJ::Init();
      }
      else if (t=="Florinskiy") {
        SEP::Diffusion::GetPitchAngleDiffusionCoefficient=SEP::Diffusion::Florinskiy::GetPitchAngleDiffusionCoefficient;
      }
      else {
        exit(__LINE__,__FILE__,"Error: the option is not decognized");
      }
      #endif
    }
    
    else if (Command == "#SEP_PITCH_ANGLE_JOKOPII_DIFFUSION") {
      cout << "PT: "  << param_list.front().second << endl;
          
      #ifdef _SEP_MODEL_ON_
      t=param_list.front().first;

      SEP::Diffusion::Jokopii1966AJ::k_ref_min=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      SEP::Diffusion::Jokopii1966AJ::k_ref_max=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
      
      t=param_list.front().first;
      SEP::Diffusion::Jokopii1966AJ::k_ref_R=atof(t.c_str())*_AU_;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      if (t=="fraction") {
        SEP::Diffusion::Jokopii1966AJ::Mode=SEP::Diffusion::Jokopii1966AJ::_fraction;        
      }
      else if (t=="awsom") {
        SEP::Diffusion::Jokopii1966AJ::Mode=SEP::Diffusion::Jokopii1966AJ::_awsom;
      }
      else exit(__LINE__,__FILE__,"Error: the option is unknown");

      t=param_list.front().first;
      SEP::Diffusion::Jokopii1966AJ::FractionValue=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      SEP::Diffusion::Jokopii1966AJ::FractionPowerIndex=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
      #endif
    }

    else if (Command == "#SEP_PARTICLE_LIMIT") {
      std::string::size_type sz;
      cout << "PT: "  << param_list.front().second << endl;
      
      #ifdef _SEP_MODEL_ON_
      t=param_list.front().first;
      SEP::MinParticleLimit=std::stoi(t.c_str(),&sz);
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      SEP::MaxParticleLimit=atof(t.c_str())*MeV2J;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
      #endif
    }
    
    
    else if (Command == "#SEP_INJECTION_FL") {
      std::string::size_type sz;
     
      #ifdef _SEP_MODEL_ON_
      t=param_list.front().first;
      SEP::FieldLine::InjectionParameters::nParticlesPerIteration=std::stoi(t.c_str(),&sz);
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      SEP::FieldLine::InjectionParameters::emin=atof(t.c_str())*MeV2J;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      t=param_list.front().first;
      SEP::FieldLine::InjectionParameters::emax=atof(t.c_str())*MeV2J;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
      
      t=param_list.front().first;
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      if (t == "beginning") {
        SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectBegginingFL;
      }
      else if (t == "shock") {
        SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectShockLocations;
      }
      else {
        exit(__LINE__,__FILE__,"Error: the option is not recognized");
      }
      #endif    
    }
    
    else if (Command == "#SEP_INJECTION_TYPE_FL") {
      cout << "PT: "  << param_list.front().second << endl;
     
      t=param_list.front().first;
      param_list.pop_front();
      
      #ifdef _SEP_MODEL_ON_
      if (t == "PowerLaw") {
        t=param_list.front().first;
        cout << "PT: "  << param_list.front().second << endl; 
        param_list.pop_front();
        SEP::FieldLine::InjectionParameters::PowerIndex=atof(t.c_str());

        t=param_list.front().first;
        cout << "PT: "  << param_list.front().second << endl;
        param_list.pop_front();
        SEP::FieldLine::InjectionParameters::InjectionEfficiency=atof(t.c_str()); 

        SEP::FieldLine::InjectionParameters::InjectionMomentumModel=SEP::FieldLine::InjectionParameters::_tenishev2005aiaa;
      }
      else if (t=="Sokolov2004AJ") {
        t=param_list.front().first;
        cout << "PT: "  << param_list.front().second << endl;
        param_list.pop_front();
        SEP::FieldLine::InjectionParameters::PowerIndex=atof(t.c_str());

        t=param_list.front().first;
        cout << "PT: "  << param_list.front().second << endl;
        param_list.pop_front();
        SEP::FieldLine::InjectionParameters::InjectionEfficiency=atof(t.c_str());

        SEP::FieldLine::InjectionParameters::InjectionMomentumModel=SEP::FieldLine::InjectionParameters::_sokolov2004aj; 
      }
      else {
        exit(__LINE__,__FILE__,"Error: the option is not recognized");
      } 
      #endif
    }
    
    else if (Command == "#SEP_SAMPLING_LOCATION_FL") {
      double r;
      int nSamplePoints=0;
      
      cout << "PT: "  << param_list.front().second << endl;
      
      t=param_list.front().first;
      nSamplePoints=atoi(t.c_str());

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
      
      #ifdef _SEP_MODEL_ON_
      for (int i=0;i<nSamplePoints;i++) {
        t=param_list.front().first;
        cout << "PT: "  << param_list.front().second << endl;
        param_list.pop_front();

        PIC::Parser::replace(t,"au","149598000.0E3");
        PIC::Parser::replace(t,"rsun","6.96345E8");

        r=PIC::Parser::Evaluate(t);
        SEP::Sampling::SamplingHeliocentricDistanceList.push_back(r);
      }
      #endif
    }

    else if (Command == "#LOCATE_SHOCK") {
      t=param_list.front().first;

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();

      if (t == "disabled") {
        AMPS2SWMF::ShockSearchMode=AMPS2SWMF::_disabled; 
      }
      else if (t=="density_variation") {
        AMPS2SWMF::ShockSearchMode=AMPS2SWMF::_density_variation;
      }
      else if (t=="density_bump") {
        AMPS2SWMF::ShockSearchMode=AMPS2SWMF::_density_bump;
      }
      else if (t=="density_ratio") {
        AMPS2SWMF::ShockSearchMode=AMPS2SWMF::_density_ratio;
      }
      else {
        exit(__LINE__,__FILE__,"Error: the option is not decognized");
      }
    }

    else if (Command == "#MAX_DISTANCE_SHOCK_LOCATOR") {
      t=param_list.front().first;
      AMPS2SWMF::MinShockSpeed=atof(t.c_str());

      PIC::Parser::replace(t,"au","149598000.0E3");
      PIC::Parser::replace(t,"rsun","6.96345E8");

      AMPS2SWMF::ShockLocationsMaxHeliocentricDistance=PIC::Parser::Evaluate(t);

      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
    }

    else if (Command == "#SHOCK_MIN_SPEED") {
      t=param_list.front().first;
      AMPS2SWMF::MinShockSpeed=atof(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
    }

    else if (Command == "#BL_POINT_IMPORT_STEP") {
      cout << "PT: "  << param_list.front().second << endl;

      t=param_list.front().first;
      AMPS2SWMF::bl_import_point_step=atoi(t.c_str()); 
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
    }

    else if (Command == "#BL_OUTPUT_STEP") {
      cout << "PT: "  << param_list.front().second << endl;

      t=param_list.front().first;
      AMPS2SWMF::bl_output_step=atoi(t.c_str());
      cout << "PT: "  << param_list.front().second << endl;
      param_list.pop_front();
    }


    else {
      if ((Command.c_str()[0]!='!')&&(StrictCommandCheck==true)) {
        cout << "PT: Can not find Command : " << Command << endl;
        exit(__LINE__,__FILE__);
      }
    }

  }


  return 0;
}

