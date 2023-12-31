//======================================================================
//$Id$
//======================================================================
//the file contains the user defined setting for compilation of the model

//#include "UserDefinition.PIC.dfn"

//the species table
/*
#define _NA_SPEC_     0
#define _NA_PLUS_SPEC_ 1
*/

#undef _NA_SPEC_
#define _NA_SPEC_ 0

#undef _NAPLUS_SPEC_
#define _NAPLUS_SPEC_ 1


//extern unsigned int _NA_SPEC_;
extern int maxLocalBackdroundDensityOffset;

//the plenat
//#define _TARGET_ _MOON_

//debugger mode
#undef _PIC_DEBUGGER_MODE_
#define _PIC_DEBUGGER_MODE_ _PIC_DEBUGGER_MODE_ON_

//define the macro for the user-defined mode's header
#define _PIC__USER_DEFINED__USER_PHYSICAL_MODEL_LIST_ "Comet.h"

//ICES
#undef _PIC_ICES_SWMF_MODE_
#define _PIC_ICES_SWMF_MODE_ _PIC_ICES_MODE_OFF_

#undef _PIC_BACKGROUND_ATMOSPHERE_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_MODE_ _PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_

#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_


//uniform time step in the computational domain
#undef _SIMULATION_TIME_STEP_MODE_
#define _SIMULATION_TIME_STEP_MODE_ _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_

//the method for removing the thermalized particles from the system
//#de fine _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__LOCAL_BACKGROUND_THERMAL_SPEED_
#define _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__ESCAPE_SPEED_

//collision cross section type
#define _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__HARD_SPHERE_

//the volume injection of the particles and the type of injection of model particles into a cell
#undef _PIC_VOLUME_PARTICLE_INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION_MODE_ _PIC_VOLUME_PARTICLE_INJECTION_MODE__OFF_

//the mode of volume injection
#undef _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_  _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE__RATE_DEPENDENT_

//the distriburion of the collision frequency in a cell
#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_


//photolytic reactions
#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_
#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_


#undef _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_
#define _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode) \
    Comet::TotalParticleAcceleration(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);


//Include Mercury model's header
//inlude the header of the COMET DUST model
#undef _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_
#define _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__ON_

//particle sampling
#define _PIC_USER_DEFING_PARTICLE_SAMPLING_(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s) \
    Comet::Sampling::SampleParticleData(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s);


//execute user defined MPI data exchange function after the syncronizatino point at the particle exchange routine
#undef _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_MODE_
#define _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_MODE_ _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_MODE__ON_

#define _PIC__USER_DEFINED__MPI_MODEL_DATA_EXCHANGE_() \
    Comet::ExchangeSurfaceAreaDensity();

//call a user defined function after outputing of the data file to clean all user's sampling buffers
#undef _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA__MODE_
#define _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA__MODE_ _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA__MODE__ON_

#define _PIC__USER_DEFINED__CLEAN_SAMPLING_DATA_() \
    Comet::Sampling::FlushSamplingDataBuffers();

//#undef _PIC_DYNAMIC_LOAD_BALANCING_MODE_
//#define _PIC_DYNAMIC_LOAD_BALANCING_MODE_ _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_

















