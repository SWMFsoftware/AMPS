/************************************************************
 * This header file contains definitions of offsets to data *
 * carried by each particle; offsets defined as macros      *
 ************************************************************/

// PROCEDURE TO ADD A NEW PARAMETER:
// 1 - add a macro of form _USE_NEW_PARAMETER_ _PIC_MODE_OFF_/_PIC_MODE_ON_
//     to the set of optional parameters
// 2 - add condition that enables/disables this parameter
// 3 - add the section that defines offset and length of this parameter;
//     NOTE: offset is the sum of lengths of ALL previously defined parameters
//           REGARDLESS of whether they are active or disabled
//=============================================================================
#ifndef _PIC_PARTICLE_DATA_MACRO_DEFINITION_
#define _PIC_PARTICLE_DATA_MACRO_DEFINITION_

// SET OF MANDATORY PARAMETERS CARRIED BY A PARTICLE
//-----------------------------------------------------------------------------
// Mandatory parameter: next particle in the stack
extern int _PIC_PARTICLE_DATA__NEXT_OFFSET_; 

// Mandatory parameter: prev particle in the stack
extern int _PIC_PARTICLE_DATA__PREV_OFFSET_;

// Mandatory parameter: species ID
extern int _PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_;

// Mandatory parameter: velocity of a particle
extern int _PIC_PARTICLE_DATA__VELOCITY_OFFSET_;

// Mandatory parameter: position of a particle
extern int _PIC_PARTICLE_DATA__POSITION_OFFSET_;

//Mabdatory parameter: keep the weight correction factor with the basic parameters
extern int _PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_;

// Total space occupied by mandatory parameters
extern int _PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_;
//-----------------------------------------------------------------------------


// SET OF OPTIONAL PARAMETERS CARRIED BY A PARTICLE
//-----------------------------------------------------------------------------
// below is the list of all parameters that can be carried by a particle;
// all are turned off by default (see picGlobal.dfn),
// user can turn parameter on via input file by enabling corresponding mode
//-----------------------------------------------------------------------------

extern int _PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_;
extern int _PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_;
extern int _PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_;
extern int _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_OFFSET_;
extern int _PIC_PARTICLE_DATA__INJECTION_FACE_OFFSET_;
extern int _PIC_PARTICLE_DATA__WEIGHT_OVER_TIME_STEP_OFFSET_;
extern int _PIC_PARTICLE_DATA__FIELD_LINE_ID_OFFSET_;
extern int _PIC_PARTICLE_DATA__FIELD_LINE_COORD_OFFSET_;

extern int _PIC_PARTICLE_DATA__FULL_DATA_LENGTH_;

void InitBasicParticleOffset();

#endif//_PIC_PARTICLE_DATA_MACRO_DEFINITION_
