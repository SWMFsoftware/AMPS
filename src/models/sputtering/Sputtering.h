// Physical model for surface sputtering
// HOW TO ADD A NEW MODEL
// Two additions are required:
// 1) add a new namespace(s):
//    - namespace SurfaceType (if needed)
//    - namespace ModelName   (within namespace SurfaceType)
//    source file Sputtering_<SurfaceType>.cpp contains yield table
//    for corresponding surface type models and velocity injection distribution
// 2) add a corresponding branch in:
//    - Sputtering::Init,
//    - Sputtering::GetSputteringYield,
//------------------------------------------------------------------------

#ifndef _SPUTTERING_ 
#define _SPUTTERING_

#include "Sputtering.dfn"
#include "Sputtering_Source.h"

#include "pic.h"


namespace Sputtering {
  // Number and list of sputtered species, populated by Init().  These are
  // declarations, not header-local storage: callers that inspect the public
  // product catalogue must see the same objects used by GetSputteringYield().
  // The former `static` declarations created a different zero-initialized
  // catalogue in every translation unit and could therefore report no
  // products even after the model had been initialized.
  extern int nSpecTarget;
  extern int *SpecTarget;
  // Dimensionless sputtered molecules per incident ion, in SpecTarget order.
  // The buffer is scratch storage owned by this module and is overwritten by
  // each GetSputteringYield() call; it is not thread-safe.
  extern double *Yield;

  // icy surface ==============================================================
  namespace Ice {

    // yield model by Martin Rubin --------------------------------------------
    namespace Rubin {

      // sputtered species implemented in the model
      const int nSpecTarget=1;
      const int SpecTarget[nSpecTarget] = {_H2O_SPEC_};

      double SputteringYield(double speed, 
			     int spec_projectile,
			     int spec_target);
    } 
    // namespace Rubin --------------------------------------------------------

    // yield model by Benjamin Teolis -----------------------------------------
    namespace Teolis {

      // sputtered species implemented in the model
      const int nSpecTarget = 2;
      const int SpecTarget[nSpecTarget] = {_O2_SPEC_, _H2O_SPEC_};

      double SputteringYield(double speed, 
			     int spec_projectile,
			     int spec_target);
    } 
    // namespace Teolis -------------------------------------------------------

    void GetSputteringSpeed(int spec, double& speed, double& weight_correction);
  }
  // namespace Ice ============================================================


  // initialization of sputtering model
  void Init();

  // 1st Master function: 
  //   input: speed and species of a projectile
  //  output: list of species and vector of corresponding yields,
  //          returns number of species sputtered
  int GetSputteringYield(double speed, int spec_projectile,
			 const int* &spec_target, double* &yield);
  
  // 2nd Master function
  //   input: species of a sputtered particle
  //  output: its speed and statistical weight correction
  void GetSputteringSpeed(int spec, double& speed, double& weight_correction);

}

#endif//_SPUTTERING_
