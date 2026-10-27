#ifndef _HARTLEY2_RADIAL_MESH_H_
#define _HARTLEY2_RADIAL_MESH_H_

#include <cstddef>

namespace Comet {
namespace RadialModel {

// All distances are in SI units. Exactly one of DomainRadialExtent and
// DomainOuterRadius must be supplied. This keeps the two common meanings of
// "domain size" explicit instead of silently choosing one.
extern double NucleusRadius;
extern double DomainRadialExtent;
extern double DomainOuterRadius;
extern double SurfaceResolution;
extern int RequestedCellNumber;

class cLogarithmicResolutionProfile {
 public:
  bool Init(double nucleusRadius,double outerRadius,double surfaceResolution,
      int requestedCellNumber,const char **errorMessage=NULL);

  double GetResolution(double radius) const;
  double GetExpansionRatio() const {return ExpansionRatio;}
  double GetOuterResolution() const;
  double GetCellWidth(int cell) const;
  double GetCellBoundary(int boundary) const;
  double GetOuterRadius() const {return OuterRadius;}
  int GetCellNumber() const {return CellNumber;}

 private:
  double InnerRadius=0.0;
  double OuterRadius=0.0;
  double InnerResolution=0.0;
  double ExpansionRatio=1.0;
  int CellNumber=0;
};

extern cLogarithmicResolutionProfile ResolutionProfile;

void InitMeshProfile();
double GetDomainOuterRadius();
double LocalResolution(double *x);
long int CountAllocatedRadialCells();
void PrintMeshDiagnostics();

}
}

#endif
