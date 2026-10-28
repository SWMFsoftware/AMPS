#include "CometModel.h"

#include <cmath>
#include <sstream>

#include "pic.h"

// Defaults preserve the historical Hartley 2 three-dimensional case.  The
// input preprocessor replaces them for a radial or future axisymmetric case.
int Comet::Model::GeometryMode=Comet::Model::Geometry::Triangulated;
double Comet::Model::HeliocentricDistanceAU=-1.0;
char Comet::Model::SurfaceMeshFileName[_MAX_STRING_LENGTH_PIC_]="hartley2.shape.nas";
double Comet::Model::SurfaceMeshLengthUnit=1.0e3;

// These arrays are the case-neutral handoff from the input parser to geometry
// backends.  Zero rates and Disabled distributions are safe defaults; a
// runnable source must be supplied explicitly by the selected comet input.
double Comet::Model::SurfaceProductionRate[PIC::nTotalSpecies]={0.0};
int Comet::Model::SurfaceDistributionMode[PIC::nTotalSpecies]={0};
double Comet::Model::SurfaceTemperature[PIC::nTotalSpecies]={0.0};
double Comet::Model::SurfaceBulkSpeed[PIC::nTotalSpecies]={0.0};
double Comet::Model::ExtendedProductionRate[PIC::nTotalSpecies]={0.0};
int Comet::Model::ExtendedDistributionMode[PIC::nTotalSpecies]={0};

namespace {
bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}
}

bool Comet::Model::Validate(std::string *error,bool requireRunnableSources) {
  if (GeometryMode<Geometry::Spherical || GeometryMode>Geometry::Triangulated)
    return Fail("unknown comet geometry mode",error);
  if (!std::isfinite(HeliocentricDistanceAU) || HeliocentricDistanceAU<=0.0)
    return Fail("the heliocentric distance must be finite and positive",error);
  if (GeometryMode==Geometry::Triangulated && SurfaceMeshFileName[0]=='\0')
    return Fail("triangulated geometry requires a surface mesh file",error);
  if (!std::isfinite(SurfaceMeshLengthUnit) || SurfaceMeshLengthUnit<=0.0)
    return Fail("the surface mesh length unit must be finite and positive",error);

  for (int species=0;species<PIC::nTotalSpecies;species++) {
    if (!std::isfinite(SurfaceProductionRate[species]) ||
        SurfaceProductionRate[species]<0.0 ||
        !std::isfinite(ExtendedProductionRate[species]) ||
        ExtendedProductionRate[species]<0.0) {
      std::ostringstream out;
      out << "species " << species << " has an invalid production rate";
      return Fail(out.str(),error);
    }

    if (SurfaceDistributionMode[species]<SurfaceDistribution::Disabled ||
        SurfaceDistributionMode[species]>SurfaceDistribution::LegacyJet ||
        ExtendedDistributionMode[species]<ExtendedDistribution::Disabled ||
        ExtendedDistributionMode[species]>ExtendedDistribution::UserDefined) {
      std::ostringstream out;
      out << "species " << species << " has an unknown source distribution";
      return Fail(out.str(),error);
    }

    if (SurfaceProductionRate[species]>0.0 &&
        SurfaceDistributionMode[species]==SurfaceDistribution::Disabled) {
      std::ostringstream out;
      out << "species " << species
          << " has a surface rate but no surface distribution";
      return Fail(out.str(),error);
    }
    if (ExtendedProductionRate[species]>0.0 &&
        ExtendedDistributionMode[species]==ExtendedDistribution::Disabled) {
      std::ostringstream out;
      out << "species " << species
          << " has an extended rate but no extended distribution";
      return Fail(out.str(),error);
    }

    if (requireRunnableSources && SurfaceProductionRate[species]>0.0 &&
        (!(SurfaceTemperature[species]>0.0) ||
         !std::isfinite(SurfaceTemperature[species]) ||
         SurfaceBulkSpeed[species]<0.0 ||
         !std::isfinite(SurfaceBulkSpeed[species]))) {
      std::ostringstream out;
      out << "species " << species
          << " requires a positive surface temperature and a non-negative bulk speed";
      return Fail(out.str(),error);
    }
  }

  return true;
}

const char *Comet::Model::GeometryName() {
  switch (GeometryMode) {
  case Geometry::Spherical: return "spherical";
  case Geometry::Axisymmetric: return "axisymmetric";
  case Geometry::Triangulated: return "triangulated";
  }

  return "unknown";
}
