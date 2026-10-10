

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "constants.h"
#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "adapters/reduced_shock_background_adapter.h"
#include "util/sep_background_runtime.h"
#include "transport_common.h"
#include "turbulence_production_adapter.h"
#include "util/sep_run_configuration.h"
#include "util/sep_initialization.h"
#include "sep.dfn"
#include "tests.h"

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "amps2swmf.h"
#endif


const double FieldLineRequestedLength=40;

//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;
double rSphere=_RADIUS_(_TARGET_);


const double xMaxDomain=_DOMAIN_SIZE_*_AU_/_RADIUS_(_SUN_);

const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100/2.5,dxMaxSphere=DebugRunMultiplier*2.0/10.0;

const double MarkNotUsedRadiusLimit=100.0;

void install_reduced_background_on_field_lines(double epochS);

namespace {

namespace fs = std::filesystem;

// Create all initialization-product parents once on rank zero and synchronize
// the result before AMPS enters its distributed Tecplot writer.  Input-deck
// paths and --initialization-output-dir therefore share exactly one filesystem
// policy, and shared filesystems never see an avoidable many-rank mkdir race.
void EnsureInitializationOutputParents(
    const SEP::Initialization::Configuration& initialization) {
  int directoryStatus = 1;
  std::string rootMessage;
  if (PIC::ThisThread == 0) {
    const std::string paths[] = {
        initialization.meshTecplotFile,
        initialization.fieldLineTecplotFile,
        initialization.dataTecplotFile};
    for (const std::string& pathText : paths) {
      const fs::path parent = fs::path(pathText).parent_path();
      if (parent.empty()) continue;
      std::error_code error;
      fs::create_directories(parent, error);
      const bool isDirectory = fs::is_directory(parent, error);
      if (error || !isDirectory) {
        directoryStatus = 0;
        rootMessage = "cannot create initialization output directory '" +
            parent.string() + "'";
        if (error) rootMessage += ": " + error.message();
        break;
      }
    }
  }
  MPI_Bcast(&directoryStatus, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (directoryStatus == 0) {
    const std::string message = PIC::ThisThread == 0
        ? rootMessage : "rank zero could not create initialization output directory";
    exit(__LINE__, __FILE__, message.c_str());
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

// AMPS writes local time step and particle weight for one DataSetNumber at a
// time. Preserve the configured file name for a single-species executable;
// mixed compiled SpeciesList tables receive deterministic sibling files so
// every initialized species remains directly inspectable.
std::string InitializationDataPath(const std::string& base, int species,
                                   int speciesCount) {
  if (speciesCount == 1) return base;
  const fs::path path(base);
  const std::string name = path.stem().string() + ".species-" +
      std::to_string(species) + path.extension().string();
  return (path.parent_path() / name).string();
}

// Return the shock radius associated with the background state that is
// currently visible to particle transport.  Keeping this translation at the
// AMPS application boundary is important: the turbulence adapter accepts SI
// metres and must not know whether the radius came from the analytic shock or
// from the prepared SWCME provider.
double CurrentShockRadiusM() {
  switch (SEP::ShockModelType) {
    case SEP::cShockModelType::Analytic1D:
      return SEP::ParticleSource::ShockWave::Tenishev2005::rShock;
    case SEP::cShockModelType::SwCme1d:
      return SEP::SW1DAdapter::ShockRadiusM();
  }

  // The enum currently has only the two cases above.  Returning NaN keeps a
  // future unsupported model fail-closed with respect to shock injection: the
  // turbulence transaction will still run, but it cannot invent a displacement
  // interval for an unknown shock representation.
  return std::numeric_limits<double>::quiet_NaN();
}

struct ShockRadiusHistory {
  bool initialized = false;
  double previousRadiusM = std::numeric_limits<double>::quiet_NaN();
};

// Install one explicit startup timestep and statistical weight on every
// allocated leaf block.  AMPS normally derives these quantities from local
// resolution and boundary rates; schema-v2 initialization instead treats both
// as reviewed campaign inputs, so re-deriving either would silently change the
// requested physical normalization.
void InstallConfiguredParticleNumerics(
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node, int species,
    double timeStepS, double particleWeight) {
  if (node == NULL) return;
  if (node->lastBranchFlag() == _BOTTOM_BRANCH_TREE_) {
    if (node->block != NULL) {
      node->block->SetLocalTimeStep(timeStepS, species);
      node->block->SetLocalParticleWeight(particleWeight, species);
    }
    return;
  }
  for (int child = 0; child < (1 << DIM); ++child)
    InstallConfiguredParticleNumerics(
        node->downNode[child], species, timeStepS, particleWeight);
}

}  // namespace

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  if (node->IsUsedInCalculationFlag==false) return 0.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
   //delete all particles that was not reflected on the surface
   //PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}



void amps_init_mesh() {
  PIC::InitMPI();

  // The first file-driven workflow intentionally defines one Parker spiral.
  // Reject compile-time domain variants explicitly instead of accepting an
  // input whose line parameters would then be ignored by a straight/FLAMPA or
  // random multi-line branch.
  if (SEP::Initialization::HasActive() &&
      (SEP::DomainType != SEP::DomainType_ParkerSpiral ||
       SEP::Domain_nTotalParkerSpirals != 1)) {
    exit(__LINE__, __FILE__,
         "file-driven initialization requires one Parker-spiral domain");
  }

  //set up the conversion factor for output of the magnetic field line length
  PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.first=1.0/_AU_;
  PIC::FieldLine::cFieldLine::OutputLengthConversionFactor.second="AU";

  // Step 5 makes field-line sampling the only srcSEP production sampling
  // path.  Observer locations remain three-dimensional through their owning
  // field lines; the removed spatial-neighborhood sampler belongs with the
  // separate Cartesian-transport application.
  SEP::Sampling::Init();

  //init the solar wind model when SEP adiabatic cooling is accounted for in 3D modeling
  if (SEP::AccountAdiabaticCoolingFlag==true) {
    SEP::SolarWind::Init();
  }

  SEP::Init();

  // The reduced provider is being coupled as a background-only authority.
  // Leaving the historical callback installed would create ordinary baseline
  // SEP particles even though the event fingerprint says particle_mode is
  // disabled.  This conditional changes no legacy path: the exact historical
  // injection function remains installed whenever no reduced event is active.
  PIC::BC::UserDefinedParticleInjectionFunction=SEP::ReducedShock::Enabled()
      ? NULL : SEP::FieldLine::InjectParticles;

  // Every srcSEP build uses field lines, so reserve the complete current and
  // previous vertex state without retaining a dead Cartesian-mode branch.
  {
      using namespace PIC::FieldLine;

      VertexAllocationManager.MagneticField=true;
      VertexAllocationManager.ElectricField=true;
      VertexAllocationManager.PlasmaVelocity=true;
      VertexAllocationManager.PlasmaDensity=true;
      VertexAllocationManager.PlasmaTemperature=true;
      VertexAllocationManager.PlasmaPressure=true;
      VertexAllocationManager.MagneticFluxFunction=true;
      VertexAllocationManager.PlasmaWaves=true;
      VertexAllocationManager.Fluence=true;
      VertexAllocationManager.ShockLocation=true;
      VertexAllocationManager.DistanceToShockLocation=true;


      VertexAllocationManager.PreviousVertexData.MagneticField=true;
      VertexAllocationManager.PreviousVertexData.ElectricField=true;
      VertexAllocationManager.PreviousVertexData.PlasmaVelocity=true;
      VertexAllocationManager.PreviousVertexData.PlasmaDensity=true;
      VertexAllocationManager.PreviousVertexData.PlasmaTemperature=true;
      VertexAllocationManager.PreviousVertexData.PlasmaPressure=true;
      VertexAllocationManager.PreviousVertexData.PlasmaWaves=true;


      PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;
   }

  PIC::Init_BeforeParser();
  SEP::RequestParticleData();

  //SetUp the alarm
//  PIC::Alarm::SetAlarm(2000);


  rnd_seed();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


  if (PIC::CPLR::SWMF::BlCouplingFlag==false) {
    //reserve data for magnetic filed
    PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
    //PIC::CPLR::DATAFILE::Offset::MagneticField.active=true;
    //PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    //PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=PIC::CPLR::DATAFILE::Offset::MagneticField.nVars*sizeof(double);
  }

  //init the Mercury model
 ////::Init_BeforeParser();
//  PIC::Init_BeforeParser();

//  ProtostellarNebula::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;
///  ProtostellarNebula::Init_AfterParser();



  //register the sphere
  {
    double sx0[3]={0.0,0.0,0.0};
    if (SEP::Initialization::HasActive()) {
      const SEP::Initialization::Configuration& initialization =
          SEP::Initialization::Active();
      sx0[0] = initialization.parkerOriginM.x;
      sx0[1] = initialization.parkerOriginM.y;
      sx0[2] = initialization.parkerOriginM.z;
      rSphere = initialization.innerRadiusM;
    }
    cInternalBoundaryConditionsDescriptor SphereDescriptor;
    cInternalSphericalData *Sphere;

    //correct radiust of the  sphere to be consistent with the location of the inner boundary of the SWMF/SC
    //use taht in case of coupling to the SWMF
    #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    if (AMPS2SWMF::Heliosphere::rMin<0.0) AMPS2SWMF::Heliosphere::rMin=1.05*_RADIUS_(_SUN_);

    if (AMPS2SWMF::Heliosphere::rMin>0.0) rSphere=AMPS2SWMF::Heliosphere::rMin;
    #endif

    //reserve memory for sampling of the surface balance of sticking species
    long int ReserveSamplingSpace[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


    cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);

    PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
    SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
    Sphere->SetSphereGeometricalParameters(sx0,rSphere);

    //set the innber bounday sphere
    SEP::InnerBoundary=Sphere;

    // Inner-boundary sphere surface mesh and surface-data Tecplot files.
    //
    // PIC::OutputDataFileDirectory is itself a char[_MAX_STRING_LENGTH_PIC_]
    // array, so "<directory>/SpheraData.dat" can need up to
    // _MAX_STRING_LENGTH_PIC_+15 bytes: the former sprintf into a buffer of
    // exactly _MAX_STRING_LENGTH_PIC_ bytes could overflow it (GCC
    // -Wformat-overflow).  The paths are built with std::string, which sizes
    // them exactly; truncation (snprintf) is avoided because it would write
    // the products to a different file.  File names and locations are
    // unchanged.  Both printers take const char* and use the name only during
    // the call, so c_str() of the temporaries is valid for the call duration.
    const std::string outputDirectory(PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceMesh((outputDirectory + "/Sphere.dat").c_str());
    Sphere->PrintSurfaceData((outputDirectory + "/SpheraData.dat").c_str(),0);


    Sphere->localResolution=SEP::Mesh::localSphericalSurfaceResolution;
    Sphere->InjectionRate=SEP::ParticleSource::InnerBoundary::sphereInjectionRate;
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=ParticleSphereInteraction;

    if (!SEP::ReducedShock::Enabled()&&
        (_DOMAIN_GEOMETRY_!=_DOMAIN_GEOMETRY_BOX_)&&
        (_SPHERICAL_SHOCK_INJECTION_!=_PIC_MODE_ON_)) {
      Sphere->InjectionBoundaryCondition=SEP::ParticleSource::InnerBoundary::sphereParticleInjection;
    }

    Sphere->PrintTitle=SEP::Sampling::OutputSurfaceDataFile::PrintTitle;
    Sphere->PrintVariableList=SEP::Sampling::OutputSurfaceDataFile::PrintVariableList;
    Sphere->PrintDataStateVector=SEP::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

    //set up the planet pointer in Mercury model
    SEP::Planet=Sphere;
    Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);
  }

  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  cout << "Init the mesh" << endl;

  int idim;
  double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

  for (idim=0;idim<DIM;idim++) {
    xmax[idim]=xMaxDomain*_RADIUS_(_SUN_);
    xmin[idim]=-xMaxDomain*_RADIUS_(_SUN_);
  }
  if (SEP::Initialization::HasActive()) {
    const SEP::Initialization::Configuration& initialization =
        SEP::Initialization::Active();
    const double origin[3] = {initialization.parkerOriginM.x,
                              initialization.parkerOriginM.y,
                              initialization.parkerOriginM.z};
    for (idim=0; idim<DIM; ++idim) {
      xmin[idim] = origin[idim] - initialization.outerRadiusM;
      xmax[idim] = origin[idim] + initialization.outerRadiusM;
    }
  }

  //generate the magneric field line
  list<SEP::cFieldLine> field_line,field_line_old,field_line_new;
  double xStart[3]={1.1,0.0,0.0};


  if (PIC::CPLR::SWMF::BlCouplingFlag==false) switch (SEP::DomainType) {
  case SEP::DomainType_ParkerSpiral:
    PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;

    PIC::FieldLine::VertexAllocationManager.PlasmaWaves=true;
    PIC::FieldLine::VertexAllocationManager.MagneticField=true;
    PIC::FieldLine::VertexAllocationManager.PlasmaVelocity=true;

    if (SEP::Domain_nTotalParkerSpirals==1) {
      if (SEP::Initialization::HasActive()) {
        const SEP::Initialization::Configuration& initialization =
            SEP::Initialization::Active();
        // Schema 1 predates the canonical SWCME normalization bridge and
        // retains its historical signed +5 nT radial normalization.  The
        // complete schema-2 contract instead consumes the canonically derived
        // signed Br(1 AU), including its configured polarity.
        double radialFieldAtOneAuT = 5.0e-9;
        if (initialization.schemaVersion >= 2) {
          radialFieldAtOneAuT =
              SEP::SW1DAdapter::GetConfigurationSummary()
                  .parker_radial_field_at_one_au_t;
        }
        const double origin[3] = {initialization.parkerOriginM.x,
                                  initialization.parkerOriginM.y,
                                  initialization.parkerOriginM.z};
        const double initial[3] = {initialization.parkerInitialPointM.x,
                                   initialization.parkerInitialPointM.y,
                                   initialization.parkerInitialPointM.z};
        SEP::ParkerSpiral::CreateFileLine(
            &field_line, origin, initial, initialization.parkerLengthM,
            initialization.parkerPointCount,
            initialization.solarWindSpeedMPerS,
            initialization.solarRotationRateRadPerS,
            radialFieldAtOneAuT);
      }
      else {
        SEP::ParkerSpiral::CreateFileLine(
            &field_line,xStart,FieldLineRequestedLength*215.0);
      }
      SEP::Mesh::ImportFieldLine(&field_line);

      PIC::FieldLine::Init();
      SEP::Mesh::InitFieldLineAMPS(&field_line);
    }
    else {
      PIC::FieldLine::Init();

      rnd_seed(10);

      for (int iline=0;iline<SEP::Domain_nTotalParkerSpirals;iline++) {
        double x0[3],r,phi;
        double phi_max=45.0*Pi/180.0;

        r=Vector3D::Length(xStart);

        phi=phi_max*rnd();
        if (rnd()<0.5) phi=-phi;

        x0[0]=r*sin(phi);
        x0[1]=r*cos(phi);
        x0[2]=0.0;


	//create a randomly located in 3D initial point
        double cos_phi_max=cos(phi_max);

	do {
          Vector3D::Distribution::Uniform(x0);
	}
        while (x0[0]<phi_max);

        for (int i=0;i<3;i++) x0[i]*=r;



        field_line.clear();

        SEP::ParkerSpiral::CreateFileLine(&field_line,x0,7*250.0);
        SEP::Mesh::ImportFieldLine(&field_line);

        SEP::Mesh::InitFieldLineAMPS(&field_line);
      }

      //init frame of references related to each segment of the field line
     // for (int iFieldLine=0; iFieldLine<PIC::FieldLine::nFieldLine; iFieldLine++) {
     //   PIC::FieldLine::FieldLinesAll[iFieldLine].InitReferenceFrame();
     // }


      if (PIC::ThisThread==0) PIC::FieldLine::Output("all-field-lines.dat",false);
    }

    break;
  case SEP::DomainType_StraitLine:
    PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;

    PIC::FieldLine::VertexAllocationManager.PlasmaWaves=true;
    PIC::FieldLine::VertexAllocationManager.MagneticField=true;
    PIC::FieldLine::VertexAllocationManager.PlasmaVelocity=true;

    field_line.clear();

    SEP::ParkerSpiral::CreateStraitFileLine(&field_line,xmin,250.0);
    SEP::Mesh::ImportFieldLine(&field_line);

    // A straight line is still a three-dimensionally embedded field line and
    // must be registered with AMPS for every supported production mover.
    SEP::Mesh::InitFieldLineAMPS(&field_line);

    break;
  case SEP::DomainType_FLAMPA_FieldLines:
    PIC::ParticleBuffer::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;

    SEP::Mesh::LoadFieldLine_flampa(&field_line_old,"FieldLineOld.in");
    SEP::Mesh::ImportFieldLine(&field_line_old);
    SEP::Mesh::PrintFieldLine(&field_line_old,"FieldLineOld.dat");

    SEP::Mesh::LoadFieldLine_flampa(&field_line_new,"FieldLineNew.in");
    SEP::Mesh::ImportFieldLine(&field_line_new);
    SEP::Mesh::PrintFieldLine(&field_line_new,"FieldLineNew.dat");


    //init the field line library in AMPS
    PIC::FieldLine::VertexAllocationManager.PlasmaWaves=true;
    PIC::FieldLine::VertexAllocationManager.MagneticField=true;
    PIC::FieldLine::VertexAllocationManager.PlasmaVelocity=true;

    PIC::FieldLine::DatumAtVertexPlasmaWaves.length=2;

    PIC::FieldLine::Init();

    SEP::Mesh::InitFieldLineAMPS(&field_line_old);
    SEP::Mesh::InitFieldLineAMPS(&field_line_new);
    break;
  default :
    exit(__LINE__,__FILE__,"Error: the domain type is not recognized");
  }


  //init domain decomposition of the field lines
  PIC::ParallelFieldLines::StaticDecompositionFieldLineLength(0.005);

  //refining the mesh along a set of magnetic field lines: use onle wher model SEP
  if (_MODEL_CASE_==_MODEL_CASE_SEP_TRANSPORT_) {
    PIC::Mesh::mesh->UserNodeSplitCriterion=SEP::Mesh::NodeSplitCriterion;
  }


  //generate only the tree
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,SEP::Mesh::localResolution);
  PIC::Mesh::mesh->memoryAllocationReport();


/*
  if (PIC::Mesh::mesh->ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh->readMeshFile("mesh.msh");
  }
*/

  PIC::Mesh::mesh->buildMesh();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh->ThisThread << endl;

  //PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();



  //mark no-used black that are far from the magnetic filed line
  list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> not_used_list;
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*>*)> MarkNotUsed;
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> GetMaxBlockRefinmentLevel;

  int cnt=0,MaxRefinmentLevel=-1;


  GetMaxBlockRefinmentLevel=[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      double x[3];

      for (int idim=0;idim<3;idim++) x[idim]=0.5*(startNode->xmin[idim]+startNode->xmax[idim]);

       if (Vector3D::Length(x)>10.0*_RADIUS_(_SUN_)) {
         if (MaxRefinmentLevel<startNode->RefinmentLevel) MaxRefinmentLevel=startNode->RefinmentLevel;
       }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        GetMaxBlockRefinmentLevel(downNode);
      }
    }
  };


  MarkNotUsed=[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,list <cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*> *not_used_list) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //if (startNode->xmin[0]<0.0) not_used_list->push_back(startNode);

      if (cnt%PIC::nTotalThreads!=PIC::ThisThread) return;

      double r,x[3],d2min=-1.0,dmax;
      int idim;

      for (idim=0;idim<3;idim++) x[idim]=0.5*(startNode->xmin[idim]+startNode->xmax[idim]);

      r=Vector3D::Length(x);

      if (r>MarkNotUsedRadiusLimit*_RADIUS_(_SUN_)) {
        list<SEP::cFieldLine>::iterator it;
        double d2,t;

        for (it=field_line.begin();it!=field_line.end();it++) {
          for (idim=0,d2=0.0;idim<3;idim++) {
            t=min(fabs(startNode->xmin[idim]-it->x[idim]),fabs(startNode->xmax[idim]-it->x[idim]));
            d2+=t*t;
          }

          if ((d2min<0.0)||(d2min>d2)) d2min=d2;
        }

        dmax=MarkNotUsedRadiusLimit*_RADIUS_(_SUN_)+MarkNotUsedRadiusLimit/200.0*(r-50.0*_RADIUS_(_SUN_));

        if (d2min>dmax*dmax) {
          not_used_list->push_back(startNode);
        }
        else {
          if (startNode->RefinmentLevel<MaxRefinmentLevel-2) not_used_list->push_back(startNode);
        }
      }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        MarkNotUsed(downNode,not_used_list);
      }
    }
  };

  GetMaxBlockRefinmentLevel(PIC::Mesh::mesh->rootTree);

  if (!SEP::Initialization::HasActive() &&
      _DOMAIN_GEOMETRY_!= _DOMAIN_GEOMETRY_BOX_)  {
    MarkNotUsed(PIC::Mesh::mesh->rootTree,&not_used_list);

    PIC::Mesh::mesh->SetTreeNodeActiveUseFlag(&not_used_list,NULL,false,NULL);
  }

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  if (SEP::Initialization::HasActive() &&
      SEP::Initialization::Active().schemaVersion >= 2) {
    const SEP::Initialization::Configuration& initialization =
        SEP::Initialization::Active();
    // outputMeshTECPLOT owns the distributed AMR serialization.  The finite
    // field line is a rank-independent ordered zone and is therefore written
    // once by rank zero after the final active-node pruning is complete.
    EnsureInitializationOutputParents(initialization);
    PIC::Mesh::mesh->outputMeshTECPLOT(
        initialization.meshTecplotFile.c_str());
    if (PIC::ThisThread == 0) {
      const SEP::Transport::Status written =
          SEP::Initialization::WriteParkerLineTecplot(
              initialization, initialization.fieldLineTecplotFile);
      if (!written.ok())
        exit(__LINE__, __FILE__, written.message.c_str());
    }
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }

  //initialize the blocks
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  //init the volume of the cells'
  PIC::Mesh::mesh->InitCellMeasure();


}

void amps_init() {



//init the PIC solver
  PIC::Init_AfterParser ();
	PIC::Mover::Init();

  // Runtime mover dispatch.  The srcSEP AMPS deck defines
  // _PIC_PARTICLE_MOVER_LEGACY_SETTINGS_=_PIC_MODE_OFF_, so generic PIC calls
  // PIC::Mover::UserDefinedParticleMover instead of the compile-time macro.
  // Install srcSEP's production dispatcher there; it runs whichever of
  // parker / fte-dmumu / fte-mfp the application parser selected
  // (--particle-mover -> SEP::Mover::SelectProductionMover, or the registry
  // default).  The deck's _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_ macro
  // still names ::SEP::ParticleMover, the same dispatcher, for the boundary
  // injection paths that invoke the macro directly, so both routes advance a
  // particle with the same selected mover.  With legacy settings ON the
  // pointer is installed but unused, preserving the historical behavior.
  PIC::Mover::SetUserDefinedParticleMover(&SEP::Mover::DispatchProductionMover);


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=SEP::Mesh::localTimeStep;
  const bool configuredNumerics = SEP::Initialization::HasActive() &&
      SEP::Initialization::Active().schemaVersion >= 2;
  // Freeze the canonical signed Br(1 AU) once for domain-cell initialization.
  // Fetching a full configuration summary inside the cell loop would copy its
  // manifest repeatedly and obscure that every cell uses one resolved model.
  double configuredRadialFieldAtOneAuT = 5.0e-9;
  if (configuredNumerics) {
    configuredRadialFieldAtOneAuT =
        SEP::SW1DAdapter::GetConfigurationSummary()
            .parker_radial_field_at_one_au_t;
  }
  if (!configuredNumerics) {
    PIC::ParticleWeightTimeStep::initTimeStep();
  }

  //create the list of mesh nodes where the injection boundary conditinos are applied
  if (_DOMAIN_GEOMETRY_==_DOMAIN_GEOMETRY_BOX_&&
      !SEP::ReducedShock::Enabled()) {
    PIC::BC::BlockInjectionBCindicatior=SEP::BoundingBoxInjection::InjectionIndicator;
    PIC::BC::userDefinedBoundingBlockInjectionFunction=SEP::BoundingBoxInjection::InjectionProcessor;
    PIC::BC::InitBoundingBoxInjectionBlockList();
  }

  //set up the particle weight
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=SEP::ParticleSource::OuterBoundary::BoundingBoxInjectionRate;
  if (configuredNumerics) {
    const SEP::Initialization::Configuration& initialization =
        SEP::Initialization::Active();
    // SpeciesList is resolved by AMPS at build time.  Initialize every slot in
    // that generated table rather than assuming that H_PLUS exists or occupies
    // a particular macro index.  The schema-v2 [species] value is explicitly a
    // common base numerical weight, so applying it to all compiled entries is
    // input semantics—not an inferred physical composition.
    for (int species = 0; species < PIC::nTotalSpecies; ++species) {
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species] =
          initialization.timeStepS;
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] =
          initialization.particleWeight;
    }
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = true;
    for (int species = 0; species < PIC::nTotalSpecies; ++species)
      InstallConfiguredParticleNumerics(
          PIC::Mesh::mesh->rootTree, species, initialization.timeStepS,
          initialization.particleWeight);
  }
  else {
    // Preserve the legacy AMPS-derived numerical policy, but apply it to the
    // complete compiled table just as the schema-v2 path does.
    for (int species = 0; species < PIC::nTotalSpecies; ++species)
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(species);
  }

  // Do not overwrite AMPS' configured base particle weight with a second,
  // hard-coded source model.  Field-line injection now computes the physical
  // source from the common area, swept volume, density, and configured
  // InjectionEfficiency, then represents it through the individual statistical
  // weight correction.  This keeps all background providers normalized alike.

  //init magnetic filed
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> InitMagneticField;

  InitMagneticField =[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) -> void {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      PIC::Mesh::cDataBlockAMR *block;

      if ((block=startNode->block)!=NULL) {
        double B[3],x[3];
        int idim,i,j,k,LocalCellNumber;
        PIC::Mesh::cDataCenterNode *cell;
        double *data;

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);

          if ((cell=block->GetCenterNode(LocalCellNumber))!=NULL) {
            cell->GetX(x);
            if (SEP::Initialization::HasActive()) {
              const SEP::Initialization::Configuration& initialization =
                  SEP::Initialization::Active();
              const double origin[3] = {initialization.parkerOriginM.x,
                                        initialization.parkerOriginM.y,
                                        initialization.parkerOriginM.z};
              SEP::ParkerSpiral::GetB(
                  B, x, origin, initialization.innerRadiusM,
                  initialization.solarWindSpeedMPerS,
                  initialization.solarRotationRateRadPerS,
                  configuredRadialFieldAtOneAuT);
            }
            else {
              SEP::ParkerSpiral::GetB(B,x);
            }

            if (cell->Measure==0.0) {
              PIC::Mesh::mesh->InitCellMeasureBlock(startNode);

              if (cell->Measure==0.0) {
                PIC::Mesh::mesh->CenterNodes.deleteElement(cell);
                startNode->block->SetCenterNode(NULL,LocalCellNumber);
                continue;
              }
            }

            data=(double*)(cell->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset);

            for (idim=0;idim<3;idim++) data[idim]=B[idim];
          }
        }
      }
    }
    else {
      int iDownNode;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

      for (iDownNode=0;iDownNode<(1<<DIM);iDownNode++) if ((downNode=startNode->downNode[iDownNode])!=NULL) {
        InitMagneticField(downNode);
      }
    }
  };

  if (_PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_) {
    InitMagneticField(PIC::Mesh::mesh->rootTree);
  }

  // Field-line data, unlike Cartesian AMR display cells, are the native
  // background consumed by srcSEP movers.  Install generation one before the
  // initialization data product is written so a preview cannot claim reduced
  // coupling while retaining the legacy Parker values on its actual line.
  install_reduced_background_on_field_lines(
      SEP::Background::SimulationTimeSeconds());

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh->ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  if (PIC::ParticleBuffer::ParticleDataBuffer==NULL) PIC::ParticleBuffer::Init(10000000);

  if (configuredNumerics &&
      SEP::Initialization::Active().schemaVersion >= 3) {
    // The native writer now sees the completed center-node Parker field, an
    // initialized particle buffer, and the already installed global/block
    // time step and statistical weight. Unlike outputMeshTECPLOT(), this
    // product contains initialized AMPS data.
    for (int species=0;species<PIC::nTotalSpecies;++species) {
      const std::string path=InitializationDataPath(
          SEP::Initialization::Active().dataTecplotFile,species,
          PIC::nTotalSpecies);
      PIC::Mesh::mesh->outputMeshDataTECPLOT(path.c_str(),species);
    }
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
//  double TimeCounter=time(NULL);
//  int LastDataOutputFileNumber=-1;


/*  //init the sampling of the particls' distribution functions: THE DECLARATION IS MOVED INTO THE INPUT FILE
  const int nSamplePoints=3;
  double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

  PIC::DistributionFunctionSample::vMin=-40.0E3;
  PIC::DistributionFunctionSample::vMax=40.0E3;
  PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

  PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);*/
}

// Stage and publish one shared-provider realization into native srcSEP vertex
// storage.  The staging vector is not an optimization: it is the transaction
// boundary that prevents a failed query near the end of a line from leaving a
// mixture of old and new epochs in the mover-visible AMPS arrays.
void install_reduced_background_on_field_lines(double epochS) {
  if(!SEP::ReducedShock::Enabled())return;
  std::string error;
  if(!SEP::ReducedShock::Prepare(epochS,&error))
    exit(__LINE__,__FILE__,error.c_str());

  const SEP::ReducedShock::EpochMetadata& metadata=
      SEP::ReducedShock::Metadata();
  const std::shared_ptr<const SEP::Background::BackgroundSnapshot> current=
      SEP::Background::SnapshotStore::Instance().Current();
  if(current&&current->provider()==SEP::Background::Provider::ReducedShock&&
      current->epoch_seconds()==metadata.epochS&&
      current->field_line_generation()==metadata.generation)return;

  struct PendingVertex {
    PIC::FieldLine::cFieldLineVertex* vertex;
    SEP::ReducedShock::AmbientSample sample;
  };
  std::vector<PendingVertex> pending;
  for(int line=0;line<PIC::FieldLine::nFieldLine;++line) {
    PIC::FieldLine::cFieldLine& fieldLine=PIC::FieldLine::FieldLinesAll[line];
    const int vertexCount=fieldLine.GetTotalSegmentNumber()+1;
    for(int index=0;index<vertexCount;++index) {
      PIC::FieldLine::cFieldLineVertex* vertex=fieldLine.GetVertex(index);
      if(vertex==NULL)exit(__LINE__,__FILE__,
          "reduced background encountered a missing field-line vertex");
      const double* x=vertex->GetX();
      PendingVertex staged;
      staged.vertex=vertex;
      if(!SEP::ReducedShock::EvaluateAmbient({{x[0],x[1],x[2]}},
          &staged.sample,&error))exit(__LINE__,__FILE__,error.c_str());
      pending.push_back(staged);
    }
  }
  if(pending.empty())exit(__LINE__,__FILE__,
      "reduced background found no native field-line vertices");

  // Publication is forbidden while a mover read phase is active.  Check that
  // contract before replacing native arrays, then preserve the outgoing state
  // in AMPS' previous-epoch datums.  At generation one the new realization is
  // copied to both slots: there is no physical epoch before event start from
  // which a derivative could legitimately be inferred.
  SEP::Background::SnapshotStore::Instance().AssertProviderMayWrite(
      SEP::Background::Provider::ReducedShock);
  const bool firstGeneration=!current;
  for(PendingVertex& item:pending) {
    PIC::FieldLine::cFieldLineVertex* vertex=item.vertex;
    double oldB[3],oldU[3],oldN=0.0,oldT=0.0,oldP=0.0;
    vertex->GetMagneticField(oldB);
    vertex->GetPlasmaVelocity(oldU);
    vertex->GetPlasmaDensity(oldN);
    vertex->GetPlasmaTemperature(oldT);
    vertex->GetPlasmaPressure(oldP);
    double newB[3]={item.sample.magneticFieldT[0],
                    item.sample.magneticFieldT[1],
                    item.sample.magneticFieldT[2]};
    double newU[3]={item.sample.velocityMPerS[0],
                    item.sample.velocityMPerS[1],
                    item.sample.velocityMPerS[2]};
    const double* previousB=firstGeneration?newB:oldB;
    const double* previousU=firstGeneration?newU:oldU;
    vertex->SetDatum(PIC::FieldLine::DatumAtVertexPrevious::
        DatumAtVertexMagneticField,const_cast<double*>(previousB));
    vertex->SetDatum(PIC::FieldLine::DatumAtVertexPrevious::
        DatumAtVertexPlasmaVelocity,const_cast<double*>(previousU));
    vertex->SetDatum(PIC::FieldLine::DatumAtVertexPrevious::
        DatumAtVertexPlasmaDensity,
        firstGeneration?item.sample.numberDensityM3:oldN);
    vertex->SetDatum(PIC::FieldLine::DatumAtVertexPrevious::
        DatumAtVertexPlasmaTemperature,
        firstGeneration?item.sample.protonTemperatureK:oldT);
    vertex->SetDatum(PIC::FieldLine::DatumAtVertexPrevious::
        DatumAtVertexPlasmaPressure,
        firstGeneration?item.sample.pressurePa:oldP);
    vertex->SetMagneticField(newB);
    vertex->SetPlasmaVelocity(newU);
    vertex->SetPlasmaDensity(item.sample.numberDensityM3);
    vertex->SetPlasmaTemperature(item.sample.protonTemperatureK);
    vertex->SetPlasmaPressure(item.sample.pressurePa);

    // Read the application-owned storage back through the same accessors used
    // by movers.  Comparing only the provider's staging buffer would prove MPI
    // determinism but not the coupling itself: an incorrect datum offset or a
    // setter wired to legacy storage could still go unnoticed.  Setters and
    // getters operate on the same double-valued native record, so this is an
    // exact representation check rather than a floating-point approximation.
    double installedB[3],installedU[3];
    double installedN=0.0,installedT=0.0,installedP=0.0;
    vertex->GetMagneticField(installedB);
    vertex->GetPlasmaVelocity(installedU);
    vertex->GetPlasmaDensity(installedN);
    vertex->GetPlasmaTemperature(installedT);
    vertex->GetPlasmaPressure(installedP);
    bool matches=installedN==item.sample.numberDensityM3&&
        installedT==item.sample.protonTemperatureK&&
        installedP==item.sample.pressurePa;
    for(int component=0;component<3;++component)
      matches=matches&&
          installedB[component]==item.sample.magneticFieldT[component]&&
          installedU[component]==item.sample.velocityMPerS[component];
    if(!matches)exit(__LINE__,__FILE__,
        "reduced background native field-line readback differs from provider");
  }

  SEP::Background::PublishModelOwnedSnapshot(
      SEP::Background::Provider::ReducedShock,metadata.epochS,
      metadata.validUntilS,
      "shared corona/SWCME reduced front and ambient on native srcSEP vertices",
      metadata.generation);

  // Every rank holds the line geometry and provider inputs.  Reread the native
  // endpoint (not the staging object) and require the complete installed state
  // to agree before a mover may consume it.  Min/max checks are stronger than
  // comparing rank-zero log text and catch a non-deterministic asset,
  // coordinate interpretation, or rank-specific storage failure.
  const PendingVertex& endpoint=pending.back();
  double endpointB[3],endpointU[3];
  double endpointN=0.0,endpointT=0.0,endpointP=0.0;
  endpoint.vertex->GetMagneticField(endpointB);
  endpoint.vertex->GetPlasmaVelocity(endpointU);
  endpoint.vertex->GetPlasmaDensity(endpointN);
  endpoint.vertex->GetPlasmaTemperature(endpointT);
  endpoint.vertex->GetPlasmaPressure(endpointP);
  const double localCheck[9]={endpointB[0],endpointB[1],endpointB[2],
      endpointU[0],endpointU[1],endpointU[2],endpointN,endpointT,endpointP};
  double minima[9],maxima[9];
  MPI_Allreduce(localCheck,minima,9,MPI_DOUBLE,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localCheck,maxima,9,MPI_DOUBLE,MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  for(int i=0;i<9;++i)if(minima[i]!=maxima[i])exit(__LINE__,__FILE__,
      "reduced native field-line state differs between MPI ranks");
  if(PIC::ThisThread==0) {
    const SEP::ReducedShock::FrontSummary front=
        SEP::ReducedShock::CurrentFrontSummary();
    std::cout<<"REDUCED_BACKGROUND epoch_s="<<metadata.epochS
             <<" generation="<<metadata.generation
             <<" vertices="<<pending.size()
             <<" apex_radius_m="<<front.apexRadiusM
             <<" accepted_area_m2="<<front.acceptedShockAreaM2
             <<" apex_shock_accepted="<<(front.apexShockAccepted?1:0)
             <<" endpoint_b_t="<<endpointB[0]<<','<<endpointB[1]<<','
             <<endpointB[2]
             <<" endpoint_u_m_s="<<endpointU[0]<<','<<endpointU[1]<<','
             <<endpointU[2]
             <<" endpoint_n_m3="<<endpointN
             <<" endpoint_t_k="<<endpointT
             <<" endpoint_p_pa="<<endpointP
             <<" event="<<metadata.eventIdentity<<std::endl;
  }
}







  //time step
// for (long int niter=0;niter<100000001;niter++) {
void amps_time_step(){

    //make the time advance
    static int LastDataOutputFileNumber=0;

    //perform test after the first coupling with the SWMF
    static bool TestCompleted=false;

    #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    if ((TestCompleted==false)&&(AMPS2SWMF::MagneticFieldLineUpdate::FirstCouplingFlag==true)) {
      TestCompleted=true;
      TestManager();
    }


    //prepopulate the magnetic filed lines with particles if needed
    if ((AMPS2SWMF::MagneticFieldLineUpdate::FirstCouplingFlag==true)&&(AMPS2SWMF::FieldLineData::ParticlePrepopulateFlag==true)) {
      AMPS2SWMF::FieldLineData::ParticlePrepopulateFlag=false;
      SEP::ParticleSource::PopulateAllFieldLines();
    }
    #endif


    static bool init_divVsw=false;

    if ((init_divVsw==false)&&(SEP::AccountAdiabaticCoolingFlag==true)) {
       init_divVsw=true;

       if (_PIC_COUPLER_MODE_!=_PIC_COUPLER_MODE__SWMF_) {
          PIC::DomainBlockDecomposition::UpdateBlockTable();
          SEP::SolarWind::SetDivSolarWindVelocity();
       }
      else {
          PIC::DomainBlockDecomposition::UpdateBlockTable();
          SEP::SolarWind::SetDivSolarWindVelocity();
      }
    }

start:

    // The application-level timestep is the sole owner of turbulence
    // evolution.  Capture the previously published shock position before the
    // particle phase and the new position after it, then pass that interval to
    // the same transaction that reduces particle-wave exchange.  A static
    // history is intentional here: amps_time_step() is the common standalone
    // and coupled entry point, so neither caller may perform a second Advance.
    // On the first call the before/after radii are identical and no artificial
    // shock source is deposited.
    static ShockRadiusHistory shockHistory;
    if (!shockHistory.initialized) {
      shockHistory.previousRadiusM = CurrentShockRadiusM();
      shockHistory.initialized = true;
    }

    // Bind the whole AMPS particle phase to one immutable background snapshot.
    // PrepareSnapshotForParticleStep() observes a new SWMF coupling epoch only
    // between steps and publishes it as a new read-only generation.  The RAII
    // guard then prevents analytic, SWCME, SWMF, or local-evolution writers from
    // replacing the authoritative state until every mover in PIC::TimeStep()
    // has returned.  Each mover independently acquires a const handle in the
    // common SEP::ParticleMover wrapper.
    try {
      SEP::Background::PrepareSnapshotForParticleStep();
      SEP::Background::ParticleReadPhase background_read =
          SEP::Background::SnapshotStore::Instance().BeginParticleRead(
              SEP::Background::SimulationTimeSeconds());
      PIC::TimeStep();
    }
    catch (const std::exception& exception) {
      // A missing/stale snapshot is a model-state error, not a condition under
      // which particles may safely continue with whichever mutable arrays happen
      // to be present.  Route the detailed reason through the existing AMPS
      // fatal-error path so all MPI ranks stop instead of diverging.
      std::cerr << "ERROR: cannot enter SEP particle step: "
                << exception.what() << std::endl;
      exit(__LINE__, __FILE__,
           "Error: invalid SEP background snapshot state");
    }

    // Coupled-library and standalone runs use this one post-particle
    // deterministic reduction.  The current radius is sampled only after the
    // particle read phase has ended, so provider publication and turbulence
    // mutation cannot overlap an immutable mover snapshot.
    const double currentShockRadiusM = CurrentShockRadiusM();
    if(!SEP::ReducedShock::Enabled()) {
      const SEP::Transport::Status turbulenceStatus =
          SEP::Turbulence::PICAdapter::Advance(
              PIC::ParticleWeightTimeStep::GlobalTimeStep[0],
              shockHistory.previousRadiusM, currentShockRadiusM);
      if (!turbulenceStatus.ok())
        exit(__LINE__, __FILE__, turbulenceStatus.message.c_str());
    }
    shockHistory.previousRadiusM = currentShockRadiusM;

    // Zero particles is a physical mode invariant, not merely an input label.
    // Check the actual allocated AMPS population after every native step and
    // reduce it across ranks so an injection on a non-root owner cannot hide.
    if(SEP::ReducedShock::Enabled()) {
      const long int localParticles=PIC::ParticleBuffer::GetAllPartNum();
      long int globalParticles=0;
      MPI_Allreduce(&localParticles,&globalParticles,1,MPI_LONG,MPI_SUM,
                    MPI_GLOBAL_COMMUNICATOR);
      if(globalParticles!=0)exit(__LINE__,__FILE__,
          "reduced background-only run created particles");
      if(PIC::ThisThread==0)std::cout<<"REDUCED_PARTICLES count=0 epoch_s="
          <<SEP::Background::SimulationTimeSeconds()<<std::endl;
    }

//    PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(50,100); //(SEP::MinParticleLimit,SEP::MaxParticleLimit);


     const SEP::Run::Configuration& run=SEP::Run::Active().get();
     PIC::ParticleSplitting::FledLine::WeightedParticleMerging(
         run.populationControl.spatialBins,
         run.populationControl.momentumBins,
         run.populationControl.pitchBins,
         run.populationControl.minimumParticlesPerCell,
         run.populationControl.maximumParticlesPerCell);
     PIC::ParticleSplitting::FledLine::WeightedParticleSplitting(
         run.populationControl.spatialBins,
         run.populationControl.momentumBins,
         run.populationControl.pitchBins,
         run.populationControl.minimumParticlesPerCell,
         run.populationControl.maximumParticlesPerCell);


     // write output file
     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
//       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }

  //check the simulation time
  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if ((SEP::FreezeSolarWindModelTime>0.0)&&(SEP::FreezeSolarWindModelTime<AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime)) {
    goto start;
  }
  #endif
}
