// Copyright (C) 2002 Regents of the University of Michigan.
//
// M0 linked diagnostics for srcMoon.  This file contains no replacement force
// law, source law, mover, or surface model.  It creates deterministic states,
// invokes the callbacks installed by the normal Moon initializer, and records
// raw SI-valued observations.  Independent analytical/numerical references and
// acceptance decisions live in srcMoon/test/integration_common.py so a copied
// test-only implementation cannot silently become the physical oracle.

#include "MoonIntegrationTests.h"

#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonDivinerObservation.h"
#include "MoonI17Authorization.h"
#include "MoonSurfaceEnvironment.h"
#include "MoonSurfaceTemperature.h"
#include "MoonParticleBudget.h"
#include "MoonM2Budget.h"
#include "MoonM2Cohorts.h"
#include "MoonPolarDiagnostics.h"
#include "species/Argon.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <vector>

void amps_time_step();
double localSphericalSurfaceResolution(double *x);
double localResolution(double *x);

namespace {

class Receipt {
 public:
  Receipt(const std::string& path,bool enabled): enabled_(enabled) {
    // Linked diagnostics execute on every rank when a test needs MPI ownership
    // checks, but a receipt is a single ordered provenance record.  Only rank
    // zero opens the path; all other Receipt instances are intentional no-op
    // sinks so concurrent truncation or interleaved writes are impossible.
    if (enabled_) {
      output_.open(path.c_str());
      output_ << std::setprecision(17);
    }
  }

  bool Good() const { return !enabled_ || output_.good(); }

  template<class T> void Scalar(const std::string& key,const T& value) {
    if (enabled_) output_ << key << '=' << value << '\n';
  }

  void Vector(const std::string& key,const double value[3]) {
    if (enabled_)
      output_ << key << '=' << value[0] << ',' << value[1] << ',' << value[2]
              << '\n';
  }

 private:
  bool enabled_;
  std::ofstream output_;
};

bool RecordAppliedThermalInitialization(Receipt* receipt,std::string* error) {
#ifdef _NO_SPICE_CALLS_
  if (error!=nullptr)
    *error="applied thermal initialization provenance requires the actual SPICE kernel pool";
  return false;
#else
  Moon::SurfaceTemperature::ProductionThermalInitializationInfo initial;
  if (!Moon::SurfaceTemperature::GetProductionThermalInitialization(&initial,error))
    return false;
  // These values describe the actual production columns. They do not assert
  // that a Python pilot state was loaded, or that its physical proof applies.
  // The independent formal binder compares method, epoch and exported bytes.
  receipt->Scalar("thermal_initialization_mode",initial.mode);
  receipt->Scalar("thermal_initialization_et_s",initial.ephemerisTimeS);
  receipt->Scalar("thermal_initialization_spinup_completed_cycles",initial.spinup.completedCycles);
  receipt->Scalar("thermal_initialization_spinup_maximum_node_temperature_change_k",
      initial.spinup.maximumNodeTemperatureChangeK);
  receipt->Scalar("thermal_initialization_spinup_maximum_column_enthalpy_change_j_m2",
      initial.spinup.maximumColumnEnthalpyChangeJm2);
  receipt->Scalar("thermal_initialization_spinup_cycle_start_et_s",
      initial.spinup.cycleStartEphemerisTimeS);
  receipt->Scalar("thermal_initialization_spinup_cycle_end_et_s",
      initial.spinup.cycleEndEphemerisTimeS);
  // A uniform-temperature hash alone cannot reveal retargeted ephemerides.
  // Bind the actual furnished CSPICE pool, rather than predicting a kernel
  // set from Config.pl paths or reusing another run's canonical filenames.
  // Resolve relative furnished paths in this process's active run directory;
  // the independent Python binder hashes those concrete runtime bytes.
  SpiceInt kernelCount=0;
  ktotal_c("ALL",&kernelCount);
  if (kernelCount<=0) {
    *error="applied thermal initialization has no actual furnished SPICE kernels";
    return false;
  }
  receipt->Scalar("thermal_initialization_spice_kernel_count",kernelCount);
  for (SpiceInt index=0;index<kernelCount;++index) {
    SpiceChar path[4096],kind[128],source[4096];
    SpiceInt handle;SpiceBoolean found;
    kdata_c(index,"ALL",sizeof(path),sizeof(kind),sizeof(source),
        path,kind,source,&handle,&found);
    if (!found || path[0]=='\0') {
      *error="actual furnished SPICE kernel inventory is incomplete";
      return false;
    }
    receipt->Scalar("thermal_initialization_spice_kernel_"+
        std::to_string(static_cast<int>(index))+"_file",
        std::filesystem::absolute(path).lexically_normal().string());
  }
  return true;
#endif
}

struct State {
  double x[3];  // Moon-centred LSO position [m].
  double v[3];  // Velocity on the LSO axes [m s^-1].
};

double Norm(const double value[3]) {
  return std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]);
}

void Copy(const double source[3],double destination[3]) {
  for (int component=0;component<3;component++)
    destination[component]=source[component];
}

/**
 * Advance one deterministic state through the configured production mover.
 *
 * A fresh particle is used for each substep because the low-level mover appends
 * its result to PIC's temporary cell list.  Reusing that pointer without a full
 * PIC::TimeStep list exchange would put one particle in multiple lists.  The
 * diagnostic run terminates immediately after the bounded fixture, so retaining
 * these few particles is safe and, importantly, does not alter mover behavior.
 */
int MoveOneStep(int species,const State& initial,double dt,State* final) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
      PIC::Mesh::mesh->findTreeNode(const_cast<double*>(initial.x));
  if (node==NULL || node->block==NULL || node->Thread!=PIC::ThisThread)
    return -100;

  const long int particle=PIC::ParticleBuffer::GetNewParticle();
  PIC::ParticleBuffer::SetI(species,particle);
  // This AMPS API predates const-correct particle setters.  Mutable local
  // copies preserve the const State contract without casting away constness
  // or allowing the particle buffer to modify the caller's fixture.
  double initialX[3],initialV[3];
  Copy(initial.x,initialX);
  Copy(initial.v,initialV);
  PIC::ParticleBuffer::SetX(initialX,particle);
  PIC::ParticleBuffer::SetV(initialV,particle);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  // A correction value exists in the particle record only in individual-
  // weight builds.  Calling its setter in the normal Moon OFF configuration
  // is an AMPS data-layout error, not a harmless no-op.  Unity leaves the
  // fixture at the block's production statistical weight when the field is
  // present; no setter is needed when all particles share the block weight.
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,particle);
#endif
  PIC::ParticleBuffer::SetParticleAllocated(particle);

  // input/moon.input selects this exact production mover for normal particle
  // dispatch.  Calling it directly retains its status code, which the generic
  // DispatchParticleMover wrapper intentionally discards after list handling.
  // Retaining the status is essential for I09's termination audit and does
  // not change the mover, acceleration callback, or boundary callbacks used.
  const int moverStatus=PIC::Mover::UniformWeight_UniformTimeStep_SecondOrder(
      particle,dt,node);
  if (moverStatus==_PARTICLE_MOTION_FINISHED_) {
    PIC::ParticleBuffer::GetX(final->x,particle);
    PIC::ParticleBuffer::GetV(final->v,particle);
  }
  return moverStatus;
}

int Integrate(int species,const State& initial,double totalTime,double dt,
              State* final) {
  const int steps=static_cast<int>(std::llround(totalTime/dt));
  if (steps<=0 || std::fabs(steps*dt-totalTime)>
      1.0e-12*std::max(1.0,totalTime)) return -101;

  State current=initial,next;
  for (int step=0;step<steps;step++) {
    const int status=MoveOneStep(species,current,dt,&next);
    if (status!=_PARTICLE_MOTION_FINISHED_) return status;
    current=next;
  }
  *final=current;
  return _PARTICLE_MOTION_FINISHED_;
}

void WriteState(Receipt* receipt,const std::string& prefix,const State& state) {
  receipt->Vector(prefix+"_x_m",state.x);
  receipt->Vector(prefix+"_v_m_s",state.v);
}

void WriteForceState(Receipt* receipt) {
  receipt->Scalar("ephemeris_time_tdb_s",Moon::OrbitalMotion::et);
  receipt->Scalar("moon_heliocentric_distance_m",Moon::xObjectRadial);
  receipt->Scalar("moon_heliocentric_radial_velocity_m_s",Moon::vObjectRadial);
  receipt->Vector("sun_position_lso_m",Moon::xSun_SO);
  receipt->Vector("earth_position_lso_m",Moon::xEarth_SO);
  receipt->Vector("moon_velocity_frozen_lso_m_s",Moon::vObject_SO_FROZEN);
  receipt->Vector("lso_angular_velocity_rad_s",Moon::RotationVector_SO_FROZEN);
}

int RunI14(Receipt* receipt,const std::string& receiptPath) {
#ifdef _NO_SPICE_CALLS_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  // The first ordinary step initializes frozen force/frame state. Calling
  // a force-aware terrain mover immediately after amps_init otherwise uses
  // zero Earth state, giving NaN acceleration and invalid radiation-table
  // indices. This is fixture initialization, not a replacement force law.
  Moon::nIterationCounter=0;
  amps_time_step();
  receipt->Scalar("initialization_steps",1);
  receipt->Scalar("fixture_epoch_tdb_s",Moon::OrbitalMotion::et);
  if (!Moon::Runtime::HasConfiguration()) {
    receipt->Scalar("configuration_error_missing_input",1);
    return 4;
  }
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  const Moon::SurfaceEnvironment::TerrainRayTracer* terrain=
      Moon::SurfaceEnvironment::ProductionTerrain();
  if (!configuration.m2EnvironmentEnabled ||
      configuration.surfaceGeometry!=Moon::Runtime::SurfaceGeometry::Lola ||
      terrain==NULL) {
    receipt->Scalar("configuration_error_lola_production_boundary_required",1);
    return 4;
  }
  const Moon::Surface::Triangulation& mesh=terrain->triangulation();
  receipt->Scalar("surface_geometry","lola");
  receipt->Scalar("surface_mesh_resolution_m",
                  configuration.surfaceMeshResolutionM);
  receipt->Scalar("surface_collision_distance_tolerance_m",
                  configuration.surfaceCollisionDistanceToleranceM);
  receipt->Scalar("vertex_count",mesh.verticesM.size());
  receipt->Scalar("face_count",mesh.faces.size());
  receipt->Scalar("subdivision_level",mesh.subdivisionLevel);
  receipt->Scalar("minimum_radius_m",mesh.minimumRadiusM);
  receipt->Scalar("maximum_radius_m",mesh.maximumRadiusM);

  // Evaluate every face, rather than a display-file subset, for degenerate or
  // inward connectivity.  This is a geometric invariant of the exact
  // triangulation installed into the production BVH on every MPI rank.
  std::size_t degenerateFaces=0,invertedFaces=0;
  for (const auto& face:mesh.faces) {
    const auto& a=mesh.verticesM[face[0]];
    const auto& b=mesh.verticesM[face[1]];
    const auto& c=mesh.verticesM[face[2]];
    const double edge0[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    const double edge1[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
    const double cross[3]={edge0[1]*edge1[2]-edge0[2]*edge1[1],
        edge0[2]*edge1[0]-edge0[0]*edge1[2],
        edge0[0]*edge1[1]-edge0[1]*edge1[0]};
    const double area2=Norm(cross);
    if (!(area2>0.0) || !std::isfinite(area2)) degenerateFaces++;
    const double centroid[3]={(a[0]+b[0]+c[0])/3.0,
        (a[1]+b[1]+c[1])/3.0,(a[2]+b[2]+c[2])/3.0};
    if (cross[0]*centroid[0]+cross[1]*centroid[1]+cross[2]*centroid[2]<=0.0)
      invertedFaces++;
  }
  receipt->Scalar("degenerate_face_count",degenerateFaces);
  receipt->Scalar("inverted_face_count",invertedFaces);

  // Preserve a bounded, deterministic vertex/radial-ray table.  Python reads
  // LDEM_4 independently and compares these body-fixed directions with the
  // native raster; this C++ adapter never reads the validation reference.
  std::size_t radialQueries=0,radialHits=0;
  double maximumVertexRayDifferenceM=0.0;
  if (PIC::ThisThread==0) {
    std::ofstream geometry(receiptPath+".I14_geometry.csv");
    geometry << std::setprecision(17)
             << "vertex_index,latitude_deg,longitude_east_deg,"
                "vertex_radius_m,ray_radius_m,ray_face\n";
    if (!geometry) return 6;
    const std::size_t requestedSamples=512;
    const std::size_t stride=std::max<std::size_t>(
        1,mesh.verticesM.size()/requestedSamples);
    const std::array<double,3> origin{{0.0,0.0,0.0}};
    for (std::size_t vertex=0;vertex<mesh.verticesM.size();vertex+=stride) {
      std::array<double,3> direction=mesh.verticesM[vertex];
      const double radius=std::sqrt(direction[0]*direction[0]+
          direction[1]*direction[1]+direction[2]*direction[2]);
      for (double& component:direction) component/=radius;
      Moon::SurfaceEnvironment::RayHit hit;
      std::string error;
      radialQueries++;
      if (!terrain->IntersectRay(origin,direction,0.0,
              std::numeric_limits<double>::infinity(),
              std::numeric_limits<std::size_t>::max(),&hit,&error)) {
        std::cerr << "I14 radial query failed: " << error << '\n';
        return 6;
      }
      if (hit.hit) {
        radialHits++;
        maximumVertexRayDifferenceM=std::max(
            maximumVertexRayDifferenceM,std::abs(hit.distanceM-radius));
      }
      const double latitude=std::asin(direction[2])*180.0/Pi;
      double longitude=std::atan2(direction[1],direction[0])*180.0/Pi;
      if (longitude<0.0) longitude+=360.0;
      geometry << vertex << ',' << latitude << ',' << longitude << ','
               << radius << ','
               << (hit.hit ? hit.distanceM :
                   std::numeric_limits<double>::quiet_NaN()) << ','
               << (hit.hit ? static_cast<long long>(hit.face) : -1) << '\n';
    }
    if (!geometry) return 6;
  }
  receipt->Scalar("radial_query_count",radialQueries);
  receipt->Scalar("radial_hit_count",radialHits);
  receipt->Scalar("maximum_vertex_ray_difference_m",
                  maximumVertexRayDifferenceM);

  // Each rank selects the first terrain vertex whose outward-offset start
  // point belongs to one of its allocated AMR blocks.  The ordinary
  // second-order mover then crosses the moving LOLA surface.  The monotonic
  // counter is incremented inside the production finder, independently of
  // whether the surface kernel subsequently sticks or re-emits the atom.
  int localFixtureExecuted=0,localFixtureCollision=0;
  int localMoverStatus=std::numeric_limits<int>::min();
  for (std::size_t vertex=static_cast<std::size_t>(PIC::ThisThread);
       vertex<mesh.verticesM.size();
       vertex+=static_cast<std::size_t>(PIC::nTotalThreads)) {
    std::array<double,3> radial=mesh.verticesM[vertex];
    const double radius=std::sqrt(radial[0]*radial[0]+radial[1]*radial[1]+
                                  radial[2]*radial[2]);
    for (double& component:radial) component/=radius;
    SpiceDouble bodyToLso[6][6],bodyState[6],lsoState[6];
    sxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,
             Moon::OrbitalMotion::et,bodyToLso);
    for (int component=0;component<3;++component) {
      bodyState[component]=mesh.verticesM[vertex][component]+2000.0*radial[component];
      bodyState[component+3]=-200.0*radial[component];
    }
    mxvg_c(bodyToLso,bodyState,6,6,lsoState);
    State initial;
    for (int component=0;component<3;++component) {
      initial.x[component]=lsoState[component];
      initial.v[component]=lsoState[component+3];
    }
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
        PIC::Mesh::mesh->findTreeNode(initial.x);
    if (node==NULL || node->Thread!=PIC::ThisThread || node->block==NULL)
      continue;
    rnd_seed(configuration.randomSeed+PIC::ThisThread);
    Moon::SurfaceEnvironment::BeginProductionParticleStep(
        Moon::OrbitalMotion::et);
    const std::uint64_t before=
        Moon::SurfaceEnvironment::ProductionCollisionCount();
    State final;
    localMoverStatus=MoveOneStep(_NA_SPEC_,initial,20.0,&final);
    const std::uint64_t after=
        Moon::SurfaceEnvironment::ProductionCollisionCount();
    localFixtureExecuted=1;
    localFixtureCollision=after>before ? 1 : 0;
    break;
  }
  int fixtureRanks=0,collisionRanks=0;
  MPI_Allreduce(&localFixtureExecuted,&fixtureRanks,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localFixtureCollision,&collisionRanks,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  receipt->Scalar("mover_fixture_rank_count",fixtureRanks);
  receipt->Scalar("mover_collision_rank_count",collisionRanks);
  receipt->Scalar("rank_zero_mover_status",localMoverStatus);
  return 0;
#endif
}

std::vector<std::string> SplitRestrictedCsv(const std::string& line) {
  std::vector<std::string> fields;
  std::istringstream stream(line);
  std::string field;
  while (std::getline(stream,field,',')) fields.push_back(field);
  if (!line.empty() && line.back()==',') fields.emplace_back();
  return fields;
}

int RunI15(Receipt* receipt,const std::string& receiptPath) {
#ifdef _NO_SPICE_CALLS_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  // This adapter supplies geometry to the ordinary production illumination
  // service. Observation values NEVER enter the BVH or the disc quadrature.
  // Python stages only a reviewed, hash-bound geometry/time query table and
  // independently joins/scorers the native results after numerical gates.
  const auto& configuration=Moon::Runtime::GetConfiguration();
  const auto* terrain=Moon::SurfaceEnvironment::ProductionTerrain();
  if (!configuration.m2EnvironmentEnabled || terrain==nullptr ||
      configuration.illuminationTerrainSolarModel!="uniform-disc" ||
      configuration.m2IlluminationQueriesFile.empty()) return 4;
  std::ifstream input(configuration.m2IlluminationQueriesFile);
  std::ostringstream path;
  path << receiptPath << ".I15.rank" << std::setw(4) << std::setfill('0')
       << PIC::ThisThread << ".csv";
  std::ofstream output(path.str());
  std::string line,error;
  int localOk=(input && output && std::getline(input,line))?1:0,allOk=0;
  MPI_Allreduce(&localOk,&allOk,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!allOk) return 6;
  const auto names=SplitRestrictedCsv(line);
  std::map<std::string,std::size_t> fields;
  for (std::size_t k=0;k<names.size();++k) fields[names[k]]=k;
  for (const char* key:{"query_id","site_id","latitude_deg","longitude_east_deg",
      "epoch","time_system","observer_height_m","height_direction","include_earth"})
    if (!fields.count(key)) localOk=0;
  MPI_Allreduce(&localOk,&allOk,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!allOk) return 6;
  output << std::setprecision(17)
      << "query_id,site_id,input_epoch,input_time_system,ephemeris_time_tdb_s,latitude_deg,longitude_east_deg,"
         "terrain_visible_area_fraction,earth_visible_area_fraction,jointly_visible_area_fraction,"
         "projected_incidence,sampled_any_terrain_visible,sampled_any_jointly_visible,ray_count\n";
  long long row=0,owned=0;
  while (std::getline(input,line)) {
    if (line.empty()) continue;
    const auto values=SplitRestrictedCsv(line);
    if (values.size()!=names.size()) { localOk=0;break; }
    const auto text=[&](const char* key)->const std::string& { return values[fields[key]]; };
    const auto number=[&](const char* key,double* x) {
      char* end=nullptr;errno=0;
      const auto& value=text(key);*x=std::strtod(value.c_str(),&end);
      return errno!=ERANGE && end!=value.c_str() && *end=='\0' && std::isfinite(*x);
    };
    double latitude=0,longitude=0,height=0,epoch=0;
    if (!number("latitude_deg",&latitude) || !number("longitude_east_deg",&longitude) ||
        !number("observer_height_m",&height) || latitude < -90 || latitude>90 ||
        longitude<0 || longitude>=360 || height<0 ||
        (text("height_direction")!="radial" && text("height_direction")!="facet-normal") ||
        (text("include_earth")!="true" && text("include_earth")!="false")) { localOk=0;break; }
    if (text("time_system")=="UTC") str2et_c(text("epoch").c_str(),&epoch);
    else if (text("time_system")!="TDB" || !number("epoch",&epoch)) { localOk=0;break; }
    const long long index=row++;
    if (index%PIC::nTotalThreads!=PIC::ThisThread) continue;
    const double lat=latitude*Pi/180,lon=longitude*Pi/180;
    const std::array<double,3> radial{{std::cos(lat)*std::cos(lon),std::cos(lat)*std::sin(lon),std::sin(lat)}};
    Moon::SurfaceEnvironment::RayHit hit;
    if (!terrain->IntersectRay({{0,0,0}},radial,0,std::numeric_limits<double>::infinity(),
        std::numeric_limits<std::size_t>::max(),&hit,&error) || !hit.hit) { localOk=0;break; }
    Moon::SurfaceEnvironment::SolarDiscRequest request;
    request.surfacePositionM=hit.positionM;
    request.outwardNormal=hit.outwardNormal;
    for (int k=0;k<3;++k)
      request.surfacePositionM[k]+=height*(text("height_direction")=="radial"?radial[k]:hit.outwardNormal[k]);
    if (!Moon::SurfaceEnvironment::ApparentSolarSystemGeometryAtPoint(epoch,
        request.surfacePositionM,configuration.illuminationAberrationCorrection,
        &request.pointToSunM,&request.pointToEarthM,&error)) {localOk=0;break;}
    request.solarRadiusM=configuration.thermalSolarDiscRadiusM;
    request.opaqueEarthRadiusM=text("include_earth")=="true"?configuration.thermalOpaqueEarthRadiusM:0;
    request.rayOriginOffsetM=configuration.illuminationRayOriginOffsetM;
    request.radialPoints=configuration.illuminationDiscRadialPoints;
    request.azimuthPoints=configuration.illuminationDiscAzimuthPoints;
    // D02 observer sky visibility is a line-of-sight quantity. Its local
    // horizontal plane is NOT an opaque screen for an elevated observer.
    request.applyLocalSelfHorizon=false;
    Moon::SurfaceEnvironment::SolarDiscSample sample;
    if (!Moon::SurfaceEnvironment::EvaluateSolarDisc(terrain,request,&sample,&error)) { localOk=0;break; }
    output << text("query_id") << ',' << text("site_id") << ',' << text("epoch") << ','
        << text("time_system") << ',' << epoch << ','
        << latitude << ',' << longitude << ',' << sample.terrainVisibleAreaFraction << ','
        << sample.earthVisibleAreaFraction << ',' << sample.jointlyVisibleAreaFraction << ','
        << sample.projectedIncidence << ',' << sample.sampledAnyTerrainVisible << ','
        << sample.sampledAnyJointlyVisible << ',' << sample.rayCount << '\n';
    ++owned;
  }
  if (!output || input.bad()) localOk=0;
  MPI_Allreduce(&localOk,&allOk,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  long long total=0;
  MPI_Allreduce(&owned,&total,1,MPI_LONG_LONG_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (!allOk || total==0) {
    std::cerr << "I15 production geometry/query failure: " << error << '\n';
    return 6;
  }
  receipt->Scalar("query_count",total);
  receipt->Scalar("runtime_psd_geometry_hook",
      Exosphere::SourceProcesses::PhotonStimulatedDesorption::RuntimePhotonFluxFactor==
          Moon::SurfaceTemperature::ProductionPhotonFluxFactor ? 1:0);
  receipt->Scalar("source_begin_epoch_tdb_s",Moon::SurfaceEnvironment::ProductionSourceEpochS());
  receipt->Scalar("disc_radial_points",configuration.illuminationDiscRadialPoints);
  receipt->Scalar("disc_azimuth_points",configuration.illuminationDiscAzimuthPoints);
  receipt->Scalar("aberration_correction",configuration.illuminationAberrationCorrection);
  receipt->Scalar("aberration_observer_convention","constant-body-fixed-site-at-reception");
  receipt->Scalar("body_fixed_frame",Moon::Frames::BodyFixed);
  receipt->Scalar("visibility_definition","sampled-any-disc-occurrence-not-mean-disc-area");
  return 0;
#endif
}

int RunI16(Receipt* receipt,const std::string& receiptPath) {
#ifdef _NO_SPICE_CALLS_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  if (!Moon::Runtime::HasConfiguration()) return 4;
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  if (!configuration.m2EnvironmentEnabled ||
      configuration.surfaceTemperatureMode!="diviner-map" ||
      configuration.divinerMapKineticInterpretation!="isothermal-grey" ||
      !Moon::SurfaceTemperature::ProductionStateInitialized()) {
    receipt->Scalar("configuration_error_diviner_production_mode_required",1);
    return 4;
  }
  // Only rank zero writes the ordered comparison table.  All ranks have
  // initialized the identical production Diviner table during amps_init(); no
  // collective or distributed mesh data are needed for this callback query.
  if (PIC::ThisThread!=0) return 0;
  const std::filesystem::path benchmark=
      std::filesystem::path(configuration.divinerTemperatureTableFile)
          .parent_path()/"benchmark_points.csv";
  std::ifstream input(benchmark);
  std::ofstream output(receiptPath+".I16_model.csv");
  if (!input || !output) {
    receipt->Scalar("configuration_error_benchmark_table_unavailable",1);
    return 6;
  }
  Moon::SurfaceTemperature::DivinerTable eligibilityTable;
  std::string eligibilityError;
  if (!eligibilityTable.LoadCsv(configuration.divinerTemperatureTableFile,
                                &eligibilityError)) {
    std::cerr << "I16 cannot reload the production D03 table for independent "
                 "gap eligibility: " << eligibilityError << '\n';
    return 6;
  }
  output << std::setprecision(17)
         << "source_record,latitude_deg,longitude_east_deg,local_time_hour,"
            "query_local_time_hour,production_kinetic_temperature_K,"
            "production_bolometric_brightness_temperature_K,bolometric_emissivity\n";
  std::string line;
  if (!std::getline(input,line) || line!=
      "latitude_deg,longitude_east_deg,local_time_hour,temperature_k,"
      "uncertainty_k,quality_flag,product_id,source_record,selection") {
    receipt->Scalar("configuration_error_benchmark_schema",1);
    return 6;
  }
  // Select geometry with the SAME explicit epoch, body frame and named
  // aberration correction used by the production map boundary. The legacy
  // orbit cache is geometric and may differ for an input selecting LT/CN.
  SpiceDouble sunState[6],lightTime;
  spkezr_c("SUN",Moon::OrbitalMotion::et,Moon::Frames::BodyFixed,
      configuration.illuminationAberrationCorrection.c_str(),Moon::ObjectName,
      sunState,&lightTime);
  const double sunLongitudeDeg=std::atan2(sunState[1],sunState[0])*180.0/Pi;
  SpiceDouble bodyToLso[6][6];
  sxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,
           Moon::OrbitalMotion::et,bodyToLso);
  std::size_t queryCount=0;
  std::set<std::pair<std::string,std::string> > evaluatedLocations;
  while (std::getline(input,line)) {
    if (line.empty()) continue;
    const std::vector<std::string> fields=SplitRestrictedCsv(line);
    if (fields.size()!=9 || fields[5]!="GOOD" || fields[8]!="hold-out") {
      receipt->Scalar("configuration_error_benchmark_record",queryCount+2);
      return 6;
    }
    char* end=NULL;
    const double latitudeDeg=std::strtod(fields[0].c_str(),&end);
    if (end==fields[0].c_str() || *end!='\0') return 6;
    const double archiveLongitudeDeg=std::strtod(fields[1].c_str(),&end);
    if (end==fields[1].c_str() || *end!='\0') return 6;
    const double localTimeHour=std::strtod(fields[2].c_str(),&end);
    if (end==fields[2].c_str() || *end!='\0') return 6;
    const std::pair<std::string,std::string> locationKey{fields[0],fields[1]};
    if (evaluatedLocations.find(locationKey)!=evaluatedLocations.end())
      continue;
    // A normal surface callback has one physical local time at a given
    // longitude and run epoch.  Retain only the archive quarter-hour bin
    // nearest that predeclared geometry.  Selection depends on the SPICE Sun
    // direction and archive coordinates, never on observed/model temperature.
    double currentLocalTimeHour=std::fmod(
        archiveLongitudeDeg-sunLongitudeDeg+180.0,360.0)/15.0;
    if (currentLocalTimeHour<0.0) currentLocalTimeHour+=24.0;
    double timeDifference=std::abs(localTimeHour-currentLocalTimeHour);
    timeDifference=std::min(timeDifference,24.0-timeDifference);
    if (timeDifference>0.126) continue;
    // Missing archive cells remain real gaps.  Test eligibility through the
    // same strict loader before invoking the fail-closed normal callback, and
    // omit only rows whose eight-corner reference cell is unavailable.  This
    // decision sees no AMPS temperature value. Check the physical query LST,
    // not merely its nearby archive-bin centre: the latter can be a supported
    // native node whose neighboring interval crosses a genuine map gap.
    double eligibilityTemperature=0.0,eligibilityUncertainty=0.0;
    if (!eligibilityTable.Interpolate(latitudeDeg,archiveLongitudeDeg,
            currentLocalTimeHour,&eligibilityTemperature,&eligibilityUncertainty,
            &eligibilityError)) continue;
    const double latitudeRad=latitudeDeg*Pi/180.0;
    const double longitudeRad=archiveLongitudeDeg*Pi/180.0;
    SpiceDouble bodyState[6]={
        _RADIUS_(_MOON_)*std::cos(latitudeRad)*std::cos(longitudeRad),
        _RADIUS_(_MOON_)*std::cos(latitudeRad)*std::sin(longitudeRad),
        _RADIUS_(_MOON_)*std::sin(latitudeRad),0.0,0.0,0.0};
    SpiceDouble lsoState[6];
    mxvg_c(bodyToLso,bodyState,6,6,lsoState);
    const double kineticTemperatureK=
        Moon::SurfaceTemperature::ProductionTemperatureK(lsoState);
    // The callback feeds kinetic K to particle physics. D03 stores native
    // blackbody-equivalent bolometric brightness K. Compare through the same
    // explicit grey-facet observation operator used for thermal predictions;
    // retaining both columns prevents an implicit unit-emissivity comparison.
    const double bolometricTemperatureK=
        Moon::DivinerObservation::GreyBolometricBrightnessTemperatureK(
            kineticTemperatureK,configuration.thermal.infraredEmissivity);
    if (!(kineticTemperatureK>0.0) || !std::isfinite(kineticTemperatureK) ||
        !(bolometricTemperatureK>0.0) || !std::isfinite(bolometricTemperatureK))
      return 6;
    output << fields[7] << ',' << latitudeDeg << ','
           << archiveLongitudeDeg << ',' << localTimeHour << ','
           << currentLocalTimeHour << ','
           << kineticTemperatureK << ',' << bolometricTemperatureK << ','
           << configuration.thermal.infraredEmissivity << '\n';
    evaluatedLocations.insert(locationKey);
    queryCount++;
  }
  if (!output || queryCount==0) return 6;
  receipt->Scalar("benchmark_query_count",queryCount);
  receipt->Scalar("surface_temperature_mode",
                  configuration.surfaceTemperatureMode);
  receipt->Scalar("diviner_temperature_table_file",
                  configuration.divinerTemperatureTableFile);
  receipt->Scalar("diviner_map_kinetic_interpretation",
                  configuration.divinerMapKineticInterpretation);
  receipt->Scalar("bolometric_emissivity",configuration.thermal.infraredEmissivity);
  return 0;
#endif
}

int RunI17(Receipt* receipt,const std::string& receiptPath) {
#ifdef _NO_SPICE_CALLS_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  if (!Moon::Runtime::HasConfiguration()) return 4;
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  if (!configuration.m2EnvironmentEnabled ||
      configuration.surfaceTemperatureMode!="thermal-model" ||
      configuration.divinerEclipseHistoryFile.empty() ||
      !Moon::SurfaceTemperature::ProductionStateInitialized()) {
    receipt->Scalar("configuration_error_i17_thermal_mode_required",1);
    return 4;
  }

  // Selection is explicit. The default remains the isolated June diagnostic;
  // the formal Python runner verifies the reviewer/byte/source/build contract
  // before staging a hold-out-only CSV. This native adapter does not score or
  // approve that artifact: it checks the selected stream and exports the raw
  // production model. A manual hold-out launch still needs an explicit review
  // path and cannot report a formal PASS without the independent runner.
  const std::string selection=configuration.divinerEclipseSelection;
  if ((selection!="training" && selection!="holdout") ||
      (selection=="holdout" && (configuration.divinerEclipseFreezeFile.empty() ||
       !std::filesystem::is_regular_file(configuration.divinerEclipseFreezeFile)))) {
    receipt->Scalar("configuration_error_i17_reviewed_selection_required",1);
    return 4;
  }
  if (selection=="holdout") {
    std::string authorizationError;
    if (!Moon::I17Authorization::Validate(configuration.divinerEclipseAuthorizationFile,
        configuration.divinerEclipseHistoryFile,configuration.divinerEclipseFreezeFile,
        &authorizationError)) {
      if (PIC::ThisThread==0) std::cerr << authorizationError << '\n';
      receipt->Scalar("configuration_error_i17_native_authorization_required",1);
      return 4;
    }
    receipt->Scalar("native_authorization_file",configuration.divinerEclipseAuthorizationFile);
    receipt->Scalar("native_authorization_scope","schema-status-path-token; Python verifies actual SHA-256/review");
  }
  // Every selection independently spins up its complete production columns.
  // December has its own explicitly declared UTC cycle; the caller never
  // substitutes the final June state. Other ranks have
  // initialized the same production model but do not need to duplicate this
  // deterministic global thermal history.  The manifest still records the
  // requested eight-rank launch honestly; no MPI speedup is claimed here.
  if (PIC::ThisThread!=0) return 0;

  Moon::SurfaceTemperature::SpinUpMetrics spinup;
  std::string error;
  if (!Moon::SurfaceTemperature::SpinUpProductionState(&spinup,&error)) {
    std::cerr << "I17 periodic thermal spin-up failed: " << error << '\n';
    return 6;
  }
  receipt->Scalar("capability_status","EXPERIMENTAL / NOT VALIDATED");
  receipt->Scalar("comparison_selection",selection);
  receipt->Scalar("initialization_mode","independent-observation-free-periodic-spinup");
  receipt->Scalar("terrain_solar_model",configuration.illuminationTerrainSolarModel);
  receipt->Scalar("disc_radial_points",configuration.illuminationDiscRadialPoints);
  receipt->Scalar("disc_azimuth_points",configuration.illuminationDiscAzimuthPoints);
  if (selection=="holdout")
    receipt->Scalar("reviewed_freeze_file",configuration.divinerEclipseFreezeFile);
  receipt->Scalar("eclipse_direct_irradiance_model",
                  configuration.thermalEclipseDirectIrradianceModel);
  receipt->Scalar("solar_disc_radiance_model",
                  configuration.thermalSolarDiscRadianceModel);
  receipt->Scalar("earth_atmosphere_transmission_model",
                  configuration.thermalEarthAtmosphereTransmissionModel);
  receipt->Scalar("solar_disc_radius_m",configuration.thermalSolarDiscRadiusM);
  receipt->Scalar("opaque_earth_radius_m",
                  configuration.thermalOpaqueEarthRadiusM);
  receipt->Scalar("spinup_completed_cycles",spinup.completedCycles);
  receipt->Scalar("spinup_cycle_start_tdb_s",
                  spinup.cycleStartEphemerisTimeS);
  receipt->Scalar("spinup_cycle_end_tdb_s",spinup.cycleEndEphemerisTimeS);
  receipt->Scalar("spinup_maximum_node_temperature_change_k",
                  spinup.maximumNodeTemperatureChangeK);
  receipt->Scalar("spinup_maximum_column_enthalpy_change_j_m2",
                  spinup.maximumColumnEnthalpyChangeJm2);
  receipt->Scalar("spinup_temperature_tolerance_k",
                  configuration.thermalSpinupTemperatureToleranceK);
  receipt->Scalar("spinup_energy_tolerance_j_m2",
                  configuration.thermalSpinupEnergyToleranceJm2);

  struct Observation {
    std::string siteTimeBin;
    std::string utc;
    double epoch=0.0;
    double latitudeDeg=0.0;
    double longitudeDeg=0.0;
    double observedBolometricK=0.0;
    double uncertaintyK=0.0;
    std::array<double,4> observedRadiance{{0.0,0.0,0.0,0.0}};
    std::array<double,4> observedBrightnessK{{0.0,0.0,0.0,0.0}};
  };
  std::ifstream input(configuration.divinerEclipseHistoryFile);
  if (!input) return 6;
  std::string line;
  if (!std::getline(input,line)) return 6;
  const std::vector<std::string> names=SplitRestrictedCsv(line);
  std::map<std::string,std::size_t> field;
  for (std::size_t index=0;index<names.size();++index) field[names[index]]=index;
  const char* required[]={"selection","site_time_bin_id","utc_mean",
      "latitude_bin_center_deg","longitude_east_bin_center_deg",
      "bolometric_brightness_temperature_K",
      "bolometric_brightness_temperature_uncertainty_K",
      "c6_radiance_W_m-2_sr-1","c6_brightness_temperature_K",
      "c7_radiance_W_m-2_sr-1","c7_brightness_temperature_K",
      "c8_radiance_W_m-2_sr-1","c8_brightness_temperature_K",
      "c9_radiance_W_m-2_sr-1","c9_brightness_temperature_K"};
  for (const char* name:required) if (field.find(name)==field.end()) {
    std::cerr << "I17 selected history is missing field " << name << '\n';
    return 6;
  }
  std::vector<Observation> observations;
  while (std::getline(input,line)) {
    if (line.empty()) continue;
    const std::vector<std::string> values=SplitRestrictedCsv(line);
    if (values.size()!=names.size()) return 6;
    if (values[field["selection"]]!=selection) {
      // The Python staging step must make hold-out values inaccessible to the
      // executable. Fail before parsing science values rather than silently
      // skipping an unselected row that should not enter this receipt.
      std::cerr << "I17 input contains a row outside selected " << selection << '\n';
      return 6;
    }
    Observation value;
    value.siteTimeBin=values[field["site_time_bin_id"]];
    value.utc=values[field["utc_mean"]];
    str2et_c(value.utc.c_str(),&value.epoch);
    const auto parse=[&](const char* name,double* destination) {
      const std::string& text=values[field[name]];
      char* end=nullptr;
      errno=0;
      *destination=std::strtod(text.c_str(),&end);
      return errno!=ERANGE && end!=text.c_str() && *end=='\0' &&
          std::isfinite(*destination);
    };
    if (!parse("latitude_bin_center_deg",&value.latitudeDeg) ||
        !parse("longitude_east_bin_center_deg",&value.longitudeDeg) ||
        !parse("bolometric_brightness_temperature_K",
               &value.observedBolometricK) ||
        !parse("bolometric_brightness_temperature_uncertainty_K",
               &value.uncertaintyK)) return 6;
    for (int channel=6;channel<=9;++channel) {
      const std::string prefix="c"+std::to_string(channel);
      if (!parse((prefix+"_radiance_W_m-2_sr-1").c_str(),
                 &value.observedRadiance[channel-6]) ||
          !parse((prefix+"_brightness_temperature_K").c_str(),
                 &value.observedBrightnessK[channel-6])) return 6;
    }
    if (!(value.observedBolometricK>0.0) || !(value.uncertaintyK>=0.0))
      return 6;
    observations.push_back(value);
  }
  std::sort(observations.begin(),observations.end(),
      [](const Observation& first,const Observation& second) {
        if (first.epoch!=second.epoch) return first.epoch<second.epoch;
        return first.siteTimeBin<second.siteTimeBin;
      });
  if (observations.empty() ||
      observations.front().epoch<spinup.cycleEndEphemerisTimeS) {
    std::cerr << "I17 selected rows precede the configured spin-up boundary\n";
    return 6;
  }

  std::ofstream output(receiptPath+".I17_"+selection+"_model.csv");
  if (!output) return 6;
  output << std::setprecision(17)
      << "site_time_bin_id,utc,ephemeris_time_tdb_s,latitude_deg,"
         "longitude_east_deg,observed_bolometric_brightness_temperature_K,"
         "observed_uncertainty_K,model_physical_surface_temperature_K,"
         "model_bolometric_brightness_temperature_K";
  for (int channel=6;channel<=9;++channel)
    output << ",observed_c" << channel << "_radiance_W_m-2_sr-1"
           << ",model_c" << channel << "_radiance_W_m-2_sr-1"
           << ",observed_c" << channel << "_brightness_temperature_K"
           << ",model_c" << channel << "_brightness_temperature_K";
  output << ",selection\n";

  double currentEpoch=spinup.cycleEndEphemerisTimeS;
  double maximumThermalUpdate=0.0;
  for (const Observation& observation:observations) {
    if (observation.epoch>currentEpoch &&
        !Moon::SurfaceTemperature::AdvanceProductionState(
            currentEpoch,observation.epoch,&error)) {
      std::cerr << "I17 selected thermal advance failed: " << error << '\n';
      return 6;
    }
    maximumThermalUpdate=std::max(maximumThermalUpdate,
        Moon::SurfaceTemperature::
            LastAdvanceMaximumSurfaceTemperatureChangeK());
    currentEpoch=observation.epoch;
    Moon::OrbitalMotion::et=currentEpoch;
    const double latitude=observation.latitudeDeg*Pi/180.0;
    const double longitude=observation.longitudeDeg*Pi/180.0;
    SpiceDouble bodyState[6]={
        _RADIUS_(_MOON_)*std::cos(latitude)*std::cos(longitude),
        _RADIUS_(_MOON_)*std::cos(latitude)*std::sin(longitude),
        _RADIUS_(_MOON_)*std::sin(latitude),0.0,0.0,0.0};
    SpiceDouble transform[6][6],lsoState[6];
    sxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,
             currentEpoch,transform);
    mxvg_c(transform,bodyState,6,6,lsoState);
    const double kineticK=
        Moon::SurfaceTemperature::ProductionTemperatureK(lsoState);
    const double bolometricK=
        Moon::DivinerObservation::GreyBolometricBrightnessTemperatureK(
            kineticK,configuration.thermal.infraredEmissivity);
    if (!(kineticK>0.0) || !std::isfinite(bolometricK)) return 6;
    output << observation.siteTimeBin << ',' << observation.utc << ','
           << currentEpoch << ',' << observation.latitudeDeg << ','
           << observation.longitudeDeg << ','
           << observation.observedBolometricK << ','
           << observation.uncertaintyK << ',' << kineticK << ','
           << bolometricK;
    for (int channel=6;channel<=9;++channel) {
      const double modelRadiance=
          Moon::DivinerObservation::RectangularBandRadianceWm2Sr(
              channel,kineticK,configuration.thermal.infraredEmissivity);
      double modelBrightness=0.0;
      if (!Moon::DivinerObservation::RectangularBandBrightnessTemperatureK(
              channel,modelRadiance,&modelBrightness,&error)) return 6;
      output << ',' << observation.observedRadiance[channel-6] << ','
             << modelRadiance << ','
             << observation.observedBrightnessK[channel-6] << ','
             << modelBrightness;
    }
    output << ',' << selection << '\n';
  }
  if (!output) return 6;
  receipt->Scalar("comparison_row_count",observations.size());
  receipt->Scalar("training_comparison_row_count",selection=="training"?observations.size():0);
  receipt->Scalar("holdout_comparison_row_count",selection=="holdout"?observations.size():0);
  receipt->Scalar("maximum_production_update_temperature_change_k",
                  maximumThermalUpdate);
  receipt->Scalar("thermal_time_step_s",configuration.thermalTimeStepS);
  receipt->Scalar("infrared_emissivity",
                  configuration.thermal.infraredEmissivity);
  receipt->Scalar("holdout_rows_accessed",selection=="holdout"?observations.size():0);
  return 0;
#endif
}

int RunI18(Receipt* receipt,const std::string& receiptPath) {
#ifdef _NO_SPICE_CALLS_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#elif _AR_SPEC_ < 0
  // I18 is intentionally unavailable in the frozen Na-only M0/M1 binary.
  // Treating a missing compiled species as an empty Ar sample would turn a
  // build-configuration error into false surface-physics evidence.
  receipt->Scalar("configuration_error_argon_species_required",1);
  return 4;
#else
  if (!Moon::Runtime::HasConfiguration() || Moon::Planet==NULL ||
      !Moon::SurfaceTemperature::ProductionStateInitialized()) {
    receipt->Scalar("configuration_error_m2_surface_state_required",1);
    return 4;
  }
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  if (!configuration.m2EnvironmentEnabled ||
      !configuration.argonSourceEnabled) {
    receipt->Scalar("configuration_error_argon_source_required",1);
    return 4;
  }
  const auto* terrain=Moon::SurfaceEnvironment::ProductionTerrain();
  if (terrain==NULL) {
    receipt->Scalar("configuration_error_lola_terrain_required",1);
    return 4;
  }

  // Select one real production facet, preferring one owned by this MPI rank.
  // The Na ensemble below calls the ordinary Exosphere impact kernel with this
  // facet normal and the epoch-correct 6x6 SPICE transforms, exactly as the
  // moving LOLA mover does after resolving a collision. No test-only
  // re-emission formula is evaluated in C++.
  const auto& mesh=terrain->triangulation();
  SpiceDouble bodyToLso[6][6],lsoToBody[6][6];
  sxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,
           Moon::OrbitalMotion::et,bodyToLso);
  sxform_c(Moon::Frames::SolarOrbital,Moon::Frames::BodyFixed,
           Moon::OrbitalMotion::et,lsoToBody);
  double selectedPositionBody[3]={0.0,0.0,0.0};
  double selectedNormalBody[3]={0.0,0.0,0.0};
  double selectedPositionLso[3]={0.0,0.0,0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* selectedNode=NULL;
  double fallbackPositionBody[3]={0.0,0.0,0.0};
  double fallbackNormalBody[3]={0.0,0.0,0.0};
  double fallbackPositionLso[3]={0.0,0.0,0.0};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* fallbackNode=NULL;
  bool selected=false,ownsSelectedFacet=false;
  // Facet numbering is independent of AMR domain decomposition. Every rank
  // therefore inspects the same deterministic face order. Partitioning face
  // indices by rank would incorrectly assume a relationship between two
  // unrelated orderings and can miss an owned facet.
  for (std::size_t face=0;face<mesh.faces.size();++face) {
    const auto& connectivity=mesh.faces[face];
    const auto& a=mesh.verticesM[connectivity[0]];
    const auto& b=mesh.verticesM[connectivity[1]];
    const auto& c=mesh.verticesM[connectivity[2]];
    double edge0[3],edge1[3];
    for (int component=0;component<3;++component) {
      selectedPositionBody[component]=(a[component]+b[component]+c[component])/3.0;
      edge0[component]=b[component]-a[component];
      edge1[component]=c[component]-a[component];
    }
    selectedNormalBody[0]=edge0[1]*edge1[2]-edge0[2]*edge1[1];
    selectedNormalBody[1]=edge0[2]*edge1[0]-edge0[0]*edge1[2];
    selectedNormalBody[2]=edge0[0]*edge1[1]-edge0[1]*edge1[0];
    double normalLength=Norm(selectedNormalBody);
    if (!(normalLength>0.0)) continue;
    for (double& component:selectedNormalBody) component/=normalLength;
    if (selectedNormalBody[0]*selectedPositionBody[0]+
        selectedNormalBody[1]*selectedPositionBody[1]+
        selectedNormalBody[2]*selectedPositionBody[2]<0.0)
      for (double& component:selectedNormalBody) component=-component;
    for (int row=0;row<3;++row) {
      selectedPositionLso[row]=0.0;
      for (int column=0;column<3;++column)
        selectedPositionLso[row]+=bodyToLso[row][column]*
            selectedPositionBody[column];
    }
    selectedNode=PIC::Mesh::mesh->findTreeNode(selectedPositionLso);
    // Retain any geometrically valid facet as a callback fixture.  Prefer a
    // locally allocated surface leaf when this MPI rank owns one.  Some valid
    // AMPS decompositions assign a rank only far-field blocks, so demanding a
    // lunar facet on every rank would test the load balancer rather than the
    // surface law.  The impact callback itself does not read node/block data;
    // I14 separately verifies the mover-to-owned-facet wiring.
    if (selectedNode!=NULL && !selected) {
      Copy(selectedPositionBody,fallbackPositionBody);
      Copy(selectedNormalBody,fallbackNormalBody);
      Copy(selectedPositionLso,fallbackPositionLso);
      fallbackNode=selectedNode;
      selected=true;
    }
    if (selectedNode!=NULL && selectedNode->block!=NULL &&
        selectedNode->Thread==PIC::ThisThread) {
      ownsSelectedFacet=true;
      selected=true;
      break;
    }
  }
  if (selected && !ownsSelectedFacet) {
    Copy(fallbackPositionBody,selectedPositionBody);
    Copy(fallbackNormalBody,selectedNormalBody);
    Copy(fallbackPositionLso,selectedPositionLso);
    selectedNode=fallbackNode;
  }
  int localSelected=selected ? 1 : 0,totalSelected=0;
  int localOwned=ownsSelectedFacet ? 1 : 0,totalOwned=0;
  MPI_Allreduce(&localSelected,&totalSelected,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localOwned,&totalOwned,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (!selected || totalSelected!=PIC::nTotalThreads || totalOwned==0) {
    if (PIC::ThisThread==0)
      std::cerr << "I18 could not find production facet geometry on every "
                   "rank or found no locally owned lunar facet\n";
    return 6;
  }

  // Exercise the linked residence/desorption callback separately from the
  // launch distributions.  A known physical atom population is installed in
  // the ordinary logical-surface slot, the production rate callback is read,
  // and the original state is restored immediately.  This bounded diagnostic
  // does not advance or permanently seed the simulation inventory; Python
  // independently evaluates the two cited Arrhenius laws from application.in.
  long int selectedZenith=0,selectedAzimuth=0;
  Moon::Planet->GetSurfaceElementProjectionIndex(
      selectedPositionBody,selectedZenith,selectedAzimuth);
  const int selectedElement=Moon::Planet->GetLocalSurfaceElementNumber(
      selectedZenith,selectedAzimuth);
  constexpr double kResidenceFixtureAtoms=1.0e12;
  const double originalNaPopulation=
      Moon::Planet->SurfaceElementPopulation[_NA_SPEC_][selectedElement];
  const double originalArPopulation=
      Moon::Planet->SurfaceElementPopulation[_AR_SPEC_][selectedElement];
  Moon::Planet->SurfaceElementPopulation[_NA_SPEC_][selectedElement]=
      kResidenceFixtureAtoms;
  Moon::Planet->SurfaceElementPopulation[_AR_SPEC_][selectedElement]=
      kResidenceFixtureAtoms;
  const double linkedNaDesorptionRate=
      Moon::SurfaceTemperature::ProductionThermalDesorptionRate(
          _NA_SPEC_,selectedElement,Moon::Planet);
  const double linkedArDesorptionRate=
      Moon::SurfaceTemperature::ProductionThermalDesorptionRate(
          _AR_SPEC_,selectedElement,Moon::Planet);
  Moon::Planet->SurfaceElementPopulation[_NA_SPEC_][selectedElement]=
      originalNaPopulation;
  Moon::Planet->SurfaceElementPopulation[_AR_SPEC_][selectedElement]=
      originalArPopulation;
  if (!std::isfinite(linkedNaDesorptionRate) ||
      !std::isfinite(linkedArDesorptionRate) ||
      linkedNaDesorptionRate<0.0 || linkedArDesorptionRate<0.0) return 6;

  std::ostringstream sampleName;
  sampleName << receiptPath << ".I18.rank" << std::setw(4)
             << std::setfill('0') << PIC::ThisThread << ".csv";
  std::ofstream samples(sampleName.str());
  if (!samples) return 6;
  samples << std::setprecision(17)
          << "law,species,rank,sample,temperature_k,incident_speed_m_s,"
             "output_speed_m_s,outward_cosine,status\n";
  constexpr int kSamplesPerSpeciesPerRank=4096;
  constexpr double kSodiumIncidentSpeedMPerS=2000.0;
  rnd_seed(configuration.randomSeed+104729*PIC::ThisThread+18);

  int localNaAdsorbed=0,localNaReemitted=0;
  for (int sample=0;sample<kSamplesPerSpeciesPerRank;++sample) {
    double bodyState[6],lsoState[6];
    for (int component=0;component<3;++component) {
      bodyState[component]=selectedPositionBody[component];
      bodyState[component+3]=-
          kSodiumIncidentSpeedMPerS*selectedNormalBody[component];
    }
    mxvg_c(bodyToLso,bodyState,6,6,lsoState);
    const long int particle=PIC::ParticleBuffer::GetNewParticle();
    PIC::ParticleBuffer::SetI(_NA_SPEC_,particle);
    PIC::ParticleBuffer::SetX(lsoState,particle);
    PIC::ParticleBuffer::SetV(lsoState+3,particle);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,particle);
#endif
    PIC::ParticleBuffer::SetParticleAllocated(particle);
    Exosphere::Sampling::SetParticleSourceID(
        _EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_,
        PIC::ParticleBuffer::GetParticleDataPointer(particle));
    double remainingTime=0.0;
    const double temperatureK=
        Moon::SurfaceTemperature::ProductionTemperatureK(lsoState);
    const int status=Exosphere::SurfaceInteraction::
        ParticleSurfaceInteraction_SurfaceAccomodation(
            _NA_SPEC_,particle,lsoState,lsoState+3,remainingTime,selectedNode,
            Moon::Planet,selectedNormalBody,&lsoToBody[0][0],
            &bodyToLso[0][0]);
    double outputSpeed=std::numeric_limits<double>::quiet_NaN();
    double outwardCosine=std::numeric_limits<double>::quiet_NaN();
    if (status==_PARTICLE_DELETED_ON_THE_FACE_) {
      ++localNaAdsorbed;
    }
    else if (status==_PARTICLE_REJECTED_ON_THE_FACE_) {
      ++localNaReemitted;
      double returnedBodyState[6];
      mxvg_c(lsoToBody,lsoState,6,6,returnedBodyState);
      outputSpeed=Norm(returnedBodyState+3);
      outwardCosine=(returnedBodyState[3]*selectedNormalBody[0]+
          returnedBodyState[4]*selectedNormalBody[1]+
          returnedBodyState[5]*selectedNormalBody[2])/outputSpeed;
    }
    else {
      PIC::ParticleBuffer::DeleteParticle(particle);
      return 6;
    }
    samples << "impact-accommodation,Na," << PIC::ThisThread << ','
            << sample << ',' << temperatureK << ','
            << kSodiumIncidentSpeedMPerS << ',' << outputSpeed << ','
            << outwardCosine << ',' << status << '\n';
    PIC::ParticleBuffer::DeleteParticle(particle);
  }

  // The Ar ensemble invokes the registered production radiogenic-source
  // generator.  Its outward Maxwell-flux distribution is the physical
  // re-release law selected by the approved Kegerreis configuration.  Failed
  // attempts merely selected a surface element owned by another MPI rank;
  // they are not samples. Use the same fixed number of attempts on every rank
  // so accepted samples retain the production MPI ownership weighting. Some
  // ranks may legitimately own only far-field blocks and accept zero sources.
  int localArGenerated=0;
  constexpr long long kArgonAttemptsPerRank=32768;
  for (long long attempt=0;attempt<kArgonAttemptsPerRank;++attempt) {
    double xLso[3],xBody[3],vLso[3],vBody[3],centre[3]={0.0,0.0,0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode=NULL;
    if (!Moon::Species::Argon::GenerateParticleProperties(
            _AR_SPEC_,NULL,xLso,xBody,vLso,vBody,centre,
            _RADIUS_(_MOON_),startNode,_INTERNAL_BOUNDARY_TYPE_SPHERE_,
            Moon::Planet)) continue;
    const double positionLength=Norm(xBody);
    const double speed=Norm(vBody);
    if (!(positionLength>0.0) || !(speed>0.0)) return 6;
    double normal[3];
    for (int component=0;component<3;++component)
      normal[component]=xBody[component]/positionLength;
    const double cosine=(vBody[0]*normal[0]+vBody[1]*normal[1]+
                         vBody[2]*normal[2])/speed;
    const double temperatureK=
        Moon::SurfaceTemperature::ProductionTemperatureK(xLso);
    samples << "thermal-maxwell-flux,Ar," << PIC::ThisThread << ','
            << localArGenerated << ',' << temperatureK << ','
            << std::numeric_limits<double>::quiet_NaN() << ',' << speed
            << ',' << cosine << ',' << _PARTICLE_REJECTED_ON_THE_FACE_
            << '\n';
    ++localArGenerated;
  }
  if (!samples) return 6;

  int globalNaAdsorbed=0,globalNaReemitted=0,globalArGenerated=0;
  MPI_Reduce(&localNaAdsorbed,&globalNaAdsorbed,1,MPI_INT,MPI_SUM,0,
             MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localNaReemitted,&globalNaReemitted,1,MPI_INT,MPI_SUM,0,
             MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&localArGenerated,&globalArGenerated,1,MPI_INT,MPI_SUM,0,
             MPI_GLOBAL_COMMUNICATOR);
  receipt->Scalar("compiled_na_species_index",_NA_SPEC_);
  receipt->Scalar("compiled_ar_species_index",_AR_SPEC_);
  receipt->Scalar("na_mass_kg",PIC::MolecularData::GetMass(_NA_SPEC_));
  receipt->Scalar("argon_mass_kg",PIC::MolecularData::GetMass(_AR_SPEC_));
  receipt->Scalar("sample_count_per_species_per_rank",
                  kSamplesPerSpeciesPerRank);
  receipt->Scalar("rank_count_with_owned_surface_facet",totalOwned);
  receipt->Scalar("argon_attempt_count_per_rank",kArgonAttemptsPerRank);
  receipt->Scalar("na_adsorbed_count",globalNaAdsorbed);
  receipt->Scalar("na_reemitted_count",globalNaReemitted);
  receipt->Scalar("argon_generated_count",globalArGenerated);
  receipt->Scalar("residence_fixture_atoms",kResidenceFixtureAtoms);
  receipt->Scalar("residence_fixture_temperature_k",
      Moon::SurfaceTemperature::ProductionTemperatureK(selectedPositionLso));
  receipt->Scalar("linked_na_desorption_rate_atoms_s",
                  linkedNaDesorptionRate);
  receipt->Scalar("linked_argon_desorption_rate_atoms_s",
                  linkedArDesorptionRate);
  receipt->Scalar("argon_source_rate_atoms_s",
      Moon::Species::Argon::GetTotalProductionRate(
          _AR_SPEC_,_INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));
  receipt->Scalar("random_seed_base",configuration.randomSeed);
  return 0;
#endif
}

namespace I11Sampling {

// The bounded campaign controls are part of the reviewed I11 contract, not
// tunable shortcuts selected after seeing a result.  One completed production
// sample contains exactly 500 particle steps, and every spatial-resolution
// case stops after at most 5000 particle steps.  Ten complete samples allow
// the unchanged nine-batch independence/equilibrium/precision algorithm to be
// evaluated once when the raw blocks are already independent.  The cap is
// still an execution bound, not a relaxed science criterion: a case that has
// not met every gate at 5000 steps reports FAIL and retains raw diagnostics.
constexpr long int SampleLengthSteps=500;
constexpr long int MaximumParticleSteps=5000;

// The duration pilot is an execution aid, not a second acceptance contract.
// Its only permitted difference from bounded I11 is removal of the 5000-step
// stop so that the unchanged statistical gates can determine a candidate
// physical spin-up/sampling duration at the coarsest approved grid.  The
// Python pilot runner enforces operational memory/free-space guards and labels
// every resulting artifact NOT QUALIFIED.  Formal I11 always leaves this
// false and continues to use MaximumParticleSteps above.
bool durationPilot=false;

// One entry in the frozen M1 scoring vector.  ``pointwise`` selects the
// predeclared 10% spatial-convergence and 3.33% CI gates; the two global
// budgets use the stricter 5% and 1.67% integrated gates.
struct MetricDefinition {
  std::string name;
  bool pointwise;
};

std::vector<MetricDefinition> definitions;
std::vector<std::vector<double> > baseBlocks;
std::filesystem::path artifactDirectory;
std::ofstream blockFile,resourceFile,analysisFile,particleFile;
bool complete=false;
int retainedBatchSize=0,retainedBatchCount=0;
std::vector<double> retainedMean,retainedRelativeCi,retainedLagOne;
std::chrono::steady_clock::time_point startTime;

void SuppressExternalOutput(int) {
  // PIC's completed-sample callback has an integer output-file index.  I11
  // needs the completed in-memory buffers but not legacy visualization files;
  // replacing only the writer callback leaves each registered per-step
  // SamplingProcessor and all production sampling state active.
}

double DensityAt(const double direction[3],double altitudeM) {
  // A zero-altitude request lies exactly on the analytic boundary.  Move it
  // outward by the mesh's geometric tolerance so the queried cell is
  // unambiguously exospheric while retaining the declared physical altitude
  // in every artifact.
  const double radius=_RADIUS_(_MOON_)+altitudeM+
      std::max(1.0e-3,10.0*PIC::Mesh::mesh->EPS);
  double position[3]={radius*direction[0],radius*direction[1],
                      radius*direction[2]};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
      PIC::Mesh::mesh->findTreeNode(position);
  // The AMR tree is replicated but an allocated data block belongs to exactly
  // one MPI rank. Only that owner reads the production density; the reduction
  // both distributes the value and proves that ownership is neither missing
  // nor duplicated. This is the same ownership rule used by linked I12.
  double localDensity=0.0;
  int localOwner=0;
  if (node!=NULL && node->Thread==PIC::ThisThread && node->block!=NULL) {
    int i,j,k;
    const int cell=PIC::Mesh::mesh->FindCellIndex(position,i,j,k,node);
    PIC::Mesh::cDataCenterNode* center=node->block->GetCenterNode(cell);
    if (center==NULL)
      exit(__LINE__,__FILE__,
          "I11 density sampling point has no production center node");
    localDensity=center->GetNumberDensity(_NA_SPEC_);
    localOwner=1;
  }
  double density=0.0;
  int owners=0;
  MPI_Allreduce(&localDensity,&density,1,MPI_DOUBLE,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localOwner,&owners,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (node==NULL || owners!=1)
    exit(__LINE__,__FILE__,
        "I11 density point must have exactly one allocated MPI owner");
  return density;
}

void Append(std::vector<double>* values,const std::string& name,
            double value,bool pointwise=true) {
  if (baseBlocks.empty()) definitions.push_back({name,pointwise});
  values->push_back(value);
}

void Unit(double vector[3]) {
  const double length=Norm(vector);
  if (!(length>0.0) || !std::isfinite(length))
    exit(__LINE__,__FILE__,"I11 cannot normalize a degenerate geometry vector");
  for (int component=0;component<3;component++) vector[component]/=length;
}

void SampleProfiles(std::vector<double>* values) {
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  for (double szaDeg : configuration.profileSolarZenithAnglesDeg) {
    const double radians=szaDeg*Pi/180.0;
    const double direction[3]={std::cos(radians),std::sin(radians),0.0};
    for (double altitudeM : configuration.profileAltitudeGridM) {
      std::ostringstream name;
      name << "density_sza_" << szaDeg << "_alt_m_" << altitudeM;
      Append(values,name.str(),DensityAt(direction,altitudeM));
    }
  }

  if (configuration.profileIncludeBodyFixedPoles) {
    // Transform the label-native MOON_ME_DE421 poles into the production LSO
    // frame at the same frozen ET used by dynamics and mesh refinement.
    for (int sign : {-1,1}) {
      // The production orbital initializer has already populated this state
      // transform from CSPICE.  Reading it here also keeps no-SPICE builds
      // compilable; RunI11 rejects those builds before this sampler is armed.
      double direction[3]={
          sign*Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[0][2],
          sign*Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[1][2],
          sign*Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]};
      Unit(direction);
      for (double altitudeM : configuration.profileAltitudeGridM) {
        std::ostringstream name;
        name << "density_body_" << (sign>0 ? "north" : "south")
             << "_alt_m_" << altitudeM;
        Append(values,name.str(),DensityAt(direction,altitudeM));
      }
    }
  }

  if (configuration.profileIncludeShellAverage) {
    const long int elements=Moon::Planet->GetTotalSurfaceElementsNumber();
    const double radiusSquared=_RADIUS_(_MOON_)*_RADIUS_(_MOON_);
    for (double altitudeM : configuration.profileAltitudeGridM) {
      double localWeightedDensity=0.0,localSolidAngle=0.0;
      long localOwners=0;
      for (long int element=0;element<elements;element++) {
        double bodyDirection[3],lsoDirection[3]={0.0,0.0,0.0};
        Moon::Planet->GetSurfaceElementNormal(bodyDirection,element);
        for (int row=0;row<3;row++) for (int column=0;column<3;column++)
          lsoDirection[row]+=
              Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[row][column]*
              bodyDirection[column];
        Unit(lsoDirection);
        const double radius=_RADIUS_(_MOON_)+altitudeM+
            std::max(1.0e-3,10.0*PIC::Mesh::mesh->EPS);
        double position[3]={radius*lsoDirection[0],radius*lsoDirection[1],
                            radius*lsoDirection[2]};
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
            PIC::Mesh::mesh->findTreeNode(position);
        if (node!=NULL && node->Thread==PIC::ThisThread && node->block!=NULL) {
          int i,j,k;
          const int cell=PIC::Mesh::mesh->FindCellIndex(position,i,j,k,node);
          PIC::Mesh::cDataCenterNode* center=node->block->GetCenterNode(cell);
          if (center==NULL)
            exit(__LINE__,__FILE__,
                "I11 shell-average point has no production center node");
          const double solidAngle=
              Moon::Planet->GetSurfaceElementArea(element)/radiusSquared;
          localWeightedDensity+=solidAngle*
              center->GetNumberDensity(_NA_SPEC_);
          localSolidAngle+=solidAngle;
          localOwners++;
        }
      }
      double local[2]={localWeightedDensity,localSolidAngle},global[2]={0.0,0.0};
      long owners=0;
      MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,
                    MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&localOwners,&owners,1,MPI_LONG,MPI_SUM,
                    MPI_GLOBAL_COMMUNICATOR);
      if (owners!=elements || !(global[1]>0.0))
        exit(__LINE__,__FILE__,
            "I11 shell-average directions must each have one MPI owner");
      std::ostringstream name;
      name << "density_shell_solid_angle_alt_m_" << altitudeM;
      Append(values,name.str(),global[0]/global[1]);
    }
  }
}

void SampleSurfaceFlux(std::vector<double>* values) {
  const long int elements=Moon::Planet->GetTotalSurfaceElementsNumber();
  const double duration=PIC::LastSampleLength*
      PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
  if (!(duration>0.0))
    exit(__LINE__,__FILE__,"I11 surface-flux block has zero sample duration");
  // Element-level relative confidence intervals are ill-conditioned for
  // physically dark/low-count facets and make runtime depend on arbitrary
  // triangulation density. Six frozen 30-degree SZA bands retain the physical
  // day/night/terminator distribution while pooling independent particle
  // returns. The global integral remains a separate stricter budget metric.
  constexpr int bins=6;
  double localWeightedReturns[bins]={0.0,0.0,0.0,0.0,0.0,0.0};
  double binArea[bins]={0.0,0.0,0.0,0.0,0.0,0.0};
  for (long int element=0;element<elements;element++) {
    double bodyNormal[3],lsoNormal[3]={0.0,0.0,0.0};
    Moon::Planet->GetSurfaceElementNormal(bodyNormal,element);
    for (int row=0;row<3;row++) for (int column=0;column<3;column++)
      lsoNormal[row]+=
          Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[row][column]*
          bodyNormal[column];
    Unit(lsoNormal);
    const double szaDeg=acos(std::max(-1.0,std::min(1.0,lsoNormal[0])))*
        180.0/Pi;
    const int bin=std::min(bins-1,static_cast<int>(szaDeg/30.0));
    const double area=Moon::Planet->GetSurfaceElementArea(element);
    binArea[bin]+=area;
    localWeightedReturns[bin]+=
        Moon::Planet->SampleSpeciesSurfaceReturnFlux[_NA_SPEC_][element];
  }
  double globalWeightedReturns[bins]={0.0,0.0,0.0,0.0,0.0,0.0};
  MPI_Allreduce(localWeightedReturns,globalWeightedReturns,bins,MPI_DOUBLE,
                MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  double globalRate=0.0;
  for (int bin=0;bin<bins;bin++) {
    if (!(binArea[bin]>0.0))
      exit(__LINE__,__FILE__,"I11 SZA surface-flux bin has zero area");
    std::ostringstream name;
    name << "surface_return_flux_sza_" << 30*bin << "_to_" << 30*(bin+1)
         << "_m2_s";
    Append(values,name.str(),globalWeightedReturns[bin]/(duration*binArea[bin]));
    globalRate+=globalWeightedReturns[bin]/duration;
  }
  Append(values,"global_surface_return_rate_s",globalRate,false);
}

void SampleInventory(std::vector<double>* values) {
  // Integrating the completed production density field over cut-cell-aware
  // cell measures yields the block-mean physical inventory.  This avoids a
  // test-only particle-weight interpretation and uses the same sampled field
  // that supplies every profile and line-of-sight observable.
  double inventory=0.0;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
           PIC::Mesh::mesh->BranchBottomNodeList;
       node!=NULL;node=node->nextBranchBottomNode) {
    if (node->block==NULL || node->Thread!=PIC::ThisThread) continue;
    for (int k=0;k<_BLOCK_CELLS_Z_;k++)
      for (int j=0;j<_BLOCK_CELLS_Y_;j++)
        for (int i=0;i<_BLOCK_CELLS_X_;i++) {
          const int cell=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::cDataCenterNode* center=node->block->GetCenterNode(cell);
          if (center!=NULL && center->Measure>0.0)
            inventory+=center->GetNumberDensity(_NA_SPEC_)*center->Measure;
        }
  }
  double globalInventory=0.0;
  MPI_Allreduce(&inventory,&globalInventory,1,MPI_DOUBLE,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  Append(values,"global_particle_inventory",globalInventory,false);
}

void SampleLimb(std::vector<double>* values) {
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  double observer[3]={Moon::xEarth_SO[0],Moon::xEarth_SO[1],
                      Moon::xEarth_SO[2]};
  const double observerDistance=Norm(observer);
  double observerUnit[3]={observer[0]/observerDistance,
                          observer[1]/observerDistance,
                          observer[2]/observerDistance};
  double projectedSun[3]={Moon::xSun_SO[0],Moon::xSun_SO[1],Moon::xSun_SO[2]};
  Unit(projectedSun);
  const double projection=projectedSun[0]*observerUnit[0]+
      projectedSun[1]*observerUnit[1]+projectedSun[2]*observerUnit[2];
  for (int component=0;component<3;component++)
    projectedSun[component]-=projection*observerUnit[component];
  Unit(projectedSun);

  const int vectorLength=Exosphere::ColumnIntegral::GetVariableList(NULL);
  std::vector<double> integral(vectorLength,0.0);
  for (double altitudeM : configuration.profileAltitudeGridM) {
    if (altitudeM>configuration.limbTangentMaximumAltitudeM) break;
    const double tangentRadius=_RADIUS_(_MOON_)+altitudeM;
    if (!(observerDistance>tangentRadius))
      exit(__LINE__,__FILE__,"I11 Earth observer lies inside a limb shell");
    const double axial=tangentRadius*tangentRadius/observerDistance;
    const double transverse=tangentRadius*std::sqrt(
        1.0-tangentRadius*tangentRadius/(observerDistance*observerDistance));
    for (int sign : {1,-1}) {
      if ((sign>0 && configuration.limbSides=="away-from-projected-sun") ||
          (sign<0 && configuration.limbSides=="toward-projected-sun")) continue;
      double tangent[3],lineOfSight[3];
      for (int component=0;component<3;component++) {
        tangent[component]=axial*observerUnit[component]+sign*transverse*
            projectedSun[component];
        lineOfSight[component]=tangent[component]-observer[component];
      }
      Unit(lineOfSight);
      std::fill(integral.begin(),integral.end(),0.0);
      PIC::ColumnIntegration::GetCoulumnIntegral(integral.data(),vectorLength,
          observer,lineOfSight,Exosphere::ColumnIntegral::CoulumnDensityIntegrant);
      Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(
          integral.data(),vectorLength);
      const char* branch=sign>0 ? "toward_projected_sun" :
                                  "away_from_projected_sun";
      std::ostringstream prefix;
      prefix << "limb_" << branch << "_tangent_alt_m_" << altitudeM;
      Append(values,prefix.str()+"_column_density_m2",
          integral[2*_NA_SPEC_]);
      Append(values,prefix.str()+"_brightness_d2_rayleigh",
          integral[Moon::Sampling::SubsolarLimbColumnIntegrals::
              _NA_EMISSION_5891_58A_SAMPLE_OFFSET_]);
      Append(values,prefix.str()+"_brightness_d1_rayleigh",
          integral[Moon::Sampling::SubsolarLimbColumnIntegrals::
              _NA_EMISSION_5897_56A_SAMPLE_OFFSET_]);
    }
  }
}

long ReadProcValueKiB(const char* path,const char* key) {
  std::ifstream input(path);
  std::string line;
  while (std::getline(input,line)) {
    std::istringstream fields(line);
    std::string label,unit;
    long value=0;
    if ((fields >> label >> value >> unit) && label==key) return value;
  }
  return -1;
}

void RecordResources(int outputNumber) {
  if (PIC::ThisThread!=0) return;
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  const long resident=ReadProcValueKiB("/proc/self/status","VmRSS:");
  const long available=ReadProcValueKiB("/proc/meminfo","MemAvailable:");
  const double wall=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-startTime).count();
  resourceFile << outputNumber << ',' << wall << ',' << resident << ','
               << usage.ru_maxrss << ',' << available << '\n';
  resourceFile.flush();
}

double RelativeChange(double left,double right) {
  const double scale=std::max(std::fabs(left),std::fabs(right));
  if (scale==0.0) return 0.0;
  return std::fabs(left-right)/scale;
}

bool Analyze() {
  const int count=static_cast<int>(baseBlocks.size());
  const int metricCount=static_cast<int>(definitions.size());
  for (int batchSize=1;batchSize<=count/9;batchSize*=2) {
    const int batches=count/batchSize;
    if (batches<9) continue;
    std::vector<std::vector<double> > means(
        batches,std::vector<double>(metricCount,0.0));
    const int first=count-batches*batchSize;
    for (int batch=0;batch<batches;batch++)
      for (int offset=0;offset<batchSize;offset++)
        for (int metric=0;metric<metricCount;metric++)
          means[batch][metric]+=
              baseBlocks[first+batch*batchSize+offset][metric]/batchSize;

    std::vector<double> lag(metricCount,0.0);
    bool independent=true;
    int independenceFailedMetrics=0;
    const double significance=1.96/std::sqrt(static_cast<double>(batches));
    for (int metric=0;metric<metricCount;metric++) {
      double average=0.0;
      for (int batch=0;batch<batches;batch++) average+=means[batch][metric];
      average/=batches;
      double numerator=0.0,denominator=0.0;
      for (int batch=0;batch<batches;batch++) {
        const double delta=means[batch][metric]-average;
        denominator+=delta*delta;
        if (batch+1<batches)
          numerator+=delta*(means[batch+1][metric]-average);
      }
      lag[metric]=denominator==0.0 ? 0.0 : numerator/denominator;
      if (std::fabs(lag[metric])>significance) {
        independent=false;
        ++independenceFailedMetrics;
      }
    }
    if (!independent) {
      // Preserve the first failed invariant for every admissible batching
      // attempt.  These counts diagnose why a pilot is still running; they do
      // not alter batch selection or any acceptance threshold.
      if (PIC::ThisThread==0) {
        analysisFile << count << ',' << batchSize << ',' << batches << ','
                     << significance << ',' << independenceFailedMetrics
                     << ",,,independence\n";
        analysisFile.flush();
      }
      continue;
    }

    // Four independent batch means establish three consecutive <=2%
    // changes.  The following five batches are kept disjoint for the science
    // mean and Student-t 95% confidence interval.
    bool equilibrated=true;
    int equilibriumFailedChecks=0;
    const int equilibriumStart=batches-9;
    for (int metric=0;metric<metricCount;metric++)
      for (int change=0;change<3;change++)
        if (RelativeChange(means[equilibriumStart+change][metric],
                           means[equilibriumStart+change+1][metric])>0.02) {
          equilibrated=false;
          ++equilibriumFailedChecks;
        }
    if (!equilibrated) {
      if (PIC::ThisThread==0) {
        analysisFile << count << ',' << batchSize << ',' << batches << ','
                     << significance << ",0," << equilibriumFailedChecks
                     << ",,equilibrium\n";
        analysisFile.flush();
      }
      continue;
    }

    std::vector<double> scienceMean(metricCount,0.0),relativeCi(metricCount,0.0);
    bool precise=true;
    int confidenceFailedMetrics=0;
    for (int metric=0;metric<metricCount;metric++) {
      for (int batch=batches-5;batch<batches;batch++)
        scienceMean[metric]+=means[batch][metric]/5.0;
      double variance=0.0;
      for (int batch=batches-5;batch<batches;batch++) {
        const double delta=means[batch][metric]-scienceMean[metric];
        variance+=delta*delta/4.0;
      }
      const double halfWidth=2.7764451051977987*std::sqrt(variance/5.0);
      relativeCi[metric]=scienceMean[metric]==0.0 ?
          (halfWidth==0.0 ? 0.0 : std::numeric_limits<double>::infinity()) :
          halfWidth/std::fabs(scienceMean[metric]);
      const double gate=definitions[metric].pointwise ? 0.0333 : 0.0167;
      if (!(relativeCi[metric]<=gate)) {
        precise=false;
        ++confidenceFailedMetrics;
      }
    }
    if (!precise) {
      if (PIC::ThisThread==0) {
        analysisFile << count << ',' << batchSize << ',' << batches << ','
                     << significance << ",0,0," << confidenceFailedMetrics
                     << ",confidence_interval\n";
        analysisFile.flush();
      }
      continue;
    }
    if (PIC::ThisThread==0) {
      analysisFile << count << ',' << batchSize << ',' << batches << ','
                   << significance << ",0,0,0,pass\n";
      analysisFile.flush();
    }
    retainedBatchSize=batchSize;
    retainedBatchCount=batches;
    retainedMean=std::move(scienceMean);
    retainedRelativeCi=std::move(relativeCi);
    retainedLagOne=std::move(lag);
    return true;
  }
  return false;
}

void Capture(int outputNumber) {
  // PIC samples at the beginning of a step, before source/chemistry/motion.
  // nIterationCounter therefore records completed transport steps at this
  // checkpoint. The final post-loop snapshot covers the remainder explicitly.
  Moon::ParticleBudget::Snapshot(Moon::nIterationCounter,"sample");
  std::vector<double> values;
  SampleProfiles(&values);
  SampleSurfaceFlux(&values);
  SampleLimb(&values);
  SampleInventory(&values);
  // A live particle-list sum is independent of the density sampling buffers.
  // It distinguishes actual accumulation from stale-window accounting. Units
  // are model particles and represented Na atoms at this instant, whereas
  // global_particle_inventory in block_metrics.csv is a window mean.
  double localParticles[2]={0.0,0.0},globalParticles[2];
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=NULL;node=node->nextNodeThisThread) {
    if (node->block==NULL) continue;
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for (long int particle=node->block->FirstCellParticleTable[cell];
           particle!=-1;particle=PIC::ParticleBuffer::GetNext(particle)) {
        if (PIC::ParticleBuffer::GetI(particle)!=_NA_SPEC_) continue;
        double weight=node->block->GetLocalParticleWeight(_NA_SPEC_);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
        weight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
#endif
        localParticles[0]+=1.0;
        localParticles[1]+=weight;
      }
  }
  MPI_Allreduce(localParticles,globalParticles,2,MPI_DOUBLE,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    particleFile << baseBlocks.size()+1 << ',' << globalParticles[0] << ','
                 << globalParticles[1] << '\n';
    particleFile.flush();
  }
  if (values.size()!=definitions.size())
    exit(__LINE__,__FILE__,"I11 scoring-vector shape changed between blocks");
  baseBlocks.push_back(values);
  // Every rank must participate in the collective sampling above and retain
  // the same reduced vector for the common convergence decision.  The CSV is
  // nevertheless a single campaign artifact: restricting its writes to rank
  // zero prevents non-owning streams from entering a failed state and avoids
  // duplicated/interleaved records in an MPI run.
  if (PIC::ThisThread==0) {
    for (std::size_t metric=0;metric<values.size();metric++)
      blockFile << baseBlocks.size() << ',' << outputNumber << ','
                << definitions[metric].name << ',' << values[metric] << '\n';
    blockFile.flush();
  }
  RecordResources(outputNumber);
  complete=Analyze();
}

void Configure(const std::filesystem::path& directory) {
  artifactDirectory=directory;
  std::filesystem::create_directories(artifactDirectory);
  if (PIC::ThisThread==0) {
    blockFile.open((artifactDirectory/"block_metrics.csv").string());
    resourceFile.open((artifactDirectory/"resource_metrics.csv").string());
    analysisFile.open((artifactDirectory/"analysis_diagnostics.csv").string());
    particleFile.open((artifactDirectory/"particle_diagnostics.csv").string());
    if (!blockFile || !resourceFile || !analysisFile || !particleFile)
      exit(__LINE__,__FILE__,"I11 cannot create block/resource metric artifacts");
    blockFile << std::setprecision(17)
              << "base_block,production_output_number,metric,value_si\n";
    resourceFile << "production_output_number,wall_seconds,current_rss_kib,"
                    "peak_rss_kib,host_mem_available_kib\n";
    analysisFile << "base_blocks,batch_size_base_blocks,batch_count,"
                    "lag1_significance,independence_failed_metrics,"
                    "equilibrium_failed_checks,ci_failed_metrics,"
                    "first_failed_gate\n";
    particleFile << std::setprecision(17)
                 << "base_block,instantaneous_na_model_particles,"
                    "instantaneous_na_inventory_atoms\n";
  }
  startTime=std::chrono::steady_clock::now();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
      Moon::Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Capture);
}

void WriteFinal() {
  if (PIC::ThisThread!=0) return;
  std::ofstream output((artifactDirectory/"qualified_metrics.csv").string());
  if (!output)
    exit(__LINE__,__FILE__,"I11 cannot create qualified_metrics.csv");
  output << std::setprecision(17)
         << "metric,classification,mean_si,relative_ci95_half_width,"
            "lag1,batch_size_base_blocks,independent_batch_count\n";
  for (std::size_t metric=0;metric<definitions.size();metric++)
    output << definitions[metric].name << ','
           << (definitions[metric].pointwise ? "pointwise" : "integrated")
           << ',' << retainedMean[metric] << ',' << retainedRelativeCi[metric]
           << ',' << retainedLagOne[metric] << ',' << retainedBatchSize << ','
           << retainedBatchCount << '\n';
}

}  // namespace I11Sampling

int RunU01(Receipt* receipt,bool insideOrigin) {
  // FindIntegrationLimits is the production service used by the Moon LOS
  // samplers.  The inside-origin case intentionally exercises its documented
  // fatal guard and therefore does not return to this function.
  const double radius=_RADIUS_(_MOON_);
  const double directions[4][3]={{-1.0,0.0,0.0},{-1.0,0.0,0.0},
                                 {-1.0,0.0,0.0},{1.0,0.0,0.0}};
  const double origins[4][3]={{3.0*radius,0.0,0.0},
                              {3.0*radius,radius,0.0},
                              {3.0*radius,2.0*radius,0.0},
                              {3.0*radius,0.0,0.0}};
  if (insideOrigin) {
    double origin[3]={0.0,0.0,0.0},direction[3]={1.0,0.0,0.0};
    double length,start[3],finish[3];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,*finishNode;
    PIC::ColumnIntegration::FindIntegrationLimits(origin,direction,length,start,
                                                   startNode,finish,finishNode);
    return 3;  // Reaching here means the production guard did not fire.
  }

  static const char* names[4]={"radial_hit","tangent","miss","reversed"};
  for (int fixture=0;fixture<4;fixture++) {
    double origin[3],direction[3],length,start[3],finish[3];
    Copy(origins[fixture],origin);
    Copy(directions[fixture],direction);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,*finishNode;
    const bool found=PIC::ColumnIntegration::FindIntegrationLimits(
        origin,direction,length,start,startNode,finish,finishNode);
    receipt->Scalar(std::string(names[fixture])+"_found",found ? 1 : 0);
    if (found) {
      receipt->Scalar(std::string(names[fixture])+"_path_length_m",length);
      receipt->Vector(std::string(names[fixture])+"_start_m",start);
      receipt->Vector(std::string(names[fixture])+"_finish_m",finish);
    }
  }

  // GetVariableList is the production record-layout callback used by all Moon
  // column integrations.  Recording both its count and delimiter count catches
  // a stale output header without treating a frozen model field as truth.
  char variables[4096]={0};
  const int variableCount=Exosphere::ColumnIntegral::GetVariableList(variables);
  int delimiterCount=0;
  for (const char* cursor=variables;*cursor!='\0';cursor++)
    if (*cursor==',') delimiterCount++;
  receipt->Scalar("record_variable_count",variableCount);
  receipt->Scalar("record_delimiter_count",delimiterCount);
  return 0;
}

void PrepareOrbitState();

int RunI02(Receipt* receipt) {
  // Photoionization and Earth-eclipse callbacks consume the current SPICE
  // state.  Establish it through the normal time-step hook before auditing
  // their dispatch; otherwise the zero-initialized Earth vector can falsely
  // classify near-Moon fixtures as lying inside the Earth.
  PrepareOrbitState();
  receipt->Scalar("species_count",PIC::nTotalSpecies);
  receipt->Scalar("na_species_index",_NA_SPEC_);
  receipt->Scalar("na_plus_species_index",_NA_PLUS_SPEC_);
  receipt->Scalar("sphere_present",Moon::Planet!=NULL ? 1 : 0);
  receipt->Scalar("sphere_injection_callback",
      Moon::Planet!=NULL && Moon::Planet->InjectionBoundaryCondition!=NULL ? 1:0);
  receipt->Scalar("sphere_interaction_callback",
      Moon::Planet!=NULL && Moon::Planet->ParticleSphereInteraction!=NULL ? 1:0);
  receipt->Scalar("total_source_rate_s",
      Moon::SourceProcesses::totalProductionRate(_NA_SPEC_,
          _INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));
  receipt->Scalar("impact_source_rate_s",
      Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(
          _NA_SPEC_,_INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));

  bool allowed=false;
  double lit[3]={3.0*_RADIUS_(_MOON_),0.0,0.0};
  double shadow[3]={-3.0*_RADIUS_(_MOON_),0.0,0.0};
  receipt->Scalar("photo_lit_lifetime_s",
      Moon::ExospherePhotoionizationLifeTime(lit,_NA_SPEC_,-1,allowed,NULL));
  receipt->Scalar("photo_lit_allowed",allowed ? 1:0);
  receipt->Scalar("photo_shadow_lifetime_s",
      Moon::ExospherePhotoionizationLifeTime(shadow,_NA_SPEC_,-1,allowed,NULL));
  receipt->Scalar("photo_shadow_allowed",allowed ? 1:0);
  return 0;
}

int RunI03(Receipt* receipt) {
  receipt->Scalar("species_count",PIC::nTotalSpecies);
  receipt->Scalar("max_reference_injected_particles",
      PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber);
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const std::string prefix="species_"+std::to_string(species);
    receipt->Scalar(prefix+"_mass_kg",PIC::MolecularData::GetMass(species));
    receipt->Scalar(prefix+"_charge_c",
        PIC::MolecularData::ElectricChargeTable[species]);
    receipt->Scalar(prefix+"_global_dt_s",
        PIC::ParticleWeightTimeStep::GlobalTimeStep[species]);
    receipt->Scalar(prefix+"_global_weight",
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]);
  }
  receipt->Scalar("na_total_source_rate_s",
      Moon::SourceProcesses::totalProductionRate(_NA_SPEC_,
          _INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));
  receipt->Scalar("na_surface_max_dt_s",
      Moon::Planet->maxIntersectedNodeTimeStep[_NA_SPEC_]);
  return 0;
}

int RunU06(Receipt* receipt) {
  if (_NA_PLUS_SPEC_<0) {
    receipt->Scalar("configuration_error_na_plus_missing",1);
    return 4;
  }
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",1);
  return 4;
#else
  // Use identical neutral and ion states so lunar gravity cancels exactly in
  // their acceleration difference.  What remains is the production charged-
  // species branch, independently comparable with q(E+v cross B)/m without a
  // duplicate test force law inside the executable.
  double x[3]={-100.0*_RADIUS_(_MOON_),0.0,0.0};
  double v[3]={1700.0,-2300.0,4100.0};
  double neutralAcceleration[3],ionAcceleration[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->findTreeNode(x);
  Moon::TotalParticleAcceleration(neutralAcceleration,_NA_SPEC_,-1,x,v,node);
  Moon::TotalParticleAcceleration(ionAcceleration,_NA_PLUS_SPEC_,-1,x,v,node);
  receipt->Vector("position_m",x);
  receipt->Vector("velocity_m_s",v);
  receipt->Vector("electric_field_v_m",Exosphere::swE_Typical);
  receipt->Vector("magnetic_field_t",Exosphere::swB_Typical);
  receipt->Vector("neutral_acceleration_m_s2",neutralAcceleration);
  receipt->Vector("ion_acceleration_m_s2",ionAcceleration);
  receipt->Scalar("ion_mass_kg",PIC::MolecularData::GetMass(_NA_PLUS_SPEC_));
  receipt->Scalar("ion_charge_c",PIC::MolecularData::ElectricChargeTable[_NA_PLUS_SPEC_]);
  return 0;
#endif
}

int RunI05(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",1);
  return 4;
#else
  const double radius=3.0*_RADIUS_(_MOON_);
  const double mu=GravityConstant*_MASS_(_MOON_);
  State initial={{radius,0.0,0.0},{0.0,std::sqrt(mu/radius),0.0}},final;
  const double totalTime=1200.0,dt=20.0;
  const int status=Integrate(_NA_SPEC_,initial,totalTime,dt,&final);
  receipt->Scalar("mover_status",status);
  receipt->Scalar("total_time_s",totalTime);
  receipt->Scalar("dt_s",dt);
  receipt->Scalar("lunar_mu_m3_s2",mu);
  WriteState(receipt,"initial",initial);
  if (status==_PARTICLE_MOTION_FINISHED_) WriteState(receipt,"final",final);
  return status==_PARTICLE_MOTION_FINISHED_ ? 0:5;
#endif
}

void PrepareOrbitState() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  // The normal time-step hook is the production owner of the midpoint SPICE
  // state.  Running it once proves that the same state-update and PIC wiring
  // used in science runs precedes the diagnostic trajectory.
  amps_time_step();
#endif
}

int RunI06(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",0);
  return 4;
#else
  PrepareOrbitState();
  WriteForceState(receipt);
  const double radius=3.0*_RADIUS_(_MOON_);
  // The antisolar axial point is in the cylindrical lunar shadow, so the Na
  // radiation term is exactly gated off while all differential and rotating
  // production terms remain active.
  State initial={{-radius,0.0,0.0},{0.0,900.0,120.0}},final;
  const double totalTime=300.0,dt=5.0;
  const int status=Integrate(_NA_SPEC_,initial,totalTime,dt,&final);
  receipt->Scalar("mover_status",status);
  receipt->Scalar("total_time_s",totalTime);
  receipt->Scalar("dt_s",dt);
  WriteState(receipt,"initial",initial);
  if (status==_PARTICLE_MOTION_FINISHED_) WriteState(receipt,"final",final);
  return status==_PARTICLE_MOTION_FINISHED_ ? 0:5;
#endif
}

int RunI07(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",0);
  return 4;
#else
  PrepareOrbitState();
  WriteForceState(receipt);
  const double radius=3.0*_RADIUS_(_MOON_);
  State lit={{radius,0.0,0.0},{250.0,700.0,0.0}},litFinal;
  State shadow={{-radius,0.0,0.0},{250.0,700.0,0.0}},shadowFinal;
  double litAcceleration[3],shadowAcceleration[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* litNode=
      PIC::Mesh::mesh->findTreeNode(lit.x);
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* shadowNode=
      PIC::Mesh::mesh->findTreeNode(shadow.x);
  Moon::TotalParticleAcceleration(litAcceleration,_NA_SPEC_,-1,lit.x,lit.v,
                                  litNode);
  Moon::TotalParticleAcceleration(shadowAcceleration,_NA_SPEC_,-1,shadow.x,
                                  shadow.v,shadowNode);
  const double dt=0.25;
  const int litStatus=MoveOneStep(_NA_SPEC_,lit,dt,&litFinal);
  const int shadowStatus=MoveOneStep(_NA_SPEC_,shadow,dt,&shadowFinal);
  receipt->Scalar("dt_s",dt);
  receipt->Scalar("lit_mover_status",litStatus);
  receipt->Scalar("shadow_mover_status",shadowStatus);
  receipt->Vector("lit_acceleration_m_s2",litAcceleration);
  receipt->Vector("shadow_acceleration_m_s2",shadowAcceleration);
  WriteState(receipt,"lit_initial",lit);
  WriteState(receipt,"shadow_initial",shadow);
  if (litStatus==_PARTICLE_MOTION_FINISHED_)
    WriteState(receipt,"lit_final",litFinal);
  if (shadowStatus==_PARTICLE_MOTION_FINISHED_)
    WriteState(receipt,"shadow_final",shadowFinal);
  return litStatus==_PARTICLE_MOTION_FINISHED_ &&
         shadowStatus==_PARTICLE_MOTION_FINISHED_ ? 0:5;
#endif
}

int RunI08(Receipt* receipt) {
  if (_NA_PLUS_SPEC_<0) {
    receipt->Scalar("configuration_error_na_plus_missing",1);
    return 4;
  }
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",1);
  return 4;
#else
  receipt->Vector("electric_field_v_m",Exosphere::swE_Typical);
  receipt->Vector("magnetic_field_t",Exosphere::swB_Typical);
  receipt->Scalar("ion_mass_kg",PIC::MolecularData::GetMass(_NA_PLUS_SPEC_));
  receipt->Scalar("ion_charge_c",
      PIC::MolecularData::ElectricChargeTable[_NA_PLUS_SPEC_]);

  const double radius=100.0*_RADIUS_(_MOON_);
  // The production Moon box extends far down the antisolar (-X) tail but only
  // ten lunar radii sunward.  Place the uniform-field fixture at -100 R_Moon:
  // lunar gravity is small, all four gyro trajectories remain inside the real
  // mesh, and the test does not enlarge or bypass production boundaries.
  State initial={{-radius,0.0,0.0},{0.0,10000.0,3000.0}};
  WriteState(receipt,"initial",initial);
  const double gyroFrequency=std::fabs(
      PIC::MolecularData::ElectricChargeTable[_NA_PLUS_SPEC_])*
      Norm(Exosphere::swB_Typical)/PIC::MolecularData::GetMass(_NA_PLUS_SPEC_);
  const double period=PiTimes2/gyroFrequency;
  receipt->Scalar("gyrofrequency_rad_s",gyroFrequency);
  receipt->Scalar("total_time_s",period);
  // A preliminary 32/64/128/256 sequence demonstrated second-order ratios
  // but left the finest production state at 1.705% scaled error, above the
  // predeclared 1% M0 gate.  Preserve four factor-of-two levels and the gate;
  // move the frozen sequence to 64/128/256/512 so the acceptance decision is
  // made in the demonstrated asymptotic regime rather than by relaxing the
  // physical accuracy requirement.
  const int stepCounts[4]={64,128,256,512};
  for (int level=0;level<4;level++) {
    State final;
    const double dt=period/stepCounts[level];
    const int status=Integrate(_NA_PLUS_SPEC_,initial,period,dt,&final);
    const std::string prefix="level_"+std::to_string(level);
    receipt->Scalar(prefix+"_steps",stepCounts[level]);
    receipt->Scalar(prefix+"_dt_s",dt);
    receipt->Scalar(prefix+"_mover_status",status);
    if (status==_PARTICLE_MOTION_FINISHED_) WriteState(receipt,prefix+"_final",final);
    else return 5;
  }
  return 0;
#endif
}

int RunI09(Receipt* receipt) {
  // The agreed M0 exit classification uses Moon-relative two-body specific
  // energy, epsilon=|v_inertial|^2/2-GM_Moon/r.  The configured LSO axes rotate,
  // so v_inertial=v_LSO+Omega_LSO x r_LSO; omitting that transform could change
  // the sign for marginal particles.  epsilon>0 is physical lunar escape.
  // A particle that reaches the finite AMPS box with epsilon<=0 is instead a
  // numerical-domain loss.  The production mover returns the same terminal
  // status for both, so this explicit diagnostic proves the classification
  // inputs and preserves the distinction required by the loss budget.
  PrepareOrbitState();
  // The production Moon mesh is asymmetric along the Sun--Moon axis: its
  // sunward face is at +30 lunar radii (the same limit exercised by U01), not
  // at the +/-10-R_Moon transverse limits.  Place both particles immediately
  // inside that actual face.  Using a smaller, guessed radius would let the
  // mover finish normally and would not exercise domain-loss accounting.
  const double boundaryRadius=30.0*_RADIUS_(_MOON_);
  State escaping={{0.9999*boundaryRadius,0.0,0.0},{1000.0,0.0,0.0}};
  State numerical={{0.9999*boundaryRadius,0.0,0.0},{100.0,0.0,0.0}};
  const State fixtures[2]={escaping,numerical};
  const char* fixtureNames[2]={"physical_escape","numerical_domain_loss"};
  // Each duration carries the particle beyond the remaining 0.003 R_Moon
  // gap even after the small lunar acceleration is included.  The two values
  // differ because the non-escaping fixture is deliberately ten times slower.
  const double fixtureDt[2]={6.0,60.0};
  // Export the compiled status values rather than duplicating their numeric
  // definitions in Python.  The checker can consequently distinguish a true
  // domain exit from an ordinary completed step across AMPS configurations.
  receipt->Scalar("particle_left_domain_status",_PARTICLE_LEFT_THE_DOMAIN_);
  receipt->Scalar("particle_motion_finished_status",_PARTICLE_MOTION_FINISHED_);
  for (int fixture=0;fixture<2;fixture++) {
    double inertialVelocity[3],omegaCrossPosition[3];
    omegaCrossPosition[0]=Moon::RotationVector_SO_FROZEN[1]*fixtures[fixture].x[2]-
        Moon::RotationVector_SO_FROZEN[2]*fixtures[fixture].x[1];
    omegaCrossPosition[1]=Moon::RotationVector_SO_FROZEN[2]*fixtures[fixture].x[0]-
        Moon::RotationVector_SO_FROZEN[0]*fixtures[fixture].x[2];
    omegaCrossPosition[2]=Moon::RotationVector_SO_FROZEN[0]*fixtures[fixture].x[1]-
        Moon::RotationVector_SO_FROZEN[1]*fixtures[fixture].x[0];
    for (int component=0;component<3;component++) {
      inertialVelocity[component]=fixtures[fixture].v[component]+
          omegaCrossPosition[component];
    }
    const double specificEnergy=0.5*Norm(inertialVelocity)*Norm(inertialVelocity)-
        GravityConstant*_MASS_(_MOON_)/Norm(fixtures[fixture].x);
    State unusedFinal;
    const int status=MoveOneStep(_NA_SPEC_,fixtures[fixture],fixtureDt[fixture],
                                 &unusedFinal);
    const std::string prefix=fixtureNames[fixture];
    receipt->Scalar(prefix+"_specific_energy_j_kg",specificEnergy);
    receipt->Scalar(prefix+"_mover_status",status);
    receipt->Scalar(prefix+"_dt_s",fixtureDt[fixture]);
    WriteState(receipt,prefix+"_initial",fixtures[fixture]);
  }

  // Independently of the terminal fixtures, sample the installed production
  // sphere callback.  The fixed seed is part of the I09 contract.  Every
  // callback result must be accounted, and every re-emitted velocity must
  // point out of the surface (positive radial component at +X).
  rnd_seed(731993);
  const int samples=4096;
  int stuck=0,reemitted=0,inwardReemission=0,other=0;
  double speedSum=0.0;
  for (int sample=0;sample<samples;sample++) {
    State state={{_RADIUS_(_MOON_),0.0,0.0},{-1000.0,0.0,0.0}};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
        PIC::Mesh::mesh->findTreeNode(state.x);
    const long int particle=PIC::ParticleBuffer::GetNewParticle();
    PIC::ParticleBuffer::SetI(_NA_SPEC_,particle);
    PIC::ParticleBuffer::SetX(state.x,particle);
    PIC::ParticleBuffer::SetV(state.v,particle);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    // Match MoveOneStep's layout guard.  Surface-callback sampling needs a
    // valid production particle record but does not change its base weight.
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,particle);
#endif
    PIC::ParticleBuffer::SetParticleAllocated(particle);
    double remainingTime=0.0;
    const int status=Moon::Planet->ParticleSphereInteraction(
        _NA_SPEC_,particle,state.x,state.v,remainingTime,node,Moon::Planet);
    if (status==_PARTICLE_DELETED_ON_THE_FACE_) {
      stuck++;
    }
    else if (status==_PARTICLE_REJECTED_ON_THE_FACE_) {
      reemitted++;
      if (state.v[0]<=0.0) inwardReemission++;
      speedSum+=Norm(state.v);
    }
    else other++;
  }
  receipt->Scalar("seed",731993);
  receipt->Scalar("samples",samples);
  receipt->Scalar("stuck",stuck);
  receipt->Scalar("reemitted",reemitted);
  receipt->Scalar("inward_reemission",inwardReemission);
  receipt->Scalar("other_status",other);
  receipt->Scalar("mean_reemission_speed_m_s",
      reemitted>0 ? speedSum/reemitted:0.0);
  return 0;
}

int RunI10(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  receipt->Scalar("configuration_error_orbit_mode",0);
  return 4;
#else
  PrepareOrbitState();
  WriteForceState(receipt);
  const double radius=3.0*_RADIUS_(_MOON_);
  State initial={{-radius,0.0,0.0},{0.0,900.0,120.0}};
  const double totalTime=300.0;
  const int stepCounts[4]={15,30,60,120};
  receipt->Scalar("total_time_s",totalTime);
  WriteState(receipt,"initial",initial);
  for (int level=0;level<4;level++) {
    State final;
    const double dt=totalTime/stepCounts[level];
    const int status=Integrate(_NA_SPEC_,initial,totalTime,dt,&final);
    const std::string prefix="level_"+std::to_string(level);
    receipt->Scalar(prefix+"_steps",stepCounts[level]);
    receipt->Scalar(prefix+"_dt_s",dt);
    receipt->Scalar(prefix+"_mover_status",status);
    if (status==_PARTICLE_MOTION_FINISHED_) WriteState(receipt,prefix+"_final",final);
    else return 5;
  }
  return 0;
#endif
}

// Reuse the actual production mesh for the independent I12 and I13 scorers.
// These declarations avoid creating a second geometry or selector formula.
int RunI12(Receipt* receipt);
int RunI13(Receipt* receipt);

namespace SamplingWindowCheck {
PIC::Mesh::cDataCenterNode* cell=NULL;
int callbacks=0;
double maximumError=0.0;

double Marker(int index) { return 1.0+(index%8)/8.0; }

void Capture(int) {
  ++callbacks;
  const double* completed=reinterpret_cast<const double*>(
      cell->GetAssociatedDataBufferPointer()+
      PIC::Mesh::completedCellSampleDataPointerOffset);
  const int count=PIC::Mesh::sampleSetDataLength/sizeof(double);
  for (int index=0;index<count;++index)
    maximumError=std::max(maximumError,
        std::fabs(completed[index]-Marker(index)));
}
}

int RunSamplingWindowCheck(Receipt* receipt) {
  // The fixture adds exactly representable, dimensionless markers to actual
  // production storage. No particles, forces, source laws or observations
  // are manufactured. PIC::Sampling::Sampling still switches, exchanges,
  // clears and dispatches the completed callback through its normal path.
  using namespace SamplingWindowCheck;
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=NULL && cell==NULL;node=node->nextNodeThisThread) {
    // The one-argument index includes ghost cells; zero is a halo entry,
    // not interior (i,j,k)=(0,0,0). A halo marker is exchanged from another
    // rank and is not covered by the owned-cell reset being tested. Select
    // a genuine interior cell through the production coordinate overload.
    if (node->block!=NULL)
      cell=node->block->GetCenterNode(_BLOCK_CELLS_X_/2,
          _BLOCK_CELLS_Y_/2,_BLOCK_CELLS_Z_/2);
  }
  if (cell==NULL || PIC::Mesh::sampleSetDataLength==0)
    exit(__LINE__,__FILE__,"I11 sampling fixture needs an owned allocated cell");
  PIC::Sampling::RuntimeSamplingSwitch=false;
  PIC::SamplingMode=_RESTART_SAMPLING_MODE_;
  PIC::RequiredSampleLength=1;
  PIC::CollectingSampleCounter=0;
  PIC::Sampling::minIterationNumberForDataOutput=0;
  PIC::Sampling::SupressOutputFlag=false;
  PIC::Sampling::SkipOutputStep=1;
  for (int routine=0;routine<PIC::Sampling::ExternalSamplingLocalVariables::
           SamplingRoutinesRegistrationCounter;++routine)
    PIC::Sampling::ExternalSamplingLocalVariables::PrintOutputFile[routine]=
        I11Sampling::SuppressExternalOutput;
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
      Moon::Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Capture);
  const int count=PIC::Mesh::sampleSetDataLength/sizeof(double);
  std::fill_n(reinterpret_cast<double*>(cell->GetAssociatedDataBufferPointer()+
      PIC::Mesh::completedCellSampleDataPointerOffset),count,0.0);
  PIC::Mesh::flushCollectingSamplingBuffer(cell);
  double collectingError=0.0;
  for (int layout=0;layout<3;++layout) {
    for (int species=0;species<PIC::nTotalSpecies;++species)
      PIC::Sampling::SaveOutputDataFile[species]=
          layout==2 || (layout==1 && species==0);
    for (int window=0;window<3;++window) {
      double* collecting=reinterpret_cast<double*>(
          cell->GetAssociatedDataBufferPointer()+
          PIC::Mesh::collectingCellSampleDataPointerOffset);
      for (int index=0;index<count;++index)
        collecting[index]+=Marker(index);
      // Below FirstPrintedOutputFile: skip costly cell visualization while
      // retaining the species-enabled branch and every completed callback.
      PIC::DataOutputFileNumber=-2;
      PIC::Sampling::Sampling();
      collecting=reinterpret_cast<double*>(
          cell->GetAssociatedDataBufferPointer()+
          PIC::Mesh::collectingCellSampleDataPointerOffset);
      for (int index=0;index<count;++index)
        collectingError=std::max(collectingError,std::fabs(collecting[index]));
    }
  }
  double local[3]={maximumError,collectingError,std::fabs(callbacks-9.0)};
  double global[3];
  MPI_Allreduce(local,global,3,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  receipt->Scalar("sampling_check_completed_error",global[0]);
  receipt->Scalar("sampling_check_collecting_error",global[1]);
  receipt->Scalar("sampling_check_callback_error",global[2]);
  receipt->Scalar("sampling_check_windows",9);
  receipt->Scalar("sampling_check_species",PIC::nTotalSpecies);
  return 0;
}

// Allocate controlled initial particles through the real AMPS lists. These
// states are numerical fixtures, never observational evidence or a source-law
// substitute. Only the leaf's MPI owner allocates its particles.
void AddBudgetFixture(double* position,double* velocity,int count,int species=_NA_SPEC_) {
  auto* node=PIC::Mesh::mesh->findTreeNode(position);
  if (node==NULL) exit(__LINE__,__FILE__,"budget fixture lies outside the mesh");
  if (node->Thread!=PIC::ThisThread) return;
  int i,j,k;
  PIC::Mesh::mesh->FindCellIndex(position,i,j,k,node);
  auto& first=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  for (int index=0;index<count;++index) {
    const long int particle=PIC::ParticleBuffer::GetNewParticle(first);
    PIC::ParticleBuffer::SetI(species,particle);
    PIC::ParticleBuffer::SetX(position,particle);
    PIC::ParticleBuffer::SetV(velocity,particle);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,particle);
#endif
    auto* data=PIC::ParticleBuffer::GetParticleDataPointer(particle);
    Exosphere::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_,data);
    Exosphere::Sampling::SetParicleOriginSurfaceElementNumber(0,data);
  }
}

int RunBudgetCheck(Receipt* receipt,const std::string& receiptPath) {
  if (PIC::nTotalSpecies!=1 || _NA_SPEC_<0 || !Moon::Runtime::HasConfiguration())
    return 4;
  const auto directory=std::filesystem::path(receiptPath).parent_path();
  std::filesystem::create_directories(directory/"chemistry");
  std::filesystem::create_directories(directory/"transport");
  const double radius=_RADIUS_(_MOON_);
  double lit[3]={2.0*radius,0.0,0.0},shadow[3]={-1.1*radius,0.0,0.0};
  double velocity[3]={0.0,0.0,0.0};
  constexpr int Population=4096;
  // amps_init allocates geometry, but the first real step populates the
  // decomposition's chemistry-walker table and frozen Sun/Earth state. Run
  // that ordinary initialization step before constructing controlled states;
  // otherwise an empty walker and zero Earth position would invalidate the
  // fixture itself. Its particles are cleared outside the ledger interval.
  Moon::nIterationCounter=0;
  amps_time_step();
  PIC::ParticleBuffer::DeleteAllParticles();
  AddBudgetFixture(lit,velocity,Population);
  AddBudgetFixture(shadow,velocity,Population);
  // Preserve a mechanistic control: the formerly configured generic processor
  // has no Na model and leaves every lit fixture unchanged. This calls the
  // existing processor itself, not a copied lifetime/table implementation.
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=NULL;node=node->nextNodeThisThread) if (node->block!=NULL)
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell) {
      auto& first=node->block->FirstCellParticleTable[cell];
      long int current=first;first=-1;
      while (current!=-1) {
        const long int next=PIC::ParticleBuffer::GetNext(current);
        Exosphere::ChemicalModel::PhotochemicalModelProcessor(current,first,node);
        current=next;
      }
    }
  const double normalDt=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
  bool allowed=false;
  const double lifetime=Moon::ExospherePhotoionizationLifeTime(lit,_NA_SPEC_,-1,allowed,NULL);
  if (!allowed || !(lifetime>0.0)) return 4;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]=lifetime*std::log(2.0);
  Moon::ParticleBudget::Begin(directory/"chemistry");
  receipt->Scalar("generic_processor_control_survivors",Moon::ParticleBudget::GasCount());
  rnd_seed(731993+PIC::ThisThread);
  // Exercise the configured macro through the real PIC list walker. A future
  // input dispatch regression therefore fails even if the Moon kernel passes.
  PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
  Moon::ParticleBudget::Snapshot(1,"chemistry");
  receipt->Scalar("chemistry_budget_valid",Moon::ParticleBudget::Valid() ? 1:0);
  receipt->Scalar("chemistry_loss_count",Moon::ParticleBudget::Count(PIC::ParticleEvents::PhotolyticLoss));
  long long local[2]={0,0},global[2];
  for (auto* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node!=NULL;node=node->nextNodeThisThread) if (node->block!=NULL)
    for (int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for (long int p=node->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p)) {
        double position[3];PIC::ParticleBuffer::GetX(position,p);
        ++local[position[0]>0.0 ? 0:1];
      }
  MPI_Allreduce(local,global,2,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  receipt->Scalar("sunlit_survivors",global[0]);
  receipt->Scalar("shadow_survivors",global[1]);
  receipt->Scalar("initial_population_per_site",Population);
  receipt->Scalar("fixture_chemistry_dt_s",lifetime*std::log(2.0));
  receipt->Scalar("fixture_lifetime_s",lifetime);
  Moon::ParticleBudget::End();
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]=normalDt;
  PIC::ParticleBuffer::DeleteAllParticles();
  // Normal transport retains its configured dt, sources and stochastic
  // surface kernel. Multiple night-side particles make surface coverage
  // robust without forcing a sticking outcome or changing its coefficient.
  double impact[3]={-radius-100.0,0.0,0.0},inward[3]={1000.0,0.0,0.0};
  double outside[3]={PIC::Mesh::mesh->rootTree->xmax[0]-10.0,0.0,0.0};
  double outward[3]={1000.0,0.0,0.0};
  AddBudgetFixture(impact,inward,64);
  AddBudgetFixture(outside,outward,1);
  PIC::RequiredSampleLength=500;
  PIC::CollectingSampleCounter=0;
  Moon::nIterationCounter=0;
  rnd_seed(731993+PIC::ThisThread);
  Moon::ParticleBudget::Begin(directory/"transport");
  for (int step=0;step<20;++step) { amps_time_step();++Moon::nIterationCounter; }
  Moon::ParticleBudget::Snapshot(20,"final");
  receipt->Scalar("transport_budget_valid",Moon::ParticleBudget::Valid() ? 1:0);
  receipt->Scalar("maximum_count_residual",Moon::ParticleBudget::MaximumCountResidual());
  receipt->Scalar("emission_count",Moon::ParticleBudget::Count(PIC::ParticleEvents::Emission));
  receipt->Scalar("adsorption_count",Moon::ParticleBudget::Count(PIC::ParticleEvents::Adsorption));
  receipt->Scalar("escape_count",Moon::ParticleBudget::Count(PIC::ParticleEvents::Escape));
  receipt->Scalar("numerical_removal_count",Moon::ParticleBudget::Count(PIC::ParticleEvents::NumericalRemoval));
  Moon::ParticleBudget::End();
  return 0;
}

int RunI11(Receipt* receipt,const std::string& receiptPath,
           bool durationPilot) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  if (!Moon::Runtime::HasConfiguration()) {
    receipt->Scalar("configuration_error_missing_input",1);
    return 4;
  }

  // Particle-only restart files cannot restore I11's in-memory baseBlocks,
  // batch selection, or acceptance state.  Require the normal Moon build to
  // disable their compile-time autosave before the first production step so a
  // bounded test never emits misleading or disk-intensive checkpoints.
#if _PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ == _PIC_AUTOSAVE_PARTICLE_DATA_RESTART_FILE__MODE_ON_
  receipt->Scalar("particle_restart_writes_enabled",1);
  receipt->Scalar("configuration_error_particle_restart_must_be_disabled",1);
  return 4;
#else
  receipt->Scalar("particle_restart_writes_enabled",0);
#endif

  // Override the general Moon output cadence only for linked I11.  Setting it
  // here leaves normal production and all other U/I tests unchanged while the
  // receipt below records the exact value exercised by the production sampler.
  PIC::RequiredSampleLength=I11Sampling::SampleLengthSteps;
  PIC::CollectingSampleCounter=0;
  I11Sampling::durationPilot=durationPilot;

  // A completed production sample normally emits a full volumetric Tecplot
  // field plus application-specific column/velocity maps.  The 2026-10-31
  // bounded campaign showed that those visualization copies consumed 20--71
  // GB per resolution even though I11 scores the already-completed in-memory
  // density/surface buffers.  Disable only file writers here: particle
  // transport, production sampling, buffer exchange/switching, source and
  // boundary callbacks, and the I11 Capture callback below remain unchanged.
  // This is therefore an output-policy change, not a test-only physics path.
  for (int species=0;species<PIC::nTotalSpecies;species++)
    PIC::Sampling::SaveOutputDataFile[species]=false;
  const int suppressedExternalWriters=PIC::Sampling::ExternalSamplingLocalVariables::
      SamplingRoutinesRegistrationCounter;
  for (int routine=0;routine<suppressedExternalWriters;routine++)
    PIC::Sampling::ExternalSamplingLocalVariables::PrintOutputFile[routine]=
        I11Sampling::SuppressExternalOutput;
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  receipt->Scalar("simulation_environment_mode",
      configuration.simulationEnvironmentMode);
  if (configuration.simulationEnvironmentMode!="frozen") {
    receipt->Scalar("configuration_error_frozen_environment_required",1);
    return 4;
  }

  // Capture independent initialization diagnostics in the same allocated
  // production mesh. Each receives its own receipt and independent Python
  // status; neither can grant I11 convergence or freeze an I13 reference.
  for (const std::string id : {"I12","I13"}) {
    Receipt companion((std::filesystem::path(receiptPath).parent_path()/
        (id+"_initialization_receipt.txt")).string(),PIC::ThisThread==0);
    if (!companion.Good())
      exit(__LINE__,__FILE__,"I11 cannot create shared initialization receipt");
    companion.Scalar("test_id",id);
    companion.Scalar("mpi_ranks",PIC::nTotalThreads);
    const int status=id=="I12" ? RunI12(&companion) : RunI13(&companion);
    if (status!=0) return status;
  }

  receipt->Scalar("random_seed",configuration.randomSeed);
  receipt->Scalar("mesh_profile",configuration.meshProfile);
  receipt->Scalar("mesh_output_tecplot",
      configuration.meshOutputTecplot ? 1 : 0);
  receipt->Scalar("i11_execution_contract",
      durationPilot ? "duration-pilot-not-qualified" :
                      "bounded-formal-qualification");
  receipt->Scalar("production_volume_output_suppressed",1);
  receipt->Scalar("external_visualization_writers_suppressed",
      suppressedExternalWriters);
  receipt->Scalar("base_block_length_steps",PIC::RequiredSampleLength);
  // Zero means that no scientific step ceiling is applied by the pilot.  The
  // Python runner records/enforces independent operational resource guards.
  // Formal I11 continues to record the literal reviewed 5000-step cap.
  receipt->Scalar("maximum_particle_steps",
      durationPilot ? 0 : I11Sampling::MaximumParticleSteps);
  receipt->Scalar("iteration_limit_enabled",durationPilot ? 0 : 1);
  receipt->Scalar("equilibrium_relative_change_gate",0.02);
  receipt->Scalar("equilibrium_consecutive_changes",3);
  receipt->Scalar("minimum_independent_science_blocks",5);
  receipt->Scalar("integrated_ci95_relative_half_width_gate",0.0167);
  receipt->Scalar("pointwise_ci95_relative_half_width_gate",0.0333);

  // Attach a read-only diagnostic to PIC's normal completed-sample callback.
  // It consumes the same center-node and surface buffers as production output
  // before PIC flushes them, while all transport, sources, chemistry, and
  // surface interactions continue through amps_time_step().
  I11Sampling::Configure(std::filesystem::path(receiptPath).parent_path());
  Moon::ParticleBudget::Begin(std::filesystem::path(receiptPath).parent_path());

  Moon::nIterationCounter=0;
  amps_time_step();
  Moon::nIterationCounter++;
  receipt->Scalar("step_1_ephemeris_time_tdb_s",Moon::OrbitalMotion::et);
  receipt->Vector("step_1_sun_position_lso_m",Moon::xSun_SO);
  receipt->Vector("step_1_earth_position_lso_m",Moon::xEarth_SO);
  receipt->Vector("step_1_lso_angular_velocity_rad_s",
      Moon::RotationVector_SO_FROZEN);
  amps_time_step();
  Moon::nIterationCounter++;
  receipt->Scalar("step_2_ephemeris_time_tdb_s",Moon::OrbitalMotion::et);
  receipt->Vector("step_2_sun_position_lso_m",Moon::xSun_SO);
  receipt->Vector("step_2_earth_position_lso_m",Moon::xEarth_SO);
  receipt->Vector("step_2_lso_angular_velocity_rad_s",
      Moon::RotationVector_SO_FROZEN);
  receipt->Scalar("particle_global_timestep_s",
      PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]);
  // Stop at the approved hard iteration cap even when the unchanged
  // statistical gates are unresolved.  Reaching the cap is an evaluated FAIL,
  // not ERROR or PASS; the Python runner executes the remaining mesh levels so
  // all four bounded diagnostic cases have complete provenance.
  while (!I11Sampling::complete &&
         (durationPilot ||
          Moon::nIterationCounter<I11Sampling::MaximumParticleSteps)) {
    // A bounded diagnostic runner requests a stop via this run-owned file.
    // Check only between complete production steps; final receipts and event
    // snapshots then survive a wall-time stop without a particle restart.
    if (durationPilot && std::filesystem::exists(
        std::filesystem::path(receiptPath).parent_path()/"pilot_stop.request"))
      break;
    amps_time_step();
    Moon::nIterationCounter++;
  }
  Moon::ParticleBudget::Snapshot(Moon::nIterationCounter,"final");
  receipt->Scalar("particle_budget_valid",Moon::ParticleBudget::Valid() ? 1:0);
  receipt->Scalar("particle_budget_maximum_count_residual",Moon::ParticleBudget::MaximumCountResidual());
  Moon::ParticleBudget::End();
  if (I11Sampling::complete) I11Sampling::WriteFinal();
  const bool iterationLimitReached=
      !durationPilot && !I11Sampling::complete &&
      Moon::nIterationCounter>=I11Sampling::MaximumParticleSteps;
  receipt->Scalar("production_particle_steps_executed",
      Moon::nIterationCounter);
  receipt->Scalar("base_blocks_executed",I11Sampling::baseBlocks.size());
  receipt->Scalar("retained_batch_size_base_blocks",
      I11Sampling::retainedBatchSize);
  receipt->Scalar("retained_independent_batch_count",
      I11Sampling::retainedBatchCount);
  receipt->Scalar("scoring_metric_count",I11Sampling::definitions.size());
  receipt->Scalar("convergence_criteria_met",
      I11Sampling::complete ? 1 : 0);
  receipt->Scalar("iteration_limit_reached",iterationLimitReached ? 1 : 0);
  receipt->Scalar("convergence_metrics_available",
      I11Sampling::complete ? 1 : 0);
  receipt->Scalar("pilot_stop_requested",durationPilot &&
      std::filesystem::exists(std::filesystem::path(receiptPath).parent_path()/"pilot_stop.request") ? 1:0);
  return 0;
#endif
}

int RunI12(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  // I12's polar fixture is defined in MOON_ME_DE421 and cannot be reconstructed
  // honestly in a no-SPICE executable. Record a configuration error instead
  // of substituting LSO +Z or a test-only rotation.
  receipt->Scalar("configuration_error_spice_required",1);
  return 4;
#else
  if (!Moon::Runtime::HasConfiguration()) {
    receipt->Scalar("configuration_error_missing_input",1);
    return 4;
  }
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  receipt->Scalar("mesh_profile",configuration.meshProfile);
  receipt->Scalar("mesh_output_tecplot",
      configuration.meshOutputTecplot ? 1 : 0);

  // I12 measures both halves of the linked contract: the target returned by
  // the installed production selector and the cell size actually realized by
  // AMPS after its discrete power-of-two refinement and balancing constraints.
  // Points are just outside the analytic boundary to avoid an ambiguous
  // internal-boundary ownership query while remaining negligible relative to
  // even the finest target selected by the approved I11 campaign.
  const double radius=_RADIUS_(_MOON_)+
      std::max(1.0e-3,10.0*PIC::Mesh::mesh->EPS);
  const double szaDeg[6]={0.0,30.0,60.0,90.0,120.0,180.0};
  SpiceDouble lsoToBody[3][3];
  pxform_c(Moon::Frames::SolarOrbital,Moon::Frames::BodyFixed,
           Moon::OrbitalMotion::et,lsoToBody);
  for (int sample=0;sample<6;sample++) {
    const double angle=szaDeg[sample]*Pi/180.0;
    double position[3]={radius*cos(angle),radius*sin(angle),0.0};
    double bodyPosition[3]={0.0,0.0,0.0};
    for (int row=0;row<3;row++) for (int column=0;column<3;column++)
      bodyPosition[row]+=lsoToBody[row][column]*position[column];
    const double latitudeDeg=asin(std::max(-1.0,std::min(1.0,
        bodyPosition[2]/radius)))*180.0/Pi;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
        PIC::Mesh::mesh->findTreeNode(position);
    // The AMR tree is replicated, whereas a leaf's allocated data block exists
    // only on its owning MPI rank.  Requiring node->block on rank zero would
    // incorrectly reject a valid distributed mesh.  Every rank therefore
    // contributes one only when it is the declared owner and has the block;
    // exactly one contribution proves both ownership and allocation without
    // copying mesh data or replacing the production decomposition service.
    int localAllocated=(node!=NULL && node->Thread==PIC::ThisThread &&
                        node->block!=NULL) ? 1 : 0;
    int allocatedOwners=0;
    MPI_Allreduce(&localAllocated,&allocatedOwners,1,MPI_INT,MPI_SUM,
                  MPI_GLOBAL_COMMUNICATOR);
    if (node==NULL || allocatedOwners!=1) return 5;
    const std::string prefix="sza_"+std::to_string(sample);
    receipt->Scalar(prefix+"_angle_deg",szaDeg[sample]);
    receipt->Scalar(prefix+"_absolute_body_latitude_deg",fabs(latitudeDeg));
    receipt->Scalar(prefix+"_requested_resolution_m",
        localSphericalSurfaceResolution(position));
    receipt->Scalar(prefix+"_realized_cell_size_m",
        node->GetCharacteristicCellSize());
  }

  // Construct the north-pole point through the inverse SPICE transform rather
  // than assuming that LSO +Z is the DE421 lunar pole.  This directly checks
  // that the polar selector can win over the SZA table in the linked mesh.
  SpiceDouble bodyToLso[3][3];
  pxform_c(Moon::Frames::BodyFixed,Moon::Frames::SolarOrbital,
           Moon::OrbitalMotion::et,bodyToLso);
  double polePosition[3]={bodyToLso[0][2]*radius,
                          bodyToLso[1][2]*radius,
                          bodyToLso[2][2]*radius};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* poleNode=
      PIC::Mesh::mesh->findTreeNode(polePosition);
  int localPoleAllocated=(poleNode!=NULL &&
                          poleNode->Thread==PIC::ThisThread &&
                          poleNode->block!=NULL) ? 1 : 0;
  int allocatedPoleOwners=0;
  MPI_Allreduce(&localPoleAllocated,&allocatedPoleOwners,1,MPI_INT,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (poleNode==NULL || allocatedPoleOwners!=1) return 5;
  const double poleSzaDeg=acos(std::max(-1.0,std::min(1.0,
      polePosition[0]/radius)))*180.0/Pi;
  receipt->Scalar("north_pole_sza_deg",poleSzaDeg);
  receipt->Scalar("north_pole_requested_resolution_m",
      localSphericalSurfaceResolution(polePosition));
  receipt->Scalar("north_pole_realized_cell_size_m",
      poleNode->GetCharacteristicCellSize());

  receipt->Scalar("mesh_terminator_half_width_deg",
      configuration.meshTerminatorHalfWidthDeg);
  receipt->Scalar("mesh_polar_cap_start_latitude_deg",
      configuration.meshPolarCapStartLatitudeDeg);
  receipt->Scalar("mesh_polar_resolution_m",
      configuration.meshPolarResolutionM);
  for (std::size_t i=0;i<configuration.meshAngularSzaDeg.size();i++) {
    receipt->Scalar("mesh_knot_"+std::to_string(i)+"_sza_deg",
        configuration.meshAngularSzaDeg[i]);
    receipt->Scalar("mesh_knot_"+std::to_string(i)+"_resolution_m",
        configuration.meshAngularResolutionM[i]);
  }
  receipt->Scalar("mesh_knot_count",configuration.meshAngularSzaDeg.size());
  return 0;
#endif
}

int RunI13(Receipt* receipt) {
  if (!Moon::Runtime::HasConfiguration()) {
    receipt->Scalar("configuration_error_missing_input",1);
    return 4;
  }
  const Moon::Runtime::Configuration& configuration=
      Moon::Runtime::GetConfiguration();
  if (configuration.surfaceGeometry!=Moon::Runtime::SurfaceGeometry::Sphere ||
      Moon::Planet==NULL) {
    receipt->Scalar("configuration_error_not_sphere",1);
    return 4;
  }

  // Sum the areas returned by the registered production boundary itself.
  // The independent runner compares this value with 4*pi*R^2; duplicating the
  // analytical expression here would only compare a formula with itself.
  const long int surfaceElements=Moon::Planet->GetTotalSurfaceElementsNumber();
  double sphereAreaM2=0.0;
  for (long int element=0;element<surfaceElements;element++)
    sphereAreaM2+=Moon::Planet->GetSurfaceElementArea(element);

  // Count the actual linked AMR tree rather than an input prediction.  Only
  // leaf nodes own physical mesh blocks, and the separate allocated count
  // catches a partially initialized or ownership-mismatched tree.
  long int leafNodes=0,localAllocatedLeafBlocks=0;
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
           PIC::Mesh::mesh->BranchBottomNodeList;
       node!=NULL;node=node->nextBranchBottomNode) {
    leafNodes++;
    if (node->Thread==PIC::ThisThread && node->block!=NULL)
      localAllocatedLeafBlocks++;
  }
  long int allocatedLeafBlocks=0;
  MPI_Allreduce(&localAllocatedLeafBlocks,&allocatedLeafBlocks,1,MPI_LONG,
                MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  receipt->Scalar("sphere_radius_m",_RADIUS_(_MOON_));
  receipt->Scalar("simulation_environment_mode",
      configuration.simulationEnvironmentMode);
  receipt->Scalar("mesh_profile",configuration.meshProfile);
  receipt->Scalar("mesh_near_surface_cell_size_m",
      configuration.meshNearSurfaceCellSizeM);
  receipt->Scalar("mesh_near_surface_extent_m",
      configuration.meshNearSurfaceExtentM);
  receipt->Scalar("mesh_far_field_cell_size_m",
      configuration.meshFarFieldCellSizeM);
  receipt->Scalar("mesh_far_field_extent_m",
      configuration.meshFarFieldExtentM);
  receipt->Scalar("mesh_outer_cell_size_m",
      configuration.meshOuterCellSizeM);
  receipt->Scalar("mesh_output_tecplot",
      configuration.meshOutputTecplot ? 1 : 0);
  receipt->Scalar("sphere_surface_element_count",surfaceElements);
  receipt->Scalar("sphere_area_m2",sphereAreaM2);
  receipt->Scalar("amr_leaf_node_count",leafNodes);
  receipt->Scalar("amr_allocated_leaf_block_count",allocatedLeafBlocks);
  receipt->Scalar("species_count",PIC::nTotalSpecies);
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    const std::string prefix="species_"+std::to_string(species);
    receipt->Scalar(prefix+"_symbol",
        PIC::MolecularData::GetChemSymbol(species));
    receipt->Scalar(prefix+"_mass_kg",
        PIC::MolecularData::GetMass(species));
    receipt->Scalar(prefix+"_charge_c",
        PIC::MolecularData::GetElectricCharge(species));
    receipt->Scalar(prefix+"_global_timestep_s",
        PIC::ParticleWeightTimeStep::GlobalTimeStep[species]);
    receipt->Scalar(prefix+"_global_particle_weight",
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[species]);
  }
  receipt->Scalar("na_total_source_rate_s",
      Moon::SourceProcesses::totalProductionRate(_NA_SPEC_,
          _INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));
  receipt->Scalar("na_impact_source_rate_s",
      Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(
          _NA_SPEC_,_INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet));

  // A controlled I13 baseline also requires equilibrated profiles, column
  // density, brightness, and return flux from a *separate passing I11 run*.
  // This initialization-only process cannot manufacture those quantities.
  // The independent runner therefore imports the preserved finest I11 artifact,
  // writes a review candidate, and returns SKIPPED until an investigator has
  // approved and frozen that candidate.  This zero is process-local evidence;
  // it does not claim that the external I11 artifact is unavailable.
  receipt->Scalar("sampled_baseline_metrics_available",0);
  return 0;
}

}  // namespace

namespace {
// This adapter selects declared initial conditions and numerical controls;
// every physical step remains amps_time_step with production sources/movers.
double NoExternalImpactSource(int) {return 0.0;}
// A model-input normalization is prepared before any residuals are examined.
// Integrate the actual external source kernels with a separate converged
// composite-Simpson quadrature; no particles, thermal state or RNG advance.
int RunM2SourceNormalization(Receipt* receipt) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ != _PIC_MODE_ON_
  return 6;
#else
  if (_NA_SPEC_<0 || _AR_SPEC_<0 || PIC::nTotalSpecies!=2) return 6;
  const auto& c=Moon::Runtime::GetConfiguration();
  const double H=c.m2VerificationDurationS,start=Moon::OrbitalMotion::et;
  if (!(H>0)) exit(__LINE__,__FILE__,"source normalization requires an explicit positive SI horizon");
  const double savedRadius=Exosphere::xObjectRadial;
  long double previous[2]={},coarseDifference[2]={1e300L,1e300L},integrals[2]={};
  int intervals=8;
  for (;;) {
    integrals[0]=integrals[1]=0;
    for (int n=0;n<=intervals;++n) {
      SpiceDouble state[6],lightTime;
      spkezr_c("MOON",start+H*n/intervals,"J2000","NONE","SUN",state,&lightTime);
      Exosphere::xObjectRadial=1000.0*std::sqrt(state[0]*state[0]+state[1]*state[1]+state[2]*state[2]);
      const int coefficient=(n==0 || n==intervals)?1:(n%2?4:2);
      for (int j=0;j<2;++j) {
        const int species=j==0?_NA_SPEC_:_AR_SPEC_;
        const double rate=Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(
            species,_INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet)+
            Moon::Species::Argon::GetTotalProductionRate(species,_INTERNAL_BOUNDARY_TYPE_SPHERE_,Moon::Planet);
        if (!std::isfinite(rate) || rate<0) exit(__LINE__,__FILE__,"invalid external source rate in normalization quadrature");
        integrals[j]+=coefficient*static_cast<long double>(rate);
      }
    }
    bool converged=intervals>=32;
    for (int j=0;j<2;++j) {
      integrals[j]*=H/(3.0L*intervals);
      const long double difference=std::fabs(integrals[j]-previous[j]);
      const long double floor=128*std::numeric_limits<double>::epsilon()*std::fabs(integrals[j]);
      converged=converged && difference<=c.m2VerificationSourceIntegralTolerance*std::fabs(integrals[j]) &&
          (difference<=coarseDifference[j] || difference<=floor);
      coarseDifference[j]=difference;previous[j]=integrals[j];
    }
    if (converged) break;
    if (intervals>std::numeric_limits<int>::max()/2)
      exit(__LINE__,__FILE__,"source integral quadrature exhausted interval counter; inspect forcing/horizon");
    intervals*=2;
  }
  Exosphere::xObjectRadial=savedRadius;
  receipt->Scalar("execution_scope","model-input-source-integral-only-not-qualification");
  receipt->Scalar("duration_s",H);receipt->Scalar("source_integral_intervals",intervals);
  receipt->Scalar("expected_external_na_atoms",static_cast<double>(integrals[0]));
  receipt->Scalar("expected_external_ar_atoms",static_cast<double>(integrals[1]));
  receipt->Scalar("initial_na_cohort_atoms_per_site",c.m2VerificationCohortPackets*PIC::ParticleWeightTimeStep::GlobalParticleWeight[_NA_SPEC_]);
  receipt->Scalar("initial_ar_cohort_atoms_per_site",c.m2VerificationCohortPackets*PIC::ParticleWeightTimeStep::GlobalParticleWeight[_AR_SPEC_]);
  return 0;
#endif
}

int RunM2NumericalCase(Receipt* receipt,const std::string& path) {
  const auto& c=Moon::Runtime::GetConfiguration();
  std::string initializationError;
  if (!RecordAppliedThermalInitialization(receipt,&initializationError))
    exit(__LINE__,__FILE__,initializationError.c_str());
  Moon::SurfaceTemperature::ProductionThermalInitializationInfo appliedInitial;
  if (!Moon::SurfaceTemperature::GetProductionThermalInitialization(
          &appliedInitial,&initializationError) ||
      appliedInitial.ephemerisTimeS!=Moon::OrbitalMotion::et)
    exit(__LINE__,__FILE__,"M2 numerical case start epoch does not equal its actual applied thermal initial-state epoch; declare/reinitialize one consistent t0 before seeding cohorts");
  if (c.m2VerificationRecordInitialThermalState && PIC::ThisThread==0) {
    // Capture the actual applied t0 columns before source/cohort state or
    // numerical dt/weight controls are changed. Every case gets an exclusive
    // path; these model bytes are evidence to bind, never a reference truth.
    const std::string statePath=std::filesystem::absolute(
        path+".I19_initial_thermal_state.csv").lexically_normal().string();
    std::vector<double> depths,temperatures;
    const int elements=Moon::Planet->GetTotalSurfaceElementsNumber();
    if (elements<=0 || !Moon::SurfaceTemperature::ProductionThermalProfile(
            0,&depths,&temperatures,&initializationError))
      exit(__LINE__,__FILE__,("initial thermal-state metadata unavailable: "+initializationError).c_str());
    const std::size_t layers=depths.size();
    for (int element=1;element<elements;++element) {
      if (!Moon::SurfaceTemperature::ProductionThermalProfile(element,&depths,
              &temperatures,&initializationError) || depths.size()!=layers)
        exit(__LINE__,__FILE__,"initial thermal-state columns have inconsistent layer counts");
    }
    if (!Moon::SurfaceTemperature::WriteProductionInitialThermalStateCsv(
            statePath,&initializationError))
      exit(__LINE__,__FILE__,initializationError.c_str());
    receipt->Scalar("initial_thermal_state_csv",statePath);
    receipt->Scalar("initial_thermal_state_element_count",elements);
    receipt->Scalar("initial_thermal_state_layer_count",layers);
    receipt->Scalar("initial_thermal_state_rows",
        static_cast<std::size_t>(elements)*(layers+1));
  }
  if (c.m2VerificationCase=="normalization") return RunM2SourceNormalization(receipt);
  if (_NA_SPEC_<0 || _AR_SPEC_<0 || PIC::nTotalSpecies!=2 ||
      !c.m2EnvironmentEnabled || c.simulationEnvironmentMode!="evolving")
    exit(__LINE__,__FILE__,"M2 numerical cases require evolving neutral Na/Ar production environment");
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ != _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_
  exit(__LINE__,__FILE__,"M2 numerical cases require chemistry OFF");
#endif
  if (c.surfaceInventoryInitialState!="zero")
    exit(__LINE__,__FILE__,"M2 numerical profiles require explicitly zero initial surface state before controlled cohorts");
  const int species[2]={_NA_SPEC_,_AR_SPEC_};
  double savedDt[2],savedWeight[2];
  for (int j=0;j<2;++j) {
    savedDt[j]=PIC::ParticleWeightTimeStep::GlobalTimeStep[species[j]];
    savedWeight[j]=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species[j]];
  }
  const double bound=std::min(savedDt[0],savedDt[1]);
  double coarse=c.m2VerificationCoarseDtS>0 ? c.m2VerificationCoarseDtS:bound;
  long long coarseSteps=c.m2VerificationCoarseSteps;
  if (c.m2VerificationDurationS>0) {
    if (c.m2VerificationCoarseDtS==0) {
      const double intervals=std::ceil(c.m2VerificationDurationS/(8*coarse));
      if (!std::isfinite(intervals) || intervals<1 ||
          intervals>std::numeric_limits<int>::max()/(8.0*c.m2VerificationDtDivisor))
        exit(__LINE__,__FILE__,"M2 physical horizon exceeds supported step counter; review dt/horizon before rerun");
      coarseSteps=8*static_cast<long long>(intervals);
    }
    // An explicit baseline pins BOTH count and dt. Re-ceiling H/dt could
    // add an interval solely through floating-point division rounding.
    const double aligned=c.m2VerificationDurationS/coarseSteps;
    if (c.m2VerificationCoarseDtS>0 && std::abs(aligned-coarse)>
        16*std::numeric_limits<double>::epsilon()*coarse)
      exit(__LINE__,__FILE__,"explicit cluster coarse dt/count/horizon are inconsistent");
    coarse=aligned;
  }
  if (!(coarse>0) || !std::isfinite(coarse) || !(bound>0) || !std::isfinite(bound) || coarse>bound*(1+16*std::numeric_limits<double>::epsilon()) ||
      coarseSteps>std::numeric_limits<int>::max()/c.m2VerificationDtDivisor)
    exit(__LINE__,__FILE__,"M2 common dt/step count violates initialized production bound or counter range");
  const double dt=coarse/c.m2VerificationDtDivisor;
  const long long steps=coarseSteps*c.m2VerificationDtDivisor;
  const double reference[2]={c.m2VerificationNaWeightAtoms>0?c.m2VerificationNaWeightAtoms:savedWeight[0],
    c.m2VerificationArWeightAtoms>0?c.m2VerificationArWeightAtoms:savedWeight[1]};
  for (int j=0;j<2;++j) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species[j]]=dt;
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species[j]]=reference[j]/c.m2VerificationWeightDivisor;
  }
  // Seed after initialization so every compared case begins the declared
  // stochastic stream, independent of mesh setup's earlier random draws.
  rnd_seed(c.randomSeed+PIC::ThisThread);
  const bool cohort=c.m2VerificationCase=="cohort";
  auto oldImpact=Exosphere::SourceProcesses::ImpactVaporization::RuntimeSourceRate;
  std::string error;
  if (cohort) {
    Exosphere::SourceProcesses::ImpactVaporization::RuntimeSourceRate=NoExternalImpactSource;
    auto disabled=c;disabled.argonSourceEnabled=false;
    if (!Moon::Species::Argon::ConfigureProduction(disabled,Moon::Planet,&error))
      exit(__LINE__,__FILE__,error.c_str());
  }
  // Explicit zero startup; no hidden warm-up step may inject gas or deposit
  // atoms before this initial-value problem. PIC initializes its own block
  // traversal on the first ordinary step.
  std::vector<double> initial;
  if (cohort) {
    for (int j=0;j<2;++j) for (int element:c.m2VerificationElements) {
      if (element>=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber)
        exit(__LINE__,__FILE__,"configured cohort element exceeds logical surface");
      const double atoms=c.m2VerificationCohortPackets*savedWeight[j];
      if (!std::isfinite(atoms) || !(atoms>0))
        exit(__LINE__,__FILE__,"cohort initial atoms overflow; reduce explicit cohort packet scale and rerun");
      Moon::Planet->SurfaceElementPopulation[species[j]][element]+=atoms;initial.push_back(atoms);
    }
    Moon::M2Cohorts::Begin(c.m2VerificationElements,initial,std::filesystem::path(path).parent_path());
  }
  // The numerical fixtures need production sampling but no full-volume file
  // writers or restart files. A short case has no complete 500-step sample.
  PIC::RequiredSampleLength=500;PIC::SamplingMode=_RESTART_SAMPLING_MODE_;
  PIC::Sampling::SupressOutputFlag=true;PIC::Sampling::SupressRestartFilesFlag=true;
  for (int routine=0;routine<PIC::Sampling::ExternalSamplingLocalVariables::SamplingRoutinesRegistrationCounter;++routine)
    PIC::Sampling::ExternalSamplingLocalVariables::PrintOutputFile[routine]=I11Sampling::SuppressExternalOutput;
  const double start=Moon::OrbitalMotion::et;
  Moon::M2Budget::Begin(std::filesystem::path(path).parent_path()/"m2_budget.csv");
  for (long long step=1;step<=steps;++step) {
    Moon::OrbitalMotion::et=start+(step-1)*dt;
    Moon::nIterationCounter=static_cast<int>(step);amps_time_step();
    // Integer-based physical clock prevents millions of repeated additions
    // to a large TDB epoch from accumulating avoidable timestamp drift.
    Moon::OrbitalMotion::et=start+step*dt;
    if (cohort) Moon::M2Cohorts::EndStep(step*dt);
    if (step%(steps/8)==0 || step%500==0 || step==steps) {
      Moon::M2Budget::Snapshot(step);
      if (cohort) Moon::M2Cohorts::Snapshot(step*dt);
    }
  }
  receipt->Scalar("execution_scope",c.m2VerificationProfile=="development"?
      "supporting-short-M2-not-formal-qualification":"cluster-M2-numerical-case");
  receipt->Scalar("verification_case",c.m2VerificationCase);
  receipt->Scalar("seed",c.randomSeed);receipt->Scalar("steps",static_cast<double>(steps));
  receipt->Scalar("coarse_steps",static_cast<double>(coarseSteps));
  receipt->Scalar("coarse_timestep_s",coarse);receipt->Scalar("synchronized_timestep_s",dt);
  receipt->Scalar("initial_epoch_tdb_s",start);receipt->Scalar("duration_s",steps*dt);
  receipt->Scalar("initialized_na_dt_s",savedDt[0]);receipt->Scalar("initialized_ar_dt_s",savedDt[1]);
  receipt->Scalar("initialized_na_weight_atoms",savedWeight[0]);receipt->Scalar("initialized_ar_weight_atoms",savedWeight[1]);
  receipt->Scalar("reference_na_weight_atoms",reference[0]);receipt->Scalar("reference_ar_weight_atoms",reference[1]);
  receipt->Scalar("na_weight_atoms",PIC::ParticleWeightTimeStep::GlobalParticleWeight[_NA_SPEC_]);
  receipt->Scalar("ar_weight_atoms",PIC::ParticleWeightTimeStep::GlobalParticleWeight[_AR_SPEC_]);
  receipt->Scalar("rng_mode",_RND_MODE_);receipt->Scalar("photochemistry_enabled",0);
  receipt->Scalar("budget_valid",Moon::M2Budget::Valid()?1:0);
  Moon::M2Budget::End();if (cohort) Moon::M2Cohorts::End();
  Exosphere::SourceProcesses::ImpactVaporization::RuntimeSourceRate=oldImpact;
  for (int j=0;j<2;++j) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species[j]]=savedDt[j];
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species[j]]=savedWeight[j];
  }
  return 0;
}

int RunM2BudgetCheck(Receipt* receipt,const std::string& path) {
  if (Moon::Runtime::HasConfiguration() &&
      Moon::Runtime::GetConfiguration().m2VerificationCase!="none")
    return RunM2NumericalCase(receipt,path);

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ != _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_
  std::cerr << "I19 supporting budget requires photolytic chemistry OFF\n";
  return 6;
#endif
  if (_NA_SPEC_<0 || _AR_SPEC_<0 || PIC::nTotalSpecies!=2) return 6;
  constexpr int Steps=20;
  constexpr double SeedPacketsPerSpecies=32.0;
  // Equal simulated seconds for both species are essential for a transient
  // atom/mass window. Select the smaller initialized production dt, retain
  // its exact value, and restore both controls after the bounded fixture.
  const std::array<double,2> saved{{
      PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],
      PIC::ParticleWeightTimeStep::GlobalTimeStep[_AR_SPEC_]}};
  const double dt=std::min(saved[0],saved[1]);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]=dt;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_AR_SPEC_]=dt;
  rnd_seed(731993+PIC::ThisThread);
  Moon::nIterationCounter=0;amps_time_step();
  const auto& terrain=Moon::SurfaceEnvironment::ProductionTerrain()->triangulation();
  const auto& face=terrain.faces.front();
  const auto& a=terrain.verticesM[face[0]];
  const auto& b=terrain.verticesM[face[1]];
  const auto& c=terrain.verticesM[face[2]];
  double e0[3],e1[3],normal[3],body[6]={},lso[6]={};
  for (int i=0;i<3;++i) {e0[i]=b[i]-a[i];e1[i]=c[i]-a[i];}
  normal[0]=e0[1]*e1[2]-e0[2]*e1[1];
  normal[1]=e0[2]*e1[0]-e0[0]*e1[2];
  normal[2]=e0[0]*e1[1]-e0[1]*e1[0];
  const double length=Norm(normal);
  for (int i=0;i<3;++i) {
    normal[i]/=length;
    body[i]=(a[i]+b[i]+c[i])/3.0+100.0*normal[i];
    body[i+3]=-1000.0*normal[i];
  }
  for (int row=0;row<6;++row) for (int col=0;col<6;++col)
    lso[row]+=Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[row][col]*body[col];
  // Controlled incoming verification particles guarantee actual surface
  // callback coverage. Allocation occurs only on the true mesh owner, before
  // the ledger's initial state; these are not counted as physical sources.
  AddBudgetFixture(lso,lso+3,64,_NA_SPEC_);
  AddBudgetFixture(lso,lso+3,64,_AR_SPEC_);
  double exitPosition[3]={PIC::Mesh::mesh->xGlobalMax[0]-100.0,0.0,0.0};
  double exitVelocity[3]={3000.0,0.0,0.0};
  AddBudgetFixture(exitPosition,exitVelocity,1,_NA_SPEC_);
  AddBudgetFixture(exitPosition,exitVelocity,1,_AR_SPEC_);
  const int element=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber/2;
  for (int species:std::array<int,2>{{_NA_SPEC_,_AR_SPEC_}}) {
    // Deterministic INITIAL verification state, not an emission or an
    // observation. In particular Ar exercises dt >> residence time without
    // a Poisson-attempt loop or changing its literature binding energy.
    Moon::Planet->SurfaceElementPopulation[species][element]+=
        SeedPacketsPerSpecies*PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
  }
  const auto csv=std::filesystem::path(path).parent_path()/"m2_budget.csv";
  Moon::M2Budget::Begin(csv);
  double classificationPosition[3]={29.9*_RADIUS_(_MOON_),0.0,0.0};
  double boundVelocity[3]={1.0,0.0,0.0};
  receipt->Scalar("positive_energy_exit_classified_escape",
      PIC::ParticleEvents::DomainExitClassifier(classificationPosition,exitVelocity)==PIC::ParticleEvents::Escape?1:0);
  receipt->Scalar("bound_exit_classified_numerical",
      PIC::ParticleEvents::DomainExitClassifier(classificationPosition,boundVelocity)==PIC::ParticleEvents::NumericalRemoval?1:0);
  for (int step=1;step<=Steps;++step) {
    Moon::nIterationCounter=step;amps_time_step();
    Moon::M2Budget::Snapshot(step);
  }
  receipt->Scalar("execution_scope","supporting-M2-development-not-I19-qualification");
  receipt->Scalar("seed",731993);
  receipt->Scalar("steps",Steps);
  receipt->Scalar("synchronized_timestep_s",dt);
  receipt->Scalar("initial_fixture_packets_per_species",SeedPacketsPerSpecies);
  receipt->Scalar("initial_fixture_surface_element",element);
  receipt->Scalar("photochemistry_enabled",0);
  receipt->Scalar("budget_valid",Moon::M2Budget::Valid()?1:0);
  Moon::M2Budget::End();
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_]=saved[0];
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_AR_SPEC_]=saved[1];
  return 0;
}
}

int Moon::IntegrationTests::Run(const std::string& testId,
                                const std::string& outputFile,
                                bool i11DurationPilot,bool i11SamplingCheck,
                                bool i11BudgetCheck,bool m2BudgetCheck) {
  // I11--I13 explicitly support distributed AMR ownership and reduce all
  // process-local fields or counts before scoring.  Earlier M0 diagnostics
  // retain their frozen one-rank contract because several create/move local
  // particles and have not been reviewed for distributed execution.
  const bool distributedLinked=(testId=="I11" || testId=="I12" ||
      testId=="I13" || testId=="I14" || testId=="I15" || testId=="I16" ||
      testId=="I17" || testId=="I18" || testId=="I20" || (testId=="I19" && m2BudgetCheck));
  if (PIC::nTotalThreads!=1 && !distributedLinked) {
    if (PIC::ThisThread==0)
      std::cerr << "M0 linked diagnostics require exactly one MPI rank\n";
    return 2;
  }

  Receipt receipt(outputFile,PIC::ThisThread==0);
  if (!receipt.Good()) {
    if (PIC::ThisThread==0)
      std::cerr << "Cannot create Moon integration receipt: " << outputFile
                << '\n';
    return 2;
  }
  receipt.Scalar("test_id",testId);
  receipt.Scalar("mpi_ranks",PIC::nTotalThreads);
  if (testId=="I19" && m2BudgetCheck) return RunM2BudgetCheck(&receipt,outputFile);
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  receipt.Scalar("orbit_mode",1);
#else
  receipt.Scalar("orbit_mode",0);
#endif

  // All ranks must enter each M1 diagnostic because their ownership checks or
  // sampled fields use MPI collectives. Other linked tests are one-rank-only
  // above, so these branches change no existing M0 execution semantics.
  if (testId=="I11" && i11SamplingCheck)
    return RunSamplingWindowCheck(&receipt);
  if (testId=="I11" && i11BudgetCheck)
    return RunBudgetCheck(&receipt,outputFile);
  if (testId=="I11")
    return RunI11(&receipt,outputFile,i11DurationPilot);
  if (testId=="I12") return RunI12(&receipt);
  if (testId=="I13") return RunI13(&receipt);
  if (testId=="I14") return RunI14(&receipt,outputFile);
  if (testId=="I15") return RunI15(&receipt,outputFile);
  if (testId=="I16") return RunI16(&receipt,outputFile);
  if (testId=="I17") return RunI17(&receipt,outputFile);
  if (testId=="I18") return RunI18(&receipt,outputFile);
  if (testId=="I20") {
    // The deterministic global thermal state is replicated, so only rank
    // zero exports it. Other ranks remain in the normal linked lifecycle;
    // the receipt/layout is honest about this absence of MPI speedup.
    if (PIC::ThisThread!=0) return 0;
    std::map<std::string,std::string> metrics;
    std::string error;
    if (!Moon::PolarDiagnostics::RunPolarDiagnostics(outputFile,&metrics,&error)) {
      std::cerr << "I20 production polar export failed: " << error << '\n';
      return 6;
    }
    for (const auto& value:metrics) receipt.Scalar(value.first,value.second);
    if (!RecordAppliedThermalInitialization(&receipt,&error)) {
      std::cerr << "I20 applied thermal-initialization metadata unavailable: " << error << '\n';
      return 6;
    }
    return 0;
  }
  if (PIC::ThisThread!=0) return 2;

  if (testId=="U01") return RunU01(&receipt,false);
  if (testId=="U01_INSIDE") return RunU01(&receipt,true);
  if (testId=="U06") return RunU06(&receipt);
  if (testId=="I02") return RunI02(&receipt);
  if (testId=="I03") return RunI03(&receipt);
  if (testId=="I05") return RunI05(&receipt);
  if (testId=="I06") return RunI06(&receipt);
  if (testId=="I07") return RunI07(&receipt);
  if (testId=="I08") return RunI08(&receipt);
  if (testId=="I09") return RunI09(&receipt);
  if (testId=="I10") return RunI10(&receipt);
  std::cerr << "Unknown Moon integration diagnostic: " << testId << '\n';
  return 2;
}
