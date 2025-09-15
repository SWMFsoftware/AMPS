

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
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "constants.h"
#include "sep.h"

#include "tests.h"

void amps_init();
void amps_init_mesh();
void amps_time_step();

enum class CMEScenario { Fast, Slow };


//the physical simulation time 
double SimulationTime=0.0;


// Configure the single model for the requested scenario (called only on change)
static inline void configure_swcme1d(CMEScenario scenario){
  static bool inited = false;

  if (inited==true) return;

  inited=true;

  if (scenario == CMEScenario::Fast){
	  SEP::sw1d.SetAmbient(400.0, 6.0, 5.0, 1.2e5)                 // Vsw[km/s], n1AU[cm^-3], B1AU[nT], T[K]
      .SetCME(1.05, 1900.0, 8e-8)                         // r0[Rs],   V0_sh[km/s],  Γ[1/km]
      .SetGeometry(0.12, 0.22)                            // sheath & ME thickness @1 AU [AU]
      .SetSmoothing(0.010, 0.020, 0.030)                  // shock/LE/TE widths @1 AU [AU]
      .SetSheathEjecta(1.25, 2.0, 1.12, 0.50, 0.80);      // rc_floor, ramp_p, Vshe_LE, fME, VME
  } else {
	  SEP::sw1d.SetAmbient(380.0, 5.0, 4.5, 1.0e5)
      .SetCME(1.05, 950.0, 3e-8)
      .SetGeometry(0.08, 0.18)
      .SetSmoothing(0.015, 0.030, 0.050)
      .SetSheathEjecta(1.15, 1.5, 1.08, 0.60, 0.90);
  }

  SEP::SW1DAdapter::EnableSheathClamp(true);      // optional stability aid

}

/**
 * Advance one global step:
 *  - builds the per-time StepState for the selected scenario
 *  - publishes it to the mover through the adapter
 *  - uses a single static Model 'sw' (no multiple instances)
 */
void advance_sw1d(double dt) { 
  static double t_since_launch=0.0;
  static int ncall=0;

  ncall++;

  // Build per-time cache and publish to mover via adapter
  swcme1d::StepState S = SEP::sw1d.prepare_step(t_since_launch);
  SEP::SW1DAdapter::SetModelAndState(&SEP::sw1d, S);

  t_since_launch+=dt;



  // One call does everything: n, V, Br, Bphi, |B|, ∇·V → Tecplot POINT file
  if ((PIC::ThisThread==0)&&(ncall%10==0)) {
    char fname[200];

    const int N = 400;
    static double r[N];
    const double rmin = 1.05*_SUN__RADIUS_, rmax = 2.00*_AU_;     // 0.2–2 AU
    for (int i=0;i<N;++i){
      double t = double(i)/(N-1);
      r[i] = rmin*std::pow(rmax/rmin, t);            // log-spacing (nice for r^-2)
    }

    sprintf(fname,"sw_profile_%i.dat",ncall); 
    SEP::sw1d.write_tecplot_radial_profile_from_r(S, r, N, fname,SimulationTime);
  }
}

int main(int argc,char **argv) {
  //      MPI_Init(&argc,&argv);


  SEP::ShockModelType=SEP::cShockModelType::SwCme1d;

  //read post-compile input file  
  if (PIC::PostCompileInputFileName!="") {
     SEP::Parser::ReadFile(PIC::PostCompileInputFileName);
  }


  //set up shock wave model 
  configure_swcme1d(CMEScenario::Fast); 

  //output parameters of the sshock 
  SEP::sw1d.write_tecplot_shock_vs_time(2.0*24.0*3600, 200, "shock_vs_time.dat");

  //setup datum to store the segment's data for the Alfven turbulence model 
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) { 
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy); 
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);

    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S_pm);


    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array);

  }

  //set up datum to store distance of a field line vertex to the location of the shock 
  PIC::FieldLine::UserDefinedfDataProcessingManager=SEP::FieldLine::CalculateVertexShockDistances; 
 
  amps_init_mesh();
  amps_init();

  //init the Alfven turbulence IC
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) SEP::AlfvenTurbulence_Kolmogorov::ModelInit::Init(); 

  if (_PIC_FIELD_LINE_MODE_==_PIC_MODE_ON_) { 
    TestManager();
  }

  int TotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) ? PIC::RequiredSampleLength+10 : 100000001;  

  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient;

  //init turbulence wave energy 
  double B0_1AU = 5.0e-9;        // 5 nT magnetic field
  double turbulence_level = 0.2; // 20% turbulence

  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);

  SEP::AlfvenTurbulence_Kolmogorov::InitializeWaveEnergyFromPhysicalParameters(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy, B0_1AU, turbulence_level, true);


  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);

  //exchenge the initial wave energy density and output in a file
  PIC::FieldLine::Parallel::MPIGatherDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,0); 

  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);

  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,1);


  PIC::FieldLine::Parallel::MPIBcastDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,0);

  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
  SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,1);


  //set background plasma density 
auto set_background_plasma_density = []() {
    // reference density at 1 AU [m⁻³]
    constexpr double n0 = 5.0e6;
    // 1 astronomical unit in meters
    constexpr double AU = 1.495978707e11;

    // Loop over all field lines
    for (int fldIdx = 0; fldIdx < PIC::FieldLine::nFieldLine; ++fldIdx) {
        auto fieldLine = &PIC::FieldLine::FieldLinesAll[fldIdx];
        if (!fieldLine) continue;

        int nSeg = fieldLine->GetTotalSegmentNumber();
        if (nSeg < 1) continue;

        // Loop over all segments in this field line
        for (int segIdx = 0; segIdx < nSeg; ++segIdx) {
            auto seg = fieldLine->GetSegment(segIdx);
            if (!seg) continue;

            // Process both end‐points (left & right) of the segment
            PIC::FieldLine::cFieldLineVertex* vertices[2] = { seg->GetBegin(), seg->GetEnd() };
            for (auto vtx : vertices) {
                if (!vtx) continue;

                // get pointer to {x,y,z} [m]
                double* X = vtx->GetX();
                // compute radial distance from Sun [m]
                double r = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);

                // scale density as n0/(r/AU)^2
                double density = n0 / ((r/AU) * (r/AU));

                // store it in the vertex
                vtx->SetPlasmaDensity(density); //SetDatum(density, PIC::FieldLine::DatumAtVertexPlasmaDensity);
            }
        }
    }
};


 set_background_plasma_density();

  // Calculate turbulence wave enregy density from wave energy integrated over the segments of the magnetic tube:
auto CalculateWaveEnergyDensity = [&]() {
    // Focus specifically on SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy
    auto& integrated_energy = SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy;
    auto& energy_density = SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity;

    int local_count = 0;

    for (int fl = 0; fl < PIC::FieldLine::nFieldLine; fl++) {
        auto& field_line = PIC::FieldLine::FieldLinesAll[fl];
        int num_segs = field_line.GetTotalSegmentNumber();

        for (int seg = 0; seg < num_segs; seg++) {
            auto* segment = field_line.GetSegment(seg);
            if (!segment) continue;

            // Only process segments belonging to this thread
            if (segment->Thread == PIC::ThisThread) {
                // Focus on getting CellIntegratedWaveEnergy values
                double* energy_data = segment->GetDatum_ptr(integrated_energy);
                double* density_data = segment->GetDatum_ptr(energy_density);

                // Use SEP::FieldLine::GetSegmentVolume for volume calculation
                double volume = SEP::FieldLine::GetSegmentVolume(segment, fl);

                if (energy_data && density_data && volume > 0.0) {
                    // Process all elements using SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy.length
                    for (int i = 0; i < SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy.length; i++) {
                        density_data[i] = energy_data[i] / volume;

                        if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                          validate_numeric(density_data[i],__LINE__,__FILE__);
			}
                    }
                    local_count++;
                }
            }
        }
    }

    // MPI operations to gather and broadcast results
    PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(energy_density);

    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);
    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,2); 

    if (PIC::ThisThread == 0) {
        std::cout << "Compact wave energy density calculation completed with MPI operations" << std::endl;
    }
};


 //calculate the wave energy density
 CalculateWaveEnergyDensity();

    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);
    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,2);


    //PIC::FieldLine::Output("fl-edge-test.dat",false); 

//hooks for calculating the magnertic tube radius and the volume of the field line segment
PIC::FieldLine::SegmentVolume=SEP::FieldLine::GetSegmentVolume;



  vector<vector<double> > DeltaE_plus, DeltaE_minus;



  //time step
  for (long int niter=0;niter<TotalIterations;niter++) {
    //SEP::InitDriftVelData();
    
    static double rsh0=SEP::ParticleSource::ShockWave::Tenishev2005::rShock; 
    
    if (niter==0) {
      switch (SEP::ShockModelType) {
      case SEP::cShockModelType::Analytic1D:
        rsh0=SEP::ParticleSource::ShockWave::Tenishev2005::rShock; 
        break;
      case SEP::cShockModelType::SwCme1d:
        rsh0=SEP::SW1DAdapter::gState.r_sh_m;
        break;
      }
    }


    amps_time_step();
    SimulationTime+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0]; 

    //advance the shock+CME model
   advance_sw1d(PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);


    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {

      // Function to increment integrated wave energy due to shock passing
      if (niter!=0) {
         double rsh1;

         switch (SEP::ShockModelType) {
         case SEP::cShockModelType::Analytic1D:
           rsh1=SEP::ParticleSource::ShockWave::Tenishev2005::rShock;
           break;
         case SEP::cShockModelType::SwCme1d:
           rsh1=SEP::SW1DAdapter::gState.r_sh_m;
           break;
         }

	 SEP::ParticleSource::ShockWave::ShockTurbulenceEnergyInjection(rsh0, rsh1, PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);

	 rsh0=rsh1;
      }

      if (SEP::ParticleMoverPtr!=SEP::ParticleMover_FocusedTransport_WaveScattering) { // in case SEP::ParticleMover_FocusedTransport_WaveScattering(), particle/turbulence coupling is already done 

      // Function to increment integrated wave energy due to shock passing
      //reduce S
      PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
      PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);

      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatum(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,0,"g+",0);
      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatum(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,0,"g-",0);

      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"g+",0);
      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,"g-",0);


      SEP::AlfvenTurbulence_Kolmogorov::AnalyzeMaxSegmentParticles(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"G_plus_streaming" ,0);
      SEP::AlfvenTurbulence_Kolmogorov::AnalyzeMaxSegmentParticles(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,"G_minus_streaming" ,0);


//      SEP::AlfvenTurbulence_Kolmogorov::SetDatumAll(0.0,SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
//      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"g+",0);

      //couple particles and turbulence  
//      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::UpdateAllSegmentsWaveEnergyWithParticleCoupling(
//		     SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
//		    SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S,
//		   PIC::ParticleWeightTimeStep::GlobalTimeStep[0]); 

      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::WaveParticleCouplingManager(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
		      PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);
      }


      //advect turbulence energy 
      SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(DeltaE_plus, DeltaE_minus,SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
			      PIC::ParticleWeightTimeStep::GlobalTimeStep[0],0.01,0.01);

      //model the effect of wave reflection 
      if (SEP::AlfvenTurbulence_Kolmogorov::Reflection::active==true) {
        double C_reflection=0.6;

        SEP::AlfvenTurbulence_Kolmogorov::Reflection::ReflectTurbulenceEnergyAllFieldLines(PIC::ParticleWeightTimeStep::GlobalTimeStep[0],C_reflection,0.0,false); 
      }

      // Configure cascade
      if (SEP::AlfvenTurbulence_Kolmogorov::Cascade::active==true) { 
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetCascadeCoefficient(0.8);             // C_nl
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetDefaultPerpendicularCorrelationLength(1.0e7); // 10,000 km
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetDefaultEffectiveArea(1.0);           // V_cell = Δs
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetElectronHeatingFraction(0.3);        // 30% to electrons
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableCrossHelicityModulation(false);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableTwoSweepIMEX(false);
  
        // Advance cascade for all field lines (ΔE arrays accumulate changes)
        // Optional: stronger physics
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableCrossHelicityModulation(true);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableTwoSweepIMEX(true);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::CascadeTurbulenceEnergyAllFieldLines(PIC::ParticleWeightTimeStep::GlobalTimeStep[0],/*enable_logging=*/ false);
      }

    
      //scatter wave energy   
      PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);


      //calculate the wave energy density 
      CalculateWaveEnergyDensity();

/*
      if ((niter+1)%2==0)  {   
         char fname[300];
	 sprintf(fname,"fl-%i-%e.dat",PIC::ThisThread,rsh0/_AU_);

	 //PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(&PIC::FieldLine::DatumAtVertexParticleWeight);
         //PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(&PIC::FieldLine::DatumAtVertexParticleCosPitchAngle);

	 PIC::FieldLine::Output(fname,false);

	 //PIC::FieldLine::Parallel::SetDatumStoredAtVertex(0.0,&PIC::FieldLine::DatumAtVertexParticleWeight);
	 //PIC::FieldLine::Parallel::SetDatumStoredAtVertex(0.0,&PIC::FieldLine::DatumAtVertexParticleCosPitchAngle);
      }
*/

    }


    //PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(10,200);
    //
    //PIC::ParticleSplitting::FledLine::WeightedParticleMerging(20,20,20,500,800); 
    //PIC::ParticleSplitting::FledLine::WeightedParticleSplitting(20,20,20,500,800);    
  }


  char fname[400];

  sprintf(fname,"%s/test_SEP.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  MPI_Finalize();
  return EXIT_SUCCESS;
}
