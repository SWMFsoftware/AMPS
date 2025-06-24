

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


int main(int argc,char **argv) {
  //      MPI_Init(&argc,&argv);


  //read post-compile input file  
  if (PIC::PostCompileInputFileName!="") {
     SEP::Parser::ReadFile(PIC::PostCompileInputFileName);
  }
  

  //setup datum to store the segment's data for the Alfven turbulence model 
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) { 
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy); 
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);

    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S_pm);
  }
	  
 
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


        if (PIC::ThisThread==2) PIC::FieldLine::Output("fl-edge-test.dat",false); 


  vector<vector<double> > DeltaE_plus, DeltaE_minus;

  //time step
  for (long int niter=0;niter<TotalIterations;niter++) {
    //SEP::InitDriftVelData();
    
    double rsh0=SEP::ParticleSource::ShockWave::Tenishev2005::rShock; 
    
    
    amps_time_step();

    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {

      // Function to increment integrated wave energy due to shock passing
      if (niter!=0) SEP::ParticleSource::ShockWave::ShockTurbulenceEnergyInjection(rsh0, SEP::ParticleSource::ShockWave::Tenishev2005::rShock, PIC::ParticleWeightTimeStep::GlobalTimeStep[0]); 


      //reduce S
      PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);

      //couple particles and turbulence  
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::UpdateAllSegmentsWaveEnergyWithParticleCoupling(
		     SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
		    SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S,
		   PIC::ParticleWeightTimeStep::GlobalTimeStep[0]); 


      //advect turbulence energy 
      SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(DeltaE_plus, DeltaE_minus,SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
			      PIC::ParticleWeightTimeStep::GlobalTimeStep[0],0.01,0.01);
    
      //scatter wave energy   
      PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);


      //calculate the wave energy density 
      CalculateWaveEnergyDensity();




      if ((PIC::ThisThread==2)&&((niter+1)%10==0))  {   
         char fname[300];
	 sprintf(fname,"fl-%e.dat",rsh0/_AU_);
	      PIC::FieldLine::Output(fname,false);
      }
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
