//the model of the Alfven turbulence. 

#include "sep.h"

bool SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag=true;

PIC::Datum::cDatumStored SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy(2,"\"W+\",\"W-\"",DATUM_DO_PRINT|DATUM_DO_DEVIDE_VOLUME_PRINT);
PIC::Datum::cDatumStored SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyGrowthRate(2,"\"dW+/dt\",\"dW-/dt\"",true);

PIC::Datum::cDatumStored SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S(SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::n_stream_intervals,"",false); 
PIC::Datum::cDatumStored SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S_pm(2,"\"S+\",\"S-\"",false);

double SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_min=log(Relativistic::Energy2Momentum(0.8*SEP::FieldLine::InjectionParameters::emin*MeV2J,_H__MASS_)); 
double SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_max=log(Relativistic::Energy2Momentum(1.2*SEP::FieldLine::InjectionParameters::emax*MeV2J,_H__MASS_));

double SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_dp_stream=
  (SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_max-SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::log_p_stream_min)/
       SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::n_stream_intervals;	

double SEP::AlfvenTurbulence_Kolmogorov::ModelInit::dB_B(double r) {
   // Example: power-law dependence on heliocentric distance
   // dB/B = dB_B_0 * (r/r_0)^alpha
   // where r is in meters, r_0 is 1 AU
        
    const double dB_B_0 = 0.1;        // Turbulence level at 1 AU
    const double r_0 = 1.496e11;      // 1 AU in meters
    const double alpha = -0.5;        // Power law index
        
    if (r <= 0.0) return 0.0;
        
     return dB_B_0 * pow(r / r_0, alpha);
}

//=========================================================================
// Main initialization function
void SEP::AlfvenTurbulence_Kolmogorov::ModelInit::Init() {      
  namespace FL = PIC::FieldLine;

  // Check if field lines are available
  if (FL::FieldLinesAll == nullptr || FL::nFieldLine == 0) {
    if (PIC::ThisThread == 0) {
      printf("SEP::AlfvenTurbulence_Kolmogorov::ModelInit: No field lines available for initialization\n");
    }

    return;
  }

  // Ensure the wave energy density datum is activated
  if (!SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy.is_active()) {
    if (PIC::ThisThread == 0) {
      printf("SEP::AlfvenTurbulence_Kolmogorov::ModelInit: Warning - CellIntegratedWaveEnergy datum not activated\n");
    }

     return;
  }

  // Permeability of free space
  const double mu_0 = 4.0 * M_PI * 1.0e-7; // H/m

  int totalSegmentsProcessed = 0;

  // Loop through all field lines
  for (int iFieldLine = 0; iFieldLine < FL::nFieldLine; iFieldLine++) {
    FL::cFieldLine& fieldLine = FL::FieldLinesAll[iFieldLine];

    // Skip uninitialized field lines
    if (!fieldLine.IsInitialized()) continue;

    // Loop through all segments in this field line
    FL::cFieldLineSegment* segment = fieldLine.GetFirstSegment();

    while (segment != nullptr) {
      if (segment->Thread!=PIC::ThisThread) {
        totalSegmentsProcessed++;

        // Move to next segment
        segment = segment->GetNext();
	continue;
      }

      // Get the vertices at the beginning and end of the segment
      FL::cFieldLineVertex* beginVertex = segment->GetBegin();
      FL::cFieldLineVertex* endVertex = segment->GetEnd();

      if (beginVertex == nullptr || endVertex == nullptr) {
        segment = segment->GetNext();
        continue;
      }

       // Get magnetic field values at both vertices
       double B_begin[3], B_end[3];
       beginVertex->GetMagneticField(B_begin);
       endVertex->GetMagneticField(B_end);

       // Calculate average magnetic field at segment midpoint
       double B_avg[3];
       double B_magnitude = 0.0;

       for (int i = 0; i < 3; i++) {
         B_avg[i] = 0.5 * (B_begin[i] + B_end[i]);
         B_magnitude += B_avg[i] * B_avg[i];
       }

       B_magnitude = sqrt(B_magnitude);

       // Get positions of vertices to calculate midpoint
       double x_begin[3], x_end[3];
       beginVertex->GetX(x_begin);
       endVertex->GetX(x_end);

       // Calculate midpoint position
       double x_mid[3];

       for (int i = 0; i < 3; i++) {
         x_mid[i] = 0.5 * (x_begin[i] + x_end[i]);
       }

       // Calculate heliocentric distance (assuming Sun at origin)
       double r = sqrt(x_mid[0]*x_mid[0] + x_mid[1]*x_mid[1] + x_mid[2]*x_mid[2]);

       // Calculate relative turbulence level dB/B
       double dB_over_B = dB_B(r);

       // Calculate wave energy density: (dB/B * B)^2 / (2 * mu_0)
       // This represents the magnetic energy density of turbulent fluctuations
       double dB = dB_over_B * B_magnitude;
       double W[2] = {0.0,0.0};

       if (B_magnitude > 0.0) {
         W[0] = (dB * dB) / (2.0 * mu_0);
         W[1]=W[0];
       }

       // Store the wave energy density in the segment
       segment->SetDatum(CellIntegratedWaveEnergy,W);

       totalSegmentsProcessed++;

       // Move to next segment
       segment = segment->GetNext();
     }

     // Report completion
     if (PIC::ThisThread == 0) {
       printf("SEP::AlfvenTurbulence_Kolmogorov::ModelInit: Initialized turbulence for %d segments across %d field lines\n",
       totalSegmentsProcessed, FL::nFieldLine);
     }
   }
}

//=========================================================================
//sampling particle particle streaming  
void SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::SampleParticleData(double s_final,double s_init,double speed,long int ptr,double dt,PIC::FieldLine::cFieldLineSegment *segment_start, int iFieldLine) {
  int pBin;

  int spec=PIC::ParticleBuffer::GetI(ptr);
  double log_p=log(Relativistic::Speed2Momentum(speed,PIC::MolecularData::GetMass(spec)));
  double e=Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec))*J2MeV;

  if ((log_p>log_p_stream_max)||(log_p<log_p_stream_min)) return;  
  pBin=static_cast<int>(std::floor((log_p-log_p_stream_min)/log_dp_stream));

  if ((pBin<0)||(pBin>=n_stream_intervals)) exit(__LINE__,__FILE__,"Error: out of range"); 

  double weight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]; 
  weight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);


  auto seg = segment_start;
  double s=s_init;

  if (!seg) return;

  // Branch A: forward motion  (s_final>s_init > 0)
  if (s_final - s > 0.0) {
    while (s_final > s && seg) {
      // Calculate fractional position within current segment
      double xi = s - std::floor(s);
        
      // Get Parker stream data array
      double* ParkerStream = seg->GetDatum_ptr(S);
      if (ParkerStream==NULL) return;
        
      // Get segment geometry
      auto v0 = seg->GetBegin();
      auto v1 = seg->GetEnd();
      double L, vol, dFlux;
        
      if (std::floor(s_final) != std::floor(s)) {
        // Cross segment boundary - travel from xi to end of segment
        L = (1.0 - xi) * seg->GetLength();
            
        // Use FULL segment volume for Parker flux normalization
        vol = SEP::FieldLine::GetSegmentVolume(seg, iFieldLine);
            
        dFlux = weight * L / (dt * vol);
        ParkerStream[pBin] += dFlux;
            
        // Move to next segment
        s += 1.0;
        s = std::floor(s);  // Ensure we're at exact segment boundary
        seg = seg->GetNext();
	if (seg==NULL) return;
      }
      else {
         // Stop within current segment
         double s_final_xi = s_final - std::floor(s_final);  // final fractional position
         L = (s_final_xi - xi) * seg->GetLength();
            
         // Use FULL segment volume for Parker flux normalization
         vol = SEP::FieldLine::GetSegmentVolume(seg, iFieldLine);
            
         dFlux = weight * L / (dt * vol);
         ParkerStream[pBin] += dFlux;
            
         // Particle stops here
         s = s_final;
         break;
       }
    }
  }


  // Branch B: Backward motion (s_final < s_init)
  if (s_final < s_init) {
    while (s_final < s && seg) {
      double L, dFlux, vol;
      double xi = s - std::floor(s);  // fractional part within current segment
        
       // Get segment data pointer
       double* ParkerStream = seg->GetDatum_ptr(S);
       if (ParkerStream==NULL) return;
        
        // Get segment endpoints for validation
        auto v0 = seg->GetBegin();
        auto v1 = seg->GetEnd();
        
       if (std::floor(s_final) != std::floor(s)) {
          // s_final and s are in different segments
          // Particle crosses segment boundary backward - travels from xi to beginning (0.0)
            
          // Distance traveled in this segment (from xi to beginning of segment)
          L = xi * seg->GetLength();
          vol = SEP::FieldLine::GetSegmentVolume(seg, iFieldLine);
            
          // Calculate Parker flux contribution (negative for backward motion)
          dFlux = -weight * L / (dt * vol);
            
          // Add flux to this segment
          ParkerStream[pBin] += dFlux;
            
          // Move to previous segment (proper two-step process)
          s = std::floor(s) - 1.0E-8;     // Step 2: Position at end of previous segment (xi = 1.0)
          seg = seg->GetPrev();
	  if (seg==NULL) return;
        }
        else {
          // Particle motion stops within the current segment
            
          // Distance traveled within current segment (backward)
          L = (s - s_final) * seg->GetLength();
          vol = SEP::FieldLine::GetSegmentVolume(seg, iFieldLine);
            
          // Calculate Parker flux contribution (negative for backward motion)
          dFlux = -weight * L / (dt * vol);
            
          // Add flux to this segment
          ParkerStream[pBin] += dFlux;
            
          // Update position to final position and stop
          s = s_final;
          break;
        }
      }
   }
}






