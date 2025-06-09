//the model of the Alfven turbulence. 

#include "sep.h"

PIC::Datum::cDatumStored SEP::AlfvenTurbulence::WaveEnergyDensity(2,"\"W+\",\"W-\"",true);

double SEP::AlfvenTurbulence::ModelInit::dB_B(double r) {
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
void SEP::AlfvenTurbulence::ModelInit::Init() {      
  namespace FL = PIC::FieldLine;

  // Check if field lines are available
  if (FL::FieldLinesAll == nullptr || FL::nFieldLine == 0) {
    if (PIC::ThisThread == 0) {
      printf("SEP::AlfvenTurbulence::ModelInit: No field lines available for initialization\n");
    }

    return;
  }

  // Ensure the wave energy density datum is activated
  if (!SEP::AlfvenTurbulence::WaveEnergyDensity.is_active()) {
    if (PIC::ThisThread == 0) {
      printf("SEP::AlfvenTurbulence::ModelInit: Warning - WaveEnergyDensity datum not activated\n");
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
       segment->SetDatum(WaveEnergyDensity,W);

       totalSegmentsProcessed++;

       // Move to next segment
       segment = segment->GetNext();
     }

     // Report completion
     if (PIC::ThisThread == 0) {
       printf("SEP::AlfvenTurbulence::ModelInit: Initialized turbulence for %d segments across %d field lines\n",
       totalSegmentsProcessed, FL::nFieldLine);
     }
   }
}
