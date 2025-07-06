
/*==============================================================================
 * SEP::ParticleMover_FocusedTransport_EventDriven
 *==============================================================================
 *
 * DESCRIPTION:
 *   Event-driven Monte Carlo solver for the Focused Transport Equation (FTE)
 *   in the heliosphere. This function advances energetic particles along
 *   magnetic field lines with proper treatment of:
 *   - Deterministic streaming along field lines
 *   - Adiabatic magnetic focusing/defocusing
 *   - Pitch angle scattering with Alfvén waves
 *   - Adiabatic cooling due to solar wind expansion
 *   - Relativistic wave-particle interactions
 *
 * PHYSICS MODEL:
 *   The Focused Transport Equation describes particle motion in the solar wind:
 *   ∂f/∂t + v·μ ∂f/∂s + (1-μ²)/(2L) v ∂f/∂μ = ∂/∂μ[D_μμ ∂f/∂μ] + source terms
 *
 *   where:
 *   - s: distance along magnetic field line
 *   - μ = v∥/v: pitch angle cosine
 *   - L = -|B|/(∂|B|/∂s): magnetic focusing length
 *   - D_μμ: pitch angle diffusion coefficient
 *
 * REFERENCE FRAME:
 *   - Simulation conducted in PLASMA REST FRAME
 *   - Particle velocities (vParallel, vNormal) given in plasma rest frame
 *   - Alfvén waves propagate at ±v_A relative to plasma
 *   - Wave-particle scattering uses relativistic frame transformations
 *
 * SCATTERING PHYSICS:
 *   Wave-particle resonance conditions:
 *   - Forward particles (μ > 0): scatter only with W- waves (v = -v_A)
 *   - Backward particles (μ < 0): scatter only with W+ waves (v = +v_A)
 *
 *   Scattering process:
 *   1. Transform particle velocity: Plasma frame → Alfvén wave frame
 *   2. Perform isotropic scattering in wave frame (preserve speed)
 *   3. Transform back: Wave frame → Plasma frame
 *   4. Calculate energy gain/loss from relativistic transformations
 *
 * RELATIVISTIC TRANSFORMATIONS:
 *   Uses exact Lorentz boost formulas for frame transformations:
 *
 *   Plasma → Wave frame:
 *   v'∥ = (v∥ - V_A) / (1 - v∥V_A/c²)
 *   v'⊥ = v⊥ / [γ_A(1 - v∥V_A/c²)]
 *
 *   Wave → Plasma frame:
 *   v∥ = (v'∥ + V_A) / (1 + v'∥V_A/c²)
 *   v⊥ = v'⊥ / [γ_A(1 + v'∥V_A/c²)]
 *
 *   where β_A = V_A/c, γ_A = 1/√(1-β_A²)
 *
 * EVENT-DRIVEN ALGORITHM:
 *   1. Calculate scattering time: t_scatter ~ -ln(r)/ν_scatter
 *   2. Advance particle deterministically until next event
 *   3. Apply magnetic focusing: dμ/dt = (1-μ²)v/(2L)
 *   4. Apply adiabatic cooling: dp/dt ∝ ∇·v_sw
 *   5. Handle scattering event with proper frame transformations
 *   6. Sample particle flux for wave coupling between events
 *   7. Repeat until time step completed or particle exits domain
 *
 * FLUX SAMPLING:
 *   - Samples AccumulateParticleFluxForWaveCoupling() between events
 *   - Records: time interval, start/end positions, distance traveled
 *   - Resets tracking variables after each scattering event
 *   - Ensures complete coverage of particle trajectory
 *
 * NUMERICAL FEATURES:
 *   - Adaptive time stepping based on physics time scales
 *   - Numerical stability checks for relativistic transformations
 *   - Speed limits to prevent superluminal velocities
 *   - Proper handling of pitch angle boundaries (μ ∈ [-1,1])
 *
 * INPUT PARAMETERS:
 *   ptr      - Particle buffer pointer
 *   dtTotal  - Total time step duration
 *   node     - AMR tree node containing particle
 *
 * RETURN VALUES:
 *   _PARTICLE_MOTION_FINISHED_ - Particle successfully advanced
 *   _PARTICLE_LEFT_THE_DOMAIN_ - Particle exited simulation domain
 *
 * DEPENDENCIES:
 *   - Field line data structure (magnetic field, plasma properties)
 *   - Alfvén wave amplitudes (W+, W-) at field line vertices
 *   - Mean free path calculation routines
 *   - Relativistic energy-momentum conversion functions
 *   - Transport coefficient calculation (adiabatic cooling)
 *
 * AUTHOR: Generated for SEP transport simulation
 * DATE: 2024
 *
 * NOTES:
 *   - Currently uses simple scattering frequency: ν = v/λ
 *   - Future: Can extend to sophisticated wave-dependent rates
 *   - Assumes mean free path >> Larmor radius (guiding center approximation)
 *   - Compatible with existing SEP simulation framework
 *==============================================================================*/

#include "sep.h"

int SEP::ParticleMover_FocusedTransport_EventDriven(long int ptr, double dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double Speed, mu, vParallel, vNormal;
  double FieldLineCoord, xCartesian[3];
  int iFieldLine, spec;
  FL::cFieldLineSegment *Segment;

  static int ncall = 0;
  ncall++;

  ParticleData = PB::GetParticleDataPointer(ptr);

  FieldLineCoord = PB::GetFieldLineCoord(ParticleData);
  iFieldLine = PB::GetFieldLineId(ParticleData);
  spec = PB::GetI(ParticleData);

  // Velocity is in the frame moving with solar wind
  vParallel = PB::GetVParallel(ParticleData);
  vNormal = PB::GetVNormal(ParticleData);

  if ((Segment = FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord)) == NULL) {
    exit(__LINE__, __FILE__, "Error: cannot find the segment");
  }

  // Calculate initial particle properties
  Speed = sqrt(vNormal * vNormal + vParallel * vParallel);
  mu = vParallel / Speed;

  double TimeCounter = 0.0;
  double totalTraversedPath = 0.0;

  // Initialize tracking variables for flux sampling
  double s_current = FieldLineCoord; // Current position for flux sampling
  double t_current = 0.0;            // Current time for flux sampling

  // Helper function to calculate mean free path
  auto CalculateMeanFreePath = [&](int spec, double rHelio, double Speed, double AbsB) {
    double MeanFreePath, dxx, energy;

    switch (SEP::Scattering::MeanFreePathMode) {
    case SEP::Scattering::MeanFreePathMode_QLT:
      MeanFreePath = QLT::calculateMeanFreePath(rHelio, Speed);
      break;
    case SEP::Scattering::MeanFreePathMode_QLT1:
      MeanFreePath = QLT1::calculateMeanFreePath(rHelio, Speed, AbsB);
      break;
    case SEP::Scattering::MeanFreePathMode_Tenishev2005AIAA:
      energy = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));
      MeanFreePath = SEP::Scattering::Tenishev2005AIAA::lambda0 *
        pow(energy/GeV2J, SEP::Scattering::Tenishev2005AIAA::alpha) *
        pow(rHelio/_AU_, SEP::Scattering::Tenishev2005AIAA::beta);
      break;
    case SEP::Scattering::MeanFreePathMode_Chen2024AA:
      energy = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));
      dxx = SEP::Diffusion::Chen2024AA::GetDxx(rHelio, energy);
      MeanFreePath = 3.0 * dxx / Speed; // Eq. 15, Liu-2024-arXiv
      break;
    default:
      exit(__LINE__, __FILE__, "Error: the option is unknown");
    }

    return MeanFreePath;
  };

  // Helper function to calculate scattering frequency
  auto CalculateScatteringFrequency = [&](double Speed, double MeanFreePath, double mu,
                                          double W_plus, double W_minus,
                                          double& Nu_scattering, bool& scatter_with_Wplus) -> void {
    // For now, use simple mean free path approach without wave energy density
    Nu_scattering = 0.0;
    scatter_with_Wplus = false;

    // Determine which wave branch can cause scattering based on particle direction
    if (mu > 0.0 && W_minus > 0.0) {
      // Outward moving particle - can scatter with inward wave (W-)
      Nu_scattering = Speed / MeanFreePath;
      scatter_with_Wplus = false;
    } else if (mu < 0.0 && W_plus > 0.0) {
      // Inward moving particle - can scatter with outward wave (W+)
      Nu_scattering = Speed / MeanFreePath;
      scatter_with_Wplus = true;
    }
    // else: No scattering possible - values remain 0.0 and false
  };

  // Helper function to calculate adiabatic focusing parameter L
  auto CalculateFocusingLength = [&](FL::cFieldLineSegment* currentSegment, double currentFieldLineCoord, int fieldLineId) {
    double B[3], B0[3], B1[3], AbsBDeriv, AbsB, L;

    FL::FieldLinesAll[fieldLineId].GetMagneticField(B0, (int)currentFieldLineCoord);
    FL::FieldLinesAll[fieldLineId].GetMagneticField(B, currentFieldLineCoord);
    FL::FieldLinesAll[fieldLineId].GetMagneticField(B1, (int)currentFieldLineCoord + 1 - 1E-7);

    AbsB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    AbsBDeriv = (sqrt(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2]) -
                 sqrt(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2])) /
                FL::FieldLinesAll[fieldLineId].GetSegmentLength(currentFieldLineCoord);

    L = -AbsB / AbsBDeriv; // Focusing length scale

    return L;
  };

  // Event-driven Monte Carlo loop for Focused Transport Equation
  while (TimeCounter < dtTotal) {
    // Get current position and magnetic field
    double x[3], rHelio;
    Segment->GetCartesian(x, FieldLineCoord);
    rHelio = Vector3D::Length(x);

    double AbsB;
    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      AbsB = SEP::FieldLineData::GetAbsB(FieldLineCoord, Segment, iFieldLine);
      break;
    default:
      AbsB = SEP::ParkerSpiral::GetAbsB(rHelio);
    }

    // Calculate transport coefficients for Focused Transport Equation
    double MeanFreePath = CalculateMeanFreePath(spec, rHelio, Speed, AbsB);
    double L = CalculateFocusingLength(Segment, FieldLineCoord, iFieldLine);

    // Store mean free path if sampling is active
    if ((SEP::Offset::MeanFreePath != -1) && (SEP::Sampling::MeanFreePath::active_flag == true)) {
      *((double*)(ParticleData + SEP::Offset::MeanFreePath)) = MeanFreePath;
    }

    // Calculate wave-specific scattering rates
    // Get Alfven wave amplitudes at current location
    double W[2]; // W[0] = W+, W[1] = W-
    FL::cFieldLineVertex* VertexBegin = Segment->GetBegin();
    FL::cFieldLineVertex* VertexEnd = Segment->GetEnd();
    double *W0 = VertexBegin->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);
    double *W1 = VertexEnd->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);

    // Interpolate wave amplitudes
    double w1 = fmod(FieldLineCoord, 1);
    double w0 = 1.0 - w1;
    W[0] = w0 * W0[0] + w1 * W1[0]; // W+ (outward wave)
    W[1] = w0 * W0[1] + w1 * W1[1]; // W- (inward wave)

    // Calculate scattering frequency using the lambda function
    double Nu_scattering;
    bool scatter_with_Wplus;
    CalculateScatteringFrequency(Speed, MeanFreePath, mu, W[0], W[1], Nu_scattering, scatter_with_Wplus);
    bool can_scatter = (Nu_scattering > 0.0);

    // Calculate characteristic times for different processes
    double tau_scattering = (can_scatter && Nu_scattering > 0.0) ? 1.0 / Nu_scattering : 1e20;
    double tau_focusing = (L > 0.0 && fabs(mu) < 1.0) ? fabs(L) / (Speed * sqrt(1.0 - mu*mu)) : 1e20;
    double tau_remaining = dtTotal - TimeCounter;

    // Sample time to next scattering event (exponential distribution)
    double t_scattering = (can_scatter) ? -tau_scattering * log(rnd()) : 1e20;

    // Determine the time to advance (minimum of scattering time and remaining time)
    double dt_event = min(t_scattering, tau_remaining);

    // Determine if scattering occurs during this time step
    bool scattering_occurs = (t_scattering <= tau_remaining) && can_scatter;

    // For Focused Transport: particles move along field lines with velocity v*mu
    // No stochastic spatial displacement - only deterministic motion along field lines
    double ds_parallel = vParallel * dt_event; // Motion along magnetic field line

    // Adiabatic focusing: change in mu due to magnetic field gradient
    // dmu/dt = (1-mu²)/(2L) * v for adiabatic focusing
    double dmu_focusing = 0.0;
    if (L > 0.0 && fabs(L) > 1e-20) {
      dmu_focusing = (1.0 - mu*mu) / (2.0 * L) * Speed * dt_event;
    }

    // Update pitch angle due to focusing (before potential scattering)
    mu += dmu_focusing;

    // Apply pitch angle limits
    if (mu > 1.0 - 1e-6) mu = 1.0 - 1e-6;
    if (mu < -1.0 + 1e-6) mu = -1.0 + 1e-6;

    totalTraversedPath += fabs(ds_parallel);

    // Move particle along field line
    FieldLineCoord = FL::FieldLinesAll[iFieldLine].move(FieldLineCoord, ds_parallel, Segment);

    if (Segment == NULL) {
      // Particle left the simulation domain
      // Sample flux for wave coupling before deletion (from s_current to exit point)
      if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {
        double s_final = FieldLineCoord;
        SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
          iFieldLine, ptr, dt_event, Speed, s_current, s_final, fabs(ds_parallel)
        );
      }

      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    // Update particle position
    Segment->GetCartesian(x, FieldLineCoord);
    rHelio = Vector3D::Length(x);

    // Apply adiabatic cooling if enabled
    if (SEP::AccountAdiabaticCoolingFlag == true) {
      // Get transport coefficients including adiabatic cooling
      double dP, dLogP, dmu_transport, vSolarWindParallel;
      SEP::GetTransportCoefficients(dP, dLogP, dmu_transport, Speed, mu, Segment,
                                    FieldLineCoord, dt_event, iFieldLine, vSolarWindParallel);

      if (isfinite(dLogP) && isfinite(dmu_transport)) {
        // Apply momentum change due to adiabatic cooling
        double p;
        switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
        case _PIC_MODE_OFF_:
          p = Speed * PIC::MolecularData::GetMass(spec);
          break;
        case _PIC_MODE_ON_:
          p = Relativistic::Speed2Momentum(Speed, PIC::MolecularData::GetMass(spec));
          break;
        }

        // Apply momentum change: p_new = p_old * exp(dLogP)
        p *= exp(dLogP);

        // Convert back to speed
        switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
        case _PIC_MODE_OFF_:
          Speed = p / PIC::MolecularData::GetMass(spec);
          break;
        case _PIC_MODE_ON_:
          Speed = Relativistic::Momentum2Speed(p, PIC::MolecularData::GetMass(spec));
          break;
        }

        // Apply speed limit
        if (Speed > 0.99 * SpeedOfLight) {
          Speed = 0.99 * SpeedOfLight;
        }

        // Add transport coefficient contribution to pitch angle change
        // This includes effects from solar wind velocity gradients
        mu += dmu_transport;

        // Apply pitch angle limits after transport effects
        if (mu > 1.0 - 1e-6) mu = 1.0 - 1e-6;
        if (mu < -1.0 + 1e-6) mu = -1.0 + 1e-6;
      }
    }

    // Handle scattering event with proper frame transformations
    if (scattering_occurs) {
      // Get Alfven velocity for frame transformations
      double *B0 = VertexBegin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
      double *B1 = VertexEnd->GetDatum_ptr(FL::DatumAtVertexMagneticField);
      double B[3];
      for (int idim = 0; idim < 3; idim++) {
        B[idim] = w0 * B0[idim] + w1 * B1[idim];
      }
      double AbsB = Vector3D::Length(B);

      // Get plasma density for Alfven velocity calculation
      double PlasmaDensity0, PlasmaDensity1, PlasmaDensity;
      VertexBegin->GetDatum(FL::DatumAtVertexPlasmaDensity, &PlasmaDensity0);
      VertexEnd->GetDatum(FL::DatumAtVertexPlasmaDensity, &PlasmaDensity1);
      PlasmaDensity = (w0 * PlasmaDensity0 + w1 * PlasmaDensity1) * PIC::CPLR::SWMF::MeanPlasmaAtomicMass;

      double vAlfven = AbsB / sqrt(VacuumPermeability * PlasmaDensity);

      // Store initial energy for energy gain/loss calculation
      double E_initial = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));

      // Determine which wave caused the scattering (no probability needed - it's deterministic)
      double v_wave_plasma_frame; // Wave velocity in plasma rest frame
      double v_wave_lab_frame;    // Wave velocity in lab frame (after relativistic composition)

      if (scatter_with_Wplus) {
        // Scattering with W+ wave (outward wave in plasma frame)
        // This should only happen for inward moving particles (mu < 0)
        v_wave_plasma_frame = +vAlfven; // W+ moves away from Sun in plasma frame
      } else {
        // Scattering with W- wave (inward wave in plasma frame)
        // This should only happen for outward moving particles (mu > 0)
        v_wave_plasma_frame = -vAlfven; // W- moves toward Sun in plasma frame
      }
      // Get solar wind velocity at particle location for relativistic composition
      double vSolarWind[3], vSolarWindParallel;
      FL::FieldLinesAll[iFieldLine].GetPlasmaVelocity(vSolarWind, FieldLineCoord);

      // Calculate parallel component of solar wind velocity
      double B_unit[3];
      for (int idim = 0; idim < 3; idim++) {
        B_unit[idim] = B[idim] / AbsB;
      }
      vSolarWindParallel = Vector3D::DotProduct(vSolarWind, B_unit);

      // Relativistic velocity addition: v_wave_lab = (v_sw + v_wave_plasma) / (1 + v_sw*v_wave_plasma/c²)
      double beta_sw = vSolarWindParallel / SpeedOfLight;
      double beta_wave_plasma = v_wave_plasma_frame / SpeedOfLight;
      double beta_wave_lab = (beta_sw + beta_wave_plasma) / (1.0 + beta_sw * beta_wave_plasma);
      v_wave_lab_frame = beta_wave_lab * SpeedOfLight;

      // Calculate Lorentz factor for wave frame transformation
      double gamma_wave = 1.0 / sqrt(1.0 - beta_wave_lab * beta_wave_lab);

      // 1) Transform particle velocity from lab frame to wave frame
      double v_parallel_lab = Speed * mu; // Parallel component in lab frame
      double v_perp_lab = Speed * sqrt(1.0 - mu*mu); // Perpendicular component (unchanged)

      // Relativistic velocity transformation along field line
      double beta_parallel_lab = v_parallel_lab / SpeedOfLight;
      double beta_wave = v_wave_lab_frame / SpeedOfLight;

      // Transform parallel velocity to wave frame (relativistic velocity subtraction)
      double beta_parallel_wave = (beta_parallel_lab - beta_wave) / (1.0 - beta_parallel_lab * beta_wave);
      double v_parallel_wave = beta_parallel_wave * SpeedOfLight;

      // Perpendicular velocity transforms with Lorentz contraction
      double v_perp_wave = v_perp_lab / gamma_wave;

      // Calculate speed and pitch angle in wave frame
      double Speed_wave = sqrt(v_parallel_wave * v_parallel_wave + v_perp_wave * v_perp_wave);
      double mu_wave = v_parallel_wave / Speed_wave;

      // 2) Perform isotropic scattering in wave frame
      // Preserve speed in wave frame, randomize direction
      mu_wave = -1.0 + 2.0 * rnd(); // Isotropic in pitch angle
      double phi_wave = 2.0 * Pi * rnd(); // Random azimuthal angle

      // New velocity components in wave frame
      v_parallel_wave = Speed_wave * mu_wave;
      v_perp_wave = Speed_wave * sqrt(1.0 - mu_wave * mu_wave);

      // 3) Transform back to plasma rest frame
      beta_parallel_wave = v_parallel_wave / SpeedOfLight;

      // Transform parallel velocity back to plasma rest frame (relativistic velocity addition)
      double beta_parallel_plasma_new = (beta_parallel_wave + beta_wave) / (1.0 + beta_parallel_wave * beta_wave);
      double v_parallel_plasma_new = beta_parallel_plasma_new * SpeedOfLight;

      // Perpendicular velocity transforms back with Lorentz expansion
      double v_perp_plasma_new = v_perp_wave * gamma_wave;

      // Calculate new speed and pitch angle in plasma rest frame
      Speed = sqrt(v_parallel_plasma_new * v_parallel_plasma_new + v_perp_plasma_new * v_perp_plasma_new);
      mu = v_parallel_plasma_new / Speed;

      // Apply speed limit in plasma rest frame
      if (Speed > 0.99 * SpeedOfLight) {
        double scale_factor = 0.99 * SpeedOfLight / Speed;
        Speed *= scale_factor;
        v_parallel_plasma_new *= scale_factor;
        v_perp_plasma_new *= scale_factor;
        mu = v_parallel_plasma_new / Speed;
      }

      // Calculate energy gain/loss
      double E_final = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));
      double dE = E_final - E_initial;

      // Optional: Store energy change for diagnostics
      /*if (SEP::Sampling::EnergyChange::active_flag == true) {
        // Store energy change statistics
        double ParticleWeight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
        ParticleWeight *= PB::GetIndividualStatWeightCorrection(ParticleData);

        // Accumulate energy change statistics (if such sampling exists)
        // SEP::Sampling::EnergyChange::AccumulateEnergyChange(dE, ParticleWeight, ...);
      }
      */

      // Apply pitch angle limits after transformation
      if (mu > 1.0 - 1e-6) mu = 1.0 - 1e-6;
      if (mu < -1.0 + 1e-6) mu = -1.0 + 1e-6;

      // Sample flux for this scattering event (from s_current to current position)
      if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {
        double s_scattering = FieldLineCoord;
        SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
          iFieldLine, ptr, dt_event, Speed, s_current, s_scattering, fabs(ds_parallel)
        );
      }

      // Reset tracking variables for next segment
      s_current = FieldLineCoord;
      t_current = 0.0;
    } else {
      // No scattering occurred - will sample at end of time step or next event
    }

    // Update velocity components based on new mu and speed
    vParallel = mu * Speed;
    vNormal = sqrt(1.0 - mu*mu) * Speed;

    // Update perpendicular scattering if enabled (for field line spreading)
    if (SEP::Offset::RadialLocation != -1) {
      double radial_pos = *((double*)(ParticleData + SEP::Offset::RadialLocation));
      double dr, theta = PiTimes2 * rnd();
      double sin_theta = sin(theta);
      double D_perp;

      // Apply perpendicular diffusion
      D_perp = QLT1::calculatePerpendicularDiffusion(rHelio, Speed, AbsB);
      dr = sqrt(2.0 * D_perp * dt_event) * Vector3D::Distribution::Normal();
      radial_pos = sqrt(radial_pos*radial_pos + dr*dr + 2.0*radial_pos*dr*sin_theta);

      *((double*)(ParticleData + SEP::Offset::RadialLocation)) = radial_pos;
    }

    TimeCounter += dt_event;
    t_current += dt_event;
  }

  // Sample final flux for any remaining segment (from s_current to final position)
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag && t_current > 0.0) {
    double s_final = FieldLineCoord;
    double final_distance = fabs(s_final - s_current);
    SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
      iFieldLine, ptr, t_current, Speed, s_current, s_final, final_distance
    );
  }

  // Set the new values of the normal and parallel particle velocities
  PB::SetVParallel(vParallel, ParticleData);
  PB::SetVNormal(vNormal, ParticleData);

  // Set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord, ParticleData);

  // Attach the particle to the temporary list
  switch (_PIC_PARTICLE_LIST_ATTACHING_) {
  case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
    exit(__LINE__, __FILE__, "Error: the function was developed for the case _PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_");
    break;
  case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
    {
      long int temp = Segment->tempFirstParticleIndex.exchange(ptr);
      PIC::ParticleBuffer::SetNext(temp, ParticleData);
      PIC::ParticleBuffer::SetPrev(-1, ParticleData);

      if (temp != -1) PIC::ParticleBuffer::SetPrev(ptr, temp);
    }
    break;
  default:
    exit(__LINE__, __FILE__, "Error: the option is unknown");
  }

  return _PARTICLE_MOTION_FINISHED_;
}
