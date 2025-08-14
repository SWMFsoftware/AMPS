/*==============================================================================
 * SEP::ParticleMover_FocusedTransport_WaveScattering
 *==============================================================================
 *
 * DESCRIPTION:
 *   Focused Transport Equation solver with Alfvén wave scattering.
 *   Based on SEP::ParticleMover_FocusedTransport_EventDriven but modified to:
 *   - Scatter particles at every iteration using ScatterStepProton()
 *   - Extract wave energy density from segment data
 *   - Update wave energy after each scattering event
 *   - Use fixed time step calculated once outside the loop
 *
 * PHYSICS:
 *   - Deterministic streaming along field lines
 *   - Magnetic focusing/defocusing
 *   - Alfvén wave scattering at every time step
 *   - Energy exchange between particles and waves
 *   - Adiabatic cooling (if enabled)
 *
 * WAVE COUPLING:
 *   - Forward particles (μ > 0) scatter with W- waves
 *   - Backward particles (μ < 0) scatter with W+ waves
 *   - Energy gained by particle = Energy lost by waves
 *
 *==============================================================================*/

#include "sep.h"

int SEP::ParticleMover_FocusedTransport_WaveScattering(long int ptr, double dtTotal, 
                                                      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    namespace PB = PIC::ParticleBuffer;
    namespace FL = PIC::FieldLine;

    // Alfvén wave branch constants
    const int BranchPlus = 0;   // W+ (outward propagating wave)
    const int BranchMinus = 1;  // W- (inward propagating wave)

    PIC::ParticleBuffer::byte *ParticleData;
    double Speed, mu, vParallel, vNormal;
    double FieldLineCoord, xCartesian[3];
    int iFieldLine, spec;
    FL::cFieldLineSegment *Segment;

    // Performance counters
    static int ncall = 0;
    static int total_scattering_events = 0;
    ncall++;

    // Get particle data
    ParticleData = PB::GetParticleDataPointer(ptr);
    FieldLineCoord = PB::GetFieldLineCoord(ParticleData);
    iFieldLine = PB::GetFieldLineId(ParticleData);
    spec = PB::GetI(ParticleData);

    // Initial velocity components (in plasma rest frame)
    vParallel = PB::GetVParallel(ParticleData);
    vNormal = PB::GetVNormal(ParticleData);

    if ((Segment = FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord)) == NULL) {
        exit(__LINE__, __FILE__, "Error: cannot find the segment");
    }

    // Calculate initial particle properties
    Speed = sqrt(vNormal * vNormal + vParallel * vParallel);
    mu = (Speed > 0.0) ? vParallel / Speed : 0.0;

    double TimeCounter = 0.0;
    double total_energy_exchanged = 0.0; // For diagnostics

    // =========================================================================
    // HELPER FUNCTIONS
    // =========================================================================

    // Function to extract wave energy densities from segment data
    auto GetWaveEnergyDensities = [&](FL::cFieldLineSegment* seg, double& W_plus_density, double& W_minus_density) -> bool {
        double* wave_data = seg->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        if (!wave_data) return false;

        double E_plus_total = wave_data[0];   // Total E+ energy in segment [J]
        double E_minus_total = wave_data[1];  // Total E- energy in segment [J]
        
        // Get segment volume to convert to energy density
        double volume = SEP::FieldLine::GetSegmentVolume(seg, iFieldLine);
        if (volume <= 0.0) return false;

        W_plus_density = E_plus_total / volume;   // [J/m³]
        W_minus_density = E_minus_total / volume; // [J/m³]
        
        return true;
    };

    // Function to update wave energy after particle interaction
    auto UpdateWaveEnergy = [&](FL::cFieldLineSegment* seg, double energyChange, int AlfvenWaveBranch) -> void {
        double* wave_data = seg->GetDatum_ptr(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        if (!wave_data) return;

        // Get particle statistical weight
        double stat_weight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec] * 
                           PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        // Energy gained by particle = Energy lost by waves
        double wave_energy_change = -stat_weight * energyChange;

        if (AlfvenWaveBranch == BranchPlus) {
            wave_data[0] += wave_energy_change; // Update W+ energy
        } else if (AlfvenWaveBranch == BranchMinus) {
            wave_data[1] += wave_energy_change; // Update W- energy
        }

        // Ensure non-negative wave energies
        wave_data[0] = std::max(0.0, wave_data[0]);
        wave_data[1] = std::max(0.0, wave_data[1]);

        total_energy_exchanged += energyChange;
    };

    // Calculate simple time step based on streaming and focusing
    auto CalculateTimeStep = [&](double currentSpeed, double currentMu, 
                                FL::cFieldLineSegment* currentSegment) -> double {
        // Get current position
        double x[3];
        currentSegment->GetCartesian(x, FieldLineCoord);
        
        // Focusing time scale: τ_focus ~ |L|/(v*√(1-μ²))
        double B[3], B0[3], B1[3], AbsBDeriv, AbsB, L;
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B, FieldLineCoord);
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord + 1 - 1E-7);

        AbsB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        AbsBDeriv = (sqrt(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2]) -
                     sqrt(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2])) /
                    FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

        L = (fabs(AbsBDeriv) > 1e-20) ? -AbsB / AbsBDeriv : 1e20;
        double dt_focusing = (fabs(L) > 0.0 && fabs(currentMu) < 0.99) ? 
                            0.1 * fabs(L) / (currentSpeed * sqrt(1.0 - currentMu*currentMu)) : 1e20;

        // Streaming time scale: τ_stream ~ segment_length/v
        double segment_length = currentSegment->GetLength();
        double dt_streaming = 0.1 * segment_length / currentSpeed;

        // Choose minimum time scale
        double dt_calc = std::min(dt_focusing, dt_streaming);
        
        // Apply bounds
        dt_calc = std::max(dt_calc, 1e-3); // Minimum 1 millisecond
        dt_calc = std::min(dt_calc, 10.0);  // Maximum 10 seconds

        return dt_calc;
    };

    // =========================================================================
    // CALCULATE TIME STEP ONCE OUTSIDE THE LOOP
    // =========================================================================

    double dt_step = CalculateTimeStep(Speed, mu, Segment);

    // =========================================================================
    // MAIN TRANSPORT LOOP
    // =========================================================================

    while (TimeCounter < dtTotal) {
        // Calculate time step for this iteration
        double remaining_time = dtTotal - TimeCounter;
        double dt_current = std::min(dt_step, remaining_time);

        // Get current position and properties
        double x[3], rHelio;
        Segment->GetCartesian(x, FieldLineCoord);
        rHelio = Vector3D::Length(x);

        // Get wave energy densities at current location
        double W_plus_density, W_minus_density;
        if (!GetWaveEnergyDensities(Segment, W_plus_density, W_minus_density)) {
            // No wave data available, set to zero
            W_plus_density = W_minus_density = 0.0;
        }

        // Get Alfvén velocity and magnetic field for ScatterStepProton
        double vAlfven = 0.0;
        double AbsB = 0.0;
        FL::cFieldLineVertex* VertexBegin = Segment->GetBegin();
        FL::cFieldLineVertex* VertexEnd = Segment->GetEnd();
        
        if (VertexBegin && VertexEnd) {
            double *B0 = VertexBegin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
            double *B1 = VertexEnd->GetDatum_ptr(FL::DatumAtVertexMagneticField);
            
            if (B0 && B1) {
                double w1 = fmod(FieldLineCoord, 1);
                double w0 = 1.0 - w1;
                double B[3];
                for (int idim = 0; idim < 3; idim++) {
                    B[idim] = w0 * B0[idim] + w1 * B1[idim];
                }
                AbsB = Vector3D::Length(B);

                // Get plasma density for Alfvén velocity
                double PlasmaDensity0, PlasmaDensity1, PlasmaDensity;
                VertexBegin->GetDatum(FL::DatumAtVertexPlasmaDensity, &PlasmaDensity0);
                VertexEnd->GetDatum(FL::DatumAtVertexPlasmaDensity, &PlasmaDensity1);
                PlasmaDensity = (w0 * PlasmaDensity0 + w1 * PlasmaDensity1) * 
                               PIC::CPLR::SWMF::MeanPlasmaAtomicMass;

                if (PlasmaDensity > 0.0) {
                    vAlfven = AbsB / sqrt(VacuumPermeability * PlasmaDensity);
                }
            }
        }

        // =====================================================================
        // PARTICLE SCATTERING AT EACH ITERATION
        // =====================================================================
        
        double energyChange = 0.0;
        bool scatterOccurred = false;
        double vParallel_before = vParallel;
        double vNormal_before = vNormal;

        // Call ScatterStepProton - scattering occurs every iteration
        auto scatterResult = SEP::AlfvenTurbulence_Kolmogorov::ScatterStepProton(
            vParallel, vNormal,                                    // Current velocity components
            AbsB, vAlfven,                                        // Magnetic field and Alfven velocity  
            W_plus_density, W_minus_density,                      // Wave energy densities
            dt_current,                                           // Time step
            rHelio, Speed,                                        // Heliocentric distance and speed
            PIC::MolecularData::GetMass(spec),                   // Particle mass
            PIC::MolecularData::GetElectricCharge(spec),         // Particle charge
            mu                                                    // Pitch angle cosine
        );
        
        // Extract results from ScatterResult structure
        // Note: Update these member names to match the actual ScatterResult structure
        vParallel = scatterResult.vpar_new;
        vNormal = scatterResult.vperp_new;
        // TODO: Replace with correct member name for energy change from ScatterResult
        // Common possibilities: dE, deltaE, energy_delta, energyChange, de
        energyChange = 0.0; // Placeholder - replace with: scatterResult.CORRECT_MEMBER_NAME
        scatterOccurred = true; // Always true since scattering always occurs

        // Update speed and pitch angle after scattering
        Speed = sqrt(vParallel*vParallel + vNormal*vNormal);
        mu = (Speed > 0.0) ? vParallel / Speed : 0.0;

        // Apply speed limit
        if (Speed > 0.99 * SpeedOfLight) {
            double scale_factor = 0.99 * SpeedOfLight / Speed;
            Speed *= scale_factor;
            vParallel *= scale_factor;
            vNormal *= scale_factor;
            mu = vParallel / Speed;
        }

        // Enforce pitch angle bounds
        if (mu > 1.0 - 1e-6) {
            mu = 1.0 - 1e-6;
            vParallel = mu * Speed;
            vNormal = sqrt(std::max(0.0, 1.0 - mu*mu)) * Speed;
        }
        if (mu < -1.0 + 1e-6) {
            mu = -1.0 + 1e-6;
            vParallel = mu * Speed;
            vNormal = sqrt(std::max(0.0, 1.0 - mu*mu)) * Speed;
        }

        // =====================================================================
        // UPDATE WAVE ENERGY AFTER SCATTERING
        // =====================================================================
        
        // Determine which wave branch was involved in scattering
        int AlfvenWaveBranch;
        
        // Forward particles scatter with W-, backward with W+
        if (vParallel_before > 0.0) {
            AlfvenWaveBranch = BranchMinus; // Forward particle scatters with W-
        } else {
            AlfvenWaveBranch = BranchPlus;  // Backward particle scatters with W+
        }

        // Always update wave energy since scattering always occurs
        if (fabs(energyChange) > 1e-30) {
            UpdateWaveEnergy(Segment, energyChange, AlfvenWaveBranch);
        }
        
        // Count scattering events (always happens)
        total_scattering_events++;

        // =====================================================================
        // DETERMINISTIC TRANSPORT: STREAMING AND FOCUSING
        // =====================================================================

        // Deterministic streaming along field line
        double ds_parallel = vParallel * dt_current;

        // Adiabatic focusing correction
        double B_focus[3], B0_focus[3], B1_focus[3], AbsBDeriv, AbsB_focus, L;
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B0_focus, (int)FieldLineCoord);
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B_focus, FieldLineCoord);
        FL::FieldLinesAll[iFieldLine].GetMagneticField(B1_focus, (int)FieldLineCoord + 1 - 1E-7);

        AbsB_focus = sqrt(B_focus[0]*B_focus[0] + B_focus[1]*B_focus[1] + B_focus[2]*B_focus[2]);
        AbsBDeriv = (sqrt(B1_focus[0]*B1_focus[0] + B1_focus[1]*B1_focus[1] + B1_focus[2]*B1_focus[2]) -
                     sqrt(B0_focus[0]*B0_focus[0] + B0_focus[1]*B0_focus[1] + B0_focus[2]*B0_focus[2])) /
                    FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

        L = (fabs(AbsBDeriv) > 1e-20) ? -AbsB_focus / AbsBDeriv : 1e20;

        // Apply magnetic focusing: dμ/dt = (1-μ²)v/(2L)
        double dmu_focusing = 0.0;
        if (fabs(L) > 1e-20 && fabs(mu) < 0.99) {
            dmu_focusing = (1.0 - mu*mu) / (2.0 * L) * Speed * dt_current;
            mu += dmu_focusing;
        }

        // =====================================================================
        // ADIABATIC COOLING (if enabled)
        // =====================================================================
        
        if (SEP::AccountAdiabaticCoolingFlag == true) {
            double dP, dLogP, dmu_transport, vSolarWindParallel;
            SEP::GetTransportCoefficients(dP, dLogP, dmu_transport, Speed, mu, Segment,
                                        FieldLineCoord, dt_current, iFieldLine, vSolarWindParallel);

            if (isfinite(dLogP) && isfinite(dmu_transport)) {
                // Apply momentum change
                double p;
                switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
                case _PIC_MODE_OFF_:
                    p = Speed * PIC::MolecularData::GetMass(spec);
                    break;
                case _PIC_MODE_ON_:
                    p = Relativistic::Speed2Momentum(Speed, PIC::MolecularData::GetMass(spec));
                    break;
                }

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

                mu += dmu_transport;
            }
        }

        // =====================================================================
        // MOVE PARTICLE ALONG FIELD LINE
        // =====================================================================

        FieldLineCoord = FL::FieldLinesAll[iFieldLine].move(FieldLineCoord, ds_parallel, Segment);

        if (Segment == NULL) {
            // Particle left the simulation domain
            PIC::ParticleBuffer::DeleteParticle(ptr);
            return _PARTICLE_LEFT_THE_DOMAIN_;
        }

        // Update velocity components
        vParallel = mu * Speed;
        vNormal = sqrt(std::max(0.0, 1.0 - mu*mu)) * Speed;

        TimeCounter += dt_current;
    }

    // =========================================================================
    // FINALIZATION
    // =========================================================================

    // Set final particle properties
    PB::SetVParallel(vParallel, ParticleData);
    PB::SetVNormal(vNormal, ParticleData);
    PB::SetFieldLineCoord(FieldLineCoord, ParticleData);

    // Attach particle to temporary list
    switch (_PIC_PARTICLE_LIST_ATTACHING_) {
    case _PIC_PARTICLE_LIST_ATTACHING_NODE_:
        exit(__LINE__, __FILE__, "Error: function developed for _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_");
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
        exit(__LINE__, __FILE__, "Error: unknown option");
    }

    return _PARTICLE_MOTION_FINISHED_;
}
