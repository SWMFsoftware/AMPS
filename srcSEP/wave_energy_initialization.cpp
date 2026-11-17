/*
================================================================================
                    WAVE ENERGY INITIALIZATION SOURCE
================================================================================

FILENAME: wave_energy_initialization.cpp

PURPOSE:
--------
Implementation of wave energy initialization functions for Alfvén wave turbulence
in solar energetic particle (SEP) transport simulations. Provides functions to 
set initial wave energy distributions across field line segments with proper 
heliospheric distance scaling and physical parameter-based initialization.

MAIN FUNCTIONS:
---------------

1. InitializeWaveEnergyInAllSegments() - Core initialization with r^-2 scaling
2. InitializeWaveEnergyFromPhysicalParameters() - Physics-based initialization  
3. CalculateTypicalWaveEnergyDensity1AU() - Convert B0 and turbulence to energy density
4. PrintWaveEnergyProfile() - Display energy distribution along field lines
5. PrintWaveEnergyInitializationSummary() - Global statistics and validation
6. ValidateWaveEnergyInitialization() - Consistency checking and validation

PHYSICS IMPLEMENTATION:
-----------------------
- Heliospheric scaling: ε(r) = ε_1AU × (1AU/r)²
- Wave energy density: ε = (δB²)/(2μ₀) where δB = turbulence_level × B₀
- Equal wave energies: E+ = E- = ε × V_segment / 2
- Segment center position: (x_left + x_right) / 2
- Distance calculation: r = √(x² + y² + z²)

TYPICAL USAGE EXAMPLES:
-----------------------

// Example 1: Basic initialization with specified energy density
void BasicWaveEnergySetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Initialize with 1 pJ/m³ at 1 AU
    double energy_density_1AU = 1.0e-12;  // [J/m³]
    InitializeWaveEnergyInAllSegments(WaveEnergy, energy_density_1AU, true);
    
    // Print profile for verification
    PrintWaveEnergyProfile(WaveEnergy);
    PrintWaveEnergyInitializationSummary(WaveEnergy);
}

// Example 2: Physics-based initialization (recommended)
void PhysicsBasedSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    double B0_1AU = 5.0e-9;        // 5 nT magnetic field
    double turbulence_level = 0.2; // 20% turbulence
    
    InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0_1AU, turbulence_level, true);
    
    // Validate initialization
    double expected_energy = CalculateTypicalWaveEnergyDensity1AU(B0_1AU, turbulence_level);
    bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected_energy, 0.01);
    
    if (!valid) {
        std::cerr << "Error: Wave energy initialization validation failed!" << std::endl;
        exit(__LINE__);
    }
}

// Example 3: Custom reference distance (Parker Solar Probe mission)
void InnerHeliosphereSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    double energy_density_at_10Rs = 1.0e-8;    // High energy near Sun [J/m³]
    double reference_distance = 10.0 * 6.96e8; // 10 solar radii [m]
    
    InitializeWaveEnergyInAllSegments(
        WaveEnergy, 
        energy_density_at_10Rs, 
        reference_distance, 
        true  // verbose
    );
    
    // Print detailed profile for inner heliosphere
    std::vector<int> inner_field_lines = {0, 1, 2, 3, 4};
    PrintWaveEnergyProfile(WaveEnergy, inner_field_lines, 50);
}

// Example 4: Complete simulation setup with validation
void CompleteSimulationSetup() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Step 1: Initialize with typical solar wind conditions
    InitializeWaveEnergyFromPhysicalParameters(
        WaveEnergy,
        WaveEnergyConstants::TYPICAL_B0_1AU,     // 5 nT
        WaveEnergyConstants::TYPICAL_TURBULENCE, // 20%
        true  // verbose output
    );
    
    // Step 2: Validate initialization
    double expected_energy = CalculateTypicalWaveEnergyDensity1AU();
    if (!ValidateWaveEnergyInitialization(WaveEnergy, expected_energy, 0.01)) {
        exit("Wave energy initialization validation failed", __LINE__, __FILE__);
    }
    
    // Step 3: Print comprehensive diagnostics
    PrintWaveEnergyInitializationSummary(WaveEnergy);
    
    // Step 4: Optional - print profile for specific field lines
    std::vector<int> sample_lines = {0, PIC::FieldLine::nFieldLine/4, PIC::FieldLine::nFieldLine/2};
    PrintWaveEnergyProfile(WaveEnergy, sample_lines, 20);
    
    std::cout << "Wave energy initialization complete. System ready for simulation." << std::endl;
}

// Example 5: Sensitivity study with different turbulence levels
void TurbulenceSensitivityStudy() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    std::vector<double> turbulence_levels = {0.1, 0.15, 0.2, 0.25, 0.3};
    double B0 = 5.0e-9;  // Fixed magnetic field
    
    for (size_t i = 0; i < turbulence_levels.size(); ++i) {
        double turb = turbulence_levels[i];
        
        std::cout << "\n=== Turbulence Level: " << turb*100 << "% ===" << std::endl;
        
        // Initialize with current turbulence level
        InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0, turb, false);
        
        // Print summary statistics
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        
        // Validate
        double expected = CalculateTypicalWaveEnergyDensity1AU(B0, turb);
        bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected, 0.01);
        std::cout << "Validation: " << (valid ? "PASSED" : "FAILED") << std::endl;
        
        // Could save results for analysis or comparison
        // SaveEnergyConfiguration(WaveEnergy, "turbulence_" + std::to_string(turb));
    }
}

// Example 6: Distance-dependent study (multi-spacecraft mission)
void DistanceDependentStudy() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Study energy scaling from 0.1 AU to 10 AU
    std::vector<double> reference_distances = {
        0.1 * WaveEnergyConstants::ONE_AU,  // Near Sun
        0.3 * WaveEnergyConstants::ONE_AU,  // Venus orbit
        1.0 * WaveEnergyConstants::ONE_AU,  // Earth orbit
        1.5 * WaveEnergyConstants::ONE_AU,  // Mars orbit
        5.2 * WaveEnergyConstants::ONE_AU   // Jupiter orbit
    };
    
    double energy_density_base = CalculateTypicalWaveEnergyDensity1AU();
    
    for (size_t i = 0; i < reference_distances.size(); ++i) {
        double ref_dist = reference_distances[i];
        double ref_AU = ref_dist / WaveEnergyConstants::ONE_AU;
        
        std::cout << "\n=== Reference Distance: " << ref_AU << " AU ===" << std::endl;
        
        // Initialize with scaled energy density
        double scaled_energy = energy_density_base;  // Energy density at reference distance
        InitializeWaveEnergyInAllSegments(WaveEnergy, scaled_energy, ref_dist, false);
        
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        
        // Print sample profile
        std::vector<int> sample_lines = {0, 1};
        PrintWaveEnergyProfile(WaveEnergy, sample_lines, 10);
    }
}

// Example 7: Integration with simulation main loop
void MainSimulationLoop() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // One-time initialization at simulation start
    static bool initialized = false;
    if (!initialized) {
        InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, 5.0e-9, 0.2, true);
        
        // Validate once at startup
        double expected = CalculateTypicalWaveEnergyDensity1AU(5.0e-9, 0.2);
        if (!ValidateWaveEnergyInitialization(WaveEnergy, expected)) {
            exit("Initial wave energy validation failed", __LINE__, __FILE__);
        }
        
        PrintWaveEnergyInitializationSummary(WaveEnergy);
        initialized = true;
    }
    
    // Main simulation timesteps would follow...
    // for (int timestep = 0; timestep < max_timesteps; ++timestep) {
    //     UpdateAllSegmentsWaveEnergyWithParticleCoupling(WaveEnergy, S_scalar, dt);
    //     // ... other physics updates ...
    // }
}

// Example 8: Error handling and robustness testing
void RobustnessTest() {
    using namespace SEP::AlfvenTurbulence_Kolmogorov;
    
    // Test with extreme parameters
    std::vector<std::pair<double, double>> test_cases = {
        {1.0e-9, 0.1},   // Low field, low turbulence
        {1.0e-8, 0.5},   // High field, high turbulence
        {5.0e-9, 0.01},  // Normal field, very low turbulence
        {5.0e-9, 0.9}    // Normal field, very high turbulence
    };
    
    for (auto& test_case : test_cases) {
        double B0 = test_case.first;
        double turb = test_case.second;
        
        std::cout << "\nTesting B0=" << B0*1e9 << " nT, turb=" << turb*100 << "%" << std::endl;
        
        try {
            InitializeWaveEnergyFromPhysicalParameters(WaveEnergy, B0, turb, false);
            
            double expected = CalculateTypicalWaveEnergyDensity1AU(B0, turb);
            bool valid = ValidateWaveEnergyInitialization(WaveEnergy, expected, 0.05); // 5% tolerance
            
            std::cout << "  Result: " << (valid ? "SUCCESS" : "WARNING") << std::endl;
            
        } catch (const std::exception& e) {
            std::cout << "  Result: ERROR - " << e.what() << std::endl;
        }
    }
}

INTEGRATION NOTES:
------------------
1. Call initialization functions once at simulation startup
2. Use ValidateWaveEnergyInitialization() to ensure proper setup
3. PrintWaveEnergyProfile() is useful for debugging field line geometry
4. PrintWaveEnergyInitializationSummary() provides quick health checks
5. For production runs, set verbose=false to reduce output
6. MPI-parallel safe - all functions handle thread ownership correctly
7. Consistent with AMPS framework data access patterns

PERFORMANCE CONSIDERATIONS:
---------------------------
- Initialization is O(N_segments) and should be fast even for large simulations
- Diagnostic functions can be expensive for large numbers of field lines
- Use selective field line printing for large simulations
- Validation should be run only during setup, not every timestep

ERROR HANDLING:
---------------
- Functions check for null pointers and invalid parameters
- MPI rank 0 handles most error reporting to avoid message flooding
- Validation function returns boolean success/failure
- Detailed error messages include segment and field line indices

AUTHOR: Generated for SEP simulation framework
DATE: 2025
VERSION: 1.0

================================================================================
*/

#include "sep.h"
#include "polarized_wave_native_adapter.h"
#include <climits>
#include <stdexcept>

namespace SEP {
namespace AlfvenTurbulence_Kolmogorov {

// ============================================================================
// CORE INITIALIZATION FUNCTIONS
// ============================================================================

void InitializeWaveEnergyInAllSegments(
    PIC::Datum::cDatumStored& WaveEnergy,
    double wave_energy_density_1AU,
    double reference_distance,
    bool verbose
) {
    SEP::Turbulence::Polarized::PICAdapter::RequireLegacyMutationAllowed("legacy radial integrated wave initialization");
    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &rank);
    
    int initialized_segments = 0;
    double total_wave_energy_initialized = 0.0;
    
    if (rank == 0 && verbose) {
        std::cout << "Initializing wave energy in all field line segments..." << std::endl;
        std::cout << "Reference wave energy density: " << wave_energy_density_1AU << " J/m³ at " 
                  << reference_distance/WaveEnergyConstants::ONE_AU << " AU" << std::endl;
    }
    
    // ========================================================================
    // LOOP THROUGH ALL FIELD LINES AND SEGMENTS
    // ========================================================================
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        // Loop through all segments in this field line
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            
            // Only process segments assigned to this MPI process
            if (!segment || segment->Thread != PIC::ThisThread) {
                continue;
            }
            
            // ================================================================
            // CALCULATE SEGMENT CENTER POSITION
            // ================================================================
            
            // Get segment vertex positions
            double* x_left = segment->GetBegin()->GetX();
            double* x_right = segment->GetEnd()->GetX();
            
            // Calculate segment center position
            double x_center[3] = {
                0.5 * (x_left[0] + x_right[0]),
                0.5 * (x_left[1] + x_right[1]),
                0.5 * (x_left[2] + x_right[2])
            };
            
            // Calculate heliospheric distance (distance from origin)
            double r_helio = CalculateHeliosphericDistance(x_center);
            
            if (r_helio <= 0.0) {
                std::cerr << "Warning: Invalid heliospheric distance (" << r_helio 
                          << ") for segment " << seg_idx << " in field line " 
                          << field_line_idx << std::endl;
                continue;
            }
            
            // ================================================================
            // CALCULATE WAVE ENERGY DENSITY WITH r^-2 SCALING
            // ================================================================
            
            // Scale wave energy density as r^-2
            double local_wave_energy_density = ApplyHeliosphericScaling(
                wave_energy_density_1AU, reference_distance, r_helio
            );
            
            // ================================================================
            // CALCULATE SEGMENT VOLUME
            // ================================================================
            
            double V_segment = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, field_line_idx);
            
            if (V_segment <= 0.0) {
                std::cerr << "Warning: Invalid segment volume (" << V_segment 
                          << ") for segment " << seg_idx << " in field line " 
                          << field_line_idx << std::endl;
                continue;
            }
            
            // ================================================================
            // CALCULATE INTEGRATED WAVE ENERGIES
            // ================================================================
            
            double E_plus, E_minus;
            ConvertDensityToIntegratedEnergy(local_wave_energy_density, V_segment, E_plus, E_minus);
            
            // ================================================================
            // STORE WAVE ENERGY DATA IN SEGMENT
            // ================================================================
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (wave_data) {
                wave_data[0] = E_plus;   // Outward wave energy [J]
                wave_data[1] = E_minus;  // Inward wave energy [J]
                
                initialized_segments++;
                total_wave_energy_initialized += (E_plus + E_minus);
                
                if (verbose && rank == 0 && initialized_segments <= 10) {
                    std::cout << "  Segment " << seg_idx << " (FL " << field_line_idx << "): "
                              << "r=" << r_helio/WaveEnergyConstants::ONE_AU << " AU, "
                              << "ε=" << local_wave_energy_density << " J/m³, "
                              << "V=" << V_segment << " m³, "
                              << "E+=" << E_plus << " J, "
                              << "E-=" << E_minus << " J, "
		              << "E+/Vol=" << E_plus/V_segment << " J/m³, "
                              << "E-/Vol=" << E_minus/V_segment << " J/m³" << std::endl;
                }
            } else {
                std::cerr << "Error: Could not access wave energy datum for segment " 
                          << seg_idx << " in field line " << field_line_idx << std::endl;
            }
        }
    }
    
    // ========================================================================
    // MPI REDUCTION FOR GLOBAL STATISTICS
    // ========================================================================
    
    int total_initialized_segments = 0;
    double global_total_wave_energy = 0.0;
    
    MPI_Allreduce(&initialized_segments, &total_initialized_segments, 1, MPI_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&total_wave_energy_initialized, &global_total_wave_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    
    if (rank == 0) {
        std::cout << "Wave energy initialization complete:" << std::endl;
        std::cout << "  Total segments initialized: " << total_initialized_segments << std::endl;
        std::cout << "  Total wave energy in system: " << global_total_wave_energy << " J" << std::endl;
        if (total_initialized_segments > 0) {
            std::cout << "  Average energy per segment: " << global_total_wave_energy / total_initialized_segments << " J" << std::endl;
        }
    }
}

// ============================================================================
// CONVENIENCE WRAPPER WITH DEFAULT 1 AU REFERENCE
// ============================================================================

void InitializeWaveEnergyInAllSegments(
    PIC::Datum::cDatumStored& WaveEnergy,
    double wave_energy_density_1AU,
    bool verbose
) {
    InitializeWaveEnergyInAllSegments(WaveEnergy, wave_energy_density_1AU, WaveEnergyConstants::ONE_AU, verbose);
}

// ============================================================================
// PHYSICS-BASED INITIALIZATION
// ============================================================================
/// ********************************************************************************************
/// InitializeWaveEnergyFromPhysicalParameters
///
/// PURPOSE
/// --------
/// Initialize the total integrated wave energies E⁺ and E⁻ in each field–line segment for the
/// SEP turbulence module. The initialization is based on:
///   - a model for the large–scale background magnetic field B(r),
///   - a prescribed turbulence–level profile δB/B along the field line,
///   - and a simple balanced–turbulence assumption W⁺ = W⁻ = ½ W_total.
///
/// This function writes, for each field–line segment,
///   wave_data[0] = E⁺  [J] = W⁺(r) * V_segment
///   wave_data[1] = E⁻  [J] = W⁻(r) * V_segment
/// where V_segment is the volume associated with that segment.
///
///
/// PHYSICS MODEL
/// -------------
/// 1. Magnetic field model:
///    - If UseScaledReferenceB == true:
///        The local magnetic field magnitude at radius r is taken from a simple heliospheric
///        scaling of a reference magnitude B0_1AU:
///
///            B_local(r) = ApplyHeliosphericScaling(B0_1AU, 1 AU, r),
///
///        which typically encodes a Parker–spiral or WKB–like radial scaling.
///    - If UseScaledReferenceB == false:
///        The local magnetic field is taken from the field–line data stored at the segment
///        vertices:
///
///            B_local = |(B_vertex_begin + B_vertex_end) / 2|,
///
///        falling back to the scaled reference model if vertex data are missing.
///
/// 2. Turbulence level (δB/B) along the field line:
///    - The turbulence level at the *beginning* and *end* of the field line are specified:
///
///            turbulence_level_beginning_fl = (δB/B)_begin
///            turbulence_level_end_fl       = (δB/B)_end
///
///    - A radial coordinate r is built from the positions of the first and last vertices:
///            r_begin = |x(begin)|, r_end = |x(end)|.
///
///    - At each segment midpoint (radius r_helio) we compute a normalized coordinate
///
///            f = (r_helio - r_begin) / (r_end - r_begin),  f ∈ [0,1],
///
///      and use a power–shaped interpolation controlled by turbulence_level_decay_power_index:
///
///         If p = turbulence_level_decay_power_index > 0:
///           w_begin = (1 - f)^p,  w_end = f^p
///           turbulence_level_local = ( (δB/B)_begin * w_begin + (δB/B)_end * w_end )
///                                   / (w_begin + w_end)
///
///         If p <= 0:
///           simple linear interpolation is used:
///              turbulence_level_local = (δB/B)_begin
///                                      + f * [ (δB/B)_end - (δB/B)_begin ].
///
///      This gives a flexible 1D profile of δB/B between the two endpoints.
///
/// 3. Wave energy density and total energy:
///    - Given B_local and turbulence_level_local = δB/B, we define
///
///            δB = turbulence_level_local * B_local
///            W_total = δB^2 / (2 μ0),
///
///      where μ0 = VacuumPermeability. This is a standard expression for magnetic–energy
///      density in Alfvenic fluctuations.
///    - We assume balanced turbulence:
///
///            W⁺ = W⁻ = ½ W_total.
///
///    - For each segment, we obtain its volume via:
///
///            V_segment = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, iLine),
///
///      and store
///
///            E⁺ = W⁺ * V_segment,
///            E⁻ = W⁻ * V_segment
///
///      into the provided WaveEnergy datum.
///
///
/// USE CASES
/// ---------
/// - Quickly setting up a *controlled background turbulence profile* for tests of SEP transport
///   (e.g., varying δB/B between the inner and outer heliosphere).
/// - Matching approximate analytical expectations for turbulence strength near the Sun and at 1 AU.
/// - Performing parameter scans over:
///     * B0_1AU,
///     * turbulence_level_beginning_fl / turbulence_level_end_fl,
///     * turbulence_level_decay_power_index,
///   to study how the initial wave field modifies SEP propagation.
///
/// NOTE
/// ----
/// This is an initialization routine only. Subsequent evolution of W⁺ and W⁻ should be governed
/// by the turbulence cascade, advection, reflection, and wave–particle coupling modules.
///
/// ********************************************************************************************
SEP::Transport::Status InitializeWaveEnergyFromPhysicalParametersChecked(
    PIC::Datum::cDatumStored& WaveEnergy,double B0_1AU,
    double beginLevel,double endLevel,double power,bool scaledReference,bool verbose) {
  namespace FL=PIC::FieldLine;
  using SEP::Transport::Status;using SEP::Transport::StatusCode;
  struct Row {double* target;int line,segment;double radiusAu,fieldT,volumeM3,level,energyJ;};
  std::vector<Row> staged;std::vector<double> physical,minimum,maximum;
  Status local=SEP::Turbulence::Polarized::PICAdapter::LegacyMutationPreflight("physical wave initializer");
  auto invalid=[&](const std::string& text) {local=Status::Error(StatusCode::InvalidArgument,text);};
  auto finiteVector=[](const double* x) {return x&&std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(x[2]);};
  auto magnitude=[](const double* x) {return std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);};
  try {
    if(local.ok()&&(WaveEnergy.length!=2||!WaveEnergy.is_active()||!std::isfinite(beginLevel)||beginLevel<0||
        !std::isfinite(endLevel)||endLevel<0||!std::isfinite(power)||
        (scaledReference&&(!std::isfinite(B0_1AU)||B0_1AU<=0))))
      invalid("physical wave initializer requires a two-value active datum, finite nonnegative amplitudes and a declared valid B model");
    for(int line=0;local.ok()&&line<FL::nFieldLine;++line) {
      auto& fieldLine=FL::FieldLinesAll[line];const int count=fieldLine.GetTotalSegmentNumber();
      if(count<=0){invalid("physical wave initializer encountered an empty field line");break;}
      auto* first=fieldLine.GetSegment(0);auto* last=fieldLine.GetSegment(count-1);
      if(!first||!last||!first->GetBegin()||!last->GetEnd()) {invalid("physical initializer radial boundary vertices are absent");break;}
      double start[3],finish[3];first->GetBegin()->GetX(start);last->GetEnd()->GetX(finish);
      const double rBegin=finiteVector(start)?magnitude(start):0,rEnd=finiteVector(finish)?magnitude(finish):0;
      if(!std::isfinite(rBegin)||rBegin<=0||!std::isfinite(rEnd)||rEnd<=rBegin) {
        invalid("physical initializer requires finite positive radially ordered line endpoints; no fictitious radial offset is allowed");break;
      }
      for(int index=0;local.ok()&&index<count;++index) {
        auto* segment=fieldLine.GetSegment(index);
        if(!segment||!segment->GetBegin()||!segment->GetEnd()||
            !std::isfinite(segment->GetLength())||segment->GetLength()<=0) {
          invalid("physical initializer segment/vertices/positive length are unavailable");break;
        }
        double x0[3],x1[3],mid[3];segment->GetBegin()->GetX(x0);segment->GetEnd()->GetX(x1);
        if(!finiteVector(x0)||!finiteVector(x1)){invalid("physical initializer vertex position is nonfinite");break;}
        double displacement[3];for(int d=0;d<3;++d)displacement[d]=x1[d]-x0[d];
        const double actualLength=magnitude(displacement);
        if(!std::isfinite(actualLength)||actualLength<=0||std::fabs(actualLength-segment->GetLength())>
            64*std::numeric_limits<double>::epsilon()*std::max(actualLength,segment->GetLength())) {
          invalid("physical initializer cached segment length disagrees with actual vertex geometry");break;
        }
        for(int d=0;d<3;++d)mid[d]=0.5*x0[d]+0.5*x1[d];
        const double radius=magnitude(mid),fraction=(radius-rBegin)/(rEnd-rBegin);
        if(!std::isfinite(radius)||radius<=0||!std::isfinite(fraction)||fraction<0||fraction>1) {
          invalid("physical initializer midpoint lies outside its declared radial interpolation domain; no clipping/extrapolation is allowed");break;
        }
        // p<=0 is the documented linear profile family, not an error fallback.
        // For p>0, evaluate the same normalized endpoint blend with scaled
        // exponentials so a large exponent cannot trigger an arithmetic mean.
        long double profile=0;
        if(power<=0)profile=static_cast<long double>(beginLevel)+fraction*(static_cast<long double>(endLevel)-beginLevel);
        else if(fraction==0)profile=beginLevel;
        else if(fraction==1)profile=endLevel;
        else {
          const long double a=static_cast<long double>(power)*std::log1p(-static_cast<long double>(fraction));
          const long double b=static_cast<long double>(power)*std::log(static_cast<long double>(fraction));
          const long double scale=std::max(a,b),left=std::exp(a-scale),right=std::exp(b-scale);
          profile=(beginLevel*left+endLevel*right)/(left+right);
        }
        const double level=static_cast<double>(profile);
        const bool positiveProfile=fraction==0?beginLevel>0:fraction==1?endLevel>0:beginLevel>0||endLevel>0;
        if(!std::isfinite(level)||level<0||(positiveProfile&&level==0)) {
          invalid("physical initializer profile is nonfinite, negative or below representable positive resolution");break;
        }
        // Validate BOTH vertex fields, including geometry's midpoint field,
        // before asking the volume provider. Invalid B cannot select its area
        // callback or the old scaled-reference fallback implicitly.
        double* beginB=segment->GetBegin()->GetDatum_ptr(FL::DatumAtVertexMagneticField);
        double* endB=segment->GetEnd()->GetDatum_ptr(FL::DatumAtVertexMagneticField);
        double middleB[3];segment->GetMagneticField(0.5,middleB);
        if(!finiteVector(beginB)||!finiteVector(endB)||!finiteVector(middleB)||
            !std::isfinite(magnitude(beginB))||!std::isfinite(magnitude(endB))||!std::isfinite(magnitude(middleB))||
            !(magnitude(beginB)>0)||!(magnitude(endB)>0)||!(magnitude(middleB)>0)) {
          invalid("physical initializer requires finite positive endpoint/midpoint B; reference-field/area substitution is prohibited");break;
        }
        double field=0;
        if(scaledReference)field=B0_1AU*std::pow(WaveEnergyConstants::ONE_AU/radius,2);
        else if(index==0)field=magnitude(beginB);
        else if(index==count-1)field=magnitude(endB);
        else {
          double average[3];for(int d=0;d<3;++d)average[d]=0.5*beginB[d]+0.5*endB[d];
          field=magnitude(average);
        }
        const double volume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment,line);
        if(!std::isfinite(field)||field<=0||!std::isfinite(volume)||volume<=0) {
          invalid("physical initializer B/flux-tube volume is nonfinite or nonpositive");break;
        }
        const long double delta=static_cast<long double>(level)*field;
        const long double branchEnergy=delta*delta*volume/(4*static_cast<long double>(VacuumPermeability));
        const double energy=static_cast<double>(branchEnergy);
        double* target=segment->GetDatum_ptr(WaveEnergy);
        if(!target||!std::isfinite(energy)||energy<0||(level>0&&energy==0)) {
          invalid("physical initializer target or finite positive magnetic energy resolution is unavailable");break;
        }
        staged.push_back({target,line,index,radius/WaveEnergyConstants::ONE_AU,field,volume,level,energy});
        physical.insert(physical.end(),{static_cast<double>(line),static_cast<double>(index),radius,field,volume,level,energy});
      }
    }
    if(local.ok()&&(staged.empty()||physical.size()>static_cast<std::size_t>(INT_MAX)))
      invalid("physical initializer has no cells or exceeds collective buffer dimensions");
    if(local.ok()){minimum.resize(physical.size());maximum.resize(physical.size());}
  } catch(const std::exception& error) {
    local=Status::Error(StatusCode::InvalidArgument,std::string("physical initializer preflight failed: ")+error.what());
  }
  int ready=local.ok()?1:0,allReady=0;
  MPI_Allreduce(&ready,&allReady,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!allReady)return local.ok()?Status::Error(StatusCode::InvalidArgument,
      "another MPI rank rejected physical wave initialization; no host values changed"):local;
  unsigned long long count=physical.size(),smallest=0,largest=0;
  MPI_Allreduce(&count,&smallest,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&count,&largest,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(smallest!=largest)return Status::Error(StatusCode::InvalidArgument,
      "MPI physical initializer cell namespaces differ; no host values changed");
  MPI_Allreduce(physical.data(),minimum.data(),static_cast<int>(count),MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(physical.data(),maximum.data(),static_cast<int>(count),MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(minimum!=maximum)return Status::Error(StatusCode::InvalidArgument,
      "MPI physical initializer geometry/B/volume/profile/energy differ; no host values changed");
  // Collective preflight succeeded. These writes cannot allocate or encounter
  // a later physical rejection; every target was validated before publication.
  for(const auto& row:staged){row.target[0]=row.energyJ;row.target[1]=row.energyJ;}
  if(verbose&&PIC::ThisThread==0) {
    std::cout<<"Physical magnetic wave initialization committed: segments="<<staged.size()
        <<"; Eplus=Eminus; no geometry/B/volume/profile fallback; deltaB endpoints="
        <<beginLevel<<","<<endLevel<<std::endl;
  }
  return Status::Ok();
}

void InitializeWaveEnergyFromPhysicalParameters(PIC::Datum::cDatumStored& WaveEnergy,
    double B0_1AU,double beginLevel,double endLevel,double power,bool scaledReference,bool verbose) {
  const auto status=InitializeWaveEnergyFromPhysicalParametersChecked(WaveEnergy,B0_1AU,
      beginLevel,endLevel,power,scaledReference,verbose);
  if(!status.ok())throw std::runtime_error(status.message);
}

// ============================================================================
// PHYSICAL PARAMETER CALCULATIONS
// ============================================================================

double CalculateTypicalWaveEnergyDensity1AU(double B0_1AU, double turbulence_level) {
    // Wave energy density: ε = (δB²)/(2μ₀) where δB = turbulence_level * B₀
    double delta_B = turbulence_level * B0_1AU;
    double wave_energy_density = (delta_B * delta_B) / (2.0 * WaveEnergyConstants::MU0);
    
    return wave_energy_density;  // [J/m³]
}

// ============================================================================
// DIAGNOSTIC AND MONITORING FUNCTIONS
// ============================================================================

void PrintWaveEnergyProfile(
    PIC::Datum::cDatumStored& WaveEnergy,
    const std::vector<int>& field_line_indices,
    int max_segments_print
) {
    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &rank);
    
    if (rank != 0) return; // Only rank 0 prints
    
    std::cout << "\n=== Wave Energy Profile ===" << std::endl;
    std::cout << "FieldLine  Segment  Distance[AU]  E+[J]        E-[J]        Total[J]     Density[J/m³]" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    
    // Determine which field lines to print
    std::vector<int> lines_to_print;
    if (field_line_indices.empty()) {
        // Print all field lines
        for (int i = 0; i < PIC::FieldLine::nFieldLine; ++i) {
            lines_to_print.push_back(i);
        }
    } else {
        lines_to_print = field_line_indices;
    }
    
    for (int field_line_idx : lines_to_print) {
        if (field_line_idx >= PIC::FieldLine::nFieldLine) continue;
        
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        int segments_printed = 0;
        
        for (int seg_idx = 0; seg_idx < num_segments && segments_printed < max_segments_print; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment) continue;
            
            // Calculate segment center distance
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
            
            // Get wave energy data
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (wave_data) {
                double E_plus = wave_data[0];
                double E_minus = wave_data[1];
                double E_total = E_plus + E_minus;
                
                double V_segment = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, field_line_idx);
                double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
                
                printf("%8d  %7d  %11.4f  %11.4e  %11.4e  %11.4e  %11.4e\n",
                       field_line_idx, seg_idx, r_AU, E_plus, E_minus, E_total, energy_density);
                
                segments_printed++;
            }
        }
        
        if (segments_printed >= max_segments_print && num_segments > max_segments_print) {
            std::cout << "    ... (" << (num_segments - max_segments_print) << " more segments)" << std::endl;
        }
    }
    
    std::cout << "=== End Profile ===" << std::endl << std::endl;
}

void PrintWaveEnergyInitializationSummary(PIC::Datum::cDatumStored& WaveEnergy) {
    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &rank);
    
    // Calculate local statistics
    int local_segment_count = 0;
    double local_total_energy = 0.0;
    double local_min_distance = 1e20;
    double local_max_distance = 0.0;
    double local_min_energy_density = 1e20;
    double local_max_energy_density = 0.0;
    
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) continue;
            
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            
            double E_plus = wave_data[0];
            double E_minus = wave_data[1];
            double E_total = E_plus + E_minus;
            
            double V_segment = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, field_line_idx);
            double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
            
            local_segment_count++;
            local_total_energy += E_total;
            local_min_distance = std::min(local_min_distance, r_helio);
            local_max_distance = std::max(local_max_distance, r_helio);
            local_min_energy_density = std::min(local_min_energy_density, energy_density);
            local_max_energy_density = std::max(local_max_energy_density, energy_density);
        }
    }
    
    // MPI reductions for global statistics
    int global_segment_count = 0;
    double global_total_energy = 0.0;
    double global_min_distance = 0.0;
    double global_max_distance = 0.0;
    double global_min_energy_density = 0.0;
    double global_max_energy_density = 0.0;
    
    MPI_Allreduce(&local_segment_count, &global_segment_count, 1, MPI_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local_total_energy, &global_total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local_min_distance, &global_min_distance, 1, MPI_DOUBLE, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local_max_distance, &global_max_distance, 1, MPI_DOUBLE, MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local_min_energy_density, &global_min_energy_density, 1, MPI_DOUBLE, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&local_max_energy_density, &global_max_energy_density, 1, MPI_DOUBLE, MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
    
    if (rank == 0) {
        std::cout << "\n=== Wave Energy Initialization Summary ===" << std::endl;
        std::cout << "Total segments: " << global_segment_count << std::endl;
        std::cout << "Total wave energy: " << std::scientific << global_total_energy << " J" << std::endl;
        
        if (global_segment_count > 0) {
            std::cout << "Average energy per segment: " << global_total_energy / global_segment_count << " J" << std::endl;
        }
        
        std::cout << "Distance range: " << global_min_distance/WaveEnergyConstants::ONE_AU 
                  << " - " << global_max_distance/WaveEnergyConstants::ONE_AU << " AU" << std::endl;
        std::cout << "Energy density range: " << global_min_energy_density 
                  << " - " << global_max_energy_density << " J/m³" << std::endl;
        std::cout << "==========================================" << std::endl << std::endl;
    }
}

bool ValidateWaveEnergyInitialization(
    PIC::Datum::cDatumStored& WaveEnergy,
    double expected_energy_1AU,
    double tolerance
) {
    int rank;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &rank);
    
    bool validation_passed = true;
    int error_count = 0;
    
    // Check segments near 1 AU for proper scaling
    for (int field_line_idx = 0; field_line_idx < PIC::FieldLine::nFieldLine; ++field_line_idx) {
        PIC::FieldLine::cFieldLine* field_line = &PIC::FieldLine::FieldLinesAll[field_line_idx];
        int num_segments = field_line->GetTotalSegmentNumber();
        
        for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
            PIC::FieldLine::cFieldLineSegment* segment = field_line->GetSegment(seg_idx);
            if (!segment || segment->Thread != PIC::ThisThread) continue;
            
            double* wave_data = segment->GetDatum_ptr(WaveEnergy);
            if (!wave_data) continue;
            
            double x_center[3];
            GetSegmentCenterPosition(segment, x_center);
            double r_helio = CalculateHeliosphericDistance(x_center);
            double r_AU = r_helio / WaveEnergyConstants::ONE_AU;
            
            double E_plus = wave_data[0];
            double E_minus = wave_data[1];
            
            // Check for negative energies
            if (E_plus < 0.0 || E_minus < 0.0) {
                if (rank == 0) {
                    std::cerr << "Error: Negative wave energy in segment " << seg_idx 
                              << " of field line " << field_line_idx << std::endl;
                }
                validation_passed = false;
                error_count++;
            }
            
            // Check E+ = E- assumption
            double energy_ratio = std::abs(E_plus - E_minus) / (E_plus + E_minus);
            if (energy_ratio > tolerance) {
                if (rank == 0) {
                    std::cerr << "Warning: E+ != E- in segment " << seg_idx 
                              << " of field line " << field_line_idx 
                              << " (ratio: " << energy_ratio << ")" << std::endl;
                }
            }
            
            // Check r^-2 scaling for segments near 1 AU
            if (r_AU > 0.8 && r_AU < 1.2) {  // Within 20% of 1 AU
                double V_segment = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, field_line_idx);
                double energy_density = ConvertIntegratedEnergyToDensity(E_plus, E_minus, V_segment);
                
                // Expected energy density with r^-2 scaling
                double expected_density = expected_energy_1AU * (1.0 / (r_AU * r_AU));
                double relative_error = std::abs(energy_density - expected_density) / expected_density;
                
                if (relative_error > tolerance) {
                    if (rank == 0) {
                        std::cerr << "Warning: Energy density scaling error in segment " << seg_idx
                                  << " (expected: " << expected_density 
                                  << ", actual: " << energy_density 
                                  << ", error: " << relative_error * 100 << "%)" << std::endl;
                    }
                }
            }
        }
    }
    
    // MPI reduction for global validation result
    int global_validation_passed = validation_passed ? 1 : 0;
    int global_validation_result = 0;
    MPI_Allreduce(&global_validation_passed, &global_validation_result, 1, MPI_INT, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
    
    if (rank == 0) {
        if (global_validation_result == 1) {
            std::cout << "Wave energy initialization validation: PASSED" << std::endl;
        } else {
            std::cout << "Wave energy initialization validation: FAILED" << std::endl;
        }
    }
    
    return (global_validation_result == 1);
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

double CalculateSegmentWaveEnergy(
    PIC::FieldLine::cFieldLineSegment* segment,
    PIC::Datum::cDatumStored& WaveEnergy
) {
    if (!segment) return 0.0;
    
    double* wave_data = segment->GetDatum_ptr(WaveEnergy);
    if (!wave_data) return 0.0;
    
    return wave_data[0] + wave_data[1];  // E+ + E-
}

void GetSegmentCenterPosition(
    PIC::FieldLine::cFieldLineSegment* segment,
    double x_center[3]
) {
    if (!segment) {
        x_center[0] = x_center[1] = x_center[2] = 0.0;
        return;
    }
    
    double* x_left = segment->GetBegin()->GetX();
    double* x_right = segment->GetEnd()->GetX();
    
    x_center[0] = 0.5 * (x_left[0] + x_right[0]);
    x_center[1] = 0.5 * (x_left[1] + x_right[1]);
    x_center[2] = 0.5 * (x_left[2] + x_right[2]);
}

} // namespace AlfvenTurbulence_Kolmogorov
} // namespace SEP
