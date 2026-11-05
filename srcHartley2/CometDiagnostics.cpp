#include "CometDiagnostics.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sys/stat.h>
#include <vector>

#include "pic.h"
#include "CometChemistry.h"
#include "CometModel.h"
#include "src/models/chemical_solver/ChemistryDiagnostics.h"

bool Comet::Diagnostics::Enabled=false;

namespace {

enum cFluxField {
  SurfaceInjection=0,
  ChemicalProductionExpected,
  ChemicalProductionRepresented,
  ChemicalProductionVariance,
  ChemicalDestruction,
  InnerBoundaryIncident,
  InnerBoundaryReemitted,
  OuterBoundaryEscape,
  FluxFieldCount
};

ChemicalSolver::cChemistryDiagnostics ChemistryCounters;
std::vector<double> SpeciesFlux;
std::vector<std::size_t> ProductBoundOffset;
std::vector<double> ProductSystematicBound;
std::vector<double> SpeciesSystematicBound;
std::uint64_t GlobalEmergencyRequests=0;
enum cEmergencyEventField {
  EmergencyEventTotal=0,
  EmergencyEventBufferUtilization,
  EmergencyEventSamplingPopulation,
  EmergencyEventAllocationBurst,
  EmergencyEventFieldCount
};
std::uint64_t GlobalEmergencyEvents[EmergencyEventFieldCount]={0,0,0,0};
bool Initialized=false;

double& Flux(int species,cFluxField field) {
  if (species<0 || species>=PIC::nTotalSpecies) {
    std::fprintf(stderr,"invalid comet diagnostic species %d\n",species);
    std::abort();
  }
  return SpeciesFlux[static_cast<std::size_t>(species)*FluxFieldCount+field];
}

void EnsureDirectory() {
  if (PIC::ThisThread==0) {
    if (mkdir("PT",0775)!=0 && errno!=EEXIST) {
      std::fprintf(stderr,"cannot create PT for comet diagnostics: %s\n",
          std::strerror(errno));
      std::abort();
    }
    if (mkdir("PT/diagnostics",0775)!=0 && errno!=EEXIST) {
      std::fprintf(stderr,
          "cannot create PT/diagnostics for comet diagnostics: %s\n",
          std::strerror(errno));
      std::abort();
    }
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void SpeciesSymbol(char *buffer,int species) {
  PIC::MolecularData::GetChemSymbol(buffer,species);
}

double Average(double accumulated,long int sampleLength) {
  return accumulated/static_cast<double>(sampleLength);
}

bool ProductionAgreement(double expected,double represented,double variance,
    double systematicBound,long int sampleLength,double &sigma,
    double &reportedBound,double &tolerance) {
  const bool systematic=Comet::Chemistry::ResidualSamplingMode==
      Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
  sigma=systematic ? std::numeric_limits<double>::quiet_NaN() :
      std::sqrt(std::max(0.0,variance))/sampleLength;
  reportedBound=systematic ? systematicBound/sampleLength : 0.0;
  tolerance=std::max(systematic ? reportedBound : 5.0*sigma,
      0.01*std::fabs(expected));
  return std::fabs(represented-expected)<=tolerance;
}

} // namespace

void Comet::Diagnostics::Initialize() {
  if (!Enabled) {
    ChemistryCounters.Configure(std::vector<std::size_t>(),false);
    SpeciesFlux.clear();
    ProductBoundOffset.clear();
    ProductSystematicBound.clear();
    SpeciesSystematicBound.clear();
    GlobalEmergencyRequests=0;
    std::fill(GlobalEmergencyEvents,
        GlobalEmergencyEvents+EmergencyEventFieldCount,0);
    Initialized=false;
    return;
  }

  const ChemicalSolver::cReactionNetwork &network=Comet::Chemistry::GetNetwork();
  std::vector<std::size_t> productCounts(network.GetReactionCount(),0);
  ProductBoundOffset.assign(network.GetReactionCount()+1,0);
  std::size_t productTermCount=0;
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
    ProductBoundOffset[reaction]=productTermCount;
    productCounts[reaction]=network.GetReaction(reaction).ProductCount;
    productTermCount+=productCounts[reaction];
  }
  ProductBoundOffset[network.GetReactionCount()]=productTermCount;
  ChemistryCounters.Configure(productCounts,true);
  ProductSystematicBound.assign(productTermCount,0.0);
  SpeciesFlux.assign(static_cast<std::size_t>(PIC::nTotalSpecies)*
      FluxFieldCount,0.0);
  SpeciesSystematicBound.assign(PIC::nTotalSpecies,0.0);
  GlobalEmergencyRequests=0;
  std::fill(GlobalEmergencyEvents,
      GlobalEmergencyEvents+EmergencyEventFieldCount,0);
  Initialized=true;
  EnsureDirectory();
}

void Comet::Diagnostics::RecordEligibleParentVisit(std::size_t reaction) {
  if (!Enabled) return;
  ChemistryCounters.RecordEligible(reaction);
}

void Comet::Diagnostics::RecordReactionSelection(std::size_t reaction,
    bool neutralSink) {
  if (!Enabled) return;
  ChemistryCounters.RecordSelection(reaction,neutralSink);
}

void Comet::Diagnostics::RecordAllocationDeferral(std::size_t reaction) {
  if (!Enabled) return;
  ChemistryCounters.RecordDeferral(reaction);
}

void Comet::Diagnostics::RecordEmergencyRequest() {
  if (!Enabled) return;
  // Allocation-burst and sampling-population requests are collective. Count
  // each request once globally, on rank zero, so the later MPI sum is
  // independent of process count. Pure buffer pressure has no separate
  // request stage and is captured by RecordEmergencyEvent instead.
  if (PIC::ThisThread==0) GlobalEmergencyRequests++;
}

void Comet::Diagnostics::RecordEmergencyEvent(bool bufferUtilization,
    bool samplingPopulation,bool allocationBurst) {
  if (!Enabled || PIC::ThisThread!=0) return;
  if (!bufferUtilization && !samplingPopulation && !allocationBurst) {
    std::fprintf(stderr,
        "comet diagnostic emergency event has no trigger\n");
    std::abort();
  }
  GlobalEmergencyEvents[EmergencyEventTotal]++;
  if (bufferUtilization)
    GlobalEmergencyEvents[EmergencyEventBufferUtilization]++;
  if (samplingPopulation)
    GlobalEmergencyEvents[EmergencyEventSamplingPopulation]++;
  if (allocationBurst)
    GlobalEmergencyEvents[EmergencyEventAllocationBurst]++;
}

void Comet::Diagnostics::RecordSystematicResidualBound(
    std::size_t reaction,std::size_t productTerm,int species,
    double carrierRate) {
  if (!Enabled) return;
  if (reaction+1>=ProductBoundOffset.size() ||
      productTerm>=ProductBoundOffset[reaction+1]-
          ProductBoundOffset[reaction] ||
      species<0 || species>=PIC::nTotalSpecies ||
      !(carrierRate>0.0) || !std::isfinite(carrierRate)) {
    std::fprintf(stderr,"invalid systematic comet residual diagnostic bound\n");
    std::abort();
  }
  ProductSystematicBound[ProductBoundOffset[reaction]+productTerm]+=
      carrierRate;
  SpeciesSystematicBound[species]+=carrierRate;
}

void Comet::Diagnostics::RecordChemistryProduct(std::size_t reaction,
    std::size_t productTerm,std::uint64_t modelMultiplicity,
    std::uint64_t allocatedSlots,std::uint64_t reusedSlots,
    double expectedNumberRate,double expectedMassRate,
    const double expectedMomentumRate[3],double expectedEnergyRate,
    double representedNumberRate,double representedMassRate,
    const double representedMomentumRate[3],double representedEnergyRate,
    double representedNumberRateVariance) {
  if (!Enabled) return;
  ChemistryCounters.RecordProduct(reaction,productTerm,modelMultiplicity,
      allocatedSlots,reusedSlots,expectedNumberRate,expectedMassRate,
      expectedMomentumRate,expectedEnergyRate,representedNumberRate,
      representedMassRate,representedMomentumRate,representedEnergyRate,
      representedNumberRateVariance);
}

void Comet::Diagnostics::RecordSurfaceInjection(int species,double rate) {
  if (!Enabled) return;
  Flux(species,SurfaceInjection)+=rate;
}

void Comet::Diagnostics::RecordChemicalProduction(int species,
    double expectedRate,double representedRate,double variance) {
  if (!Enabled) return;
  Flux(species,ChemicalProductionExpected)+=expectedRate;
  Flux(species,ChemicalProductionRepresented)+=representedRate;
  Flux(species,ChemicalProductionVariance)+=variance;
}

void Comet::Diagnostics::RecordChemicalDestruction(int species,double rate) {
  if (!Enabled) return;
  Flux(species,ChemicalDestruction)+=rate;
}

void Comet::Diagnostics::RecordInnerBoundary(int species,double rate) {
  if (!Enabled) return;
  // Every incident neutral is immediately thermalized and re-emitted by the
  // radial boundary callback.  Keep gross inward and outward contributions
  // distinct even though this boundary's net particle flux is exactly zero.
  Flux(species,InnerBoundaryIncident)+=rate;
  Flux(species,InnerBoundaryReemitted)+=rate;
}

void Comet::Diagnostics::RecordOuterBoundary(int species,double rate) {
  if (!Enabled) return;
  Flux(species,OuterBoundaryEscape)+=rate;
}

void Comet::Diagnostics::WriteSamplingReport(long int completedIteration,
    long int outputIndex,long int sampleLength,
    const std::vector<double> &relativeL1,bool hasPreviousProfile) {
  if (!Enabled) return;
  if (!Initialized || sampleLength<=0 ||
      relativeL1.size()!=static_cast<std::size_t>(PIC::nTotalSpecies)) {
    std::fprintf(stderr,"invalid comet diagnostic sampling report state\n");
    std::abort();
  }

  // Convert every owner block's active systematic accumulator into one
  // conservative fixed-carrier-epoch error bound before reducing diagnostic
  // tables. This does not reset the physical phase.
  Comet::Chemistry::CloseResidualDiagnosticWindow();

  const std::vector<std::uint64_t> localCounts=
      ChemistryCounters.CountValues();
  const std::vector<double> localMoments=ChemistryCounters.MomentValues();
  std::vector<std::uint64_t> globalCounts(localCounts.size(),0);
  std::vector<double> globalMoments(localMoments.size(),0.0);
  std::vector<double> globalFlux(SpeciesFlux.size(),0.0);
  std::vector<double> globalProductBounds(ProductSystematicBound.size(),0.0);
  std::vector<double> globalSpeciesBounds(SpeciesSystematicBound.size(),0.0);
  std::uint64_t globalEmergencies=0;
  std::uint64_t globalEmergencyEvents[EmergencyEventFieldCount]={0,0,0,0};
  MPI_Reduce(localCounts.data(),globalCounts.data(),
      static_cast<int>(localCounts.size()),MPI_UINT64_T,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMoments.data(),globalMoments.data(),
      static_cast<int>(localMoments.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(SpeciesFlux.data(),globalFlux.data(),
      static_cast<int>(SpeciesFlux.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(ProductSystematicBound.data(),globalProductBounds.data(),
      static_cast<int>(ProductSystematicBound.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(SpeciesSystematicBound.data(),globalSpeciesBounds.data(),
      static_cast<int>(SpeciesSystematicBound.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&GlobalEmergencyRequests,&globalEmergencies,1,MPI_UINT64_T,
      MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(GlobalEmergencyEvents,globalEmergencyEvents,
      EmergencyEventFieldCount,MPI_UINT64_T,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    char csvName[256],jsonName[256];
    std::snprintf(csvName,sizeof(csvName),
        "PT/diagnostics/comet_stage0.out=%ld.csv",outputIndex);
    std::snprintf(jsonName,sizeof(jsonName),
        "PT/diagnostics/comet_stage0.out=%ld.json",outputIndex);
    FILE *csv=std::fopen(csvName,"w");
    FILE *json=std::fopen(jsonName,"w");
    if (csv==NULL || json==NULL) {
      std::fprintf(stderr,"cannot open comet diagnostic report output\n");
      std::abort();
    }

    std::fprintf(csv,
        "record_type,iteration,output_index,sampling_length,reaction_index,reaction_id,product_term,species,stoichiometric_coefficient,eligible_parent_visits,reaction_selections,selected_channels,parent_destructions,neutral_sinks,allocation_deferrals,emergency_requests,multiplicity_zero,multiplicity_one,multiplicity_many,allocated_slots,reused_slots,model_carriers,expected_number_rate_s-1,represented_number_rate_s-1,represented_number_sigma_s-1,production_tolerance_s-1,production_agreement,expected_mass_rate_kg_s-1,represented_mass_rate_kg_s-1,expected_momentum_x_N,expected_momentum_y_N,expected_momentum_z_N,represented_momentum_x_N,represented_momentum_y_N,represented_momentum_z_N,expected_energy_rate_W,represented_energy_rate_W,configured_surface_source_s-1,sampled_surface_source_s-1,chemical_expected_production_s-1,chemical_represented_production_s-1,chemical_destruction_s-1,inner_incident_s-1,inner_reemitted_s-1,outer_escape_s-1,net_expected_flux_s-1,net_represented_flux_s-1,relative_l1_change,systematic_error_bound_s-1,uncertainty_method\n");
    std::fprintf(csv,
        "metadata,%ld,%ld,%ld,,,,,,,,,,,,%llu,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n",
        completedIteration,outputIndex,sampleLength,
        static_cast<unsigned long long>(globalEmergencies));

    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    std::size_t countOffset=0,momentOffset=0;
    std::fprintf(json,
        "{\n  \"schema\": \"amps.comet.stage0-diagnostics.v1\",\n"
        "  \"iteration\": %ld,\n  \"output_index\": %ld,\n"
        "  \"sampling_length\": %ld,\n"
        "  \"random_seed\": %llu,\n"
        "  \"random_seed_configured\": %s,\n"
        "  \"mpi_ranks\": %d,\n"
        "  \"residual_sampling\": \"%s\",\n"
        "  \"emergency_requests\": %llu,\n"
        "  \"emergency_events\": {\"total\": %llu, "
        "\"buffer_utilization\": %llu, \"sampling_population\": %llu, "
        "\"allocation_burst\": %llu},\n  \"reactions\": [\n",
        completedIteration,outputIndex,sampleLength,
        static_cast<unsigned long long>(
            Comet::Model::RandomSeedConfigured ?
            Comet::Model::RandomSeed : 0),
        Comet::Model::RandomSeedConfigured ? "true" : "false",
        PIC::nTotalThreads,
        Comet::Chemistry::ResidualSamplingMode==
            Comet::Chemistry::ResidualSampling::RandomPhaseSystematic ?
            "SYSTEMATIC" : "BERNOULLI",
        static_cast<unsigned long long>(globalEmergencies),
        static_cast<unsigned long long>(
            globalEmergencyEvents[EmergencyEventTotal]),
        static_cast<unsigned long long>(
            globalEmergencyEvents[EmergencyEventBufferUtilization]),
        static_cast<unsigned long long>(
            globalEmergencyEvents[EmergencyEventSamplingPopulation]),
        static_cast<unsigned long long>(
            globalEmergencyEvents[EmergencyEventAllocationBurst]));
    for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++) {
      const ChemicalSolver::cReaction &entry=network.GetReaction(reaction);
      const std::uint64_t eligible=globalCounts[countOffset++];
      const std::uint64_t selections=globalCounts[countOffset++];
      const std::uint64_t channels=globalCounts[countOffset++];
      const std::uint64_t destructions=globalCounts[countOffset++];
      const std::uint64_t sinks=globalCounts[countOffset++];
      const std::uint64_t deferrals=globalCounts[countOffset++];
      const std::uint64_t reactionEmergencies=globalCounts[countOffset++];
      std::fprintf(csv,
          "reaction,%ld,%ld,%ld,%zu,%s,,,,%llu,%llu,%llu,%llu,%llu,%llu,%llu,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n",
          completedIteration,outputIndex,sampleLength,reaction,entry.Id,
          static_cast<unsigned long long>(eligible),
          static_cast<unsigned long long>(selections),
          static_cast<unsigned long long>(channels),
          static_cast<unsigned long long>(destructions),
          static_cast<unsigned long long>(sinks),
          static_cast<unsigned long long>(deferrals),
          static_cast<unsigned long long>(reactionEmergencies));
      std::fprintf(json,
          "    {\"index\": %zu, \"id\": \"%s\", \"eligible_parent_visits\": %llu, "
          "\"reaction_selections\": %llu, \"selected_channels\": %llu, "
          "\"parent_destructions\": %llu, \"neutral_sinks\": %llu, "
          "\"allocation_deferrals\": %llu, \"products\": [",
          reaction,entry.Id,static_cast<unsigned long long>(eligible),
          static_cast<unsigned long long>(selections),
          static_cast<unsigned long long>(channels),
          static_cast<unsigned long long>(destructions),
          static_cast<unsigned long long>(sinks),
          static_cast<unsigned long long>(deferrals));

      for (std::size_t product=0;product<entry.ProductCount;product++) {
        const ChemicalSolver::cStoichiometricTerm &term=
            network.GetProductTerm(reaction,product);
        char symbol[_MAX_STRING_LENGTH_PIC_];
        SpeciesSymbol(symbol,term.Species);
        const std::uint64_t multiplicityZero=globalCounts[countOffset++];
        const std::uint64_t multiplicityOne=globalCounts[countOffset++];
        const std::uint64_t multiplicityMany=globalCounts[countOffset++];
        const std::uint64_t allocated=globalCounts[countOffset++];
        const std::uint64_t reused=globalCounts[countOffset++];
        const std::uint64_t carriers=globalCounts[countOffset++];
        const double expectedNumber=Average(globalMoments[momentOffset++],sampleLength);
        const double expectedMass=Average(globalMoments[momentOffset++],sampleLength);
        const double expectedMomentum[3]={
            Average(globalMoments[momentOffset++],sampleLength),
            Average(globalMoments[momentOffset++],sampleLength),
            Average(globalMoments[momentOffset++],sampleLength)};
        const double expectedEnergy=Average(globalMoments[momentOffset++],sampleLength);
        const double representedNumber=Average(globalMoments[momentOffset++],sampleLength);
        const double representedMass=Average(globalMoments[momentOffset++],sampleLength);
        const double representedMomentum[3]={
            Average(globalMoments[momentOffset++],sampleLength),
            Average(globalMoments[momentOffset++],sampleLength),
            Average(globalMoments[momentOffset++],sampleLength)};
        const double representedEnergy=Average(globalMoments[momentOffset++],sampleLength);
        const double variance=globalMoments[momentOffset++];
        double sigma=0.0,systematicBound=0.0,tolerance=0.0;
        const bool agreement=ProductionAgreement(expectedNumber,
            representedNumber,variance,
            globalProductBounds[ProductBoundOffset[reaction]+product],
            sampleLength,sigma,systematicBound,tolerance);
        const bool systematic=Comet::Chemistry::ResidualSamplingMode==
            Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
        char sigmaText[64]="",boundText[64]="";
        if (systematic)
          std::snprintf(boundText,sizeof(boundText),"%.17e",systematicBound);
        else
          std::snprintf(sigmaText,sizeof(sigmaText),"%.17e",sigma);
        std::fprintf(csv,
            "product,%ld,%ld,%ld,%zu,%s,%zu,%s,%.17e,,,,,,,,%llu,%llu,%llu,%llu,%llu,%llu,%.17e,%.17e,%s,%.17e,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,,,,,,,,,,,,%s,%s\n",
            completedIteration,outputIndex,sampleLength,reaction,entry.Id,
            product,symbol,term.Coefficient,
            static_cast<unsigned long long>(multiplicityZero),
            static_cast<unsigned long long>(multiplicityOne),
            static_cast<unsigned long long>(multiplicityMany),
            static_cast<unsigned long long>(allocated),
            static_cast<unsigned long long>(reused),
            static_cast<unsigned long long>(carriers),expectedNumber,
            representedNumber,sigmaText,tolerance,
            agreement ? "true" : "false",
            expectedMass,representedMass,expectedMomentum[0],
            expectedMomentum[1],expectedMomentum[2],representedMomentum[0],
            representedMomentum[1],representedMomentum[2],expectedEnergy,
            representedEnergy,boundText,systematic ?
                "fixed-carrier-epoch-bound" :
                "independent-bernoulli-analytic");
        std::fprintf(json,
            "%s{\"term\": %zu, \"species\": \"%s\", \"coefficient\": %.17g, "
            "\"multiplicity_zero\": %llu, \"multiplicity_one\": %llu, "
            "\"multiplicity_many\": %llu, \"allocated_slots\": %llu, "
            "\"reused_slots\": %llu, \"model_carriers\": %llu, "
            "\"expected_number_rate_s-1\": %.17e, "
            "\"represented_number_rate_s-1\": %.17e, ",
            product==0 ? "" : ", ",product,symbol,term.Coefficient,
            static_cast<unsigned long long>(multiplicityZero),
            static_cast<unsigned long long>(multiplicityOne),
            static_cast<unsigned long long>(multiplicityMany),
            static_cast<unsigned long long>(allocated),
            static_cast<unsigned long long>(reused),
            static_cast<unsigned long long>(carriers),expectedNumber,
            representedNumber);
        if (Comet::Chemistry::ResidualSamplingMode==
            Comet::Chemistry::ResidualSampling::RandomPhaseSystematic)
          std::fprintf(json,
              "\"represented_number_sigma_s-1\": null, "
              "\"systematic_error_bound_s-1\": %.17e, "
              "\"uncertainty_method\": \"fixed-carrier-epoch-bound\", ",
              systematicBound);
        else
          std::fprintf(json,
              "\"represented_number_sigma_s-1\": %.17e, "
              "\"systematic_error_bound_s-1\": null, "
              "\"uncertainty_method\": \"independent-bernoulli-analytic\", ",
              sigma);
        std::fprintf(json,
            "\"production_tolerance_s-1\": %.17e, "
            "\"production_agreement\": %s}",tolerance,
            agreement ? "true" : "false");
      }
      std::fprintf(json,"]}%s\n",
          reaction+1==network.GetReactionCount() ? "" : ",");
    }

    std::fprintf(json,"  ],\n  \"species_flux\": [\n");
    for (int species=0;species<PIC::nTotalSpecies;species++) {
      char symbol[_MAX_STRING_LENGTH_PIC_];
      SpeciesSymbol(symbol,species);
      const std::size_t base=static_cast<std::size_t>(species)*FluxFieldCount;
      const double surface=Average(globalFlux[base+SurfaceInjection],sampleLength);
      const double expectedProduction=Average(
          globalFlux[base+ChemicalProductionExpected],sampleLength);
      const double representedProduction=Average(
          globalFlux[base+ChemicalProductionRepresented],sampleLength);
      const double destruction=Average(globalFlux[base+ChemicalDestruction],sampleLength);
      const double innerIncident=Average(globalFlux[base+InnerBoundaryIncident],sampleLength);
      const double innerReemitted=Average(globalFlux[base+InnerBoundaryReemitted],sampleLength);
      const double outerEscape=Average(globalFlux[base+OuterBoundaryEscape],sampleLength);
      const double netExpected=surface+expectedProduction+innerReemitted-
          destruction-innerIncident-outerEscape;
      const double netRepresented=surface+representedProduction+innerReemitted-
          destruction-innerIncident-outerEscape;
      double sigma=0.0,systematicBound=0.0,tolerance=0.0;
      const bool agreement=ProductionAgreement(expectedProduction,
          representedProduction,globalFlux[base+ChemicalProductionVariance],
          globalSpeciesBounds[species],sampleLength,sigma,systematicBound,
          tolerance);
      const bool systematic=Comet::Chemistry::ResidualSamplingMode==
          Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
      char sigmaText[64]="",boundText[64]="";
      if (systematic)
        std::snprintf(boundText,sizeof(boundText),"%.17e",systematicBound);
      else
        std::snprintf(sigmaText,sizeof(sigmaText),"%.17e",sigma);
      std::fprintf(csv,
          "species_flux,%ld,%ld,%ld,,,,%s,,,,,,,,,,,,,,,%.17e,%.17e,%s,%.17e,%s,,,,,,,,,,,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,",
          completedIteration,outputIndex,sampleLength,symbol,
          expectedProduction,representedProduction,sigmaText,tolerance,
          agreement ? "true" : "false",
          Comet::Model::SurfaceProductionRate[species],surface,
          expectedProduction,representedProduction,destruction,innerIncident,
          innerReemitted,outerEscape,netExpected,netRepresented);
      if (hasPreviousProfile) std::fprintf(csv,"%.17e",relativeL1[species]);
      std::fprintf(csv,",%s,%s\n",boundText,systematic ?
          "fixed-carrier-epoch-bound" :
          "independent-bernoulli-analytic");
      std::fprintf(json,
          "    {\"species\": \"%s\", \"configured_surface_source_s-1\": %.17e, "
          "\"sampled_surface_source_s-1\": %.17e, "
          "\"chemical_expected_production_s-1\": %.17e, "
          "\"chemical_represented_production_s-1\": %.17e, ",
          symbol,Comet::Model::SurfaceProductionRate[species],surface,
          expectedProduction,representedProduction);
      if (Comet::Chemistry::ResidualSamplingMode==
          Comet::Chemistry::ResidualSampling::RandomPhaseSystematic)
        std::fprintf(json,
          "\"chemical_production_sigma_s-1\": null, "
          "\"systematic_error_bound_s-1\": %.17e, "
          "\"uncertainty_method\": \"fixed-carrier-epoch-bound\", ",
          systematicBound);
      else
        std::fprintf(json,
          "\"chemical_production_sigma_s-1\": %.17e, "
          "\"systematic_error_bound_s-1\": null, "
          "\"uncertainty_method\": \"independent-bernoulli-analytic\", ",
          sigma);
      std::fprintf(json,
          "\"production_tolerance_s-1\": %.17e, "
          "\"production_agreement\": %s, "
          "\"chemical_destruction_s-1\": %.17e, "
          "\"inner_incident_s-1\": %.17e, \"inner_reemitted_s-1\": %.17e, "
          "\"outer_escape_s-1\": %.17e, \"net_expected_flux_s-1\": %.17e, "
          "\"net_represented_flux_s-1\": %.17e, \"relative_l1_change\": ",
          tolerance,agreement ? "true" : "false",destruction,innerIncident,
          innerReemitted,outerEscape,netExpected,netRepresented);
      if (hasPreviousProfile) std::fprintf(json,"%.17e",relativeL1[species]);
      else std::fprintf(json,"null");
      std::fprintf(json,"}%s\n",species+1==PIC::nTotalSpecies ? "" : ",");
    }
    std::fprintf(json,"  ]\n}\n");
    std::fclose(csv);
    std::fclose(json);
    std::fprintf(PIC::DiagnospticMessageStream,
        "$PREFIX: comet diagnostics wrote output=%ld iteration=%ld sampling-length=%ld emergencies=%llu\n",
        outputIndex,completedIteration,sampleLength,
        static_cast<unsigned long long>(globalEmergencies));
    std::fflush(PIC::DiagnospticMessageStream);
  }

  ChemistryCounters.Reset();
  std::fill(SpeciesFlux.begin(),SpeciesFlux.end(),0.0);
  std::fill(ProductSystematicBound.begin(),ProductSystematicBound.end(),0.0);
  std::fill(SpeciesSystematicBound.begin(),SpeciesSystematicBound.end(),0.0);
  GlobalEmergencyRequests=0;
  std::fill(GlobalEmergencyEvents,
      GlobalEmergencyEvents+EmergencyEventFieldCount,0);
}
