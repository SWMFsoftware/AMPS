#include "CometDiagnostics.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
std::uint64_t GlobalEmergencyRequests=0;
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
    long int sampleLength,double &sigma,double &tolerance) {
  sigma=std::sqrt(std::max(0.0,variance))/sampleLength;
  tolerance=std::max(5.0*sigma,0.01*std::fabs(expected));
  return std::fabs(represented-expected)<=tolerance;
}

} // namespace

void Comet::Diagnostics::Initialize() {
  if (!Enabled) {
    ChemistryCounters.Configure(std::vector<std::size_t>(),false);
    SpeciesFlux.clear();
    GlobalEmergencyRequests=0;
    Initialized=false;
    return;
  }

  const ChemicalSolver::cReactionNetwork &network=Comet::Chemistry::GetNetwork();
  std::vector<std::size_t> productCounts(network.GetReactionCount(),0);
  for (std::size_t reaction=0;reaction<network.GetReactionCount();reaction++)
    productCounts[reaction]=network.GetReaction(reaction).ProductCount;
  ChemistryCounters.Configure(productCounts,true);
  SpeciesFlux.assign(static_cast<std::size_t>(PIC::nTotalSpecies)*
      FluxFieldCount,0.0);
  GlobalEmergencyRequests=0;
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
  // RequestAllocationBurstEmergency is collective.  Count it once globally,
  // on rank zero, so the later MPI sum is independent of process count.
  if (PIC::ThisThread==0) GlobalEmergencyRequests++;
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

  const std::vector<std::uint64_t> localCounts=
      ChemistryCounters.CountValues();
  const std::vector<double> localMoments=ChemistryCounters.MomentValues();
  std::vector<std::uint64_t> globalCounts(localCounts.size(),0);
  std::vector<double> globalMoments(localMoments.size(),0.0);
  std::vector<double> globalFlux(SpeciesFlux.size(),0.0);
  std::uint64_t globalEmergencies=0;
  MPI_Reduce(localCounts.data(),globalCounts.data(),
      static_cast<int>(localCounts.size()),MPI_UINT64_T,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localMoments.data(),globalMoments.data(),
      static_cast<int>(localMoments.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(SpeciesFlux.data(),globalFlux.data(),
      static_cast<int>(SpeciesFlux.size()),MPI_DOUBLE,MPI_SUM,0,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(&GlobalEmergencyRequests,&globalEmergencies,1,MPI_UINT64_T,
      MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

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
        "record_type,iteration,output_index,sampling_length,reaction_index,reaction_id,product_term,species,stoichiometric_coefficient,eligible_parent_visits,reaction_selections,selected_channels,parent_destructions,neutral_sinks,allocation_deferrals,emergency_requests,multiplicity_zero,multiplicity_one,multiplicity_many,allocated_slots,reused_slots,model_carriers,expected_number_rate_s-1,represented_number_rate_s-1,represented_number_sigma_s-1,production_tolerance_s-1,production_agreement,expected_mass_rate_kg_s-1,represented_mass_rate_kg_s-1,expected_momentum_x_N,expected_momentum_y_N,expected_momentum_z_N,represented_momentum_x_N,represented_momentum_y_N,represented_momentum_z_N,expected_energy_rate_W,represented_energy_rate_W,configured_surface_source_s-1,sampled_surface_source_s-1,chemical_expected_production_s-1,chemical_represented_production_s-1,chemical_destruction_s-1,inner_incident_s-1,inner_reemitted_s-1,outer_escape_s-1,net_expected_flux_s-1,net_represented_flux_s-1,relative_l1_change\n");
    std::fprintf(csv,
        "metadata,%ld,%ld,%ld,,,,,,,,,,,,%llu,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n",
        completedIteration,outputIndex,sampleLength,
        static_cast<unsigned long long>(globalEmergencies));

    const ChemicalSolver::cReactionNetwork &network=
        Comet::Chemistry::GetNetwork();
    std::size_t countOffset=0,momentOffset=0;
    std::fprintf(json,
        "{\n  \"schema\": \"amps.comet.stage0-diagnostics.v1\",\n"
        "  \"iteration\": %ld,\n  \"output_index\": %ld,\n"
        "  \"sampling_length\": %ld,\n"
        "  \"emergency_requests\": %llu,\n  \"reactions\": [\n",
        completedIteration,outputIndex,sampleLength,
        static_cast<unsigned long long>(globalEmergencies));
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
          "reaction,%ld,%ld,%ld,%zu,%s,,,,%llu,%llu,%llu,%llu,%llu,%llu,%llu,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n",
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
        double sigma=0.0,tolerance=0.0;
        const bool agreement=ProductionAgreement(expectedNumber,
            representedNumber,variance,sampleLength,sigma,tolerance);
        std::fprintf(csv,
            "product,%ld,%ld,%ld,%zu,%s,%zu,%s,%.17e,,,,,,,,%llu,%llu,%llu,%llu,%llu,%llu,%.17e,%.17e,%.17e,%.17e,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,,,,,,,,,,,\n",
            completedIteration,outputIndex,sampleLength,reaction,entry.Id,
            product,symbol,term.Coefficient,
            static_cast<unsigned long long>(multiplicityZero),
            static_cast<unsigned long long>(multiplicityOne),
            static_cast<unsigned long long>(multiplicityMany),
            static_cast<unsigned long long>(allocated),
            static_cast<unsigned long long>(reused),
            static_cast<unsigned long long>(carriers),expectedNumber,
            representedNumber,sigma,tolerance,agreement ? "true" : "false",
            expectedMass,representedMass,expectedMomentum[0],
            expectedMomentum[1],expectedMomentum[2],representedMomentum[0],
            representedMomentum[1],representedMomentum[2],expectedEnergy,
            representedEnergy);
        std::fprintf(json,
            "%s{\"term\": %zu, \"species\": \"%s\", \"coefficient\": %.17g, "
            "\"multiplicity_zero\": %llu, \"multiplicity_one\": %llu, "
            "\"multiplicity_many\": %llu, \"allocated_slots\": %llu, "
            "\"reused_slots\": %llu, \"model_carriers\": %llu, "
            "\"expected_number_rate_s-1\": %.17e, "
            "\"represented_number_rate_s-1\": %.17e, "
            "\"represented_number_sigma_s-1\": %.17e, "
            "\"production_tolerance_s-1\": %.17e, "
            "\"production_agreement\": %s}",
            product==0 ? "" : ", ",product,symbol,term.Coefficient,
            static_cast<unsigned long long>(multiplicityZero),
            static_cast<unsigned long long>(multiplicityOne),
            static_cast<unsigned long long>(multiplicityMany),
            static_cast<unsigned long long>(allocated),
            static_cast<unsigned long long>(reused),
            static_cast<unsigned long long>(carriers),expectedNumber,
            representedNumber,sigma,tolerance,agreement ? "true" : "false");
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
      double sigma=0.0,tolerance=0.0;
      const bool agreement=ProductionAgreement(expectedProduction,
          representedProduction,globalFlux[base+ChemicalProductionVariance],
          sampleLength,sigma,tolerance);
      std::fprintf(csv,
          "species_flux,%ld,%ld,%ld,,,,%s,,,,,,,,,,,,,,,%.17e,%.17e,%.17e,%.17e,%s,,,,,,,,,,,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,",
          completedIteration,outputIndex,sampleLength,symbol,
          expectedProduction,representedProduction,sigma,tolerance,
          agreement ? "true" : "false",
          Comet::Model::SurfaceProductionRate[species],surface,
          expectedProduction,representedProduction,destruction,innerIncident,
          innerReemitted,outerEscape,netExpected,netRepresented);
      if (hasPreviousProfile) std::fprintf(csv,"%.17e",relativeL1[species]);
      std::fprintf(csv,"\n");
      std::fprintf(json,
          "    {\"species\": \"%s\", \"configured_surface_source_s-1\": %.17e, "
          "\"sampled_surface_source_s-1\": %.17e, "
          "\"chemical_expected_production_s-1\": %.17e, "
          "\"chemical_represented_production_s-1\": %.17e, "
          "\"chemical_production_sigma_s-1\": %.17e, "
          "\"production_tolerance_s-1\": %.17e, "
          "\"production_agreement\": %s, "
          "\"chemical_destruction_s-1\": %.17e, "
          "\"inner_incident_s-1\": %.17e, \"inner_reemitted_s-1\": %.17e, "
          "\"outer_escape_s-1\": %.17e, \"net_expected_flux_s-1\": %.17e, "
          "\"net_represented_flux_s-1\": %.17e, \"relative_l1_change\": ",
          symbol,Comet::Model::SurfaceProductionRate[species],surface,
          expectedProduction,representedProduction,sigma,tolerance,
          agreement ? "true" : "false",destruction,innerIncident,
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
  GlobalEmergencyRequests=0;
}
