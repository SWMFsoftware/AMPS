#ifndef _CHEMICAL_SOLVER_CHEMISTRY_DIAGNOSTICS_H_
#define _CHEMICAL_SOLVER_CHEMISTRY_DIAGNOSTICS_H_

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ChemicalSolver {

// ChemistryDiagnostics is deliberately independent of PIC, MPI, and comet
// geometry.  The chemistry executor records local events here; a caller may
// then reduce CountValues() with an unsigned-integer sum and MomentValues()
// with a floating-point sum.  Keeping MPI outside this class makes the event
// accounting unit-testable and permits other AMPS applications to use the
// same definitions.
class cChemistryDiagnostics {
 public:
  enum ReactionCountField {
    EligibleParentVisits=0,
    PhysicalReactionSelections,
    SelectedChannels,
    ParentDestructions,
    NeutralSinks,
    AllocationDeferrals,
    EmergencyRequests,
    ReactionCountFieldNumber
  };

  enum ProductCountField {
    ProductMultiplicityZero=0,
    ProductMultiplicityOne,
    ProductMultiplicityMany,
    AllocatedSlots,
    ReusedSlots,
    ModelProductCarriers,
    ProductCountFieldNumber
  };

  enum MomentField {
    ExpectedNumber=0,
    ExpectedMass,
    ExpectedMomentumX,
    ExpectedMomentumY,
    ExpectedMomentumZ,
    ExpectedEnergy,
    RepresentedNumber,
    RepresentedMass,
    RepresentedMomentumX,
    RepresentedMomentumY,
    RepresentedMomentumZ,
    RepresentedEnergy,
    RepresentedNumberVariance,
    MomentFieldNumber
  };

  struct cProductRecord {
    std::uint64_t MultiplicityZero=0;
    std::uint64_t MultiplicityOne=0;
    std::uint64_t MultiplicityMany=0;
    std::uint64_t Allocated=0;
    std::uint64_t Reused=0;
    double ExpectedNumber=0.0;
    double ExpectedMass=0.0;
    double ExpectedMomentum[3]={0.0,0.0,0.0};
    double ExpectedEnergy=0.0;
    double Number=0.0;
    double Mass=0.0;
    double Momentum[3]={0.0,0.0,0.0};
    double Energy=0.0;
    // Sum of independent residual split/roulette variances in physical
    // number-rate squared.  After a window average over N outer iterations,
    // sqrt(value)/N is the one-sigma uncertainty of represented production.
    double NumberVariance=0.0;
  };

  struct cReactionRecord {
    std::uint64_t Eligible=0;
    std::uint64_t Reactions=0;
    std::uint64_t Channels=0;
    std::uint64_t Destructions=0;
    std::uint64_t Sinks=0;
    std::uint64_t Deferrals=0;
    std::uint64_t Emergencies=0;
    std::vector<cProductRecord> Products;
  };

  cChemistryDiagnostics() : Enabled(false) {}

  // ProductCounts has one entry per reaction and gives the number of product
  // terms, not the number of distinct species.  Consequently repeated
  // products remain separately auditable and can later be reconciled with
  // the exact generated reaction-table term that produced them.
  void Configure(const std::vector<std::size_t> &productCounts,bool enabled) {
    Enabled=enabled;
    Records.clear();
    if (!Enabled) return;
    Records.resize(productCounts.size());
    for (std::size_t reaction=0;reaction<productCounts.size();reaction++)
      Records[reaction].Products.resize(productCounts[reaction]);
  }

  bool IsEnabled() const { return Enabled; }
  std::size_t GetReactionCount() const { return Records.size(); }

  const cReactionRecord& GetReaction(std::size_t reaction) const {
    CheckReaction(reaction);
    return Records[reaction];
  }

  void RecordEligible(std::size_t reaction,std::uint64_t count=1) {
    if (!Enabled) return;
    CheckReaction(reaction);
    Records[reaction].Eligible+=count;
  }

  void RecordSelection(std::size_t reaction,bool neutralSink) {
    if (!Enabled) return;
    CheckReaction(reaction);
    cReactionRecord &record=Records[reaction];
    record.Reactions++;
    record.Channels++;
    record.Destructions++;
    if (neutralSink) record.Sinks++;
  }

  // RecordProduct describes one stoichiometric product term for one selected
  // physical reaction event.  multiplicity is the number of model carriers
  // created by residual split/roulette.  represented* values are physical
  // rates carried by those model particles after all weight/timestep factors
  // have been applied, so their reduction is meaningful across blocks.
  void RecordProduct(std::size_t reaction,std::size_t product,
      std::uint64_t multiplicity,std::uint64_t allocated,
      std::uint64_t reused,double expectedNumber,double expectedMass,
      const double expectedMomentum[3],double expectedEnergy,
      double representedNumber,
      double representedMass,const double representedMomentum[3],
      double representedEnergy,double representedNumberVariance) {
    if (!Enabled) return;
    CheckProduct(reaction,product);
    if (allocated+reused!=multiplicity)
      throw std::invalid_argument(
          "chemistry diagnostic allocation/reuse does not equal multiplicity");
    cProductRecord &record=Records[reaction].Products[product];
    if (multiplicity==0) record.MultiplicityZero++;
    else if (multiplicity==1) record.MultiplicityOne++;
    else record.MultiplicityMany++;
    record.Allocated+=allocated;
    record.Reused+=reused;
    record.ExpectedNumber+=expectedNumber;
    record.ExpectedMass+=expectedMass;
    for (int component=0;component<3;component++)
      record.ExpectedMomentum[component]+=expectedMomentum[component];
    record.ExpectedEnergy+=expectedEnergy;
    record.Number+=representedNumber;
    record.Mass+=representedMass;
    for (int component=0;component<3;component++)
      record.Momentum[component]+=representedMomentum[component];
    record.Energy+=representedEnergy;
    record.NumberVariance+=representedNumberVariance;
  }

  void RecordDeferral(std::size_t reaction,std::uint64_t count=1) {
    if (!Enabled) return;
    CheckReaction(reaction);
    Records[reaction].Deferrals+=count;
  }

  void RecordEmergency(std::size_t reaction,std::uint64_t count=1) {
    if (!Enabled) return;
    CheckReaction(reaction);
    Records[reaction].Emergencies+=count;
  }

  // Clear one completed sampling window while preserving the generated
  // reaction/product shape and the enabled state.  Production diagnostics
  // call this only after collective reduction and durable output.
  void Reset() {
    if (!Enabled) return;
    for (std::size_t reaction=0;reaction<Records.size();reaction++) {
      const std::size_t productCount=Records[reaction].Products.size();
      Records[reaction]=cReactionRecord();
      Records[reaction].Products.resize(productCount);
    }
  }

  // The flattened layout is stable and documented by the two enums.  Each
  // reaction begins with its seven reaction-wide counters, followed by six
  // counters for each product term.  Disabled diagnostics intentionally
  // serialize to empty vectors: they add no output fields.
  std::vector<std::uint64_t> CountValues() const {
    std::vector<std::uint64_t> values;
    if (!Enabled) return values;
    for (std::size_t reaction=0;reaction<Records.size();reaction++) {
      const cReactionRecord &entry=Records[reaction];
      values.push_back(entry.Eligible);
      values.push_back(entry.Reactions);
      values.push_back(entry.Channels);
      values.push_back(entry.Destructions);
      values.push_back(entry.Sinks);
      values.push_back(entry.Deferrals);
      values.push_back(entry.Emergencies);
      for (std::size_t product=0;product<entry.Products.size();product++) {
        const cProductRecord &term=entry.Products[product];
        values.push_back(term.MultiplicityZero);
        values.push_back(term.MultiplicityOne);
        values.push_back(term.MultiplicityMany);
        values.push_back(term.Allocated);
        values.push_back(term.Reused);
        // The final value is the total number of model product carriers.  It
        // is retained independently of allocation/reuse to catch accounting
        // errors when a transaction reuses its parent slot.
        values.push_back(term.Allocated+term.Reused);
      }
    }
    return values;
  }

  std::vector<double> MomentValues() const {
    std::vector<double> values;
    if (!Enabled) return values;
    for (std::size_t reaction=0;reaction<Records.size();reaction++)
      for (std::size_t product=0;
           product<Records[reaction].Products.size();product++) {
        const cProductRecord &term=Records[reaction].Products[product];
        values.push_back(term.ExpectedNumber);
        values.push_back(term.ExpectedMass);
        values.push_back(term.ExpectedMomentum[0]);
        values.push_back(term.ExpectedMomentum[1]);
        values.push_back(term.ExpectedMomentum[2]);
        values.push_back(term.ExpectedEnergy);
        values.push_back(term.Number);
        values.push_back(term.Mass);
        values.push_back(term.Momentum[0]);
        values.push_back(term.Momentum[1]);
        values.push_back(term.Momentum[2]);
        values.push_back(term.Energy);
        values.push_back(term.NumberVariance);
      }
    return values;
  }

 private:
  bool Enabled;
  std::vector<cReactionRecord> Records;

  void CheckReaction(std::size_t reaction) const {
    if (reaction>=Records.size())
      throw std::out_of_range("chemistry diagnostic reaction index");
  }

  void CheckProduct(std::size_t reaction,std::size_t product) const {
    CheckReaction(reaction);
    if (product>=Records[reaction].Products.size())
      throw std::out_of_range("chemistry diagnostic product index");
  }
};

} // namespace ChemicalSolver

#endif
