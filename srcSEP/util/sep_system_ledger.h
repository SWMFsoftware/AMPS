#ifndef SEP_UTIL_SEP_SYSTEM_LEDGER_H
#define SEP_UTIL_SEP_SYSTEM_LEDGER_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace SystemVerification {

// Each component follows an explicit sign convention: positive values enter
// the modeled particle-plus-wave system and negative values leave it.  Charge
// [C], energy [J], and parallel momentum [kg m/s] remain species-summed only
// after per-species adapters have recorded their own transaction.
struct ConservedQuantities {
  double number = 0.0;
  double chargeC = 0.0;
  double energyJ = 0.0;
  double parallelMomentumKgMPerS = 0.0;
};

struct Transaction {
  std::string seam;
  ConservedQuantities change;
  std::string provenance;
};

struct GlobalLedger {
  ConservedQuantities initial;
  ConservedQuantities final;
  std::vector<Transaction> transactions;
  ConservedQuantities residual;
};

Transport::Status CloseLedger(GlobalLedger* ledger, double relativeTolerance);

struct ParticleIdentityState {
  std::uint64_t stableId = 0;
  std::uint64_t lineageGeneration = 0;
  std::uint64_t rngEventCounter = 0;
  double residualHazard = 0.0;
};

struct CampaignCheckpoint {
  std::string configurationFingerprint;
  std::uint64_t fieldLineGeneration = 0;
  double epochS = 0.0;
  std::vector<ParticleIdentityState> particles;
  GlobalLedger ledger;
  std::string turbulenceStateHash;
  std::string outputManifestChecksum;
};

// The checkpoint text is canonical and ends with an FNV checksum over every
// preceding byte.  This protects residual hazards, RNG counters, identities,
// ledgers, and remap generation from partial or duplicated seam transactions.
Transport::Status SerializeCheckpoint(const CampaignCheckpoint& checkpoint,
                                      std::string* text);
Transport::Status DeserializeCheckpoint(const std::string& text,
                                        CampaignCheckpoint* checkpoint);
std::uint64_t CheckpointHash(const CampaignCheckpoint& checkpoint);

}  // namespace SystemVerification
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SYSTEM_LEDGER_H
