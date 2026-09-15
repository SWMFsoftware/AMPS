#include "sep_system_ledger.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace SEP {
namespace SystemVerification {
namespace {

Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                  message);
}

std::uint64_t HashBytes(const std::string& text) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (std::size_t i = 0; i < text.size(); ++i) {
    hash ^= static_cast<unsigned char>(text[i]);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

void Add(ConservedQuantities* target, const ConservedQuantities& value) {
  target->number += value.number;
  target->chargeC += value.chargeC;
  target->energyJ += value.energyJ;
  target->parallelMomentumKgMPerS += value.parallelMomentumKgMPerS;
}

bool Finite(const ConservedQuantities& q) {
  return std::isfinite(q.number) && std::isfinite(q.chargeC) &&
         std::isfinite(q.energyJ) && std::isfinite(q.parallelMomentumKgMPerS);
}

void WriteQuantities(std::ostream& out, const ConservedQuantities& q) {
  out << q.number << ' ' << q.chargeC << ' ' << q.energyJ << ' '
      << q.parallelMomentumKgMPerS << '\n';
}

bool ReadQuantities(std::istream& in, ConservedQuantities* q) {
  return static_cast<bool>(in >> q->number >> q->chargeC >> q->energyJ >>
                           q->parallelMomentumKgMPerS);
}

std::string EncodeToken(const std::string& value) {
  static const char digits[] = "0123456789abcdef";
  std::string encoded("x");
  for (std::size_t i = 0; i < value.size(); ++i) {
    const unsigned char byte = static_cast<unsigned char>(value[i]);
    encoded.push_back(digits[byte >> 4]);
    encoded.push_back(digits[byte & 15U]);
  }
  return encoded;
}

bool DecodeToken(const std::string& encoded, std::string* value) {
  if (!value || encoded.empty() || encoded[0] != 'x' || encoded.size() % 2 != 1)
    return false;
  value->clear();
  for (std::size_t i = 1; i < encoded.size(); i += 2) {
    const std::size_t hi = std::string("0123456789abcdef").find(encoded[i]);
    const std::size_t lo = std::string("0123456789abcdef").find(encoded[i + 1]);
    if (hi == std::string::npos || lo == std::string::npos) return false;
    value->push_back(static_cast<char>((hi << 4) | lo));
  }
  return true;
}

std::string CanonicalWithoutChecksum(const CampaignCheckpoint& c) {
  std::ostringstream out;
  out << "SEP_SYSTEM_CHECKPOINT 1\n" << std::setprecision(17)
      << EncodeToken(c.configurationFingerprint) << ' '
      << c.fieldLineGeneration << ' ' << c.epochS << '\n'
      << EncodeToken(c.turbulenceStateHash) << ' '
      << EncodeToken(c.outputManifestChecksum) << '\n';
  WriteQuantities(out, c.ledger.initial);
  WriteQuantities(out, c.ledger.final);
  out << c.ledger.transactions.size() << '\n';
  for (std::size_t i = 0; i < c.ledger.transactions.size(); ++i) {
    out << EncodeToken(c.ledger.transactions[i].seam) << ' '
        << EncodeToken(c.ledger.transactions[i].provenance) << '\n';
    WriteQuantities(out, c.ledger.transactions[i].change);
  }
  out << c.particles.size() << '\n';
  for (std::size_t i = 0; i < c.particles.size(); ++i)
    out << c.particles[i].stableId << ' '
        << c.particles[i].lineageGeneration << ' '
        << c.particles[i].rngEventCounter << ' '
        << c.particles[i].residualHazard << '\n';
  return out.str();
}

}  // namespace

Transport::Status CloseLedger(GlobalLedger* ledger, double tolerance) {
  if (!ledger || !(tolerance >= 0.0) || !std::isfinite(tolerance) ||
      !Finite(ledger->initial) || !Finite(ledger->final))
    return Error("global ledger or closure tolerance is invalid");
  ConservedQuantities expected = ledger->initial;
  for (std::size_t i = 0; i < ledger->transactions.size(); ++i) {
    if (ledger->transactions[i].seam.empty() ||
        ledger->transactions[i].provenance.empty() ||
        !Finite(ledger->transactions[i].change))
      return Error("global ledger contains an invalid seam transaction");
    Add(&expected, ledger->transactions[i].change);
  }
  ledger->residual.number = ledger->final.number - expected.number;
  ledger->residual.chargeC = ledger->final.chargeC - expected.chargeC;
  ledger->residual.energyJ = ledger->final.energyJ - expected.energyJ;
  ledger->residual.parallelMomentumKgMPerS =
      ledger->final.parallelMomentumKgMPerS - expected.parallelMomentumKgMPerS;
  const double residuals[] = {ledger->residual.number, ledger->residual.chargeC,
      ledger->residual.energyJ, ledger->residual.parallelMomentumKgMPerS};
  const double finals[] = {ledger->final.number, ledger->final.chargeC,
      ledger->final.energyJ, ledger->final.parallelMomentumKgMPerS};
  for (std::size_t i = 0; i < 4; ++i) {
    const double scale = std::max(1.0, std::fabs(finals[i]));
    if (std::fabs(residuals[i]) > tolerance * scale)
      return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                      "global system ledger does not close");
  }
  return Transport::Status::Ok();
}

Transport::Status SerializeCheckpoint(const CampaignCheckpoint& checkpoint,
                                      std::string* text) {
  if (!text || checkpoint.configurationFingerprint.empty() ||
      checkpoint.turbulenceStateHash.empty() ||
      checkpoint.outputManifestChecksum.empty() ||
      !std::isfinite(checkpoint.epochS))
    return Error("system checkpoint metadata is incomplete");
  const std::string payload = CanonicalWithoutChecksum(checkpoint);
  std::ostringstream out;
  out << payload << "CHECKSUM " << std::hex << HashBytes(payload) << '\n';
  *text = out.str();
  return Transport::Status::Ok();
}

Transport::Status DeserializeCheckpoint(const std::string& text,
                                        CampaignCheckpoint* checkpoint) {
  if (!checkpoint) return Error("system checkpoint destination is null");
  const std::size_t checksumPosition = text.rfind("CHECKSUM ");
  if (checksumPosition == std::string::npos)
    return Error("system checkpoint has no completion checksum");
  const std::string payload = text.substr(0, checksumPosition);
  std::istringstream checksumInput(text.substr(checksumPosition));
  std::string checksumLabel;
  std::uint64_t expectedHash = 0;
  if (!(checksumInput >> checksumLabel >> std::hex >> expectedHash) ||
      checksumLabel != "CHECKSUM" || expectedHash != HashBytes(payload))
    return Error("system checkpoint checksum mismatch");

  CampaignCheckpoint parsed;
  std::istringstream in(payload);
  std::string magic, fingerprint, turbulence, output;
  int version = 0;
  if (!(in >> magic >> version) || magic != "SEP_SYSTEM_CHECKPOINT" || version != 1 ||
      !(in >> fingerprint >> parsed.fieldLineGeneration >> parsed.epochS >>
        turbulence >> output) ||
      !DecodeToken(fingerprint, &parsed.configurationFingerprint) ||
      !DecodeToken(turbulence, &parsed.turbulenceStateHash) ||
      !DecodeToken(output, &parsed.outputManifestChecksum) ||
      !ReadQuantities(in, &parsed.ledger.initial) ||
      !ReadQuantities(in, &parsed.ledger.final))
    return Error("system checkpoint header is malformed");
  std::size_t transactions = 0;
  if (!(in >> transactions)) return Error("system transaction count is malformed");
  parsed.ledger.transactions.resize(transactions);
  for (std::size_t i = 0; i < transactions; ++i) {
    std::string seam, provenance;
    if (!(in >> seam >> provenance) ||
        !DecodeToken(seam, &parsed.ledger.transactions[i].seam) ||
        !DecodeToken(provenance, &parsed.ledger.transactions[i].provenance) ||
        !ReadQuantities(in, &parsed.ledger.transactions[i].change))
      return Error("system transaction record is malformed");
  }
  std::size_t particles = 0;
  if (!(in >> particles)) return Error("system particle count is malformed");
  parsed.particles.resize(particles);
  for (std::size_t i = 0; i < particles; ++i) {
    ParticleIdentityState& particle = parsed.particles[i];
    if (!(in >> particle.stableId >> particle.lineageGeneration >>
          particle.rngEventCounter >> particle.residualHazard) ||
        particle.stableId == 0 || !(particle.residualHazard >= 0.0) ||
        !std::isfinite(particle.residualHazard))
      return Error("system particle identity record is malformed");
  }
  std::string trailing;
  if (in >> trailing) return Error("system checkpoint contains trailing payload");
  *checkpoint = parsed;
  return Transport::Status::Ok();
}

std::uint64_t CheckpointHash(const CampaignCheckpoint& checkpoint) {
  std::string text;
  return SerializeCheckpoint(checkpoint, &text).ok() ? HashBytes(text) : 0;
}

}  // namespace SystemVerification
}  // namespace SEP
