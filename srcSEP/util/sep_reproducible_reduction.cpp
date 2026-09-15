#include "sep_reproducible_reduction.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace Reproducibility {
namespace {

bool KeyLess(const Contribution& a, const Contribution& b) {
  if (a.key.schema != b.key.schema) return a.key.schema < b.key.schema;
  if (a.key.source != b.key.source) return a.key.source < b.key.source;
  if (a.key.fieldLine != b.key.fieldLine)
    return a.key.fieldLine < b.key.fieldLine;
  if (a.key.segment != b.key.segment) return a.key.segment < b.key.segment;
  if (a.key.branch != b.key.branch) return a.key.branch < b.key.branch;
  if (a.key.spectralBin != b.key.spectralBin)
    return a.key.spectralBin < b.key.spectralBin;
  if (a.key.species != b.key.species) return a.key.species < b.key.species;
  if (a.key.particle != b.key.particle) return a.key.particle < b.key.particle;
  if (a.key.step != b.key.step) return a.key.step < b.key.step;
  if (a.key.event != b.key.event) return a.key.event < b.key.event;
  if (a.key.interval != b.key.interval) return a.key.interval < b.key.interval;
  return a.key.purpose < b.key.purpose;
}

bool SameKey(const Contribution& a, const Contribution& b) {
  return !KeyLess(a, b) && !KeyLess(b, a);
}

bool SameSegment(const Contribution& value, const SegmentAccumulator& segment) {
  return value.key.fieldLine == segment.fieldLine &&
         value.key.segment == segment.segment &&
         value.key.branch == segment.branch &&
         value.key.spectralBin == segment.spectralBin;
}

Transport::Status ReduceFlat(std::vector<Contribution> values,
                             std::vector<SegmentAccumulator>* segments) {
  if (!segments)
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "reduction output pointer is null");
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (!std::isfinite(values[i].waveEnergyJ) ||
        !std::isfinite(values[i].streaming))
      return Transport::Status::Error(
          Transport::StatusCode::InvalidArgument,
          "non-finite particle contribution reached deterministic reduction");
  }
  std::sort(values.begin(), values.end(), KeyLess);
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (SameKey(values[i - 1], values[i]))
      return Transport::Status::Error(
          Transport::StatusCode::InvalidArgument,
          "duplicate physical contribution key reached deterministic reduction");
  }
  segments->clear();
  std::size_t begin = 0;
  while (begin < values.size()) {
    SegmentAccumulator result;
    result.fieldLine = values[begin].key.fieldLine;
    result.segment = values[begin].key.segment;
    result.branch = values[begin].key.branch;
    result.spectralBin = values[begin].key.spectralBin;
    long double energy = 0.0L;
    long double streaming = 0.0L;
    std::uint64_t count = 0;
    std::size_t end = begin;
    while (end < values.size() && SameSegment(values[end], result)) {
      energy += static_cast<long double>(values[end].waveEnergyJ);
      streaming += static_cast<long double>(values[end].streaming);
      if (std::numeric_limits<std::uint64_t>::max() - count <
          values[end].resonantCount)
        return Transport::Status::Error(
            Transport::StatusCode::InvalidArgument,
            "resonant event counter overflow in deterministic reduction");
      count += values[end].resonantCount;
      ++end;
    }
    result.waveEnergyJ = static_cast<double>(energy);
    result.streaming = static_cast<double>(streaming);
    result.resonantCount = count;
    segments->push_back(result);
    begin = end;
  }
  return Transport::Status::Ok();
}

void HashByte(std::uint64_t* hash, unsigned char value) {
  *hash ^= value;
  *hash *= 1099511628211ULL;
}

void HashU64(std::uint64_t* hash, std::uint64_t value) {
  // An explicit little-endian byte order makes the evidence hash independent
  // of host endianness and avoids hashing compiler-dependent struct padding.
  for (unsigned i = 0; i < 8; ++i)
    HashByte(hash, static_cast<unsigned char>((value >> (8 * i)) & 0xffU));
}

void HashDouble(std::uint64_t* hash, double value) {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value), "unexpected double width");
  std::memcpy(&bits, &value, sizeof(bits));
  HashU64(hash, bits);
}

}  // namespace

void ThreadLocalBuffer::Add(const Contribution& value) {
  values_.push_back(value);
}

void ThreadLocalBuffer::Clear() {
  values_.clear();
}

Transport::Status CanonicalReduction(
    const std::vector<ThreadLocalBuffer>& workers,
    std::vector<SegmentAccumulator>* segments) {
  std::size_t count = 0;
  for (std::size_t worker = 0; worker < workers.size(); ++worker)
    count += workers[worker].values().size();
  std::vector<Contribution> flat;
  flat.reserve(count);
  for (std::size_t worker = 0; worker < workers.size(); ++worker)
    flat.insert(flat.end(), workers[worker].values().begin(),
                workers[worker].values().end());
  return ReduceFlat(flat, segments);
}

Transport::Status CanonicalPartitionReduction(
    const std::vector<std::vector<Contribution> >& partitions,
    std::vector<SegmentAccumulator>* segments) {
  std::vector<Contribution> flat;
  std::size_t count = 0;
  for (std::size_t rank = 0; rank < partitions.size(); ++rank)
    count += partitions[rank].size();
  flat.reserve(count);
  // MPI implementations may return gathered rank blocks in rank order.  The
  // subsequent complete-key sort erases that decomposition-dependent order.
  for (std::size_t rank = 0; rank < partitions.size(); ++rank)
    flat.insert(flat.end(), partitions[rank].begin(), partitions[rank].end());
  return ReduceFlat(flat, segments);
}

Transport::KeyedRandomStream MakeRandomStream(std::uint64_t campaignSeed,
                                               std::uint64_t particleId,
                                               std::uint64_t step,
                                               std::uint64_t purpose) {
  return Transport::KeyedRandomStream(campaignSeed, particleId, purpose, step);
}

void AtomicCounters::AddWarning(std::uint64_t count) {
  warnings_.fetch_add(count, std::memory_order_relaxed);
}

void AtomicCounters::AddLimiterActivation(std::uint64_t count) {
  limiterActivations_.fetch_add(count, std::memory_order_relaxed);
}

std::uint64_t AtomicCounters::warnings() const {
  return warnings_.load(std::memory_order_relaxed);
}

std::uint64_t AtomicCounters::limiterActivations() const {
  return limiterActivations_.load(std::memory_order_relaxed);
}

void AtomicCounters::Reset() {
  warnings_.store(0, std::memory_order_relaxed);
  limiterActivations_.store(0, std::memory_order_relaxed);
}

std::uint64_t EvidenceHash(const std::vector<SegmentAccumulator>& segments) {
  std::uint64_t hash = 1469598103934665603ULL;
  HashU64(&hash, static_cast<std::uint64_t>(segments.size()));
  for (std::size_t i = 0; i < segments.size(); ++i) {
    HashU64(&hash, segments[i].fieldLine);
    HashU64(&hash, segments[i].segment);
    HashU64(&hash, segments[i].branch);
    HashU64(&hash, segments[i].spectralBin);
    HashDouble(&hash, segments[i].waveEnergyJ);
    HashDouble(&hash, segments[i].streaming);
    HashU64(&hash, segments[i].resonantCount);
  }
  return hash;
}

std::string EvidenceHashHex(const std::vector<SegmentAccumulator>& segments) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16)
      << EvidenceHash(segments);
  return out.str();
}

}  // namespace Reproducibility
}  // namespace SEP
