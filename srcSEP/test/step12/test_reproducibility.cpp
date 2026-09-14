#include "sep_reproducible_reduction.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

namespace {

using SEP::Reproducibility::Contribution;
using SEP::Reproducibility::SegmentAccumulator;
using SEP::Reproducibility::ThreadLocalBuffer;

Contribution Make(std::uint64_t particle, double energy) {
  Contribution value;
  value.key.fieldLine = particle % 2;
  value.key.segment = particle % 3;
  value.key.branch = particle % 2;
  value.key.particle = particle;
  value.key.step = 7;
  value.key.purpose = 11;
  value.waveEnergyJ = energy;
  value.streaming = 0.5 * energy;
  value.resonantCount = 1;
  return value;
}

void Report(const char* id, bool passed, int* failures) {
  std::cout << (passed ? "PASS " : "FAIL ") << id << '\n';
  if (!passed) ++*failures;
}

}  // namespace

int main() {
  int failures = 0;
  std::vector<Contribution> input;
  for (std::uint64_t i = 0; i < 200; ++i)
    input.push_back(Make(i, i % 5 == 0 ? 1.0e12 : 1.0e-6));

  std::vector<ThreadLocalBuffer> oneWorker(1);
  for (std::size_t i = 0; i < input.size(); ++i) oneWorker[0].Add(input[i]);
  std::vector<SegmentAccumulator> reference;
  const bool referenceOk =
      SEP::Reproducibility::CanonicalReduction(oneWorker, &reference).ok();

  std::vector<ThreadLocalBuffer> manyWorkers(7);
  for (std::size_t i = 0; i < input.size(); ++i)
    manyWorkers[(13 * i + 3) % manyWorkers.size()].Add(input[i]);
  std::vector<SegmentAccumulator> scheduled;
  const bool scheduledOk =
      SEP::Reproducibility::CanonicalReduction(manyWorkers, &scheduled).ok();
  Report("PAR01 single-writer reduction", referenceOk && !reference.empty(),
         &failures);
  Report("PAR02 scheduler independence", scheduledOk &&
         SEP::Reproducibility::EvidenceHash(reference) ==
         SEP::Reproducibility::EvidenceHash(scheduled), &failures);

  std::vector<std::vector<Contribution> > partitions(5);
  for (std::size_t i = 0; i < input.size(); ++i)
    partitions[(17 * i + 1) % partitions.size()].push_back(input[i]);
  std::vector<SegmentAccumulator> mpi;
  const bool mpiOk = SEP::Reproducibility::CanonicalPartitionReduction(
      partitions, &mpi).ok();
  Report("PAR03 MPI decomposition independence", mpiOk &&
         SEP::Reproducibility::EvidenceHash(reference) ==
         SEP::Reproducibility::EvidenceHash(mpi), &failures);

  SEP::Transport::KeyedRandomStream moverA =
      SEP::Reproducibility::MakeRandomStream(5, 42, 9, 100);
  SEP::Transport::KeyedRandomStream moverB =
      SEP::Reproducibility::MakeRandomStream(5, 42, 9, 100);
  SEP::Transport::KeyedRandomStream diagnostic =
      SEP::Reproducibility::MakeRandomStream(5, 42, 9, 200);
  const double firstA = moverA.UniformOpen01();
  (void)diagnostic.UniformOpen01();
  const double firstB = moverB.UniformOpen01();
  Report("PAR04 RNG stream independence",
         std::memcmp(&firstA, &firstB, sizeof(double)) == 0, &failures);

  SEP::Reproducibility::AtomicCounters counters;
  counters.AddWarning(2);
  counters.AddLimiterActivation(3);
  const std::string hash = SEP::Reproducibility::EvidenceHashHex(reference);
  const bool evidence = hash.size() == 16 && counters.warnings() == 2 &&
                        counters.limiterActivations() == 3;
  counters.Reset();
  Report("PAR05 reproducibility evidence", evidence &&
         counters.warnings() == 0 && counters.limiterActivations() == 0,
         &failures);

  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
