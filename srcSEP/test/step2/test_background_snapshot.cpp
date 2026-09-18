#include "sep_background_snapshot.h"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

namespace {

int failures = 0;

void Check(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

template <typename ExceptionType>
void CheckThrows(const std::function<void()>& operation,
                 const std::string& message) {
  try {
    operation();
    Check(false, message + " (no exception)");
  }
  catch (const ExceptionType&) {
    // Expected failure path.
  }
  catch (...) {
    Check(false, message + " (wrong exception type)");
  }
}

SEP::Background::BackgroundSnapshot AnalyticSnapshot(
    double epoch, double valid_until, std::uint64_t generation = 1) {
  return SEP::Background::BackgroundSnapshot(
      SEP::Background::Provider::Analytic,
      SEP::Background::Ownership::ModelOwned, epoch, epoch, valid_until,
      generation, SEP::Background::FingerprintConfiguration("analytic-test"),
      "deterministic Step 2 unit fixture");
}

void TestConstructionAndImmutability() {
  using SEP::Background::BackgroundSnapshot;
  using SEP::Background::Ownership;
  using SEP::Background::Provider;

  // Const data members plus a deleted assignment operator make accidental
  // post-publication metadata mutation a compile-time failure.
  static_assert(!std::is_copy_assignable<BackgroundSnapshot>::value,
                "background snapshots must not be assignable");

  const BackgroundSnapshot snapshot = AnalyticSnapshot(10.0, 20.0, 7);
  Check(snapshot.provider() == Provider::Analytic,
        "provider identity is preserved");
  Check(snapshot.ownership() == Ownership::ModelOwned,
        "ownership identity is preserved");
  Check(snapshot.epoch_seconds() == 10.0,
        "epoch is preserved in seconds");
  Check(snapshot.field_line_generation() == 7,
        "field-line generation is preserved");
  Check(snapshot.Covers(10.0) && snapshot.Covers(20.0) &&
            !snapshot.Covers(20.01),
        "closed validity interval is enforced");

  CheckThrows<std::invalid_argument>([]() {
    BackgroundSnapshot invalid(
        Provider::Swmf, Ownership::ModelOwned, 0.0, 0.0, 1.0, 1,
        "fingerprint", "invalid mutable SWMF fixture");
    (void)invalid;
  }, "SWMF snapshots reject mutable ownership");

  CheckThrows<std::invalid_argument>([]() {
    BackgroundSnapshot invalid(
        Provider::Analytic, Ownership::ModelOwned, 2.0, 0.0, 1.0, 1,
        "fingerprint", "epoch outside validity fixture");
    (void)invalid;
  }, "epoch outside validity is rejected");

  CheckThrows<std::invalid_argument>([]() {
    BackgroundSnapshot invalid(
        Provider::Analytic, Ownership::ModelOwned,
        std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0, 1,
        "fingerprint", "non-finite epoch fixture");
    (void)invalid;
  }, "non-finite epoch is rejected");
}

void TestPublicationAndReadPhase() {
  using namespace SEP::Background;
  SnapshotStore& store = SnapshotStore::Instance();
  store.ResetForTests();

  store.Publish(AnalyticSnapshot(0.0, 10.0));
  {
    ParticleReadPhase phase = store.BeginParticleRead(5.0);
    const BackgroundSnapshot& mover_view = store.AcquireForMover();
    Check(mover_view.epoch_seconds() == phase.snapshot().epoch_seconds(),
          "mover and enclosing read phase consume one snapshot");

    CheckThrows<std::logic_error>([&store]() {
      store.Publish(AnalyticSnapshot(10.0, 20.0));
    }, "publication is blocked while particles read the background");

    CheckThrows<std::logic_error>([&store]() {
      store.AssertProviderMayWrite(Provider::Analytic);
    }, "provider backing storage cannot change while particles read it");
  }

  store.AssertProviderMayWrite(Provider::Analytic);

  CheckThrows<std::logic_error>([&store]() {
    store.AssertProviderMayWrite(Provider::Swcme);
  }, "a provider cannot mutate another provider's backing state");

  store.Publish(AnalyticSnapshot(10.0, 20.0));
  Check(store.Current()->epoch_seconds() == 10.0,
        "next epoch publishes after the read phase ends");

  CheckThrows<std::logic_error>([&store]() {
    ParticleReadPhase phase = store.BeginParticleRead(21.0);
    (void)phase;
  }, "particle time outside validity is rejected");

  CheckThrows<std::logic_error>([&store]() {
    (void)store.AcquireForMover();
  }, "mover entry outside a global read phase is rejected");
}

void TestProviderIsolationAndHandoff() {
  using namespace SEP::Background;
  SnapshotStore& store = SnapshotStore::Instance();
  store.ResetForTests();

  const std::string swmf_fingerprint =
      FingerprintConfiguration("provider=swmf;fixture=handoff");
  store.Publish(BackgroundSnapshot(
      Provider::Swmf, Ownership::ImportedReadOnly, 100.0, 100.0,
      std::numeric_limits<double>::infinity(), 3, swmf_fingerprint,
      "unit-test SWMF import"));

  CheckThrows<std::logic_error>([&store, &swmf_fingerprint]() {
    store.Publish(BackgroundSnapshot(
        Provider::LocalEvolution, Ownership::HandoffCopy, 100.0, 100.0,
        110.0, 4, swmf_fingerprint, "unannounced local copy"));
  }, "one provider cannot overwrite another provider's state");

  store.PublishHandoff(BackgroundSnapshot(
      Provider::LocalEvolution, Ownership::HandoffCopy, 100.0, 100.0,
      110.0, 4, swmf_fingerprint,
      "copied once from SWMF generation 3 at t=100 s"),
      Provider::Swmf);
  Check(store.Current()->provider() == Provider::LocalEvolution &&
            store.Current()->field_line_generation() == 4,
        "explicit SWMF handoff records new owner and generation");

  CheckThrows<std::logic_error>([&store, &swmf_fingerprint]() {
    store.Publish(BackgroundSnapshot(
        Provider::LocalEvolution, Ownership::HandoffCopy, 99.0, 99.0,
        100.0, 4, swmf_fingerprint, "time-regressing local update"));
  }, "background epoch cannot move backwards");
}

void TestConcurrentMoverViews() {
  using namespace SEP::Background;
  SnapshotStore& store = SnapshotStore::Instance();
  store.ResetForTests();
  store.Publish(AnalyticSnapshot(0.0, 1.0, 11));

  ParticleReadPhase phase = store.BeginParticleRead(0.5);
  std::atomic<int> matching_views(0);
  std::vector<std::thread> workers;

  // Multiple scheduler threads must see the exact same immutable generation;
  // this is the concurrency condition used by the common production mover
  // wrapper during an OpenMP/threaded PIC particle phase.
  for (int i = 0; i < 16; ++i) {
    workers.push_back(std::thread([&store, &matching_views]() {
      const BackgroundSnapshot& view = store.AcquireForMover();
      if (view.field_line_generation() == 11 &&
          view.epoch_seconds() == 0.0) {
        ++matching_views;
      }
    }));
  }

  for (std::size_t i = 0; i < workers.size(); ++i) workers[i].join();
  Check(matching_views == 16,
        "all concurrent movers acquire the same snapshot generation");
}

void TestFingerprintStability() {
  using SEP::Background::FingerprintConfiguration;
  const std::string first =
      FingerprintConfiguration("provider=analytic;domain=parker");
  const std::string second =
      FingerprintConfiguration("provider=analytic;domain=parker");
  const std::string changed =
      FingerprintConfiguration("provider=swmf;domain=parker");

  Check(first == second && first.size() == 16,
        "configuration fingerprint is deterministic and fixed-width");
  Check(first != changed,
        "provider-affecting configuration changes the fingerprint");
}

}  // namespace

int main() {
  TestConstructionAndImmutability();
  TestPublicationAndReadPhase();
  TestProviderIsolationAndHandoff();
  TestConcurrentMoverViews();
  TestFingerprintStability();

  if (failures != 0) {
    std::cerr << "Step 2 background-snapshot tests: " << failures
              << " failure(s)\n";
    return EXIT_FAILURE;
  }

  std::cout << "Step 2 background-snapshot tests: PASS\n";
  return EXIT_SUCCESS;
}
