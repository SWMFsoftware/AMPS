#include "test_harness.h"

#include "specfunc.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// This fixture exercises the production inline generator and distribution
// algorithms without calling rnd_seed(), whose rank-derived default is an MPI
// integration contract.  Directly resetting the same production seed state
// supplies deterministic one-process streams and keeps the statistical oracle
// independent of MPI startup behavior.
int ThisThread = 0;
int TotalThreadsNumber = 1;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;
unsigned long int RandomNumberGenerator::rndLastSeed = 0;
unsigned long int* RandomNumberGenerator::rndLastSeedArray = NULL;
thread_local std::mt19937 RandomNumberGenerator::gen(42);
thread_local std::uniform_real_distribution<double>
    RandomNumberGenerator::dist(0.0, 1.0);
thread_local bool RandomNumberGenerator::NormalSpareAvailable = false;
thread_local double RandomNumberGenerator::NormalSpare = 0.0;

long int nint(double value) {
  const long int lower = static_cast<long int>(value);
  return value - lower > 0.5 ? lower + 1 : lower;
}

void PrintErrorLog(const char*) {}
void PrintErrorLog(long int, const char*, const char*) {}

void exit(long int line, const char* file, const char* message) {
  std::string text = AMPS::Test::Location(file, static_cast<int>(line));
  text += message == NULL ? "AMPS trap exit" : message;
  throw AMPS::Test::Failure(text);
}

namespace {

void ResetSeed(unsigned long seed) {
  RandomNumberGenerator::rndLastSeed = seed;
  RandomNumberGenerator::NormalSpareAvailable = false;
  RandomNumberGenerator::NormalSpare = 0.0;
}

struct Moments {
  double mean;
  double variance;
  double third_central;
  double fourth_central;
};

Moments GetMoments(const std::vector<double>& values) {
  long double sum = 0.0L;
  for (std::size_t i = 0; i < values.size(); ++i) sum += values[i];
  const long double mean = sum / values.size();
  long double second = 0.0L, third = 0.0L, fourth = 0.0L;
  for (std::size_t i = 0; i < values.size(); ++i) {
    const long double d = values[i] - mean;
    second += d * d;
    third += d * d * d;
    fourth += d * d * d * d;
  }
  Moments result;
  result.mean = static_cast<double>(mean);
  result.variance = static_cast<double>(second / values.size());
  result.third_central = static_cast<double>(third / values.size());
  result.fourth_central = static_cast<double>(fourth / values.size());
  return result;
}

template <class Callable>
void ExpectInvalidArgument(Callable callable) {
  bool rejected = false;
  try {
    callable();
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  AMPS_EXPECT_TRUE(rejected);
}

}  // namespace

AMPS_TEST_CASE("GENERAL-RNG-REPRO-001",
               "fixed seeds reproduce an independently tabulated recurrence") {
  // The constants are the exact integer recurrence
  // s_(n+1)=(48828125*s_n) mod 2^31, divided by 2^31, starting at 17.
  // They were evaluated in integer arithmetic, not by calling rnd twice.
  const double expected[8] = {
      0.38653524825349450, 0.41862766118720174,
      0.76890633488073945, 0.63284860597923398,
      0.83882978418841958, 0.55607517482712865,
      0.14585589105263352, 0.68030437128618360};
  cRndSeedContainer stream;
  stream.Seed = 17;
  for (int i = 0; i < 8; ++i) AMPS_EXPECT_NEAR(rnd(&stream), expected[i], 0.0);
  ResetSeed(17);
  for (int i = 0; i < 8; ++i) AMPS_EXPECT_NEAR(rnd(), expected[i], 0.0);
}

AMPS_TEST_CASE("GENERAL-RNG-ENDPOINT-001",
               "uniform samples strictly exclude both endpoints for edge seeds") {
  const unsigned long seeds[4] = {0UL, 1UL, 17UL, 2147483647UL};
  for (int stream_index = 0; stream_index < 4; ++stream_index) {
    cRndSeedContainer stream;
    stream.Seed = seeds[stream_index];
    for (int i = 0; i < 50000; ++i) {
      const double value = rnd(&stream);
      AMPS_EXPECT_TRUE(value > 0.0);
      AMPS_EXPECT_TRUE(value < 1.0);
    }
  }
}

AMPS_TEST_CASE("GENERAL-RNG-STREAMS-001",
               "explicit streams remain reproducible and independent when interleaved") {
  cRndSeedContainer reference_a, reference_b, interleaved_a, interleaved_b;
  reference_a.Seed = interleaved_a.Seed = 101;
  reference_b.Seed = interleaved_b.Seed = 202;
  std::vector<double> sequence_a(1000), sequence_b(1000);
  for (int i = 0; i < 1000; ++i) sequence_a[i] = rnd(&reference_a);
  for (int i = 0; i < 1000; ++i) sequence_b[i] = rnd(&reference_b);
  int equal_pairs = 0;
  for (int i = 0; i < 1000; ++i) {
    const double a = rnd(&interleaved_a);
    const double b = rnd(&interleaved_b);
    AMPS_EXPECT_NEAR(a, sequence_a[i], 0.0);
    AMPS_EXPECT_NEAR(b, sequence_b[i], 0.0);
    if (a == b) ++equal_pairs;
  }
  AMPS_EXPECT_EQ(equal_pairs, 0);
}

AMPS_TEST_CASE("GENERAL-RNG-UNIFORM-001",
               "uniform stream moments and bin occupancies match U(0,1)") {
  const int sample_count = 200000;
  const int bin_count = 20;
  int bins[bin_count] = {};
  std::vector<double> values;
  values.reserve(sample_count);
  cRndSeedContainer stream;
  stream.Seed = 314159;
  for (int i = 0; i < sample_count; ++i) {
    const double value = rnd(&stream);
    values.push_back(value);
    int bin = static_cast<int>(value * bin_count);
    if (bin == bin_count) bin = bin_count - 1;
    ++bins[bin];
  }
  const Moments moments = GetMoments(values);
  AMPS_EXPECT_NEAR(moments.mean, 0.5, 0.0030);
  AMPS_EXPECT_NEAR(moments.variance, 1.0 / 12.0, 0.0008);
  const double expected = static_cast<double>(sample_count) / bin_count;
  double chi_square = 0.0;
  for (int i = 0; i < bin_count; ++i) {
    const double delta = bins[i] - expected;
    chi_square += delta * delta / expected;
  }
  // For 19 degrees of freedom, 60 is far above the 1-1e-5 quantile. This
  // loose deterministic gate catches gross recurrence/binning defects without
  // treating ordinary finite-sample noise as a regression.
  AMPS_EXPECT_TRUE(chi_square < 60.0);
}

AMPS_TEST_CASE("GENERAL-RNG-NORMAL-001",
               "Box-Muller Gaussian samples recover N(0,1) moments") {
  const int sample_count = 200000;
  ResetSeed(271828);
  std::vector<double> values(sample_count);
  Vector3D::Distribution::Normal(values.data(), sample_count);
  const Moments moments = GetMoments(values);
  AMPS_EXPECT_NEAR(moments.mean, 0.0, 0.015);
  AMPS_EXPECT_NEAR(moments.variance, 1.0, 0.03);
  AMPS_EXPECT_NEAR(moments.third_central, 0.0, 0.08);
  AMPS_EXPECT_NEAR(moments.fourth_central, 3.0, 0.18);

  // The public bulk API treats an empty null buffer as a no-op, but rejects
  // invalid counts and nonempty null buffers before any pointer arithmetic.
  Vector3D::Distribution::Normal(static_cast<double*>(NULL), 0);
  double sentinel = 17.0;
  ExpectInvalidArgument([&sentinel]() {
    Vector3D::Distribution::Normal(&sentinel, -1);
  });
  ExpectInvalidArgument([]() {
    Vector3D::Distribution::Normal(static_cast<double*>(NULL), 1);
  });
  AMPS_EXPECT_NEAR(sentinel, 17.0, 0.0);
}

AMPS_TEST_CASE("GENERAL-RNG-GAMMA-001",
               "gamma sampler recovers moments and rejects invalid shape/rate") {
  struct Fixture { double shape; double rate; unsigned long seed; };
  const Fixture fixtures[4] = {
      {0.5, 1.0, 4001}, {1.0, 1.0, 4002},
      {2.0, 1.0, 4003}, {5.0, 2.0, 4004}};
  const int sample_count = 120000;
  for (int fixture = 0; fixture < 4; ++fixture) {
    ResetSeed(fixtures[fixture].seed);
    std::vector<double> values(sample_count);
    for (int i = 0; i < sample_count; ++i) {
      values[i] = Vector3D::Distribution::Gamma(
          fixtures[fixture].shape, fixtures[fixture].rate);
      AMPS_EXPECT_TRUE(values[i] > 0.0);
    }
    const Moments moments = GetMoments(values);
    const double exact_mean = fixtures[fixture].shape / fixtures[fixture].rate;
    const double exact_variance = fixtures[fixture].shape /
        (fixtures[fixture].rate * fixtures[fixture].rate);
    AMPS_EXPECT_NEAR(moments.mean, exact_mean,
                     0.025 * std::max(1.0, exact_mean));
    AMPS_EXPECT_NEAR(moments.variance, exact_variance,
                     0.07 * std::max(1.0, exact_variance));
  }
  const double invalid_shapes[4] = {
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};
  for (int i = 0; i < 4; ++i) {
    ExpectInvalidArgument([&invalid_shapes, i]() {
      (void)Vector3D::Distribution::Gamma(invalid_shapes[i]);
    });
  }

  const double invalid_rates[4] = {
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};
  for (int i = 0; i < 4; ++i) {
    ExpectInvalidArgument([&invalid_rates, i]() {
      (void)Vector3D::Distribution::Gamma(2.0, invalid_rates[i]);
    });
  }
}

AMPS_TEST_CASE("GENERAL-RNG-DIRECTION-001",
               "uniform spherical directions have unit radius and isotropic moments") {
  const int sample_count = 150000;
  ResetSeed(161803);
  long double mean[3] = {};
  long double second[3] = {};
  for (int i = 0; i < sample_count; ++i) {
    double direction[3];
    Vector3D::Distribution::Uniform(direction);
    AMPS_EXPECT_NEAR(Vector3D::Length(direction), 1.0, 2.0e-14);
    for (int d = 0; d < 3; ++d) {
      mean[d] += direction[d];
      second[d] += direction[d] * direction[d];
    }
  }
  for (int d = 0; d < 3; ++d) {
    AMPS_EXPECT_NEAR(static_cast<double>(mean[d] / sample_count), 0.0, 0.01);
    AMPS_EXPECT_NEAR(static_cast<double>(second[d] / sample_count), 1.0 / 3.0,
                     0.01);
  }
}
