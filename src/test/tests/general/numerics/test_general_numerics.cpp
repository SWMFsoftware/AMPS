#define _GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_ 2
#define _GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_ 4

#include "test_harness.h"

#include "array_1d.h"
#include "array_2d.h"
#include "array_3d.h"
#include "array_4d.h"
#include "array_5d.h"
#include "quadrature.h"
#include "specfunc.h"
#include "stack_managed.h"

#include <cmath>
#include <climits>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <thread>
#include <array>

// The general headers route invalid input through AMPS' trap-style exit.
// Replacing only that terminal action with an exception preserves the exact
// branch and message while allowing a dependency-light case to assert the
// bounds policy without aborting the remaining registry.
int ThisThread = 0;
int TotalThreadsNumber = 1;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;

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

int polynomial_degree = 0;

double Polynomial(double* x) {
  return std::pow(x[0], polynomial_degree);
}

double ProductPolynomial(double* x) {
  return x[0] * x[0] * x[1] * x[1] * x[1];
}

double Exponential(double* x) {
  return std::exp(x[0]);
}

double EndpointSquareRoot(double* x) {
  return std::sqrt(1.0 + x[0]);
}

double ExactPowerIntegral(double lower, double upper, int degree) {
  // Antiderivative x^(p+1)/(p+1), evaluated independently of the production
  // quadrature nodes and weights.
  const int exponent = degree + 1;
  return (std::pow(upper, exponent) - std::pow(lower, exponent)) / exponent;
}

}  // namespace

AMPS_TEST_CASE("GENERAL-NUM-ARRAY-EMPTY-001",
               "default arrays expose a stable unallocated empty state") {
  array_1d<double> a1;
  array_2d<double> a2;
  array_3d<double> a3;
  array_4d<double> a4;
  array_5d<double> a5;
  AMPS_EXPECT_FALSE(a1.IsAllocated());
  AMPS_EXPECT_EQ(a1.GetSize(), 0);
  AMPS_EXPECT_NEAR(a1.sum(), 0.0, 0.0);
  AMPS_EXPECT_FALSE(a2.IsAllocated());
  AMPS_EXPECT_EQ(a2.GetElementNumber(), 0);
  AMPS_EXPECT_FALSE(a3.IsAllocated());
  AMPS_EXPECT_EQ(a3.size(), 0);
  AMPS_EXPECT_EQ(a4.size(), 0);
  AMPS_EXPECT_EQ(a5.size(), 0);
}

AMPS_TEST_CASE("GENERAL-NUM-ARRAY-INDEX-001",
               "asymmetric 1-D through 5-D arrays use documented linear layouts") {
  // Distinct prime-like extents make every stride observable; symmetric cubes
  // could allow two swapped dimensions to pass accidentally.
  array_1d<int> a1(7);
  for (int i = 0; i < 7; ++i) a1(i) = 10 + i;
  AMPS_EXPECT_EQ(a1.get_data_ptr()[6], 16);

  array_2d<int> a2(2, 5);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j) a2(i, j) = 100 * i + j;
  AMPS_EXPECT_EQ(a2.GetBufferPointer()[1 * 5 + 4], 104);
  AMPS_EXPECT_TRUE(a2(1) == a2.GetBufferPointer() + 5);

  array_3d<int> a3(2, 3, 5);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 5; ++k) a3(i, j, k) = 100 * i + 10 * j + k;
  AMPS_EXPECT_EQ(a3(1, 2)[4], 124);
  AMPS_EXPECT_EQ(a3.size(0), 2);
  AMPS_EXPECT_EQ(a3.size(1), 3);
  AMPS_EXPECT_EQ(a3.size(2), 5);

  array_4d<int> a4(2, 3, 4, 5);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 4; ++k)
        for (int m = 0; m < 5; ++m)
          a4(i, j, k, m) = 1000 * i + 100 * j + 10 * k + m;
  AMPS_EXPECT_EQ(a4.size(), 120);
  AMPS_EXPECT_EQ(*a4.GetPtr(1, 2, 3, 4), 1234);
  AMPS_EXPECT_EQ(a4.GetPtr(1, 2, 3, 4) - a4.GetPtr(),
                 1 + 2 * 2 + 3 * 2 * 3 + 4 * 2 * 3 * 4);
  array_4d<int> a4_offset(2, 3, 4, 5);
  a4_offset = 1;
  array_4d<int> a4_difference = a4 - a4_offset;
  AMPS_EXPECT_EQ(a4_difference.size(), 120);
  AMPS_EXPECT_EQ(a4_difference(1, 2, 3, 4), 1233);

  array_5d<int> a5(2, 3, 2, 4, 3);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
        for (int m = 0; m < 4; ++m)
          for (int n = 0; n < 3; ++n)
            a5(i, j, k, m, n) = 10000 * i + 1000 * j + 100 * k + 10 * m + n;
  AMPS_EXPECT_EQ(a5.size(), 144);
  AMPS_EXPECT_EQ(*a5.GetPtr(1, 2, 1, 3, 2), 12132);
  for (int d = 0; d < 5; ++d) {
    const int expected[5] = {2, 3, 2, 4, 3};
    AMPS_EXPECT_EQ(a5.size(d), expected[d]);
  }
}

AMPS_TEST_CASE("GENERAL-NUM-ARRAY-OWNERSHIP-001",
               "owning copies, resizing assignment, moves, and external aliases do not overlap") {
  array_1d<int> original(3);
  original(0) = 4; original(1) = 5; original(2) = 6;
  array_1d<int> copied(original);
  copied(0) = 40;
  AMPS_EXPECT_EQ(original(0), 4);
  array_1d<int> resized(1);
  resized = original;
  AMPS_EXPECT_EQ(resized.GetSize(), 3);
  array_1d<int> moved(std::move(resized));
  AMPS_EXPECT_EQ(moved(2), 6);
  AMPS_EXPECT_EQ(resized.GetSize(), 0);

  array_2d<int> source2(2, 5);
  source2 = 3;
  array_2d<int> copied2(source2);
  copied2(1, 4) = 9;
  AMPS_EXPECT_EQ(source2(1, 4), 3);
  array_2d<int> assigned2(1, 1);
  assigned2 = source2;
  AMPS_EXPECT_EQ(assigned2.GetElementNumber(), 10);
  array_2d<int> moved2(std::move(copied2));
  AMPS_EXPECT_EQ(moved2(1, 4), 9);
  AMPS_EXPECT_EQ(copied2.GetElementNumber(), 0);

  array_3d<int> source3(2, 3, 5);
  source3 = 4;
  array_3d<int> copied3(source3);
  copied3(1, 2, 4) = 10;
  AMPS_EXPECT_EQ(source3(1, 2, 4), 4);
  array_3d<int> assigned3(1, 1, 1);
  assigned3 = source3;
  AMPS_EXPECT_EQ(assigned3.size(), 30);
  array_3d<int> moved3(std::move(copied3));
  AMPS_EXPECT_EQ(moved3(1, 2, 4), 10);
  AMPS_EXPECT_EQ(copied3.size(), 0);

  // The 2-D/3-D manager constructors are explicitly non-owning views. Their
  // destructor must leave caller storage alive while indexed writes alias it.
  int external2[6] = {0, 1, 2, 3, 4, 5};
  {
    array_2d<int> view(external2, 2, 3);
    view(1, 2) = 51;
  }
  AMPS_EXPECT_EQ(external2[5], 51);
  int external3[12] = {};
  {
    array_3d<int> view(external3, 2, 2, 3);
    view(1, 1, 2) = 112;
  }
  AMPS_EXPECT_EQ(external3[11], 112);

  int external4[16] = {};
  {
    array_4d<int> view(external4, 2, 2, 2, 2);
    view(1, 1, 1, 1) = 41;
  }
  AMPS_EXPECT_EQ(external4[15], 41);
  int external5[32] = {};
  {
    array_5d<int> view(external5, 2, 2, 2, 2, 2);
    view(1, 1, 1, 1, 1) = 51;
  }
  AMPS_EXPECT_EQ(external5[31], 51);

  // Multi-dimensional assignment is supported for equal extents. Separate
  // destination buffers prove that assignment copies values rather than
  // transferring or aliasing ownership.
  array_4d<int> source4(2, 3, 2, 2), target4(2, 3, 2, 2);
  source4 = 7;
  target4 = source4;
  target4(1, 2, 1, 1) = 19;
  AMPS_EXPECT_EQ(source4(1, 2, 1, 1), 7);
  array_4d<int> copied4(source4);
  copied4(1, 2, 1, 1) = 29;
  AMPS_EXPECT_EQ(source4(1, 2, 1, 1), 7);
  array_4d<int> moved4(std::move(copied4));
  AMPS_EXPECT_EQ(moved4(1, 2, 1, 1), 29);
  AMPS_EXPECT_EQ(copied4.size(), 0);
  array_5d<int> source5(2, 2, 2, 2, 3), target5(2, 2, 2, 2, 3);
  source5 = 11;
  target5 = source5;
  target5(1, 1, 1, 1, 2) = 23;
  AMPS_EXPECT_EQ(source5(1, 1, 1, 1, 2), 11);
  array_5d<int> copied5(source5);
  copied5(1, 1, 1, 1, 2) = 31;
  AMPS_EXPECT_EQ(source5(1, 1, 1, 1, 2), 11);
  array_5d<int> moved5(std::move(copied5));
  AMPS_EXPECT_EQ(moved5(1, 1, 1, 1, 2), 31);
  AMPS_EXPECT_EQ(copied5.size(), 0);
}

AMPS_TEST_CASE("GENERAL-NUM-ARRAY-RESIZE-001",
               "every array rank resizes asymmetrically and returns to canonical empty state") {
  array_1d<int> a1(2);
  a1.resize(7);
  AMPS_EXPECT_EQ(a1.GetSize(), 7);
  a1.resize(0);
  AMPS_EXPECT_EQ(a1.GetSize(), 0);
  AMPS_EXPECT_FALSE(a1.IsAllocated());

  array_2d<int> a2(1, 1);
  a2.resize(2, 5);
  AMPS_EXPECT_EQ(a2.GetDim0(), 2);
  AMPS_EXPECT_EQ(a2.GetDim1(), 5);
  a2.resize(0, 0);
  AMPS_EXPECT_EQ(a2.GetElementNumber(), 0);
  AMPS_EXPECT_FALSE(a2.IsAllocated());

  array_3d<int> a3(1, 1, 1);
  a3.resize(2, 3, 5);
  AMPS_EXPECT_EQ(a3.size(), 30);
  a3.resize(0, 0, 0);
  AMPS_EXPECT_EQ(a3.size(), 0);
  AMPS_EXPECT_FALSE(a3.IsAllocated());

  array_4d<int> a4(1, 1, 1, 1);
  a4.resize(2, 3, 4, 5);
  AMPS_EXPECT_EQ(a4.size(), 120);
  AMPS_EXPECT_EQ(a4.size(3), 5);
  a4.resize(0, 0, 0, 0);
  AMPS_EXPECT_EQ(a4.size(), 0);

  array_5d<int> a5(1, 1, 1, 1, 1);
  a5.resize(2, 3, 2, 4, 3);
  AMPS_EXPECT_EQ(a5.size(), 144);
  a5.resize(0, 0, 0, 0, 0);
  AMPS_EXPECT_EQ(a5.size(), 0);
}

AMPS_TEST_CASE("GENERAL-NUM-ARRAY-BOUNDS-001",
               "array ranks reject invalid indices and overflowing extent products") {
  int bounds_trapped = 0;
  array_1d<int> a1(2);
  array_2d<int> a2(2, 3);
  array_3d<int> a3(2, 3, 4);
  array_4d<int> a4(2, 3, 4, 5);
  array_5d<int> a5(2, 3, 2, 4, 3);
  try { (void)a1(2); } catch (const AMPS::Test::Failure&) { ++bounds_trapped; }
  try { (void)a2(0, 3); } catch (const AMPS::Test::Failure&) { ++bounds_trapped; }
  try { (void)a3(-1, 0, 0); } catch (const AMPS::Test::Failure&) { ++bounds_trapped; }
  try { (void)a4(0, 0, 4, 0); } catch (const AMPS::Test::Failure&) { ++bounds_trapped; }
  try { (void)a5(0, 0, 0, 0, 3); } catch (const AMPS::Test::Failure&) { ++bounds_trapped; }
  AMPS_EXPECT_EQ(bounds_trapped, 5);

  // Every probe is rejected before allocation.  INT_MAX*2 used to overflow
  // as signed int in the allocation expression and could create a much
  // smaller buffer than the indexing strides advertised.
  int extent_trapped = 0;
  int expected_extent_traps = 4;
  array_2d<int> oversized2;
  array_3d<int> oversized3;
  array_4d<int> oversized4;
  array_5d<int> oversized5;
  try { oversized2.resize(INT_MAX, 2); }
  catch (const AMPS::Test::Failure&) { ++extent_trapped; }
  try { oversized3.resize(INT_MAX, 2, 1); }
  catch (const AMPS::Test::Failure&) { ++extent_trapped; }
  try { oversized4.resize(INT_MAX, 2, 1, 1); }
  catch (const AMPS::Test::Failure&) { ++extent_trapped; }
  try { oversized5.resize(INT_MAX, 2, 1, 1, 1); }
  catch (const AMPS::Test::Failure&) { ++extent_trapped; }

  // array_1d accepts long for historical compatibility but stores an int
  // extent.  Exercise the narrowing guard where long is wider than int.
  const unsigned long beyond_int = static_cast<unsigned long>(INT_MAX) + 1UL;
  if (beyond_int <= static_cast<unsigned long>(LONG_MAX)) {
    ++expected_extent_traps;
    array_1d<int> oversized1;
    try { oversized1.resize(static_cast<long>(beyond_int)); }
    catch (const AMPS::Test::Failure&) { ++extent_trapped; }
  }
  AMPS_EXPECT_EQ(extent_trapped, expected_extent_traps);
}

AMPS_TEST_CASE("GENERAL-NUM-MANAGED-STACK-001",
               "managed stack grows, recycles deleted entries, resets, and releases storage") {
  cStackManaged<int> stack;
  std::vector<int*> pointers;
  for (int i = 0; i < 11; ++i) {
    int* item = stack.newElement();
    *item = 100 + i;
    pointers.push_back(item);
    AMPS_EXPECT_EQ(stack.GetEntryCountingNumber(item), i);
    AMPS_EXPECT_TRUE(stack.GetEntryPointer(i) == item);
  }
  AMPS_EXPECT_EQ(stack.usedElements(), 11);
  AMPS_EXPECT_EQ(stack.totalSize(), 12);
  AMPS_EXPECT_TRUE(stack.getAllocatedMemory() > 0);
  stack.deleteElement(pointers[4]);
  int* recycled = stack.newElement();
  AMPS_EXPECT_TRUE(recycled == pointers[4]);
  stack.resetStack();
  AMPS_EXPECT_EQ(stack.usedElements(), 0);
  AMPS_EXPECT_TRUE(stack.newElement() == stack.GetEntryPointer(0));
  stack.clear();
  AMPS_EXPECT_EQ(stack.totalSize(), 0);
  AMPS_EXPECT_EQ(stack.usedElements(), 0);
  AMPS_EXPECT_EQ(stack.getAllocatedMemory(), 0);
}

AMPS_TEST_CASE("GENERAL-NUM-SPECFUNC-001",
               "integer and vector helpers satisfy independent algebraic identities") {
  AMPS_EXPECT_EQ(GetGcd(48, 18), 6);
  AMPS_EXPECT_EQ(GetGcd(-48, 18), 6);
  AMPS_EXPECT_EQ(GetGcd(0, 27), 27);
  AMPS_EXPECT_EQ(GetGcd(INT_MIN, 2), 2);
  int gcd_overflow = 0;
  try { (void)GetGcd(INT_MIN, 0); }
  catch (const std::overflow_error&) { ++gcd_overflow; }
  AMPS_EXPECT_EQ(gcd_overflow, 1);
  double a[3] = {1.0, 2.0, 3.0};
  double b[3] = {-2.0, 4.0, 1.0};
  double cross[3];
  Vector3D::CrossProduct(cross, a, b);
  AMPS_EXPECT_NEAR(cross[0], -10.0, 1.0e-14);
  AMPS_EXPECT_NEAR(cross[1], -7.0, 1.0e-14);
  AMPS_EXPECT_NEAR(cross[2], 8.0, 1.0e-14);
  AMPS_EXPECT_NEAR(Vector3D::DotProduct(cross, a), 0.0, 1.0e-14);
  AMPS_EXPECT_NEAR(Vector3D::DotProduct(cross, b), 0.0, 1.0e-14);
  double v[3] = {3.0, 4.0, 0.0};
  AMPS_EXPECT_NEAR(Vector3D::Normalize(v, 2.0), 5.0, 1.0e-14);
  AMPS_EXPECT_NEAR(Vector3D::Length(v), 2.0, 1.0e-14);
}

AMPS_TEST_CASE("GENERAL-NUM-QUADRATURE-POLY-001",
               "Gauss-Legendre rules integrate their full polynomial degree exactly") {
  double lower[1] = {-2.0};
  double upper[1] = {3.0};
  for (int order = 2; order <= 8; order += 2) {
    for (polynomial_degree = 0; polynomial_degree <= 2 * order - 1;
         ++polynomial_degree) {
      const double value = Quadrature::Gauss::Cube::GaussLegendre(
          1, order, Polynomial, lower, upper);
      const double exact = ExactPowerIntegral(lower[0], upper[0], polynomial_degree);
      const double tolerance = 2.0e-11 * std::max(1.0, std::fabs(exact));
      AMPS_EXPECT_NEAR(value, exact, tolerance);
    }
  }
  double lower2[2] = {-1.0, -2.0};
  double upper2[2] = {2.0, 1.0};
  const double value2 = Quadrature::Gauss::Cube::GaussLegendre(
      2, 4, ProductPolynomial, lower2, upper2);
  const double exact2 = ExactPowerIntegral(-1.0, 2.0, 2) *
                        ExactPowerIntegral(-2.0, 1.0, 3);
  AMPS_EXPECT_NEAR(value2, exact2, 2.0e-12);

  // Order five selects the maintained odd-order table path; 21 and 23 are
  // deliberately absent from the static table and exercise runtime Newton
  // generation.  A Gauss rule of order n must integrate all degrees through
  // 2n-1, independently of how its nodes were obtained.
  double symmetric_lower[1] = {-1.0};
  double symmetric_upper[1] = {1.0};
  const int additional_orders[3] = {5, 21, 23};
  for (int order : additional_orders) {
    for (polynomial_degree = 0; polynomial_degree <= 2 * order - 1;
         ++polynomial_degree) {
      const double value = Quadrature::Gauss::Cube::GaussLegendre(
          1, order, Polynomial, symmetric_lower, symmetric_upper);
      const double exact = (polynomial_degree & 1)
          ? 0.0 : 2.0 / static_cast<double>(polynomial_degree + 1);
      AMPS_EXPECT_NEAR(value, exact, 2.0e-9);
    }

    const int half_length = (order + 1) >> 1;
    std::vector<double> nodes(half_length), weights(half_length);
    Quadrature::Gauss::Cube::Generate_GaussLegendreTable(
        order, nodes.data(), weights.data(), 1.0e-13);
    double full_weight_sum = (order & 1) ? weights[0] : 0.0;
    for (int i = 0; i < half_length; ++i) {
      AMPS_EXPECT_TRUE(std::isfinite(nodes[i]));
      AMPS_EXPECT_TRUE(std::isfinite(weights[i]));
      AMPS_EXPECT_TRUE(nodes[i] >= 0.0 && nodes[i] < 1.0);
      AMPS_EXPECT_TRUE(weights[i] > 0.0);
      if (i > 0) AMPS_EXPECT_TRUE(nodes[i] > nodes[i - 1]);
      if (!(order & 1) || i > 0) full_weight_sum += 2.0 * weights[i];
    }
    if (order & 1) AMPS_EXPECT_NEAR(nodes[0], 0.0, 2.0e-14);
    AMPS_EXPECT_NEAR(full_weight_sum, 2.0, 2.0e-11);
  }
}

AMPS_TEST_CASE("GENERAL-NUM-QUADRATURE-CONVERGENCE-001",
               "quadrature error decreases for smooth and endpoint-sensitive functions") {
  double lower[1] = {-1.0};
  double upper[1] = {1.0};
  const int orders[5] = {2, 4, 8, 16, 32};
  const double exact_exp = std::exp(1.0) - std::exp(-1.0);
  const double exact_sqrt = std::pow(2.0, 2.5) / 3.0;
  double previous_exp = 1.0;
  double previous_sqrt = 1.0;
  for (int i = 0; i < 5; ++i) {
    const double exp_value = Quadrature::Gauss::Cube::GaussLegendre(
        1, orders[i], Exponential, lower, upper);
    const double sqrt_value = Quadrature::Gauss::Cube::GaussLegendre(
        1, orders[i], EndpointSquareRoot, lower, upper);
    const double exp_error = std::fabs(exp_value - exact_exp);
    const double sqrt_error = std::fabs(sqrt_value - exact_sqrt);
    // Once the error reaches double-precision roundoff, a higher rule may be
    // equal or a few ulps worse; require convergence until that floor and a
    // tight absolute bound thereafter.
    if (previous_exp > 1.0e-13) AMPS_EXPECT_TRUE(exp_error < previous_exp);
    else AMPS_EXPECT_TRUE(exp_error < 2.0e-13);
    AMPS_EXPECT_TRUE(sqrt_error < previous_sqrt);
    previous_exp = exp_error;
    previous_sqrt = sqrt_error;
  }
  AMPS_EXPECT_TRUE(previous_exp < 1.0e-13);
  AMPS_EXPECT_TRUE(previous_sqrt < 2.0e-5);
}

AMPS_TEST_CASE("GENERAL-NUM-QUADRATURE-DOMAIN-001",
               "quadrature accepts oriented/zero domains and rejects invalid input") {
  double lower[1] = {-1.0};
  double upper[1] = {1.0};
  int rejected = 0;
  try {
    (void)Quadrature::Gauss::Cube::GaussLegendre(
        0, 4, Exponential, lower, upper);
  } catch (const std::invalid_argument&) { ++rejected; }
  try {
    (void)Quadrature::Gauss::Cube::GaussLegendre(
        1, 0, Exponential, lower, upper);
  } catch (const std::invalid_argument&) { ++rejected; }
  try {
    (void)Quadrature::Gauss::Cube::GaussLegendre(
        1, 4, NULL, lower, upper);
  } catch (const std::invalid_argument&) { ++rejected; }
  double nonfinite_lower[1] = {
      std::numeric_limits<double>::quiet_NaN()};
  try {
    (void)Quadrature::Gauss::Cube::GaussLegendre(
        1, 4, Exponential, nonfinite_lower, upper);
  } catch (const std::invalid_argument&) { ++rejected; }
  double nonfinite_upper[1] = {
      std::numeric_limits<double>::infinity()};
  try {
    (void)Quadrature::Gauss::Cube::GaussLegendre(
        1, 4, Exponential, lower, nonfinite_upper);
  } catch (const std::invalid_argument&) { ++rejected; }
  AMPS_EXPECT_EQ(rejected, 5);

  const double forward = Quadrature::Gauss::Cube::GaussLegendre(
      1, 8, Exponential, lower, upper);
  const double reversed = Quadrature::Gauss::Cube::GaussLegendre(
      1, 8, Exponential, upper, lower);
  AMPS_EXPECT_NEAR(reversed, -forward, 2.0e-13);
  AMPS_EXPECT_NEAR(Quadrature::Gauss::Cube::GaussLegendre(
      1, 8, Exponential, lower, lower), 0.0, 0.0);
}

namespace {
// Independent bitwise construction of each byte's polynomial contribution.
// No production lookup table or slicing code is used. All operations retain
// unsigned-long width, including the historical nonmasked LP64 upper bits.
unsigned long ReferenceSignature(const unsigned char* data,std::size_t length,unsigned long state=0) {
  for (std::size_t i=0;i<length;i++) {
    unsigned long remainder=(((state>>24)^data[i])&255UL)<<24;
    for (int bit=0;bit<8;bit++)
      remainder=(remainder&0x80000000UL) ? (remainder<<1)^0x04c11db7UL:remainder<<1;
    state=(state<<8)^remainder;
  }
  return state;
}
}

AMPS_TEST_CASE("PERF-CRC-STREAM-001","optimized CRC matches independent full-width bit recurrence on unaligned/all-byte/long streams") {
  std::vector<unsigned char> data(131079);
  unsigned long random=0x12345678UL;
  for (std::size_t i=0;i<data.size();i++) {
    random=random*1664525UL+1013904223UL;
    data[i]=i<256 ? static_cast<unsigned char>(i):static_cast<unsigned char>(random>>16);
  }
  for (int offset=0;offset<4;offset++) {
    CRC32 actual;
    actual.crc_accum=~0UL;
    actual.add(data.data()+offset,data.size()-offset);
    AMPS_EXPECT_EQ(actual.checksum(),ReferenceSignature(data.data()+offset,data.size()-offset,~0UL));
  }
}

AMPS_TEST_CASE("PERF-CRC-SPLIT-001","chunk boundaries, scalar templates, clear and empty updates preserve the exact stream") {
  std::array<unsigned char,517> data{};
  for (std::size_t i=0;i<data.size();i++) data[i]=static_cast<unsigned char>(17*i+3);
  for (int chunk=1;chunk<=19;chunk++) {
    CRC32 actual;
    for (std::size_t offset=0;offset<data.size();offset+=std::min<std::size_t>(chunk,data.size()-offset))
      actual.add(data.data()+offset,std::min<std::size_t>(chunk,data.size()-offset));
    AMPS_EXPECT_EQ(actual.checksum(),ReferenceSignature(data.data(),data.size()));
    actual.clear(); actual.add(data.data(),0);
    AMPS_EXPECT_EQ(actual.checksum(),0UL);
  }
  unsigned long words[4]={0UL,~0UL,17UL,913UL};
  CRC32 scalar;
  for (auto word:words) scalar.add(word);
  AMPS_EXPECT_EQ(scalar.checksum(),ReferenceSignature(reinterpret_cast<unsigned char*>(words),sizeof(words)));
}

AMPS_TEST_CASE("PERF-CRC-COPY-001","copyable CRC objects retain the base table and independent update state") {
  unsigned char a[7]={0,255,7,129,33,1,2},b[5]={3,4,5,6,7};
  CRC32 original; original.add(a,7);
  CRC32 copy=original; copy.add(b,5);
  AMPS_EXPECT_EQ(original.checksum(),ReferenceSignature(a,7));
  AMPS_EXPECT_EQ(copy.checksum(),ReferenceSignature(b,5,ReferenceSignature(a,7)));
}

AMPS_TEST_CASE("PERF-CRC-THREAD-001","shared immutable tables initialize safely for independent host threads") {
  std::array<unsigned long,8> results{};
  std::vector<std::thread> threads;
  for (int id=0;id<8;id++) threads.emplace_back([&,id] {
    unsigned char bytes[17]; for (int j=0;j<17;j++) bytes[j]=static_cast<unsigned char>(id*31+j);
    CRC32 actual; actual.add(bytes,17); results[id]=actual.checksum();
  });
  for (auto& thread:threads) thread.join();
  for (int id=0;id<8;id++) {
    unsigned char bytes[17]; for (int j=0;j<17;j++) bytes[j]=static_cast<unsigned char>(id*31+j);
    AMPS_EXPECT_EQ(results[id],ReferenceSignature(bytes,17));
  }
}
