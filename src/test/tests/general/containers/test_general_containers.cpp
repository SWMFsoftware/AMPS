#define _GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_ 2
#define _GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_ 4

#include "test_harness.h"

#include "array_1d.h"
#include "bitArray.h"
#include "indexer_4d.h"
#include "linked_container.h"
#include "stack.h"
#include "string_func.h"

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// The tested headers use the AMPS trap-style exit signature.  Unit tests do
// not exercise fatal paths here, but defining the symbol makes this
// dependency-light executable independent of the production error logger.
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

AMPS_TEST_CASE("GENERAL-ARRAY-001", "array_1d arithmetic and vector identities") {
  // The hand-computed cross product is (-10,-7,8).  Orthogonality supplies two
  // independent identities, while the later mutation/resizing operations pin
  // owning value semantics and expose shallow-copy/double-free regressions.
  array_1d<double> a(3), b(3);
  a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;
  b(0) = -2.0; b(1) = 4.0; b(2) = 1.0;
  AMPS_EXPECT_EQ(a.GetSize(), 3);
  AMPS_EXPECT_NEAR(a.sum(), 6.0, 1.0e-14);
  AMPS_EXPECT_NEAR(dot_product(a, b), 9.0, 1.0e-14);
  array_1d<double> cross = cross_product(a, b);
  AMPS_EXPECT_NEAR(cross(0), -10.0, 1.0e-14);
  AMPS_EXPECT_NEAR(cross(1), -7.0, 1.0e-14);
  AMPS_EXPECT_NEAR(cross(2), 8.0, 1.0e-14);
  AMPS_EXPECT_NEAR(dot_product(cross, a), 0.0, 1.0e-14);
  AMPS_EXPECT_NEAR(dot_product(cross, b), 0.0, 1.0e-14);
  array_1d<double> copied(cross);
  cross(0) = 100.0;
  AMPS_EXPECT_NEAR(copied(0), -10.0, 1.0e-14);
  array_1d<double> assigned(1);
  assigned = copied;
  AMPS_EXPECT_EQ(assigned.GetSize(), 3);
  AMPS_EXPECT_NEAR(assigned(2), 8.0, 1.0e-14);
  array_1d<double> moved(std::move(assigned));
  AMPS_EXPECT_EQ(moved.GetSize(), 3);
  AMPS_EXPECT_EQ(assigned.GetSize(), 0);
  a.normalize();
  AMPS_EXPECT_NEAR(a.abs(), 1.0, 1.0e-14);
}

AMPS_TEST_CASE("GENERAL-BITARRAY-001", "bit-array word boundaries and range-rejecting setters") {
  // Indices 31/32 straddle a 32-bit storage word; 69 exercises the partial
  // final word.  Invalid setters must reject input without touching storage.
  cBitArray bits(70);
  AMPS_EXPECT_FALSE(bits.Get(0));
  AMPS_EXPECT_TRUE(bits.SetTrue(0));
  AMPS_EXPECT_TRUE(bits.SetTrue(31));
  AMPS_EXPECT_TRUE(bits.SetTrue(32));
  AMPS_EXPECT_TRUE(bits.SetTrue(69));
  AMPS_EXPECT_TRUE(bits.Get(0));
  AMPS_EXPECT_TRUE(bits.Get(31));
  AMPS_EXPECT_TRUE(bits.Get(32));
  AMPS_EXPECT_TRUE(bits.Get(69));
  AMPS_EXPECT_TRUE(bits.SetFalse(32));
  AMPS_EXPECT_FALSE(bits.Get(32));
  AMPS_EXPECT_FALSE(bits.SetTrue(-1));
  AMPS_EXPECT_FALSE(bits.SetFalse(70));
}

AMPS_TEST_CASE("GENERAL-STACK-001", "stack growth, pointer numbering, deletion, and reuse") {
  // Test constants force two four-element chunks.  Deleting the middle entry
  // verifies free-list reuse without relying on allocator address accidents.
  cStack<int> stack;
  std::vector<int*> pointers;
  for (int i = 0; i < 7; ++i) {
    int* item = stack.newElement();
    *item = 10 + i;
    pointers.push_back(item);
    AMPS_EXPECT_EQ(stack.GetEntryCountingNumber(item), i);
    AMPS_EXPECT_TRUE(stack.GetEntryPointer(i) == item);
  }
  AMPS_EXPECT_EQ(stack.usedElements(), 7);
  AMPS_EXPECT_EQ(stack.totalSize(), 8);
  stack.deleteElement(pointers[3]);
  AMPS_EXPECT_EQ(stack.usedElements(), 6);
  int* recycled = stack.newElement();
  AMPS_EXPECT_TRUE(recycled == pointers[3]);
  AMPS_EXPECT_EQ(stack.usedElements(), 7);
}

AMPS_TEST_CASE("GENERAL-LINKED-001", "linked-container front/back insertion and iteration order") {
  // The mixed insertion sequence has an unambiguous expected traversal 1,2,3;
  // this catches broken head/tail updates as well as iterator termination.
  cLinkedContainer<int> values;
  int two = 2, one = 1, three = 3;
  values.push_back(two);
  values.push_front(one);
  values.push_back(three);
  std::vector<int> observed;
  for (cLinkedContainer<int>::iterator it = values.begin(); it != values.end(); ++it) {
    observed.push_back(*it);
  }
  AMPS_EXPECT_EQ(observed.size(), static_cast<std::size_t>(3));
  AMPS_EXPECT_EQ(observed[0], 1);
  AMPS_EXPECT_EQ(observed[1], 2);
  AMPS_EXPECT_EQ(observed[2], 3);
}

AMPS_TEST_CASE("GENERAL-STRING-001", "string replacement, trimming, and first-word extraction") {
  // Exercise first-only and all-occurrence replacement separately, then make
  // CutFirstWord prove that it updates both the result and remaining input.
  std::string value = "  alpha beta alpha \t";
  trim(value);
  AMPS_EXPECT_EQ(value, std::string("alpha beta alpha"));
  AMPS_EXPECT_TRUE(FindAndReplaceFirst(value, "alpha", "A"));
  AMPS_EXPECT_EQ(value, std::string("A beta alpha"));
  AMPS_EXPECT_TRUE(FindAndReplaceAll(value, "alpha", "B"));
  AMPS_EXPECT_EQ(value, std::string("A beta B"));
  std::string word;
  AMPS_EXPECT_TRUE(CutFirstWord(word, value));
  AMPS_EXPECT_EQ(word, std::string("A"));
  AMPS_EXPECT_EQ(value, std::string("beta B"));
}

AMPS_TEST_CASE("GENERAL-INDEXER4D-001", "sparse 4-D indexing, reverse lookup, and reset") {
  // Capacity is the independent product 3*5*2*3=90.  Assigned indices are
  // insertion-order labels, not dense-coordinate offsets; reverse lookup and
  // reset verify both sides of that mapping contract.
  cIndexer4D indexer;
  int minimum[4] = {-1, -2, 0, 0};
  int maximum[4] = {1, 2, 1, 2};
  indexer.init(minimum, maximum);
  AMPS_EXPECT_TRUE(indexer.isInitialized());
  AMPS_EXPECT_EQ(indexer.capacity(), 90);
  AMPS_EXPECT_EQ(indexer.get_idx(0, 0, 0, 0), 0);
  AMPS_EXPECT_EQ(indexer.get_idx(1, -2, 1, 2), 1);
  AMPS_EXPECT_EQ(indexer.get_idx(0, 0, 0, 0), 0);
  AMPS_EXPECT_EQ(indexer.get_index_count(), 2);
  AMPS_EXPECT_TRUE(indexer.is_assigned(1, -2, 1, 2));
  AMPS_EXPECT_FALSE(indexer.is_assigned(-1, 1, 1, 1));
  int i = 0, j = 0, k = 0, m = 0;
  AMPS_EXPECT_TRUE(indexer.find_idx(1, i, j, k, m));
  AMPS_EXPECT_EQ(i, 1);
  AMPS_EXPECT_EQ(j, -2);
  AMPS_EXPECT_EQ(k, 1);
  AMPS_EXPECT_EQ(m, 2);
  indexer.reset();
  AMPS_EXPECT_EQ(indexer.get_index_count(), 0);
  AMPS_EXPECT_EQ(indexer.get_idx(-1, 2, 1, 1), 0);
}
