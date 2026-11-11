#include "test_harness.h"
#include "ordered_sparse_matvec.h"
#include <vector>

namespace {
struct Operator {
  std::vector<std::vector<double>> coefficient;
  std::vector<std::vector<int>> column;
  std::vector<double> x;
  void Build(amps_ordered_sparse::View& view) {
    view.Build(coefficient.size(),[&](std::size_t row) {return int(coefficient[row].size());},
        [&](std::size_t row,std::size_t i) {
          return amps_ordered_sparse::Entry{coefficient[row][i],&x[column[row][i]]};
        });
  }
  // Independent ordered dense-entry oracle, using indices rather than the
  // production pointer/offset arrays. Repeated and remote-like columns are
  // legal; row order is deliberately significant for cancellation.
  double Reference(std::size_t row) const {
    double answer=0;
    for (std::size_t i=0;i<column[row].size();i++) answer+=coefficient[row][i]*x[column[row][i]];
    return answer;
  }
};
}

AMPS_TEST_CASE("PERF-SPMV-ORDER-001","irregular and cancellation-sensitive rows preserve scalar accumulation order") {
  Operator op{{{1e16,1,-1e16},{},{2,-3,5}},{{0,1,2},{},{3,0,3}},{1,1,1,7}};
  amps_ordered_sparse::View view; op.Build(view);
  AMPS_EXPECT_NEAR(view.Row(0),0.0,0.0);
  for (int row=0;row<3;row++) AMPS_EXPECT_NEAR(view.Row(row),op.Reference(row),0.0);
}

AMPS_TEST_CASE("PERF-SPMV-REFRESH-001","new solves refresh coefficients/storage while iterative RHS values remain live") {
  Operator op{{{2,3},{7}},{{0,1},{2}},{1,4,6}};
  amps_ordered_sparse::View view; op.Build(view);
  for (int i=0;i<3;i++) op.x[i]+=2;
  for (int row=0;row<2;row++) AMPS_EXPECT_NEAR(view.Row(row),op.Reference(row),0.0);
  op.coefficient={{-3,4,1},{},{9,-7}}; op.column={{2,1,0},{},{3,1}};
  op.x={8,7,6,5}; op.Build(view);
  for (int row=0;row<3;row++) AMPS_EXPECT_NEAR(view.Row(row),op.Reference(row),0.0);
}

AMPS_TEST_CASE("PERF-SPMV-EMPTY-INVALID-001","empty snapshots and invalid rows, lengths or RHS pointers are rejected") {
  amps_ordered_sparse::View view;
  view.Build(0,[](std::size_t) {return 0;},[](std::size_t,std::size_t) {
      return amps_ordered_sparse::Entry{0,nullptr};});
  bool rejected=false;
  try {view.Row(0);} catch (const std::out_of_range&) {rejected=true;}
  AMPS_EXPECT_TRUE(rejected);
  rejected=false;
  try {view.Build(1,[](std::size_t) {return -1;},[](std::size_t,std::size_t) {
      return amps_ordered_sparse::Entry{0,nullptr};});}
  catch (const std::invalid_argument&) {rejected=true;}
  AMPS_EXPECT_TRUE(rejected);
  rejected=false;
  try {view.Build(1,[](std::size_t) {return 1;},[](std::size_t,std::size_t) {
      return amps_ordered_sparse::Entry{0,nullptr};});}
  catch (const std::invalid_argument&) {rejected=true;}
  AMPS_EXPECT_TRUE(rejected);
}
