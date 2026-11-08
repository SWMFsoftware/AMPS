#include "test_harness.h"
#include "boundary_row.h"

#include <cmath>

using PIC::FieldSolver::Electromagnetic::ECSIM::BoundaryRow::Policy;
using PIC::FieldSolver::Electromagnetic::ECSIM::BoundaryRow::Resolve;
using PIC::FieldSolver::Electromagnetic::ECSIM::BoundaryRow::Status;

AMPS_TEST_CASE("ECSIM-BC-PERIODIC-001",
               "periodic faces create no physical increment constraint row") {
  const Policy row=Resolve(true,1,0);
  AMPS_EXPECT_TRUE(row.status==Status::InteriorOrPeriodic);
  AMPS_EXPECT_EQ(row.entries,0);
}

AMPS_TEST_CASE("ECSIM-BC-DIRICHLET-001",
               "fixed-field boundary imposes delta-E equals zero") {
  const Policy row=Resolve(false,1,0);
  AMPS_EXPECT_TRUE(row.status==Status::Dirichlet);
  AMPS_EXPECT_EQ(row.entries,1);
  AMPS_EXPECT_NEAR(row.diagonal,1.0,0.0);
  AMPS_EXPECT_NEAR(row.donor_weight,0.0,0.0);
}

AMPS_TEST_CASE("ECSIM-BC-NEUMANN-001",
               "zero-normal-gradient row averages one to three inward donors") {
  for (int donors=1;donors<=3;++donors) {
    const Policy row=Resolve(false,0,donors);
    AMPS_EXPECT_TRUE(row.status==Status::Neumann);
    AMPS_EXPECT_EQ(row.entries,1+donors);
    AMPS_EXPECT_NEAR(row.diagonal,1.0,0.0);
    AMPS_EXPECT_NEAR(row.donor_weight,-1.0/donors,0.0);
    // A constant delta-E must be in the null space of the homogeneous
    // Neumann row, including multi-face edges and corners.
    AMPS_EXPECT_NEAR(row.diagonal+donors*row.donor_weight,0.0,2.0e-16);
  }
  AMPS_EXPECT_TRUE(Resolve(false,0,0).status==Status::Invalid);
  AMPS_EXPECT_TRUE(Resolve(false,0,4).status==Status::Invalid);
}

AMPS_TEST_CASE("ECSIM-BC-CONDUCTOR-UNSUPPORTED-001",
               "perfect-conductor closure is typed unsupported until component-aware physics exists") {
  const Policy row=Resolve(false,2,1);
  AMPS_EXPECT_TRUE(row.status==Status::Unsupported);
  AMPS_EXPECT_EQ(row.entries,0);
}
