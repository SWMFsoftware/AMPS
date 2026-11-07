// Two-rank regression for array_3d::gather_double.  The legacy signature has
// no destination large enough for rank_count*size() elements, so multi-rank
// gathering is an explicitly unsupported contract until a result-bearing API
// is designed.  This test proves rejection occurs before MPI_Gather can write
// past the fixed-size root buffer; MPI_COMM_SELF remains the supported identity
// case and must preserve the array values.

#include "array_3d.h"

#include <cstdio>
#include <stdexcept>
#include <string>

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

void exit(long int, const char*, const char* message) {
  throw std::runtime_error(message == NULL ? "AMPS trap exit" : message);
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);

  int rank=0,nrank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nrank);

  array_3d<double> values(1,1,2);
  values(0,0,0)=10.0+rank;
  values(0,0,1)=20.0+rank;

  bool self_identity=true;
  try {
    values.gather_double(0,MPI_COMM_SELF);
  }
  catch (...) {
    self_identity=false;
  }
  self_identity=self_identity && values(0,0,0)==10.0+rank &&
      values(0,0,1)==20.0+rank;

  bool multirank_rejected=false;
  try {
    values.gather_double(0,MPI_COMM_WORLD);
  }
  catch (const std::runtime_error& error) {
    multirank_rejected=std::string(error.what()).find(
        "explicit destination")!=std::string::npos;
  }

  const int local_ok=(nrank==2 && self_identity && multirank_rejected) ? 1 : 0;
  int global_ok=0;
  MPI_Allreduce(&local_ok,&global_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if (rank==0) {
    std::printf("ARRAY_3D_MPI_GATHER_%s\n",global_ok ? "PASS" : "FAIL");
    std::fflush(stdout);
  }

  MPI_Finalize();
  return global_ok ? 0 : 1;
}
