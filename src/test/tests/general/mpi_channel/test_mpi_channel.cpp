// Two-rank behavioral tests for src/general/mpichannel.h.
//
// CMPI_channel serializes trivially copyable values into a byte buffer and
// transfers that buffer through MPI. These tests deliberately use independent
// hand-written payloads: they do not use a second AMPS serializer as an oracle.
// Rank zero sends and rank one validates, then MPI_Allreduce makes the process
// exit status identical on both ranks so the global runner cannot miss a
// receiver-only failure.

#include "mpichannel.h"
#include "rnd.h"
#include "specfunc.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

int ThisThread = 0;
int TotalThreadsNumber = 1;
int deviceThisThread = 0;
int deviceTotalThreadsNumber = 1;
int cudaThreadLimitMallocHeapSize = 0;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;

// rnd.cpp uses the AMPS fatal symbol for an invalid compile-time selector.
// The tested selector is valid; this definition keeps the isolated executable
// link-complete and still fails collectively if that supposedly unreachable
// branch is ever entered.
void exit(long int, const char*, const char*) {
  MPI_Abort(MPI_GLOBAL_COMMUNICATOR, 97);
  std::abort();
}

namespace {

struct Record {
  int identifier;
  double position[3];
  unsigned char active;
};

bool scalar_case() {
  CMPI_channel channel(256);
  if (ThisThread == 0) {
    channel.openSend(1);
    channel.send<int>(-173);
    channel.send<double>(6.25);
    channel.send<char>('x');
    channel.send<double>(-2.5);
    channel.closeSend();
    return true;
  }

  channel.openRecv(0);
  const int integer_value = channel.recv<int>(0);
  const double real_value = channel.recv<double>(0);
  const char marker = channel.recv<char>(0);
  const double pointer_value = channel.recvPointer<double>(1, 0)[0];
  channel.closeRecv(0);
  return integer_value == -173 && std::fabs(real_value - 6.25) < 1.0e-14 &&
         marker == 'x' && std::fabs(pointer_value + 2.5) < 1.0e-14;
}

bool vector_case() {
  double expected[5] = {-3.0, -0.5, 0.0, 2.5, 9.0};
  // Exact-capacity payload is intentional.  It guards the boundary where an
  // older implementation emitted a zero-byte frame before the real vector.
  CMPI_channel channel(sizeof(expected));
  if (ThisThread == 0) {
    channel.openSend(1);
    channel.send(expected, 5);
    channel.closeSend();
    return true;
  }

  double observed[5] = {};
  channel.openRecv(0);
  channel.recv(observed, 5, 0);
  channel.closeRecv(0);
  return std::memcmp(expected, observed, sizeof(expected)) == 0;
}

bool structured_case() {
  Record expected = {42, {-1.25, 7.5, 0.125}, 1};
  CMPI_channel channel(256);
  if (ThisThread == 0) {
    channel.openSend(1);
    channel.send(expected);
    channel.closeSend();
    return true;
  }

  channel.openRecv(0);
  const Record observed = channel.recv<Record>(0);
  channel.closeRecv(0);
  return observed.identifier == expected.identifier &&
         observed.active == expected.active &&
         std::memcmp(observed.position, expected.position,
                     sizeof(expected.position)) == 0;
}

bool zero_length_case() {
  // A zero-length operation is a stream no-op. Placing it between two real
  // values proves that it neither emits an MPI frame nor advances the byte
  // cursor and therefore cannot desynchronize the following payload.
  CMPI_channel channel(256);
  int ignored = 99;
  if (ThisThread == 0) {
    channel.openSend(1);
    channel.send<int>(11);
    channel.send(&ignored, 0);
    channel.send<int>(29);
    channel.closeSend();
    return true;
  }

  int untouched = -1;
  channel.openRecv(0);
  const int first = channel.recv<int>(0);
  channel.recv(&untouched, 0, 0);
  const int second = channel.recv<int>(0);
  channel.closeRecv(0);
  return first == 11 && second == 29 && untouched == -1;
}

bool rank_rng_case() {
  // Production rnd_seed(-1) maps MPI rank r to recurrence seed r+1.  Compute
  // the expected first value through the explicit-stream overload so this
  // check compares two documented entry points rather than hard-coding a
  // machine-dependent launcher order.
  rnd_seed(-1);
  cRndSeedContainer oracle;
  oracle.Seed = static_cast<unsigned long>(ThisThread + 1);
  const bool uniform_matches = rnd() == rnd(&oracle);

  // One normal draw leaves the second Box--Muller variate cached.  Reseeding
  // must clear that derived state as well as the uniform recurrence; otherwise
  // the first normal after a restart depends on how many normals preceded it.
  rnd_seed(-1);
  const double first_normal = Vector3D::Distribution::Normal();
  rnd_seed(-1);
  const double replayed_normal = Vector3D::Distribution::Normal();
  return uniform_matches && first_normal == replayed_normal;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &TotalThreadsNumber);

  // A raw sender deliberately violates the channel's two-message frame
  // contract: the advertised byte count is larger than the actual payload.
  // Rank one must reject the discrepancy after MPI_Get_count rather than
  // treating uninitialized buffer bytes as serialized data.  This mode is
  // launched as its own expected-failure MPI job by the Python adapter.
  if (argc == 3 && std::string(argv[1]) == "--fatal-probe" &&
      std::string(argv[2]) == "frame-payload-mismatch") {
    if (TotalThreadsNumber != 2) {
      MPI_Finalize();
      return 2;
    }
    if (ThisThread == 0) {
      long int advertised = 4 * static_cast<long int>(sizeof(int));
      int payload = 17;
      MPI_Send(&advertised, 1, MPI_LONG, 1, 0, MPI_GLOBAL_COMMUNICATOR);
      MPI_Send(&payload, 1, MPI_INT, 1, 0, MPI_GLOBAL_COMMUNICATOR);
      MPI_Finalize();
      return 0;
    }
    CMPI_channel channel(64);
    channel.openRecv(0);
    int values[4] = {};
    channel.recv(values, 4, 0);
    MPI_Finalize();
    return 0;  // The invalid frame must terminate before this point.
  }

  std::string identifier;
  if (argc == 3 && std::string(argv[1]) == "--case") identifier = argv[2];
  bool local_pass = TotalThreadsNumber == 2;
  if (local_pass) {
    if (identifier == "GENERAL-MPI-CHANNEL-SCALAR-001") local_pass = scalar_case();
    else if (identifier == "GENERAL-MPI-CHANNEL-VECTOR-001") local_pass = vector_case();
    else if (identifier == "GENERAL-MPI-CHANNEL-STRUCT-001") local_pass = structured_case();
    else if (identifier == "GENERAL-MPI-CHANNEL-ZERO-001") local_pass = zero_length_case();
    else if (identifier == "GENERAL-RNG-RANK-STREAM-001") local_pass = rank_rng_case();
    else local_pass = false;
  }

  int local = local_pass ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (ThisThread == 0) {
    std::printf("[%s] %s %s\n", identifier.c_str(),
                global == 1 ? "PASS" : "FAIL",
                TotalThreadsNumber == 2 ? "two-rank channel contract" :
                                           "test requires exactly two ranks");
    std::fflush(stdout);
  }
  MPI_Finalize();
  return global == 1 ? 0 : 1;
}
