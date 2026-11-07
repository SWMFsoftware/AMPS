// Process-level rejection tests for CMPI_channel's byte-count and frame
// contracts.  This executable is intentionally compiled without MPI_ON: each
// requested violation must terminate locally and nonzero without risking a
// blocked peer in an MPI launcher.  The two-rank executable separately proves
// successful wire behavior.

#include "mpichannel.h"

#include <climits>
#include <string>

int main(int argc,char** argv) {
  if (argc!=2) return 2;
  const std::string probe=argv[1];

  if (probe=="negative-count") {
    CMPI_channel channel(32);
    channel.openSend(0);
    int value=7;
    channel.send(&value,-1);
  }
  else if (probe=="oversized-send") {
    CMPI_channel channel(4);
    channel.openSend(0);
    int values[2]={1,2};
    channel.send(values,2);
  }
  else if (probe=="missing-frame") {
    CMPI_channel channel(32);
    channel.openRecv(0);
    int value=0;
    channel.recv(value,0);
  }
  else if (probe=="count-overflow") {
    CMPI_channel channel(32);
    channel.openSend(0);
    double value=0.0;
    channel.send(&value,INT_MAX);
  }
  else return 2;

  // Every recognized probe is a contract violation and must have terminated.
  return 0;
}
