// Independent payload/order tests for the production transport helper. The
// oracle is an explicit integer byte law, not a second invocation of the
// serializer. Every receiving rank checks its complete payload; MPI_Allreduce
// publishes a rank-local failure instead of accepting a rank-zero-only PASS.
#define OMPI_SKIP_MPICXX 1
#include "parallel/ordered_periodic_exchange.h"
#include <array>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using amps_periodic_exchange::Pair;

unsigned char Value(int source,int pair,int epoch,int byte) {
  return static_cast<unsigned char>((29*source+11*pair+7*epoch+byte)%251);
}

bool Payload(amps_periodic_exchange::Workspace& workspace,int rank,int ranks,
    MPI_Comm comm,int epoch,std::size_t bytes,bool local,bool reverse,std::size_t chunk) {
  std::vector<Pair> pairs;
  for (int i=0;i<3*ranks;i++) pairs.push_back({i%ranks,local ? i%ranks : (i+1)%ranks});
  if (reverse) std::reverse(pairs.begin(),pairs.end());
  std::vector<std::vector<char>> received(pairs.size(),std::vector<char>(bytes,0));
  const auto stats=workspace.Exchange(pairs,bytes,rank,ranks,comm,
      [&](std::size_t i,char* data) {
        for (std::size_t j=0;j<bytes;j++) data[j]=Value(rank,int(i),epoch,int(j));
      },
      [&](std::size_t i,char* data) {std::memcpy(received[i].data(),data,bytes);},true,71,chunk);
  bool ok=true;
  std::size_t expected_send=0,expected_receive=0;
  for (std::size_t i=0;i<pairs.size();i++) {
    const auto pair=pairs[i];
    if (pair.destination==rank) for (std::size_t j=0;j<bytes;j++)
      ok=ok && static_cast<unsigned char>(received[i][j])==Value(pair.source,int(i),epoch,int(j));
    if (pair.source==rank && pair.destination!=rank) expected_send+=bytes;
    if (pair.destination==rank && pair.source!=rank) expected_receive+=bytes;
  }
  return ok && stats.send_bytes==expected_send && stats.receive_bytes==expected_receive;
}

bool Run(const std::string& id,int rank,int ranks) {
  amps_periodic_exchange::Workspace workspace;
  if (id=="PERF-PERIODIC-LOCAL-001") return Payload(workspace,rank,ranks,MPI_COMM_WORLD,1,37,true,false,INT_MAX);
  if (id=="PERF-PERIODIC-REMOTE-001") return Payload(workspace,rank,ranks,MPI_COMM_WORLD,2,83,false,false,INT_MAX);
  if (id=="PERF-PERIODIC-LAYOUT-001") {
    // Reuse the same workspace after record size, source contents and pair
    // order change. No mesh-generation assumption can rescue stale offsets.
    bool ok=true;
    for (int epoch=0;epoch<4;epoch++)
      ok=Payload(workspace,rank,ranks,MPI_COMM_WORLD,epoch,epoch%2 ? 131:7,false,epoch%2,INT_MAX) && ok;
    return ok;
  }
  if (id=="PERF-PERIODIC-CHUNK-001") return Payload(workspace,rank,ranks,MPI_COMM_WORLD,4,79,false,true,17);
  if (id=="PERF-PERIODIC-ALIAS-001") {
    // A remote write becomes the next source. Pack-ahead would transport the
    // stale 99; the ordered reference must propagate 17 through both links.
    const int partner=ranks>1 ? 1:0;
    std::array<unsigned char,3> storage{{17,99,0}};
    std::vector<Pair> pairs{{0,partner},{partner,0}};
    workspace.Exchange(pairs,1,rank,ranks,MPI_COMM_WORLD,
        [&](std::size_t i,char* p) {*p=storage[i];},
        [&](std::size_t i,char* p) {storage[i+1]=static_cast<unsigned char>(*p);},false,72);
    return (rank!=0 || storage[2]==17) && (rank!=partner || storage[1]==17);
  }
  if (id=="PERF-PERIODIC-EMPTY-001") {
    int calls=0;
    auto pack=[&](std::size_t,char*) {calls++;};
    auto stats=workspace.Exchange({},15,rank,ranks,MPI_COMM_WORLD,pack,pack,true);
    auto zero=workspace.Exchange({{rank,rank}},0,rank,ranks,MPI_COMM_WORLD,pack,pack,true);
    bool rejected=false;
    try {workspace.Exchange({{-1,0}},1,rank,ranks,MPI_COMM_WORLD,pack,pack,true);}
    catch (const std::invalid_argument&) {rejected=true;}
    return calls==0 && stats.send_bytes==0 && zero.receive_bytes==0 && rejected;
  }
  if (id=="PERF-PERIODIC-COMMUNICATOR-001") {
    // Reversed rank keys and split groups differ from world ownership; a
    // hidden MPI_COMM_WORLD use in the helper would fail or deadlock here.
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD,rank%2,-rank,&comm);
    int local_rank,local_size;
    MPI_Comm_rank(comm,&local_rank); MPI_Comm_size(comm,&local_size);
    const bool ok=Payload(workspace,local_rank,local_size,comm,5,91,false,false,23);
    MPI_Comm_free(&comm);
    return ok;
  }
  return false;
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  int rank,ranks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  const char* ids[]={"PERF-PERIODIC-LOCAL-001","PERF-PERIODIC-REMOTE-001",
      "PERF-PERIODIC-LAYOUT-001","PERF-PERIODIC-CHUNK-001",
      "PERF-PERIODIC-ALIAS-001","PERF-PERIODIC-EMPTY-001",
      "PERF-PERIODIC-COMMUNICATOR-001"};
  int failures=0;
  for (const char* id:ids) {
    int failed=0,global=0;
    try {failed=Run(id,rank,ranks) ? 0:1;} catch (...) {failed=1;}
    MPI_Allreduce(&failed,&global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (rank==0) std::printf("[%s] %s ranks=%d all receiving ranks checked\n",id,global ? "FAIL":"PASS",ranks);
    failures+=global;
  }
  MPI_Finalize();
  return failures ? 1:0;
}
