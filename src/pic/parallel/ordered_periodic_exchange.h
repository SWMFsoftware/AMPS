#ifndef AMPS_ORDERED_PERIODIC_EXCHANGE_H
#define AMPS_ORDERED_PERIODIC_EXCHANGE_H

// Transport-only periodic copy planner. The caller supplies a replicated,
// ordered pair list and exact byte serializers; this header knows no PIC
// records, mesh types, physical fields or dimensional units. RecordBytes is
// an opaque byte count, never a count of doubles or particles.
//
// snapshot_safe=true is a physical data-dependency contract: packing any source
// must be independent of every destination write, and callbacks must implement
// pure copies. Receives/posts may then overlap, but destination writes still
// occur in the original pair order. A caller that cannot prove independence
// uses the pairwise reference path. No collective is introduced: ranks with
// no pair participation can return without waiting on unrelated ranks.
#include <mpi.h>
#include <algorithm>
#include <climits>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace amps_periodic_exchange {
struct Pair { int source, destination; };
struct Statistics {
  std::size_t send_bytes=0, receive_bytes=0, send_messages=0, receive_messages=0, local_pairs=0;
};

class Workspace {
  // Only byte capacity is retained across calls. All pair offsets/owners are
  // rebuilt from the current list, so repartition, table replacement/reordering,
  // record-layout changes and communicator-size changes cannot reuse a stale
  // communication plan or pointer. No field value survives as current state.
  std::vector<std::vector<char>> send_, receive_;
  std::vector<char> local_;
  std::vector<std::size_t> send_offset_, receive_offset_;
  std::vector<MPI_Request> requests_;

  static void Check(int code) {
    if (code!=MPI_SUCCESS) throw std::runtime_error("periodic exchange MPI failure");
  }
  static void Add(std::size_t& size,std::size_t bytes) {
    if (bytes>std::numeric_limits<std::size_t>::max()-size)
      throw std::overflow_error("periodic exchange byte count overflow");
    size+=bytes;
  }
public:
  template<class Pack,class Unpack>
  Statistics Exchange(const std::vector<Pair>& pairs,std::size_t bytes,
      int rank,int ranks,MPI_Comm comm,Pack pack,Unpack unpack,
      bool snapshot_safe,int tag=0,std::size_t max_chunk=INT_MAX) {
    if (ranks<1 || rank<0 || rank>=ranks || max_chunk==0 || max_chunk>INT_MAX)
      throw std::invalid_argument("invalid periodic exchange rank/chunk");
    for (const auto& pair:pairs)
      if (pair.source<0 || pair.source>=ranks || pair.destination<0 || pair.destination>=ranks)
        throw std::invalid_argument("invalid periodic pair owner");
    Statistics stats;
    if (pairs.empty() || bytes==0) return stats;
    local_.resize(bytes);

    if (!snapshot_safe) {
      // This is the original streaming semantics, including pack/unpack order
      // for aliasing sources or stateful user callbacks. Fragmentation changes
      // only message sizing; MPI's same-peer/tag ordering preserves byte order.
      for (std::size_t i=0;i<pairs.size();i++) {
        const auto pair=pairs[i];
        if (pair.source==rank) pack(i,local_.data());
        if (pair.source==pair.destination) {
          if (pair.source==rank) { unpack(i,local_.data()); stats.local_pairs++; }
        }
        else for (std::size_t offset=0;offset<bytes;offset+=std::min(max_chunk,bytes-offset)) {
          const int count=static_cast<int>(std::min(max_chunk,bytes-offset));
          if (pair.source==rank) {
            Check(MPI_Send(local_.data()+offset,count,MPI_BYTE,pair.destination,tag,comm));
            stats.send_bytes+=count; stats.send_messages++;
          }
          if (pair.destination==rank) {
            Check(MPI_Recv(local_.data()+offset,count,MPI_BYTE,pair.source,tag,comm,MPI_STATUS_IGNORE));
            stats.receive_bytes+=count; stats.receive_messages++;
          }
        }
        if (pair.destination==rank && pair.source!=rank) unpack(i,local_.data());
      }
      return stats;
    }

    send_.resize(ranks); receive_.resize(ranks);
    std::vector<std::size_t> send_size(ranks,0), receive_size(ranks,0);
    send_offset_.resize(pairs.size()); receive_offset_.resize(pairs.size());
    for (std::size_t i=0;i<pairs.size();i++) {
      const auto pair=pairs[i];
      if (pair.source==pair.destination) continue;
      if (pair.source==rank) { send_offset_[i]=send_size[pair.destination]; Add(send_size[pair.destination],bytes); }
      if (pair.destination==rank) { receive_offset_[i]=receive_size[pair.source]; Add(receive_size[pair.source],bytes); }
    }
    requests_.clear();
    // All receives are posted before packing/sending. Each peer's records use
    // global pair order on both ends. No tag per block or receive-any ordering
    // is needed, and each byte count remains representable by the MPI C API.
    for (int peer=0;peer<ranks;peer++) {
      send_[peer].resize(send_size[peer]); receive_[peer].resize(receive_size[peer]);
      for (std::size_t offset=0;offset<receive_size[peer];offset+=std::min(max_chunk,receive_size[peer]-offset)) {
        MPI_Request request;
        const int count=static_cast<int>(std::min(max_chunk,receive_size[peer]-offset));
        Check(MPI_Irecv(receive_[peer].data()+offset,count,MPI_BYTE,peer,tag,comm,&request));
        requests_.push_back(request); stats.receive_messages++; stats.receive_bytes+=count;
      }
    }
    for (std::size_t i=0;i<pairs.size();i++) {
      const auto pair=pairs[i];
      if (pair.source==rank && pair.destination!=rank)
        pack(i,send_[pair.destination].data()+send_offset_[i]);
    }
    for (int peer=0;peer<ranks;peer++)
      for (std::size_t offset=0;offset<send_size[peer];offset+=std::min(max_chunk,send_size[peer]-offset)) {
        MPI_Request request;
        const int count=static_cast<int>(std::min(max_chunk,send_size[peer]-offset));
        Check(MPI_Isend(send_[peer].data()+offset,count,MPI_BYTE,peer,tag,comm,&request));
        requests_.push_back(request); stats.send_messages++; stats.send_bytes+=count;
      }
    if (requests_.size()>INT_MAX) throw std::overflow_error("too many periodic MPI requests");
    if (!requests_.empty()) Check(MPI_Waitall(static_cast<int>(requests_.size()),requests_.data(),MPI_STATUSES_IGNORE));
    for (std::size_t i=0;i<pairs.size();i++) {
      const auto pair=pairs[i];
      if (pair.destination!=rank) continue;
      if (pair.source==rank) {
        pack(i,local_.data()); unpack(i,local_.data()); stats.local_pairs++;
      }
      else unpack(i,receive_[pair.source].data()+receive_offset_[i]);
    }
    return stats;
  }
};
} // namespace amps_periodic_exchange
#endif
