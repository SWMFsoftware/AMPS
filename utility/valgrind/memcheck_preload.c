/* Preload library for valgrind memcheck runs of AMPS (run_memcheck.sh builds and preloads it).
 *
 * 1. MPI send buffers.  Memcheck tracks definedness within one process only.  Bytes that a rank
 *    sends through shared memory arrive as "defined" in the receiving process, so uninitialised
 *    data that travel between ranks are not reported anywhere.  The PMPI wrappers below call
 *    VALGRIND_CHECK_MEM_IS_DEFINED on the send buffer of the common point-to-point and collective
 *    calls.  Only contiguous data are checked (a datatype whose size equals its true extent:
 *    MPI_BYTE, MPI_DOUBLE, a contiguous type of those, ...); a derived type with holes is skipped.
 * 2. fwrite().  Memcheck reports undefined bytes of a file only when stdio flushes its buffer,
 *    with the stack of the flush and once per flush site.  The fwrite() wrapper checks the bytes
 *    at the call, so every fwrite() of AMPS that writes undefined bytes is reported at its own line.
 *
 * An undefined byte gives "Uninitialised byte(s) found during client check request" with the stack
 * of the AMPS caller and, with --track-origins=yes, the allocation site.  Without valgrind the
 * checks are no-ops.
 *
 * Build: mpicc -shared -fPIC -O1 -g memcheck_preload.c -o libmemcheck_preload.so -ldl
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <mpi.h>
#include <stddef.h>
#include <valgrind/memcheck.h>

static void check(const void *buf, MPI_Aint offset, int count, MPI_Datatype type) {
  int size;
  MPI_Aint lb, extent;

  if (buf == NULL || buf == MPI_IN_PLACE || buf == MPI_BOTTOM || count <= 0) return;
  if (type == MPI_DATATYPE_NULL) return;
  if (PMPI_Type_size(type, &size) != MPI_SUCCESS || size <= 0) return;
  if (PMPI_Type_get_true_extent(type, &lb, &extent) != MPI_SUCCESS) return;
  if ((MPI_Aint)size != extent) return; /* a type with holes: not checked */

  (void)VALGRIND_CHECK_MEM_IS_DEFINED((const char *)buf + offset * extent + lb, (size_t)count * (size_t)size);
}

static int comm_rank(MPI_Comm comm) {
  int r = -1;
  PMPI_Comm_rank(comm, &r);
  return r;
}

static int comm_size(MPI_Comm comm) {
  int s = 0, inter = 0;
  PMPI_Comm_test_inter(comm, &inter);
  if (inter) PMPI_Comm_remote_size(comm, &s);
  else PMPI_Comm_size(comm, &s);
  return s;
}

int MPI_Send(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  check(buf, 0, count, type);
  return PMPI_Send(buf, count, type, dest, tag, comm);
}

int MPI_Ssend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  check(buf, 0, count, type);
  return PMPI_Ssend(buf, count, type, dest, tag, comm);
}

int MPI_Rsend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  check(buf, 0, count, type);
  return PMPI_Rsend(buf, count, type, dest, tag, comm);
}

int MPI_Bsend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm) {
  check(buf, 0, count, type);
  return PMPI_Bsend(buf, count, type, dest, tag, comm);
}

int MPI_Isend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm, MPI_Request *req) {
  check(buf, 0, count, type);
  return PMPI_Isend(buf, count, type, dest, tag, comm, req);
}

int MPI_Issend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm, MPI_Request *req) {
  check(buf, 0, count, type);
  return PMPI_Issend(buf, count, type, dest, tag, comm, req);
}

int MPI_Sendrecv(const void *sbuf, int scount, MPI_Datatype stype, int dest, int stag,
                 void *rbuf, int rcount, MPI_Datatype rtype, int source, int rtag, MPI_Comm comm, MPI_Status *status) {
  check(sbuf, 0, scount, stype);
  return PMPI_Sendrecv(sbuf, scount, stype, dest, stag, rbuf, rcount, rtype, source, rtag, comm, status);
}

int MPI_Bcast(void *buf, int count, MPI_Datatype type, int root, MPI_Comm comm) {
  if (comm_rank(comm) == root) check(buf, 0, count, type);
  return PMPI_Bcast(buf, count, type, root, comm);
}

int MPI_Reduce(const void *sbuf, void *rbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
  check(sbuf == MPI_IN_PLACE ? rbuf : sbuf, 0, count, type);
  return PMPI_Reduce(sbuf, rbuf, count, type, op, root, comm);
}

int MPI_Allreduce(const void *sbuf, void *rbuf, int count, MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
  check(sbuf == MPI_IN_PLACE ? rbuf : sbuf, 0, count, type);
  return PMPI_Allreduce(sbuf, rbuf, count, type, op, comm);
}

int MPI_Scan(const void *sbuf, void *rbuf, int count, MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
  check(sbuf == MPI_IN_PLACE ? rbuf : sbuf, 0, count, type);
  return PMPI_Scan(sbuf, rbuf, count, type, op, comm);
}

int MPI_Exscan(const void *sbuf, void *rbuf, int count, MPI_Datatype type, MPI_Op op, MPI_Comm comm) {
  check(sbuf == MPI_IN_PLACE ? rbuf : sbuf, 0, count, type);
  return PMPI_Exscan(sbuf, rbuf, count, type, op, comm);
}

int MPI_Gather(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype,
               int root, MPI_Comm comm) {
  check(sbuf, 0, scount, stype);
  return PMPI_Gather(sbuf, scount, stype, rbuf, rcount, rtype, root, comm);
}

int MPI_Gatherv(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, const int rcounts[], const int displs[],
                MPI_Datatype rtype, int root, MPI_Comm comm) {
  check(sbuf, 0, scount, stype);
  return PMPI_Gatherv(sbuf, scount, stype, rbuf, rcounts, displs, rtype, root, comm);
}

int MPI_Allgather(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype,
                  MPI_Comm comm) {
  check(sbuf, 0, scount, stype);
  return PMPI_Allgather(sbuf, scount, stype, rbuf, rcount, rtype, comm);
}

int MPI_Allgatherv(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, const int rcounts[],
                   const int displs[], MPI_Datatype rtype, MPI_Comm comm) {
  check(sbuf, 0, scount, stype);
  return PMPI_Allgatherv(sbuf, scount, stype, rbuf, rcounts, displs, rtype, comm);
}

int MPI_Scatter(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype,
                int root, MPI_Comm comm) {
  if (comm_rank(comm) == root) check(sbuf, 0, scount * comm_size(comm), stype);
  return PMPI_Scatter(sbuf, scount, stype, rbuf, rcount, rtype, root, comm);
}

int MPI_Scatterv(const void *sbuf, const int scounts[], const int displs[], MPI_Datatype stype, void *rbuf,
                 int rcount, MPI_Datatype rtype, int root, MPI_Comm comm) {
  if (comm_rank(comm) == root) {
    int i, n = comm_size(comm);
    for (i = 0; i < n; i++) check(sbuf, displs[i], scounts[i], stype);
  }
  return PMPI_Scatterv(sbuf, scounts, displs, stype, rbuf, rcount, rtype, root, comm);
}

int MPI_Alltoall(const void *sbuf, int scount, MPI_Datatype stype, void *rbuf, int rcount, MPI_Datatype rtype,
                 MPI_Comm comm) {
  check(sbuf, 0, scount * comm_size(comm), stype);
  return PMPI_Alltoall(sbuf, scount, stype, rbuf, rcount, rtype, comm);
}

int MPI_Alltoallv(const void *sbuf, const int scounts[], const int sdispls[], MPI_Datatype stype, void *rbuf,
                  const int rcounts[], const int rdispls[], MPI_Datatype rtype, MPI_Comm comm) {
  if (sbuf != MPI_IN_PLACE) {
    int i, n = comm_size(comm);
    for (i = 0; i < n; i++) check(sbuf, sdispls[i], scounts[i], stype);
  }
  return PMPI_Alltoallv(sbuf, scounts, sdispls, stype, rbuf, rcounts, rdispls, rtype, comm);
}

size_t fwrite(const void *ptr, size_t size, size_t n, FILE *stream) {
  static size_t (*real_fwrite)(const void *, size_t, size_t, FILE *) = NULL;

  if (real_fwrite == NULL) real_fwrite = (size_t (*)(const void *, size_t, size_t, FILE *))dlsym(RTLD_NEXT, "fwrite");
  if (ptr != NULL && size != 0 && n != 0) (void)VALGRIND_CHECK_MEM_IS_DEFINED(ptr, size * n);
  return real_fwrite(ptr, size, n, stream);
}
