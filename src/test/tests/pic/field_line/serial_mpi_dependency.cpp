#include <mpi.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Test-only, single-process dependency adapter for the real native container.
// This is NOT an MPI implementation/qualification result. It permits testing
// allocator, particle lists and sampling arithmetic in an environment where
// MPI startup sockets are unavailable. Every supported collective is exactly
// its one-rank identity operation; an unexpected type/root/count aborts rather
// than manufacturing a successful communication result.
namespace {
std::size_t Width(MPI_Datatype type) {
  if(type==MPI_DOUBLE)return sizeof(double);
  if(type==MPI_INT)return sizeof(int);
  if(type==MPI_LONG)return sizeof(long);
  if(type==MPI_UNSIGNED_LONG_LONG)return sizeof(unsigned long long);
  if(type==MPI_UNSIGNED_LONG)return sizeof(unsigned long);
  if(type==MPI_CHAR||type==MPI_BYTE||type==MPI_UNSIGNED_CHAR)return 1;
  std::fprintf(stderr,"unsupported datatype in explicit serial MPI dependency adapter\n");
  std::abort();
}
void Copy(const void* from,void* to,int count,MPI_Datatype type) {
  if(count<0)std::abort();
  if(count&&from!=MPI_IN_PLACE)std::memcpy(to,from,static_cast<std::size_t>(count)*Width(type));
}
void Root(int root){if(root!=0)std::abort();}
}
extern "C" {
int MPI_Init(int*,char***) {
  std::puts("SCOPE: serial MPI dependency adapter; no native MPI communication qualification");
  return MPI_SUCCESS;
}
int MPI_Finalize(){return MPI_SUCCESS;}
int MPI_Initialized(int* flag){*flag=1;return MPI_SUCCESS;}
int MPI_Finalized(int* flag){*flag=0;return MPI_SUCCESS;}
int MPI_Comm_rank(MPI_Comm,int* rank){*rank=0;return MPI_SUCCESS;}
int MPI_Comm_size(MPI_Comm,int* size){*size=1;return MPI_SUCCESS;}
int MPI_Comm_test_inter(MPI_Comm,int* flag){*flag=0;return MPI_SUCCESS;}
int MPI_Comm_dup(MPI_Comm comm,MPI_Comm* result){*result=comm;return MPI_SUCCESS;}
int MPI_Comm_free(MPI_Comm* comm){*comm=MPI_COMM_NULL;return MPI_SUCCESS;}
int MPI_Barrier(MPI_Comm){return MPI_SUCCESS;}
int MPI_Bcast(void*,int count,MPI_Datatype type,int root,MPI_Comm) {
  Root(root);if(count<0)std::abort();if(count)Width(type);return MPI_SUCCESS;
}
int MPI_Allreduce(const void* from,void* to,int count,MPI_Datatype type,MPI_Op,MPI_Comm) {
  Copy(from,to,count,type);return MPI_SUCCESS;
}
int MPI_Reduce(const void* from,void* to,int count,MPI_Datatype type,MPI_Op,int root,MPI_Comm) {
  Root(root);Copy(from,to,count,type);return MPI_SUCCESS;
}
int MPI_Allgather(const void* from,int count,MPI_Datatype type,void* to,int targetCount,MPI_Datatype targetType,MPI_Comm) {
  if(count<0||targetCount<0||static_cast<std::size_t>(count)*Width(type)!=static_cast<std::size_t>(targetCount)*Width(targetType))std::abort();
  Copy(from,to,count,type);return MPI_SUCCESS;
}
double MPI_Wtime() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
int MPI_Abort(MPI_Comm,int code){std::exit(code?code:1);}
}
