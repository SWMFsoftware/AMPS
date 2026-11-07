//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================================
//$Id$
//====================================================================

#ifndef MPI_CHANNEL
#define MPI_CHANNEL

#include "global.h"

// This header is instantiated directly by coupling and test translation units;
// do not rely on a caller's incidental include order for the byte-copy,
// allocation, diagnostic, or timestamp facilities used below.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <new>
#include <type_traits>
#include <cstdint>
#include <vector>

using namespace std;

#define MPI_CHANNEL_DEFAULT_BUFFER_SIZE 100000

#define _MPI_CHANNEL_MODE_BCAST_    0
#define _MPI_CHANNEL_MODE_SENDRECV_ 1

class CMPI_channel {
private:
  long int max_MPIbuffer_size;

  char *sendBuffer,**recvBuffer;
  long int sendptr,*recvptr;
  
  long int *RecvDataLength;
  int TotalThreadsNumber,BcastThread;

  int ChannelMode;

  void ValidateThreadIndex(int thread) {
    if (thread<0 || thread>=TotalThreadsNumber)
      exit(__LINE__,"Error: CMPI_channel thread index is out of range");
  }

  template<class T> long int CheckedByteLength(int count) {
    static_assert(std::is_trivially_copyable<T>::value,
                  "CMPI_channel payloads must be trivially copyable");
    if (count<0)
      exit(__LINE__,"Error: CMPI_channel element count cannot be negative");
    const std::size_t elements=static_cast<std::size_t>(count);
    if (elements>static_cast<std::size_t>(
          std::numeric_limits<long int>::max())/sizeof(T))
      exit(__LINE__,"Error: CMPI_channel byte count overflow");
    return static_cast<long int>(elements*sizeof(T));
  }

  void ReceiveFrame(int thread) {
    RecvDataLength[thread]=0;
#ifdef MPI_ON
    if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
      MPI_Status status;
      MPI_Recv(RecvDataLength+thread,1,MPI_LONG,thread,0,
               MPI_GLOBAL_COMMUNICATOR,&status);
    }
    else {
      MPI_Bcast(RecvDataLength+thread,1,MPI_LONG,thread,
                MPI_GLOBAL_COMMUNICATOR);
    }
#endif

    // A frame length is untrusted wire data. Validate it before supplying it
    // as MPI's receive count or using it in pointer arithmetic.
    if (RecvDataLength[thread]<=0 ||
        RecvDataLength[thread]>max_MPIbuffer_size ||
        RecvDataLength[thread]>std::numeric_limits<int>::max())
      exit(__LINE__,"Error: CMPI_channel received an invalid frame length");

#ifdef MPI_ON
    if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
      MPI_Status status;
      MPI_Recv(recvBuffer[thread],static_cast<int>(RecvDataLength[thread]),
               MPI_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
      int received=0;
      MPI_Get_count(&status,MPI_CHAR,&received);
      if (received!=RecvDataLength[thread])
        exit(__LINE__,"Error: CMPI_channel frame header/payload length mismatch");
    }
    else {
      MPI_Bcast(recvBuffer[thread],static_cast<int>(RecvDataLength[thread]),
                MPI_CHAR,thread,MPI_GLOBAL_COMMUNICATOR);
    }
#endif
    recvptr[thread]=0;
  }

  char* ReceiveBytes(long int length,int thread) {
    ValidateThreadIndex(thread);
    if (recvBuffer==NULL || recvBuffer[thread]==NULL)
      exit(__LINE__,"Error: CMPI_channel receive buffer is not open");
    if (length<=0 || length>max_MPIbuffer_size)
      exit(__LINE__,"Error: CMPI_channel requested receive length is invalid");

    if (recvptr[thread]>=RecvDataLength[thread]) ReceiveFrame(thread);
    if (recvptr[thread]<0 || recvptr[thread]>RecvDataLength[thread] ||
        length>RecvDataLength[thread]-recvptr[thread])
      exit(__LINE__,"Error: CMPI_channel receive request exceeds the current frame");

    char* result=recvBuffer[thread]+recvptr[thread];
    recvptr[thread]+=length;
    return result;
  }

public:
  int ThisThread,sendThread;


  void exit(long int nline,const char* msg=NULL) {
    char str[2000];
    FILE* errorlog=fopen("$ERRORLOG","a+");

    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);

    const char* diagnostic=(msg==NULL) ? "CMPI_channel fatal error" : msg;
    snprintf(str,sizeof(str)," exit: line=%ld, file=%s, message=%s\n",
             nline,__FILE__,diagnostic);

    if (errorlog!=NULL) {
      fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
      fprintf(errorlog,"%s\n\n",str);
    }

    printf("$PREFIX:Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
    printf("$PREFIX:%s\n\n",str);

    if (errorlog!=NULL) fclose(errorlog);
    // A channel contract violation is fatal.  Returning success made MPI
    // launchers and the nightly runner accept truncated transfers as PASS.
    ::exit(EXIT_FAILURE);
  }

public:



  void RedirectSendBuffer(int newSendThread) {
    if (sendThread==-1) exit(__LINE__,"Error: need to init first");
    ValidateThreadIndex(newSendThread);
    flush();
    sendThread=newSendThread;
  }

  void RedirectRecvBuffer(int newRecvThread,int oldRecvThread=-1) {

    ValidateThreadIndex(newRecvThread);
    if (oldRecvThread!=-1) ValidateThreadIndex(oldRecvThread);

    if (recvBuffer[newRecvThread]!=NULL) return;

//    if (recvBuffer[newRecvThread]!=NULL) exit(__LINE__,"Error: the recv direction is already initialized");

    if (oldRecvThread==-1) { //search for an allocated buffer
      for (int i=0;i<TotalThreadsNumber;i++) if (recvBuffer[i]!=NULL) {
        oldRecvThread=i;
        break;
      }

      if (oldRecvThread==-1) exit(__LINE__,"Error: need to allocate the recv buffer first");
    }

    if ((recvBuffer[oldRecvThread]==NULL)||(recvptr[oldRecvThread]!=RecvDataLength[oldRecvThread])) exit(__LINE__,"Error: cannot redirect the recv buffer");

    recvBuffer[newRecvThread]=recvBuffer[oldRecvThread],recvptr[newRecvThread]=0,RecvDataLength[newRecvThread]=0;
    recvBuffer[oldRecvThread]=NULL,recvptr[oldRecvThread]=0,RecvDataLength[oldRecvThread]=0;
  }

  void init(long int buffersize) {
    max_MPIbuffer_size=buffersize;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    #ifdef MPI_ON
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&TotalThreadsNumber);
    #else
    ThisThread=0,TotalThreadsNumber=1;
    #endif

    if (buffersize<=0 || buffersize>std::numeric_limits<int>::max())
      exit(__LINE__,"Error: CMPI_channel buffer size must fit a positive MPI count");

    try {
      recvBuffer=new char* [TotalThreadsNumber];
      recvptr=new long int[TotalThreadsNumber];
      RecvDataLength=new long int[TotalThreadsNumber];
    }
    catch (const bad_alloc&) {exit(__LINE__,"Error: CMPI_channel::CMPI_channel(long) cannot allocate buffer");}
 
    for (int thread=0;thread<TotalThreadsNumber;thread++) recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

  CMPI_channel(long int buffersize) {
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_,BcastThread=-1;

    init(buffersize);
  }

  CMPI_channel() {
    max_MPIbuffer_size=0;
    sendBuffer=NULL,sendptr=0,recvptr=NULL;
    sendThread=-1;

    ThisThread=0,TotalThreadsNumber=0;

    recvBuffer=NULL;
    recvptr=NULL;
    RecvDataLength=NULL;

    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_,BcastThread=-1;
  };

  void remove() {
    int thread;

    if (sendBuffer!=NULL) {
      delete [] sendBuffer;
      sendBuffer=NULL;
    }

    if (recvBuffer!=NULL) {
      for (thread=0;thread<TotalThreadsNumber;thread++) if (recvBuffer[thread]!=NULL) delete [] recvBuffer[thread];
      delete [] recvBuffer;
      delete [] recvptr;
      delete [] RecvDataLength;

      recvBuffer=NULL,recvptr=NULL,RecvDataLength=NULL;
    }
  }; 

  ~CMPI_channel() {
     flush();
     remove();
  } 


  void openSend(int send) {
    if (max_MPIbuffer_size==0) init(MPI_CHANNEL_DEFAULT_BUFFER_SIZE);
    ValidateThreadIndex(send);
    if (sendBuffer!=NULL) exit(__LINE__,"CMPI_channel::openSend reinitialization of the send buffer"); 

    try {
      sendBuffer=new char[max_MPIbuffer_size];
    }
    catch (const bad_alloc&) {
      char msg[200]; 

      sprintf(msg,"Error: CMPI_channel::openSend cannot allocate send buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    sendptr=0,sendThread=send; 
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_;
  };

  void flush() {
    if (sendptr!=0) {

#ifdef MPI_ON
      if (ChannelMode==_MPI_CHANNEL_MODE_SENDRECV_) {
        MPI_Send(&sendptr,1,MPI_LONG,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Send(sendBuffer,sendptr,MPI_CHAR,sendThread,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else {
        MPI_Bcast(&sendptr,1,MPI_LONG,sendThread,MPI_GLOBAL_COMMUNICATOR);
        MPI_Bcast(sendBuffer,sendptr,MPI_CHAR,sendThread,MPI_GLOBAL_COMMUNICATOR);
      }
#endif

      sendptr=0;
    }
  };


  void closeSend() {
    if (sendBuffer!=NULL) {
      flush();

      delete [] sendBuffer;

      sendBuffer=NULL;
      sendThread=-1,sendptr=0;
    }
  };

  template<class T> inline void send(T* data,int nsend) {
    const long int length=CheckedByteLength<T>(nsend);

    if (length==0) return;

    if (data==NULL) exit(__LINE__,"Error: CMPI_channel cannot send from a null buffer");
    if (sendBuffer==NULL || sendThread<0)
      exit(__LINE__,"Error: CMPI_channel send buffer is not open");

    if (length>max_MPIbuffer_size) {
      exit(__LINE__,"Error: the length of the data vector exeeds that of the internal data buffer. Need to increase the size of the data buffer used when initialize the channel obsect.");
    }

    // Equality is safe: the current payload exactly fills the allocated byte
    // buffer.  Flushing on equality when sendptr==0 used to emit an empty MPI
    // frame before the real frame and desynchronized the receiver.
    if (sendptr<0 || sendptr>max_MPIbuffer_size ||
        length>max_MPIbuffer_size-sendptr) flush();
     
    memcpy(sendBuffer+sendptr,data,length);
    sendptr+=length;
  }


  template<class T> inline void send(T data) {
    send(&data,1);
  }


  void openRecv(int thread) { 
    if (max_MPIbuffer_size==0) init(MPI_CHANNEL_DEFAULT_BUFFER_SIZE);
    ValidateThreadIndex(thread);
    if (recvBuffer[thread]!=NULL) exit(__LINE__,"CMPI_channel::openRecv reinitialization of the recv buffer");   

    try {
      recvBuffer[thread]=new char[max_MPIbuffer_size];
    }
    catch (const bad_alloc&) {
      char msg[200];

      sprintf(msg,"Error: CMPI_channel::openRecv cannot allocate recv buffer, need %ld bytes\n", max_MPIbuffer_size);
      exit(__LINE__,msg);
    }

    recvptr[thread]=0;
    ChannelMode=_MPI_CHANNEL_MODE_SENDRECV_;
  };

  void openRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) openRecv(thread);
  };
 

  void closeRecv(int thread) {
    ValidateThreadIndex(thread);
    if (recvBuffer[thread]==NULL)
      exit(__LINE__,"Error: CMPI_channel receive buffer is not open");
    if (recvptr[thread]!=RecvDataLength[thread]) {
      char msg[200];

      sprintf(msg,"Error: closeRecv: (recvptr[thread]!=RecvDataLength[thread]) for thread=%i\n",thread);
      exit(__LINE__,msg);
    }

    delete [] recvBuffer[thread];
    recvBuffer[thread]=NULL,recvptr[thread]=0,RecvDataLength[thread]=0; 
  };

  void closeRecvAll() {
    int thread;

    for (thread=0;thread<TotalThreadsNumber;thread++) if (thread!=ThisThread) closeRecv(thread);
  };

  template<class T> inline T* recvPointer(int nrecv,int thread) {
    const long int length=CheckedByteLength<T>(nrecv);

    if (length==0) return NULL;
    static_assert(std::is_default_constructible<T>::value,
                  "CMPI_channel recvPointer payloads must be default constructible");

    // The packed wire stream has byte alignment: for example, an int followed
    // by a double puts the double at offset four.  Returning a T* directly into
    // recvBuffer would therefore expose an unaligned pointer and invoke
    // undefined behavior on dereference.  Copy into per-thread typed storage
    // so the public pointer API remains usable by its legacy double-vector
    // callers.  As with the historic frame view, callers must consume the
    // result before a subsequent recvPointer<T>() call on the same thread.
    static thread_local std::vector<T> alignedStorage;
    alignedStorage.resize(static_cast<std::size_t>(nrecv));
    std::memcpy(alignedStorage.data(),ReceiveBytes(length,thread),
                static_cast<std::size_t>(length));
    return alignedStorage.data();
  }

  template<class T> inline void recv(T& data,int thread) {
    const long int length=CheckedByteLength<T>(1);
    std::memcpy(&data,ReceiveBytes(length,thread),static_cast<std::size_t>(length));
  }

  template<class T> inline void recv(T* data,int nrecv,int thread) {
    const long int length=CheckedByteLength<T>(nrecv);
    if (length==0) return;
    if (data==NULL) exit(__LINE__,"Error: CMPI_channel cannot receive into a null buffer");
    std::memcpy(data,ReceiveBytes(length,thread),static_cast<std::size_t>(length));
  } 

  template<class T> inline T recv(int thread) {
    T value;
    recv(value,thread);
    return value;
  }

  //open/close procedures in the Bcast mode
  void openBcast(int thread) {
    if (ThisThread==thread) openSend(thread);
    else openRecv(thread);

    ChannelMode=_MPI_CHANNEL_MODE_BCAST_,BcastThread=thread;
  }

  void closeBcast() {
    if (ChannelMode!=_MPI_CHANNEL_MODE_BCAST_) exit(__LINE__,"wrong option");

    if (ThisThread==BcastThread) closeSend();
    else closeRecv(BcastThread);

    BcastThread=-1;
  }


};

#endif 
