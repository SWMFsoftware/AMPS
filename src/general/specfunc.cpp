//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf


#ifdef MPI_ON
#include "mpi.h"
#endif

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <pthread.h>
#include <unistd.h>


#include <iostream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "rnd.h"
#include "specfunc.h"


int ThisThread;
int TotalThreadsNumber;
unsigned int ExitErrorCode=0;

long int nint(double a)
{
   long int n;
   n=(long int) a;
   if (a-n>0.5) n++;
   return n;
}


//===================================================
/*
double erf(double s) { 
  double b,d,c,t;

  b=fabs(s);
  if (b>4.0)
    d=1.0;
  else { 
    c=exp(-b*b);
    t=1.0/(1.0+0.3275911*b);
    d=1.0-(0.254829592*t-0.284496736*t*t+1.421413741*t*t*t-
      1.453152027*t*t*t*t+1.061405429*t*t*t*t*t)*c;
  }

  if (s<0.0) d=-d;
  return d;
}  


//===================================================
double gam(double x) {
  double a,y;

  a=1.0;
  y=x;

  if (y<1.0) 
    a/=y;
  else {
    y--;
    while (y>1.0) {
      a*=y;
      y--;
    }
  } 

  return a*(1.0-0.5748646*y+0.9512363*y*y-0.6998588*y*y*y+
    0.4245549*y*y*y*y-0.1010678*y*y*y*y*y); 
}
*/
//===================================================
void PrintErrorLog(const char* message) {
  FILE* errorlog=fopen("$ERRORLOG","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"%s\n\n",message);

  fclose(errorlog);
}

//use: PrintErrorLog(__LINE__,__FILE__, "mesage")
void PrintErrorLog(long int nline, const char* fname, const char* message) {
  FILE* errorlog=fopen("$ERRORLOG","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"file=%s, line=%ld\n",fname,nline);
  fprintf(errorlog,"%s\n\n",message);

#if _STDOUT_ERRORLOG_MODE_ == _STDOUT_ERRORLOG_MODE__ON_
  printf("$PREFIX:Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  printf("$PREFIX:file=%s, line=%ld\n",fname,nline);
  printf("$PREFIX:%s\n\n",message);
#endif

  fclose(errorlog);
}

//===================================================
//use: exit(__LINE__,__FILE__, "mesage")
_TARGET_HOST_ _TARGET_DEVICE_
void exit(long int nline, const char* fname, const char* msg) {
  int t1,t2;

  UnpackExitErrorCode(t1,t2); 

  if (msg==NULL) {
    printf("$PREFIX: exit: line=%ld, file=%s (error code=%i.%i)\n",nline,fname,t1,t2);
  }
  else {
    printf("$PREFIX: exit: line=%ld, file=%s, message=%s (error code=%i.%i)\n",nline,fname,msg,t1,t2);
  }

  #ifndef __CUDA_ARCH__  
  char str[1000];
  PrintErrorLog(str);

  switch (_GENERIC_EXIT_FUNCTION_MODE_) {
  case  _GENERIC_EXIT_FUNCTION__MPI_ABORT_: 
    MPI_Abort(MPI_COMM_WORLD,t2);
    break;
  default:
    exit(0);
  }
  #endif

  exit(0);
}

void PrintLineMark(long int nline ,char* fname ,char* msg) {
#ifdef MPI_ON
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#endif

  if (ThisThread==0) {
    if (msg==NULL) printf("$PREFIX:linemark: line=%ld, file=%s\n",nline,fname);
    else printf("$PREFIX:linemark: line=%ld, file=%s, message=%s\n",nline,fname,msg);
  }
}

//=============================================================
//print a message into the debugger stream
void Debugger::SaveDataIntoStream(void* data,int length,const char* msg) {

  struct cStreamBuffer {
    int CollCounter;
    unsigned long CheckSum;
    char CallPoint[200];
  };

  const int CheckSumBufferLength=500;
  static cStreamBuffer StreamBuffer[CheckSumBufferLength];

  static int BufferPointer=0;
  static int CallCounter=0;
  static CRC32 CheckSum;

  //create a new copy of the dbuffer stream file at the first call of this function
  if (CallCounter==0) {
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",ThisThread);
    fout=fopen(fn,"w");
    fclose(fout);
  }

  //increment the call counter
  CallCounter++;

  //update the check sum
  CheckSum.add((char*)data,length);

  //save the checksum in the buffer
  StreamBuffer[BufferPointer].CollCounter=CallCounter;
  StreamBuffer[BufferPointer].CheckSum=CheckSum.checksum();
  sprintf(StreamBuffer[BufferPointer].CallPoint,"%s",msg);
  BufferPointer++;

  if (BufferPointer>=CheckSumBufferLength) {
    //save the accumulated checksums into a file
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",ThisThread);
    fout=fopen(fn,"a");

    for (int i=0;i<BufferPointer;i++) fprintf(fout,"%i: 0x%lx\t%s\n",StreamBuffer[i].CollCounter,StreamBuffer[i].CheckSum,StreamBuffer[i].CallPoint);

    BufferPointer=0;
    fclose(fout);
  }
}

template <class T>
void Debugger::SaveDataIntoStream(T data,const char* msg) {
  Debugger::SaveDataIntoStream(&data,sizeof(T),msg);
}

//=================================================================
//functionnality to operate POSIX threads
std::vector<pthread_t> Thread::ThreadsTable;

void Thread::CreateThread(void* (*buildMesh_OneLevelRefinment_External)(void*),void* Data) {
  pthread_t th;

  pthread_create(&th, NULL, buildMesh_OneLevelRefinment_External,Data);
  ThreadsTable.push_back(th);
}

void Thread::JoinThreads() {
  for (auto th : ThreadsTable) pthread_join(th, NULL);

  ThreadsTable.clear();
}


void Thread::Sync::Spinlock::AcquireLock(cLockData* lock) {
  while(lock->flag.test_and_set(std::memory_order_acquire)) {
    sched_yield();
  }
}


void Thread::Sync::Spinlock::ReleaseLock(cLockData* lock) {
  lock->flag.clear(std::memory_order_release);
}

void Thread::Sync::SpinlockBarrier::Init(cSpinlockBarrier* barrier,int ntot) {
  barrier->nTotalBarrierThreads=ntot;
  barrier->State=0;
  barrier->counter=0;

  barrier->enter_flag=true;
  barrier->exit_flag=false;


  barrier->lock_enter.clear(std::memory_order_release);
  barrier->lock_exit.test_and_set(std::memory_order_acquire);
}


void Thread::Sync::SpinlockBarrier::Wait(cSpinlockBarrier* barrier) {
  while (barrier->enter_flag==false) {
    sched_yield();
  }

  barrier->counter++;

  if (barrier->lock_enter.test_and_set(std::memory_order_acquire)==false) {
    //this is the first thread that enters the barrier
    while (barrier->counter!=barrier->nTotalBarrierThreads) {
      sched_yield();
    }

    //close enterence into the barrier
    barrier->enter_flag=false;

    //open exit from the barrier
    barrier->exit_flag=true;

    while (barrier->counter!=1) {
      sched_yield();
    }

    //open enterence into the barrier
    barrier->lock_enter.clear(std::memory_order_release);
    barrier->counter=0;
    barrier->exit_flag=false;
    barrier->enter_flag=true;
  }
  else {
    while (barrier->exit_flag==false) {
      sched_yield();
    }

    barrier->counter--;
  }
}



