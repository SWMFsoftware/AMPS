//functions to run mutiple AMPS simulations and compare the resuilts runtime

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <semaphore.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <iostream>
#include <ctime>


#include "pic.h" 

char PIC::Debugger::ConcurrentDebug::Key[200];
sem_t *PIC::Debugger::ConcurrentDebug::sem_id; 
PIC::Debugger::ConcurrentDebug::cData *PIC::Debugger::ConcurrentDebug::data_ptr;


void PIC::Debugger::ConcurrentDebug::GenerateKey() {
  char base[200];

  if (PIC::ThisThread==0) {
    cout << "Enter base name:";
    cin.getline(base,200);
  }

  MPI_Bcast(base,200,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR); 

  sprintf(Key,"%s.amps_rank=%i",base,PIC::ThisThread);

  FILE *fkey=fopen(Key,"w");
  fclose(fkey);
}

void PIC::Debugger::ConcurrentDebug::RemoveKeyFile() {
  remove(Key);
}

void PIC::Debugger::ConcurrentDebug::InitSharedMomery() {
  key_t ShmKey;
  int ShmID;
  int size=sizeof(cData);

  ShmKey=ftok(Key,'a');
  ShmID=shmget(ShmKey,size,IPC_CREAT);

  if (ShmID<0) {
    exit(__LINE__,__FILE__,"Error: cannot allocate shared memory");
  }

  data_ptr=(cData*)shmat(ShmID,NULL,0);


  data_ptr->clear();


//  if (sem_wait(sem_id)<0) {
//    exit(__LINE__,__FILE__,"Error: sem_wait fail");
//  }

}

void PIC::Debugger::ConcurrentDebug::InitSemaphore() {
  sem_id=sem_open(Key,O_CREAT,0600,0);

  if (sem_id==SEM_FAILED) {
    perror("Child: [sem_open] failed\n");
    exit(0);
  }
}

void PIC::Debugger::ConcurrentDebug::Trap() {
  while (true);
} 

void PIC::Debugger::ConcurrentDebug::NewEntry(cData* d,int nline,char const *fname) {

  //wait semaphore
  sem_wait(sem_id);

  //save data
  *data_ptr=*d;
  data_ptr->nline=nline;
  sprintf(data_ptr->fname,"fname=%s",fname); 

  //post semaphore
  sem_post(sem_id);
}
  
  


