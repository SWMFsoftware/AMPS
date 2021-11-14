
//class can beused to cleate a log file of a funning process vy forking a logger 
//fork() is used to create a process the conduct looking; the logger process is terminated after the partent process is terminated 
//data exchange between the logger and the parent process is dome though shared memory 

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

#include "rnd.h"

#ifndef _LOGGER_H_
#define _LOGGER_H_ 

template <class T>
class cLogger {
public:
   static const int InFunctionDataTableLength=20;
   int thread_mpi_rank;

   //in-function location
   class cInFunctionDataEl {
   public:
     int nline;
     clock_t time;
     T data;
     time_t CallTimeValue;
   }; 


   //unique data file used to reserve shared memory and init the semaphore 
   const static int fname_length=15;
   char fname[fname_length];
 
   cLogger() {
     thread_mpi_rank=0;
   }

   class cInFunctionData {
   public:
     cInFunctionDataEl InFunctionDataTable[InFunctionDataTableLength]; 
     int DataTableIndex;

     void NewEntry(int nline,cInFunctionDataEl* InData) {
       if (DataTableIndex!=0) DataTableIndex++;

       if (DataTableIndex==InFunctionDataTableLength) {
         for (int i=1;i<DataTableIndex;i++) InFunctionDataTable[i-1]=InFunctionDataTable[i];

         DataTableIndex=InFunctionDataTableLength-1;
       }

       InFunctionDataTable[DataTableIndex].data=*InData;
       InFunctionDataTable[DataTableIndex].nline=nline;
       InFunctionDataTable[DataTableIndex].time=clock()/CLOCKS_PER_SEC;    
       InFunctionDataTable[DataTableIndex].CallTimeValue=time(NULL);
     }

     void ResetIndex() {DataTableIndex=0;}

     cInFunctionData() {
       ResetIndex();
     }
   };


   //function call table 
   class cFunctionCallTable {
   public:
     cInFunctionDataEl DataTable[InFunctionDataTableLength];
     int DataTableIndex,PrintParameter;
     char fname[200];

     bool TimedFunctionExecution;
     clock_t start_time,last_access_time,time_limit;
     time_t CallTimeValue;

     void ResetIndex() {DataTableIndex=0;}

     cFunctionCallTable() {
       ResetIndex(); 

       PrintParameter=-1;
     }

     void NewEntry(int nline,T* InData) {
       if (DataTableIndex==InFunctionDataTableLength) {
         for (int i=1;i<DataTableIndex;i++) DataTable[i-1]=DataTable[i];
     
         DataTableIndex=InFunctionDataTableLength-1;
       }

       DataTable[DataTableIndex].data=*InData;
       DataTable[DataTableIndex].nline=nline;
       DataTable[DataTableIndex].time=clock()/CLOCKS_PER_SEC;
       DataTable[DataTableIndex].CallTimeValue=time(NULL);

       DataTableIndex++;
     }

     void PrintLog(FILE* fout) {
       for (int i=0;i<DataTableIndex;i++) {
         tm *ct=localtime(&DataTable[i].CallTimeValue);

         fprintf(fout,"line=%d: %i/%i %i:%i:%i\n",DataTable[i].nline,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec); 

         DataTable[i].data.PrintLog(PrintParameter,fout);
       }
     } 
   };   


   class cLoggerData {
   public:
     pid_t parent_pid;
     cFunctionCallTable FunctionCallTable[InFunctionDataTableLength];
     int FunctionCallTableIndex;

     cLoggerData() {
       FunctionCallTableIndex=-1;
     }
   };

   cLoggerData *data_ptr;

   void add_data_point(int nline,T* InData) {
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].NewEntry(nline,InData); 
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].last_access_time=clock()/CLOCKS_PER_SEC;
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].CallTimeValue=time(NULL);
   }

   void func_enter(int nline,const char* fname,T* InData,int iPrintParameter,double time=-1.0){
     data_ptr->FunctionCallTableIndex++;

     //if InFunctionDataTableIndex reached ith max value -> move all 
     if (data_ptr->FunctionCallTableIndex==InFunctionDataTableLength) { 
       for (int i=1;i<data_ptr->FunctionCallTableIndex;i++) {
         data_ptr->FunctionCallTable[i-1]=data_ptr->FunctionCallTable[i];
       } 

       data_ptr->FunctionCallTableIndex=InFunctionDataTableLength-1;
     }

     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].start_time=clock()/CLOCKS_PER_SEC;
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].TimedFunctionExecution=(time>0.0) ? true : false;
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].time_limit=time;
     std::time(&data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].CallTimeValue);

     sprintf(data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].fname,"%s",fname);

     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].NewEntry(nline,InData);  
   }

   void func_exit() {
     data_ptr->FunctionCallTableIndex--;
   }

   void InitLogger(int InMpiRank) {

     thread_mpi_rank=InMpiRank; 

    //generate unique file name
//    const int fname_length=15;
//    char fname[fname_length];

    do {
      for (int i=0;i<fname_length-1;i++) {
        int c_min,c_max;
        double r=rnd();

        if (r<0.3) {
          c_min=48,c_max=58;
        }
        else if (r<0.6) {
          c_min=65,c_max=91;
        }
        else {
          c_min=97,c_max=123;
        }    

        fname[i]=c_min+(int)(rnd()*(c_max-c_min));
      } 

      fname[fname_length-1]=0;

      struct stat buffer;
        
      if (stat(fname,&buffer)!=0) {
        //the file does not exist
        FILE *fout=fopen(fname,"w");
           
        fclose(fout);
        break;
      }
    }
    while (true);

    //fork process
    pid_t pid=fork();

    if (pid==0) {
      //this is the child. The child will be the logger 
      key_t ShmKey;
      int ShmID;

      //set the semaphore and wait for the parent to finish initialization of the data
      sem_t *sem_id=sem_open(fname,O_CREAT,0600,0);

      if (sem_id==SEM_FAILED) {
        perror("Child: [sem_open] failed\n");
        exit(0);
      }

      //wait 
      if (sem_wait(sem_id)<0) {
        perror("Child: [sem_wait] fail\n");
      }

      //remove the semaphore
      if (sem_close(sem_id)!=0) {
        perror("Child: [sem_close] failed\n");
        return;
      }

      if (sem_unlink(fname)<0) {
        perror("Child: [sem_unlink] failed\n");
        return;
      }

      //allocated shared memory
      ShmKey=ftok(fname,'a');
      ShmID=shmget(ShmKey,sizeof(cLoggerData),IPC_CREAT);

      if (ShmID<0) {
        exit(0);
      }

      data_ptr=(cLoggerData*)shmat(ShmID,NULL,0);
      Server();

    }
    else {
      //this the the parent process
      key_t ShmKey;
      int ShmID;

      sem_t *sem_id=sem_open(fname,O_CREAT,0600,0);

      if (sem_id==SEM_FAILED) {
        perror("Parent: [sem_open] failed\n");
        return;
      }

      int size=sizeof(cLoggerData);

      //get the shared memory 
      ShmKey=ftok(fname,'a');
      ShmID=shmget(ShmKey,size,IPC_CREAT|0666);

      if (ShmID<0) {
        exit(0);
      }

      data_ptr=(cLoggerData*)shmat(ShmID,NULL,0);
      data_ptr->parent_pid=getpid();
      data_ptr->FunctionCallTableIndex=-1;

      //release semaphore 
      if (sem_post(sem_id)<0) {
        perror("Parent: [semp_post] failed\n");
        return;
      }

    }
 
  }

  void InitLogger() {
   InitLogger(0);
  }

   void Server() {
     while (true) {
       //scan through all functins to check timing
       for (int i=0;i<=data_ptr->FunctionCallTableIndex;i++) {
         if (data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].TimedFunctionExecution==true) {
           if (clock()/CLOCKS_PER_SEC-data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].start_time>data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].time_limit) {
             PrintLog();

             exit(0);
           }
         }
       }           


       //verify that the parent process is still alive
       if (0!=kill(data_ptr->parent_pid,0)) {
         //the parent process was exited 
         printf("The parent process (%i) was terminated: output the log\n",data_ptr->parent_pid);
         PrintLog();


         //remove the key file and exit
        char cmd[1000];

        sprintf(cmd,"rm -f %s", fname);
        system(cmd);
         
         exit(0);
       }

       sleep(0.1);
     }  
   }


   void PrintLog() {
     char fname[200];
     FILE *fout;

     sprintf(fname,"logger.rank=%i.log",thread_mpi_rank); 
     fout=fopen(fname,"w");

     for (int i=0;i<=data_ptr->FunctionCallTableIndex;i++) {
       fprintf(fout,"%i: function=%s\n",i,data_ptr->FunctionCallTable[i].fname);

       data_ptr->FunctionCallTable[i].PrintLog(fout);  
     }

     fclose(fout);
   } 
};




#endif
