
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
#include <cerrno>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>

#include "rnd.h"

#ifndef _LOGGER_H_
#define _LOGGER_H_

/*
 * Bounded synchronous diagnostics for normal core/module logging.
 *
 * The historical cLogger below is a watchdog/call-stack recorder built around
 * fork, System-V shared memory, and named semaphores. Its lifetime semantics
 * remain unchanged for existing PIC debugging. cFileLogger supplies the
 * separate everyday contract: atomic in-process lines, explicit durable
 * flush, bounded rotation, and a fatal path visible as failure to launchers.
 *
 * Timestamps are UTC with one-second resolution. Before a complete line would
 * exceed rotation_bytes, the active file becomes `.1` and older generations
 * move toward `.N`. A single oversized line is kept intact. std::mutex
 * serializes threads in one process; MPI ranks must use rank-qualified paths.
 */
namespace AMPS {
namespace Diagnostics {

class cFileLogger {
 private:
  std::string path_;
  std::size_t rotation_bytes_;
  unsigned int backup_count_;
  FILE* stream_;
  std::mutex mutex_;

  cFileLogger(const cFileLogger&);
  cFileLogger& operator=(const cFileLogger&);

  static std::string Generation(const std::string& path,
                                unsigned int generation) {
    return path + "." + std::to_string(generation);
  }

  void Open(const char* mode) {
    stream_ = fopen(path_.c_str(), mode);
    if (stream_ == NULL)
      throw std::runtime_error("cannot open diagnostic log '" + path_ + "'");
  }

  long int CurrentSize() {
    if (fseek(stream_, 0, SEEK_END) != 0)
      throw std::runtime_error("cannot seek diagnostic log '" + path_ + "'");
    const long int size = ftell(stream_);
    if (size < 0)
      throw std::runtime_error("cannot size diagnostic log '" + path_ + "'");
    return size;
  }

  void Rotate() {
    if (fflush(stream_) != 0)
      throw std::runtime_error("cannot flush diagnostic log before rotation");
    fclose(stream_);
    stream_ = NULL;
    if (backup_count_ == 0) {
      Open("w");
      return;
    }
    (void)remove(Generation(path_, backup_count_).c_str());
    for (unsigned int generation = backup_count_; generation > 1;
         --generation) {
      const std::string old_name = Generation(path_, generation - 1);
      const std::string new_name = Generation(path_, generation);
      if (rename(old_name.c_str(), new_name.c_str()) != 0 && errno != ENOENT)
        throw std::runtime_error("cannot rotate diagnostic log generation");
    }
    if (rename(path_.c_str(), Generation(path_, 1).c_str()) != 0 &&
        errno != ENOENT)
      throw std::runtime_error("cannot rotate active diagnostic log");
    Open("w");
  }

  static std::string TimestampUTC() {
    const time_t now = time(NULL);
    struct tm broken_down;
#if defined(_WIN32)
    if (gmtime_s(&broken_down, &now) != 0)
      throw std::runtime_error("cannot construct UTC diagnostic timestamp");
#else
    if (gmtime_r(&now, &broken_down) == NULL)
      throw std::runtime_error("cannot construct UTC diagnostic timestamp");
#endif
    char buffer[32];
    if (strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ",
                 &broken_down) == 0)
      throw std::runtime_error("cannot format UTC diagnostic timestamp");
    return std::string(buffer);
  }

 public:
  explicit cFileLogger(const std::string& path,
                       std::size_t rotation_bytes = 0,
                       unsigned int backup_count = 0)
      : path_(path), rotation_bytes_(rotation_bytes),
        backup_count_(backup_count), stream_(NULL) {
    if (path_.empty())
      throw std::invalid_argument("diagnostic log path must not be empty");
    Open("a");
  }

  ~cFileLogger() {
    if (stream_ != NULL) {
      (void)fflush(stream_);
      fclose(stream_);
      stream_ = NULL;
    }
  }

  void Log(const std::string& level, const std::string& message) {
    if (level.empty())
      throw std::invalid_argument("diagnostic log level must not be empty");
    const std::string line = TimestampUTC() + " [" + level + "] " +
                             message + "\n";
    std::lock_guard<std::mutex> guard(mutex_);
    const long int current_size = CurrentSize();
    if (rotation_bytes_ > 0 && current_size > 0 &&
        static_cast<std::size_t>(current_size) + line.size() >
            rotation_bytes_)
      Rotate();
    if (fwrite(line.data(), 1, line.size(), stream_) != line.size())
      throw std::runtime_error("cannot write diagnostic log '" + path_ + "'");
  }

  void Flush() {
    std::lock_guard<std::mutex> guard(mutex_);
    if (fflush(stream_) != 0)
      throw std::runtime_error("cannot flush diagnostic log '" + path_ + "'");
#if !defined(_WIN32)
    // fflush transfers libc data to the kernel; fsync is the stronger receipt
    // needed before deliberate fatal termination.
    if (fsync(fileno(stream_)) != 0)
      throw std::runtime_error("cannot synchronize diagnostic log '" + path_ + "'");
#endif
  }

  static void Fatal(const std::string& path, const std::string& message,
                    int status = EXIT_FAILURE) {
    cFileLogger logger(path);
    logger.Log("FATAL", message);
    logger.Flush();
    // Normalize an accidental zero so a fatal diagnostic never reports PASS.
    std::exit(status == 0 ? EXIT_FAILURE : status);
  }
};

}  // namespace Diagnostics
}  // namespace AMPS

template <class T>
class cLogger {
public:
   static const int InFunctionDataTableLength=20;
   int thread_mpi_rank;

   //semaphore to sync Server and the process that forked the logger
   sem_t *sem_id;

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
 
   //remove the semaphore
   void remove_semaphore() {
     if (sem_close(sem_id)!=0) {
       perror("Child: [sem_close] failed\n");
       return;
     }
      
     if (sem_unlink(fname)<0) {
       perror("Child: [sem_unlink] failed\n");
       return;
     }
   }

   //remove key file
   void remove_key_file() {
     char cmd[1000];

     sprintf(cmd,"rm -f %s", fname);
     system(cmd);
   }
      
   cLogger() {
     thread_mpi_rank=0;
     sem_id=NULL;
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
     double start_time,time_limit;
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
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].CallTimeValue=time(NULL);
   }

   void func_enter(int nline,const char* fname,T* InData,int iPrintParameter,double time=-1.0){
     //wait semaphore
     sem_wait(sem_id);

     //update the function table
     data_ptr->FunctionCallTableIndex++;

     //if InFunctionDataTableIndex reached ith max value -> move all 
     if (data_ptr->FunctionCallTableIndex==InFunctionDataTableLength) { 
       for (int i=1;i<data_ptr->FunctionCallTableIndex;i++) {
         data_ptr->FunctionCallTable[i-1]=data_ptr->FunctionCallTable[i];
       } 

       data_ptr->FunctionCallTableIndex=InFunctionDataTableLength-1;
     }


     timespec begin; 
     clock_gettime(CLOCK_REALTIME, &begin);

     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].start_time=begin.tv_sec+1e-9*begin.tv_nsec;
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].TimedFunctionExecution=(time>0.0) ? true : false;
     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].time_limit=time;
     std::time(&data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].CallTimeValue);

     sprintf(data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].fname,"%s",fname);

     data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].NewEntry(nline,InData);  

     //post the semaphore
     sem_post(sem_id); 
   }

   void func_exit() {
     sem_wait(sem_id);

     data_ptr->FunctionCallTableIndex--;

     sem_post(sem_id);
   }

   void InitLogger(int InMpiRank,bool fork_daemon=false) {
     thread_mpi_rank=InMpiRank; 

     //remove the terminate signal file 
     remove("AMPS_LOGGER_TERMINATE");

    //generate unique file name
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
    pid_t pid=-1;

    if (fork_daemon==true) fork();

    if (pid==0) {
      //this is the child. The child will be the logger 
      key_t ShmKey;
      int ShmID;

      //cleate second fork
      pid=fork();

      if (pid<0) {
        perror("Unable to second fork\n");
      }
      else if (pid>0) {
        exit(0);
      } 

      //int i, fd;

      /* child continues here */
      /* now close all extra fds */
//      for (i=getdtablesize()-1; i>=0; --i) close(i);

      /* change tty */
//      fd = open("/dev/tty", O_RDWR);
//      ioctl(fd, TIOCNOTTY, 0);
//      close(fd);
//      chdir("/");
      umask(022); /* set a default for dumb programs */

      setpgid(0,0);  /* set the process group */
//      fd=open("/dev/null", O_RDWR); /* stdin */
//      dup(fd); /* stdout */
//      dup(fd); /* stderr */

      //set the semaphore and wait for the parent to finish initialization of the data
      sem_id=sem_open(fname,O_CREAT,0600,0);

      if (sem_id==SEM_FAILED) {
        perror("Child: [sem_open] failed\n");
        exit(0);
      }

      //wait 
      if (sem_wait(sem_id)<0) {
        perror("Child: [sem_wait] fail\n");
      }

      //allocated shared memory
      ShmKey=ftok(fname,'a');
      ShmID=shmget(ShmKey,sizeof(cLoggerData),IPC_CREAT);

      if (ShmID<0) {
        exit(0);
      }

      remove_key_file();

      setsid();
      sem_post(sem_id);

      data_ptr=(cLoggerData*)shmat(ShmID,NULL,0);
      Server();
      exit(0);
    }
    else {
      //this the the parent process
      key_t ShmKey;
      int ShmID;

      sem_id=sem_open(fname,O_CREAT,0600,0);

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

    //in case when logger is not forked, start the logger with system()
    if (fork_daemon==false) {
    //  char cmd[200];


        setsid();

      pid_t pid=fork();

    signal(SIGCHLD, SIG_IGN);
    signal(SIGHUP, SIG_IGN);


      if (pid==0) {
        setsid();
        remove(fname);

       pid=fork();

       if (pid>0) exit(EXIT_SUCCESS);

        if (pid==0) {

setsid();
    signal(SIGCHLD, SIG_IGN);
    signal(SIGHUP, SIG_IGN);

      umask(0); /* set a default for dumb programs */
      setpgid(0,0);  /* set the process group */

    /* Close all open file descriptors */
    int x;
    for (x = sysconf(_SC_OPEN_MAX); x>=0; x--)
    {
 //       close (x);
    }

//
        printf("Logger started [mpi rank=%i], pid=%i,parent pid=%i\n",thread_mpi_rank,getpid(),getppid());

//        sprintf(cmd,"nohup logger -name %s &",fname);
//        system(cmd);


        Server();
}

//sleep(1);
        exit(EXIT_SUCCESS);
      }
    }

  }

  void InitLogger() {
   InitLogger(0);
  }

   void Server() {
     while (true) {
       //scan through all functins to check timing
       sem_wait(sem_id); 
       
       for (int i=0;i<=data_ptr->FunctionCallTableIndex;i++) {
         if (data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].TimedFunctionExecution==true) {
           double t_now;
           timespec now;     

           clock_gettime(CLOCK_REALTIME, &now);
           t_now=now.tv_sec+1e-9*now.tv_nsec;

           if (t_now-data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].start_time>data_ptr->FunctionCallTable[data_ptr->FunctionCallTableIndex].time_limit) {
             char msg[200];

             sprintf(msg,"The parent process (%i) was terminated: function %s [FunctionCallTableIndex=%i, mpi rank=%i] ran out of time",
                data_ptr->parent_pid,data_ptr->FunctionCallTable[i].fname,i,thread_mpi_rank);

             printf("%s\n",msg);

             PrintLog(msg);


             //output the termination code file and waite a bit before killing the parent process 
             FILE *fsignal=fopen("AMPS_LOGGER_TERMINATE","w");
             fprintf(fsignal,"Oops....");
             fclose(fsignal);

             sleep(10);

             kill(data_ptr->parent_pid,9); ////SIGINT);

             remove_semaphore();
             return;
           }
         }
       }           

       sem_post(sem_id);

       //verify that the terminate signal file is present
       struct stat buffer;

       if (stat("AMPS_LOGGER_TERMINATE",&buffer)==0) {
         //the file does not exist
         
         printf("The parent process (%i) was terminated: output the log\n",data_ptr->parent_pid);
         PrintLog();
         remove_semaphore();
         return;         
       } 

       //verify that the parent process is still alive
       if (0!=kill(data_ptr->parent_pid,0)) {
         //the parent process was exited 
         printf("The parent process (%i) was terminated: output the log\n",data_ptr->parent_pid);
         PrintLog();

         //output the termination code file and waite a bit before killing the parent process 
         FILE *fsignal=fopen("AMPS_LOGGER_TERMINATE","w");
         fprintf(fsignal,"Oops....");
         fclose(fsignal);
         
         sleep(10);

         remove_semaphore();
         return;
       }

       sleep(1);
     }  
   }


   void PrintLog(const char* msg=NULL) {
     char fname[200];
     FILE *fout;

     sprintf(fname,"logger.rank=%i.log",thread_mpi_rank); 
     fout=fopen(fname,"w");

     if (msg!=NULL) fprintf(fout,"%s\n",msg);

     for (int i=0;i<=data_ptr->FunctionCallTableIndex;i++) {
       fprintf(fout,"%i: function=%s\n",i,data_ptr->FunctionCallTable[i].fname);

       data_ptr->FunctionCallTable[i].PrintLog(fout);  
     }

     fclose(fout);
   } 
};




#endif
