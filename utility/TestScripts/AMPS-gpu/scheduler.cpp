/*
 * main.cpp
 *
 *  Created on: Feb 4, 2020
 *      Author: vtenishe
 */





#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <map>
#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <iostream>   
#include <thread>

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include <sys/time.h>
#include <sys/resource.h>



const int status_default=0;
const int status_compiled=1;
const int status_launched=2;
const int status_completed=3;



class cJobTableElement {
public:
  int status;
  char cmd[500];

  cJobTableElement() {
    status=status_default;
  }
};



int main(int argc, char* argv[]) {
  int i;
  bool status_pgi_test=false;
  bool status_intel_test=false;
  bool status_gcc_test=false;
  char test_dir[500];
  int nTestRoutineThreads=1;
  int nTotalCompilerCases=0;

  const int nParallelTestExecutionThreads=3;

  /*
  parse command line parameters:
  argument line parameter: -pgi -gcc -intel -path [directory that containes installations of AMPS for 
  each compiler, e.g. path/Intel/AMPS, path/GNU/AMPS, path/PGI/AMPS, etc] -threads [the number of threads for executing each set of tests]
  */

  for (i=1;i<argc;i++) {
    char *endptr;

    if (strcmp("-pgi",argv[i])==0) {
      status_pgi_test=true;
      nTotalCompilerCases++;
    }
    else if (strcmp("-gcc",argv[i])==0) {
      status_gcc_test=true;
      nTotalCompilerCases++;
    }
    else if (strcmp("-intel",argv[i])==0) {
      status_intel_test=true;
      nTotalCompilerCases++;
    }
    else if (strcmp("-path",argv[i])==0) {
      if (++i>=argc) break;
      sprintf(test_dir,"%s",argv[i]);
    }
    else if (strcmp("-threads",argv[i])==0) {
      if (++i>=argc) break;
      nTestRoutineThreads=(int)strtol(argv[i],&endptr,10);
    }
    else {
      printf("Option %s is not found\n",argv[i]);
      exit(0);
    }

  }


  //populate the job table
  int index=0;
  int index_pgi_min=-1,index_pgi_max=-1;
  int index_gcc_min=-1,index_gcc_max=-1;
  int index_intel_min=-1,index_intel_max=-1;

  cJobTableElement* JobTable=new cJobTableElement[nTotalCompilerCases*nTestRoutineThreads];

  if (status_pgi_test==true) {
    index_pgi_min=index;

    for (i=1;i<=nTestRoutineThreads;i++) {
      sprintf(JobTable[index].cmd,"cd %s/PGI/AMPS; make TESTMPIRUN4=\"mpirun -np 4\" MPIRUN=\"mpirun -np 4\" test_run_thread%i   >& test_amps.run.thread%i.log",test_dir,i,i);
      index_pgi_max=index++;
    }
  }

  if (status_gcc_test==true) {
    index_gcc_min=index;

    for (i=1;i<=nTestRoutineThreads;i++) {
      sprintf(JobTable[index].cmd,"cd %s/GNU/AMPS; make TESTMPIRUN4=\"mpirun -np 4\" MPIRUN=\"mpirun -np 4\" test_run_thread%i   >& test_amps.run.thread%i.log",test_dir,i,i);
      index_gcc_max=index++;
    }
  }

  if (status_intel_test==true) {
    index_intel_min=index;

    for (i=1;i<=nTestRoutineThreads;i++) {
      sprintf(JobTable[index].cmd,"cd %s/Intel/AMPS; make TESTMPIRUN4=\"mpirun -np 4\" MPIRUN=\"mpirun -np 4\" test_run_thread%i   >& test_amps.run.thread%i.log",test_dir,i,i);
      index_intel_max=index++;
    }
  }

  //launch compiling threads
  std::thread CompilingThreadTable[nTotalCompilerCases];
  index=0;

  class cCompileThread {
  public:
    char cmd[500];
    int job_index_min,job_index_max;
    cJobTableElement *JobTable;

    void run() {
      system(cmd);

      for (int i=job_index_min;i<=job_index_max;i++) JobTable[i].status=status_compiled;
    }
  } CompileThreadTable[3];

  if (status_pgi_test==true) {
    sprintf(CompileThreadTable[index].cmd,"cd %s; PGI/AMPS/utility/TestScripts/AMPS-gpu/CompilePGI.sh",test_dir);
    CompileThreadTable[index].job_index_min=index_pgi_min,CompileThreadTable[index].job_index_max=index_pgi_max;
    CompileThreadTable[index].JobTable=JobTable;

    CompilingThreadTable[index]=std::thread(&cCompileThread::run,CompileThreadTable+index);
    index++;
  }

  if (status_gcc_test==true) {
    sprintf(CompileThreadTable[index].cmd,"cd %s; GNU/AMPS/utility/TestScripts/AMPS-gpu/CompileGNU.sh",test_dir);
    CompileThreadTable[index].job_index_min=index_gcc_min,CompileThreadTable[index].job_index_max=index_gcc_max;
    CompileThreadTable[index].JobTable=JobTable;

    CompilingThreadTable[index]=std::thread(&cCompileThread::run,CompileThreadTable+index);
    index++;
  }

  if (status_intel_test==true) {
    sprintf(CompileThreadTable[index].cmd,"cd %s; Intel/AMPS/utility/TestScripts/AMPS-gpu/CompileIntel.sh",test_dir);
    CompileThreadTable[index].job_index_min=index_intel_min,CompileThreadTable[index].job_index_max=index_intel_max;
    CompileThreadTable[index].JobTable=JobTable;


    CompilingThreadTable[index]=std::thread(&cCompileThread::run,CompileThreadTable+index);
    index++;
  }

  bool status_all_launched=false;
  int id_job_execute;

  class cExecutionThread {
  public:
    cJobTableElement *job;
    std::thread thread;
    bool status_complete;
    int call_cnt;

    cExecutionThread() {
      status_complete=true;
      job=NULL;
      call_cnt=0;
    }

    void run() {
      system(job->cmd);
      job->status=status_completed;
      status_complete=true;
    }

  } ExecutionThreadTable[nParallelTestExecutionThreads];


  do {
    status_all_launched=true;
    id_job_execute=-1;

    //find a new job to execute
    for (int i=0;i<nTotalCompilerCases*nTestRoutineThreads;i++) {
      if (JobTable[i].status!=status_launched) {
        status_all_launched=false;
      }

      if (JobTable[i].status==status_compiled) {
        status_all_launched=false;

        id_job_execute=i;
        JobTable[i].status=status_launched;

        break;
      }
    }

    if (id_job_execute==-1) { 
      //no jobs for execution are found
      sleep(1);
      continue;
    }

    //waite until any of the execution threads is available and then lounch the job
    bool thread_found=false;
    int id_execution_thread=-1;

    do {
      for (int i=0;i<nParallelTestExecutionThreads;i++) {
        if (ExecutionThreadTable[i].status_complete==true) {
          //found thread to execute the test
          thread_found=true;
          ExecutionThreadTable[i].status_complete=false;
          ExecutionThreadTable[i].job=JobTable+id_job_execute;

	  //join threads used before
	  if (ExecutionThreadTable[i].call_cnt!=0) ExecutionThreadTable[i].thread.join();

	  ExecutionThreadTable[i].call_cnt++;
          ExecutionThreadTable[i].thread=std::thread(&cExecutionThread::run,ExecutionThreadTable+i);
	  break;
        }
      }

      sleep(1);
    }
    while (thread_found==false);

  }
  while (status_all_launched==false);


  //wait for completing all test execution
  bool execution_completed;

  do {
    execution_completed=true;

    for (int i=0;i<nParallelTestExecutionThreads;i++) {
      if (JobTable[i].status!=status_completed) {
        execution_completed=false;
        sleep(1);
        break;
      }
    }
  }
  while (execution_completed==false);

  //join compiling threads
  for (int i=0;i<nTotalCompilerCases;i++) CompilingThreadTable[i].join(); 

  return 0;
}




















