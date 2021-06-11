//$Id$
//the namespace contains function that reads the AMPS' section of the SWMF PARAM.in file

#ifndef _AMPS2SWMF_
#define _AMPS2SWMF_

#include <iomanip>
#include <iostream>
#include <sstream>
#include <climits>

using namespace std;

#include "pic.h"

//AMPS component code 
#define _AMPS_SWMF_UNDEFINED_ -1
#define _AMPS_SWMF_PT_ 0
#define _AMPS_SWMF_PC_ 1 

namespace AMPS2SWMF {
  //the name of hte component within the SWMF (PT/PC)
  extern char ComponentName[10];
  extern int ComponentID;

  //parameters of the current SWMF session
  extern int iSession;
  extern double swmfTimeSimulation;
  extern bool swmfTimeAccurate;

  //amps_init_flag
  extern bool amps_init_flag;

  //init amps mesh flag
  extern bool amps_init_mesh_flag;

  //the counter of the field line update events since beginning of a new session
  extern int FieldLineUpdateCounter; 

  //the table containing the field line segment indexes where the CME shock is currently localted
  class cShockData {
  public:
    int iSegmentShock;
    double ShockSpeed,DownStreamDensity;

    cShockData() {
      iSegmentShock=-1;
      ShockSpeed=0.0,DownStreamDensity=0.0;
    }
  };

  extern cShockData *ShockData;

//  extern int *iShockWaveSegmentTable;
  

  //AMPS execution timer 
  extern PIC::Debugger::cTimer ExecutionTimer; 

  //hook that AMPS applications can use so a user-defined function is called at the end of the SWMF simulation
  typedef void (*fUserFinalizeSimulation)(); 
  extern fUserFinalizeSimulation UserFinalizeSimulation;

  //the namespace containds variables used in heliosphere simulations 
  namespace Heliosphere {
    extern double rMin;
  }

  //the location of the Earth as calcualted with the SWMF. Used for heliophysics modeling  
  extern double xEarthHgi[3]; 

  namespace PARAMIN {
    int read_paramin(stringstream *param);


    //reading of the swmf PARAM.in file
    inline void char_to_stringstream(char* chararray, int nlines, int linelength,std::stringstream  *ss){
      int i,j;

      //add the end of the line symbols
      for (i=0;i<nlines;i+=nlines) chararray[linelength*(i+1)-1]='\n';

      //remove the "empty" lines
      for (i=0;i<nlines;i++) {
        bool flag=false;

        //check is the line is "empty";
        for (j=0;j<linelength-1;j++) if (chararray[i*linelength+j]!=' ') {
          flag=true;
          break;
        }

        //if the line is "empty" -> move the next lines
        if (flag==false) {
          if (i!=nlines-1) {
            for (;i<nlines-1;i++) for (j=0;j<linelength;j++) chararray[i*linelength+j]=chararray[(i+1)*linelength+j];
            --nlines,i=-1;
            continue;
          }
          else nlines--;
        }
      }


      ss->str("");
      ss->write(chararray,nlines*linelength);
    }

    template <class T>
    void read_var(std::stringstream *ss, std::string description, T *var){
      *ss >> *var;
      ss->ignore(INT_MAX, '\n');


      if (PIC::ThisThread==0) std::cout<<std::left<<std::setw(50)<<*var<<description<<std::endl;
    }

    inline void read_var(std::stringstream *ss, std::string description, bool *var){
      std::string text;
      *ss >> text;
      ss->ignore(INT_MAX, '\n');
      *var = false;
      if(text == "T") *var = true;


      if (PIC::ThisThread==0) std::cout<<std::left<<std::setw(50)<<text<<description<<std::endl;
    }

    inline void get_next_command(std::stringstream *ss, std::string *id){
      ss->ignore(INT_MAX, '#');
      ss->unget();
      *ss >> *id;


      if (PIC::ThisThread==0) std::cout<<"\n"<<*id<<std::endl;
    }
  }
}




#endif
