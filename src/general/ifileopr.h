//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef IFILEOPR 
#define IFILEOPR


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstddef>
#include <cstring>
#include <string>

#include "specfunc.h"

#define init_str_maxlength 10000

namespace CiFileOperationsDetail {
  /*
   * CiFileOperations::fname is a long-standing public fixed-size field.  Keep
   * that representation for source/ABI compatibility, but never use sprintf
   * (or a silently truncating snprintf) to populate it.  A truncated path can
   * name a different file, so an input which does not fit is an explicit fatal
   * configuration error rather than a recoverable shortening operation.
   *
   * The function returns false only for debugger/test configurations in which
   * the AMPS exit(line,file,message) trap is intercepted and allowed to
   * return.  Production exit handlers terminate at the call below.  Keeping
   * the return value makes callers safe in both execution modes: no file I/O
   * is attempted after a rejected name.
   */
  template <std::size_t DestinationSize>
  inline bool StoreFileName(char (&destination)[DestinationSize],const char *source) {
    if (source==NULL) {
      exit(__LINE__,__FILE__,
          "CiFileOperations received a null input-file name");
      return false;
    }

    const std::size_t sourceLength=std::strlen(source);

    // Reserve one byte for the terminating NUL required by fopen and by
    // legacy clients that inspect the public fname character array.
    if (sourceLength>=DestinationSize) {
      const std::string message=
          "Input-file name requires "+std::to_string(sourceLength)+
          " characters, but CiFileOperations::fname can store at most "+
          std::to_string(DestinationSize-1)+
          "; refusing to truncate the path";

      exit(__LINE__,__FILE__,message.c_str());
      return false;
    }

    // The bounds check above proves that the complete name and its NUL fit.
    // memmove also remains correct if a legacy caller passes fname itself.
    std::memmove(destination,source,sourceLength+1);
    return true;
  }
}

class CiFileOperations {
public:
  FILE* fd;
  char fname[1000],init_str[init_str_maxlength]; 
  long int line;

  CiFileOperations() {
    fd=NULL;
    fname[0]='\0';
    line=-1;
  };

  FILE* openfile(const char* ifile) {
    /*
     * Store the exact path before opening it.  StoreFileName validates the
     * fixed legacy field and reports an overlength name without either a
     * buffer overflow or ambiguous truncation.  The explicit return protects
     * debugger builds whose AMPS exit trap may be intercepted and resumed.
     */
    if (CiFileOperationsDetail::StoreFileName(fname,ifile)==false) {
      fd=NULL;
      return NULL;
    }

    line=0;
    if ((fd=fopen(fname,"r"))==NULL) {
      // Capture errno immediately: constructing the dynamic diagnostic may
      // call library routines that are permitted to change it.  std::string
      // removes the second fixed-buffer overflow reported by fortified libc
      // when a long (but otherwise valid) path cannot be opened.
      const int openError=errno;
      const std::string message=
          "Cannot open input file '"+std::string(fname)+"': "+
          std::strerror(openError);

      exit(__LINE__,__FILE__,message.c_str());

      // AMPS exit normally terminates.  Return NULL if a debugger intercepts
      // the trap so that the caller cannot mistake the failed open for a
      // usable stream.
      return NULL;
    } 
     
    return fd;
  }; 

  void closefile() {
    fclose(fd);
  };

  void rewindfile() {
    rewind(fd);
    line=0;
  }

  void moveLineBack() {
    int c;

    fseek(fd,-1,SEEK_CUR);

    do {
      fseek(fd,-1,SEEK_CUR);
      c=getc(fd);
      fseek(fd,-1,SEEK_CUR);
    }
    while ((c!='\0')&&(c!='\n')&&(c!='\r'));

    c=getc(fd);
    line--;
  };


  void setfile(FILE* input_fd,long int input_line,char* InputFile) {
    // Apply the same exact, non-truncating filename contract used by
    // openfile().  Validate first so an intercepted failure leaves the
    // existing FILE pointer and line number unchanged.
    if (CiFileOperationsDetail::StoreFileName(fname,InputFile)==false) return;

    fd=input_fd;
    line=input_line;
  };

  long int& CurrentLine() {
    return line;
  }; 

  //Separators:' ', ',', '=', ';', ':', '(', ')', '[', ']' 
  void CutInputStr(char* dest, char* src) {
    int i,j;

    if (src[0]!='"') {
      for (i=0;(src[i]!='\0')&&(src[i]!=' ')&&
        (src[i]!=',')&&(src[i]!='=')&&(src[i]!=';')&&(src[i]!=':')&&
        (src[i]!='(')&&(src[i]!=')')&&(src[i]!='[')&&(src[i]!=']');i++) dest[i]=src[i];

      dest[i]='\0';
    }
    else {
      for (i=1,j=0;(src[i]!='\0')&&(src[i]!='"');i++,j++) dest[j]=src[i];  
      dest[j]='\0';
      ++i;
    }


    for (;(src[i]!='\0')&&((src[i]==' ')||
      (src[i]==',')||(src[i]=='=')||(src[i]==';')||(src[i]==':')||
      (src[i]=='(')||(src[i]==')')||(src[i]=='[')||(src[i]==']'));i++);

    if (src[i]=='\0')
      src[0]='\0';
    else {
      for (j=0;src[j+i]!='\0';j++) src[j]=src[j+i];
      src[j]='\0';
    }
  }; 

  bool eof() {
    return (!feof(fd)) ? false : true;
  };

  bool GetInputStr(char* str,long int n, bool ConvertToUpperCase=true){
    int i,j;

    str[0]='\0';

    if (!feof(fd)) do {
      line++;
      if (fgets(str,n,fd)==NULL) {
        //error();
        str[0]='\0',init_str[0]='\0';
        return false;
      }

      for (i=0;(str[i]!='\0')&&(str[i]!='\n')&&(str[i]==' ');i++);
      for (j=0;(str[i+j]!='\0')&&(str[i+j]!='\n');j++) str[j]=str[i+j];
      str[j]='\0';

      for (i=0;str[i]!='\0';i++) if ((str[i]=='!')||(str[i]=='\r')) str[i]='\0';
    } while ((str[0]=='\0')&&(!feof(fd)));

    if (str[0]=='\0') {
      str[0]='\0',init_str[0]='\0';   
      return false; 
    }

    for(i=0;str[i]!='\0';i++) {
//      if (str[i]=='"') str[i]=' ';

      
      if (i<init_str_maxlength-1) init_str[i]=str[i];
      if ((str[i]>='a')&&(str[i]<='z')&&ConvertToUpperCase)
	str[i]=(char)(str[i]-(char)32);
      if (str[i]=='\t') str[i]=' ';
    }
    
    init_str[(i<init_str_maxlength-1) ? i : init_str_maxlength-1 ]='\0';
    return true;
  }; 

  void error(const char *msg=NULL) {
    exit(__LINE__,__FILE__,msg);
  }; 
};

#endif
