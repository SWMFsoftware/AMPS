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

#define init_str_maxlength 10000

class CiFileOperations {
public:
  FILE* fd;
  char fname[100],init_str[init_str_maxlength]; 
  long int line;

  CiFileOperations() {
    line=-1;
  };

  FILE* openfile(const char* ifile) {
    sprintf(fname,"%s",ifile); 
    line=0;
    if ((fd=fopen(ifile,"r"))==NULL) {
      printf ("FILE %s is not found\n",ifile);
      exit(0);
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
    fd=input_fd;
    line=input_line;
    sprintf(fname,"%s",InputFile);
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
    printf("$PREFIX:Error in file %s, line=%ld\n",fname,line);
    printf("$PREFIX:%s\n",init_str);
    if (msg!=NULL) printf("$PREFIX:Error message: %s\n",msg);
    exit(0);
  }; 
};

#endif
