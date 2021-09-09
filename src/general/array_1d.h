//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef ARRAY_1R 
#define ARRAY_1R

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "specfunc.h"

using namespace std;

template<class T>

class array_1d {
protected:
  T* data;
  int size;

public:

  //===================================================
  _TARGET_HOST_ _TARGET_DEVICE_
    void init(long int n) {
    if (size!=0) {
      printf("Error: initialization of allocated array_1d object\n");
      exit(__LINE__,__FILE__);
    }

    data=new T [n];
    size=n;
  }



_TARGET_HOST_ _TARGET_DEVICE_
  array_1d() {
    data=NULL,size=0;
  }

_TARGET_HOST_ _TARGET_DEVICE_
  array_1d(int n) {
    data=NULL,size=0;
    init(n);
  }

_TARGET_HOST_ _TARGET_DEVICE_
 ~array_1d() {
   if (data!=NULL) delete [] data;

   data=NULL,size=0;
 }

_TARGET_HOST_ _TARGET_DEVICE_
 bool IsAllocated() {return (data!=NULL);}

//===================================================
_TARGET_HOST_ _TARGET_DEVICE_
  int GetSize() const {
    return size;
  }

_TARGET_HOST_ _TARGET_DEVICE_
T* get_data_ptr() {return data;}
//===================================================
  T sum() const {
    T f = (T)0;

    for (long int i=0; i<size;i++) f=f+data[i];

    return f;
  }
//===================================================
  T abs() const {
    T f = (T)0;

    for (int i=0; i<size;i++) f=f+data[i]*data[i];

    return (T)sqrt(f);
  }

//===================================================
  void normalize() {
    T m=abs();

    for (int i=0;i<size;i++) data[i]/=m;
  }

//===================================================
  void clear() {
    for (int i=0;i<size;i++) data[i]=0.0;
  }

//===================================================

_TARGET_HOST_ _TARGET_DEVICE_
  T operator () (int i) const {
    if ((i<0)||(i>=size)) exit(__LINE__,__FILE__,"Error: out of range");

    return data[i];
  }

_TARGET_HOST_ _TARGET_DEVICE_
  T & operator () (int i) {
    if ((i<0)||(i>=size)) exit(__LINE__,__FILE__,"Error: out of range");

    return data[i];
  }

//===================================================
  array_1d<T>& operator = (const array_1d<T>& v) {
    for(int i=0;i<size;i++) data[i]=v.data[i];

    return *this;
  }

//===================================================
  array_1d<T>& operator = (T f) {
    for (long int i=0; i<size;i++) data[i]=f;

    return *this;
  }
//===================================================
  friend array_1d<T> operator + (const array_1d<T> &v1,const array_1d<T> &v2) {
    if (v1.size!=v2.size) {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    array_1d<T> v3(v1.size);
    for (int i=0;i<v1.size;i++) v3.data[i]=v1.data[i]+v2.data[i];

    return v3;
  }
//===================================================

  friend array_1d<T> operator - (const array_1d<T> &v1,const array_1d<T> &v2) {
    if (v1.size!=v2.size) {
      printf("Error: subtract two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    array_1d<T> v3(v1.size);
    for (int i=0;i<v1.size;i++) v3.data[i]=v1.data[i]-v2.data[i];

    return v3;
  }

//===================================================
  friend array_1d<T> operator * (const array_1d<T> &v1, const T t) {
    array_1d<T> v3(v1.size);
    for (long int i=0;i<v1.size;i++) v3.data[i]=t*v1.data[i];
    return v3;
  }

//===================================================
  friend array_1d<T> operator / (const array_1d<T> &v1, const T t) {
    if (t == 0) {
      printf("Error: divide vector by 0.\n");
      exit(__LINE__,__FILE__);
    }
    array_1d<T> v3(v1.size);
    for (int i=0;i<v1.size;i++) v3.data[i]=v1.data[i]/t;
    return v3;
  }

//===================================================
  friend array_1d<T>& operator += (array_1d<T> &v1,const array_1d<T> &v2) {
    if (v1.size!=v2.size) {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    for (int i=0;i<v1.size;i++) v1.data[i]+=v2.data[i];
    return v1;
  }
//===================================================

  friend array_1d<T>& operator -= (array_1d<T> &v1,const array_1d<T> &v2) {
    if (v1.size!=v2.size) {
      printf("Error: subtract two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    for (int i=0;i<v1.size;i++) v1.data[i]-=v2.data[i];
    return v1;
  }
//===================================================

  friend array_1d<T>& operator *= (array_1d<T> &v1,const T t) {
    for (int i=0;i<v1.size;i++) v1.data[i]*=t;
    return v1;
  }
//===================================================

  friend array_1d<T>& operator /= (array_1d<T> &v1,const T t) {
    if (t == 0) {
      printf("Error: divide vector by 0.\n");
      exit(__LINE__,__FILE__);
    }

    for (int i=0;i<v1.size;i++) v1.data[i]/=t;
    return v1;
  }
//===================================================
  friend T dot_product (const array_1d<T> &v1,const array_1d<T> &v2) {
    if (v1.size!=v2.size) {
      printf("Error: dot product of two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    T d=0;
    for (int i=0;i<v1.size;i++) d+=v1.data[i]*v2.data[i];
    return d;
  }
//===================================================

  friend array_1d<T> cross_product (const array_1d<T> &v1,
  const array_1d<T> &v2) {
    if ((v1.size!=3)||(v2.size!=3)) {
      printf("Error: cross product of non-3D vectors.\n");
      exit(__LINE__,__FILE__);
    }

    array_1d<T> v3(3);
    v3.data[0]=v1.data[1]*v2.data[2]-v1.data[2]*v2.data[1];
    v3.data[1]=v1.data[2]*v2.data[0]-v1.data[0]*v2.data[2];
    v3.data[2]=v1.data[0]*v2.data[1]-v1.data[1]*v2.data[0];

    return v3;
  }
//===================================================

  friend double mix_product (const array_1d<T> &v1, const array_1d<T> &v2,
  const array_1d<T> &v3) {
    double res;

    if ((v1.size!=3)||(v2.size!=3)||(v3.size!=3)) {
      printf("Error: mix product of non-3D vectors.\n");
      exit(__LINE__,__FILE__);
    }

    res=v1.data[0]*(v2.data[1]*v3.data[2]-v3.data[1]*v2.data[2]);
    res-=v1.data[1]*(v2.data[0]*v3.data[2]-v3.data[0]*v2.data[2]);
    res+=v1.data[2]*(v2.data[0]*v3.data[1]-v3.data[0]*v2.data[1]);

    return res;
  }
//===================================================

  friend void printf (const array_1d<T> &v) {
    if (v.size<4) {
      for (long int i=0;i<v.size;i++)
        if (v.data[i]<0) printf("  %e",(double)v.data[i]);
        else printf("   %e",(double)v.data[i]);
      printf("\n");
    }
    else for (long int i=0;i<v.size;i++) printf(" %e\n",(double)v.data[i]);
  }

  friend void fprintf (FILE* fout,const array_1d<T> &v) {
    if (v.size<4) {
      for (long int i=0;i<v.size;i++)
        if (v.data[i]<0) fprintf(fout,"  %e",(double)v.data[i]);
        else fprintf(fout,"   %e",(double)v.data[i]);
      fprintf(fout,"\n");
    }
    else for (long int i=0;i<v.size;i++) fprintf(fout," %e\n",(double)v.data[i]);
  }
//===================================================
};

#endif
