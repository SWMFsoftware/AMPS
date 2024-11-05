//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef ARRAY_2D
#define ARRAY_2D

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "specfunc.h"

template <class T>
class array_2d {
protected:
  T* data;
  int size_dim0,size_dim1;
  bool locally_allocated_data_buiffer;

public:


//===================================================

  void init(int n0,int n1) {
    if ((n0<=0)||(n1<=0)) {
      printf("Error: allocation of array_2d object\n");
      printf("with negative number of elemens\n");
      exit(__LINE__,__FILE__);
    }

    if (size_dim0!=0) {
      printf("Error: initialization of allocated of array_2d object\n");
      exit(__LINE__,__FILE__);
    }

   data=new T[n0*n1];
   locally_allocated_data_buiffer=true;

   size_dim0=n0;
   size_dim1=n1;
  };

  int GetDim0() {return size_dim0;}
  int GetDim1() {return size_dim1;}
  int GetElementNumber() {return size_dim0*size_dim1;}
  T* GetBufferPointer() {return data;};

  void Deallocate() {
    if ((data!=NULL)&&(locally_allocated_data_buiffer==true)) delete [] data;

    data=NULL;
    size_dim0=0;
    size_dim1=0;
    locally_allocated_data_buiffer=false;
  }

//===================================================
  array_2d() {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    locally_allocated_data_buiffer=false;
  };

 ~array_2d() { 
   Deallocate();
  };

  array_2d(int n0,int n1) {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    locally_allocated_data_buiffer=false;

    init(n0,n1);
  };

  //the object does not allocate its own data buffer but acts as a manager providing indexed access to anothe data buffer
  array_2d(T* t,int n0,int n1) {
    locally_allocated_data_buiffer=false;

    data=t;
    size_dim0=n0;
    size_dim1=n1;
  }

  bool IsAllocated() {return (data!=NULL);}

//===================================================
  inline T operator () (int i0,int i1) const {
    if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)) exit(__LINE__,__FILE__,"Error: out of range");

    return data[i0*size_dim1+i1];
  };

  inline T* operator () (int i0) const {
    if ((i0<0)||(i0>=size_dim0)) exit(__LINE__,__FILE__,"Error: out of range");

    return data+i0*size_dim1;
  };
//===================================================
  inline T & operator () (int i0,int i1) {
    if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)) exit(__LINE__,__FILE__,"Error: out of range");

    return data[i0*size_dim1+i1];
  };

//===================================================
  inline array_2d<T>& operator = (const array_2d<T>& v) {
    int i,imax;
  
    imax=size_dim0*size_dim1;
    for (i=0;i<imax;i++) data[i]=v.data[i];
 
    return *this;
  };

//===================================================
  array_2d<T>& operator = (T f) {
    int i,imax;

    imax=size_dim0*size_dim1;
    for (i=0;i<imax;i++) data[i]=f;
 
    return *this;
  };

//===================================================
  friend array_2d<T> operator + (const array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      printf("Error: add two array_2d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]+v2.data[i]; 
   
    return v3;
  };

//===================================================
  friend array_2d<T> operator - (const array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      printf("Error: add two array_2d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]-v2.data[i];

    return v3;
  };

//===================================================
  friend array_2d<T> operator * (const array_2d<T> &v1, const T t) {
    int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=t*v1.data[i];

    return v3;
  };

//===================================================
  friend array_2d<T> operator / (const array_2d<T> &v1, const T t) {
    if (t == 0) {
      printf("Error: divide vector by 0.\n");
      exit(__LINE__,__FILE__);
    }
    int i,imax;
    array_2d<T> v3(v1.size_dim0,v1.size_dim1);

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]/t;

    return v3;
  };

//===================================================
  friend array_2d<T>& operator += (array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) { 
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]+=v2.data[i];

    return v1;
  };

//===================================================
  friend array_2d<T>& operator -= (array_2d<T> &v1,const array_2d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)) {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]-=v2.data[i];

    return v1;
  };

//===================================================
  friend array_2d<T>& operator *= (array_2d<T> &v1,const T t) {
    int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]*=t;

    return v1;
  };

//===================================================
  friend array_2d<T>& operator /= (array_2d<T> &v1,const T t) {
    if (t == 0) {
      printf("Error: divide array_2d<T> by 0.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.size_dim0*v1.size_dim1;
    for (i=0;i<imax;i++) v1.data[i]/=t;

    return v1;
  };

  //===================================================
    //get pointer to an element of the array
    T* GetPtr(int i0,int i1) {
      if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)) {
        printf("i0=%i, i1=%i; size_dim0=%i, size_dim1=%i\n",i0,i1,size_dim0,size_dim1);

        exit(__LINE__,__FILE__,"Error: out of range");
      }

      return data+i0*size_dim1+i1;
    }

    T* GetPtr() {
      return data;
    }

    int size() {
      return size_dim0*size_dim1;
    }

   void find_nan() {
    int i,imax=size_dim0*size_dim1;

    for (i=0;i<imax;i++) {
      if (isfinite(data[i])==false) {
        exit(__LINE__,__FILE__,"Error: a non-finite number is found");
      }
    }
  }

  void gather_double(int gather_rank,MPI_Comm &comm) { 
    double *temp;
    int rank;

    MPI_Comm_rank(comm,&rank);
    if (rank==gather_rank) temp=new double[size()];

    MPI_Gather(data,size(),MPI_DOUBLE,temp,size(),MPI_DOUBLE,gather_rank,comm); 

    if (rank==gather_rank) {
      memcpy(data,temp,size()*sizeof(double));
      delete [] temp;
    }
  }
		  

};

#endif
