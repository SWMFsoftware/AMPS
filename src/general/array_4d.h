//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#ifndef ARRAY_4D
#define ARRAY_4D

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <new>
#include <stdint.h>
#include "specfunc.h"

template <class T>
class array_4d {
protected:
  T* data;
  int size_dim0,size_dim1,size_dim2,size_dim3;
  int ndim0_ndim1,ndim0_ndim1_ndim2;
  int element_count;
  bool locally_allocated_data_buffer;

  // Array extents and MPI element counts are represented by int throughout
  // this legacy interface.  Validate in size_t before evaluating any signed
  // products, so an oversized shape cannot wrap to an undersized allocation.
  static int CheckedElementCount(int n0,int n1,int n2,int n3) {
    if (n0<=0 || n1<=0 || n2<=0 || n3<=0)
      exit(__LINE__,__FILE__,"Error: nonpositive array_4d extent");

    size_t count=1;
    const int extent[4]={n0,n1,n2,n3};
    for (int i=0;i<4;i++) {
      const size_t next=static_cast<size_t>(extent[i]);
      if (count>static_cast<size_t>(INT_MAX)/next)
        exit(__LINE__,__FILE__,"Error: array_4d extent product overflows");
      count*=next;
    }
    if (count>SIZE_MAX/sizeof(T))
      exit(__LINE__,__FILE__,"Error: array_4d allocation size overflows");

    return static_cast<int>(count);
  }

public:

  array_4d() {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    size_dim2=0;
    size_dim3=0;
    ndim0_ndim1=ndim0_ndim1_ndim2=0;
    element_count=0;
    locally_allocated_data_buffer=false;
  };

//===================================================
 ~array_4d() { Deallocate(); };

  void Deallocate() {
    if (data!=NULL && locally_allocated_data_buffer) delete [] data;
    data=NULL;
    size_dim0=size_dim1=size_dim2=size_dim3=0;
    ndim0_ndim1=ndim0_ndim1_ndim2=0;
    element_count=0;
    locally_allocated_data_buffer=false;
  }



//===================================================
  void init(int n0,int n1,int n2,int n3) {
    const int count=CheckedElementCount(n0,n1,n2,n3);

    if (data!=NULL || size_dim0!=0 || size_dim1!=0 || size_dim2!=0 || size_dim3!=0)
      exit(__LINE__,__FILE__,"Error: initialization of allocated array_4d object");

    try {
      data=new T[count];
    }
    catch (const std::bad_alloc&) {
      printf("Memory Error: array_4d() cannot allocate %zu bytes\n",
          static_cast<size_t>(count)*sizeof(T));
      exit(__LINE__,__FILE__);
    }


    size_dim0=n0;
    size_dim1=n1;
    size_dim2=n2;
    size_dim3=n3;

    ndim0_ndim1=n0*n1;
    ndim0_ndim1_ndim2=n0*n1*n2;
    element_count=count;
    locally_allocated_data_buffer=true;
  };
  
  array_4d(int n0,int n1,int n2,int n3) {
    data=NULL;
    size_dim0=0;
    size_dim1=0;
    size_dim2=0;
    size_dim3=0;
    ndim0_ndim1=ndim0_ndim1_ndim2=0;
    element_count=0;
    locally_allocated_data_buffer=false;

    init(n0,n1,n2,n3);
  }

  array_4d(T* external,int n0,int n1,int n2,int n3) {
    if (external==NULL)
      exit(__LINE__,__FILE__,"Error: invalid external array_4d buffer");
    const int count=CheckedElementCount(n0,n1,n2,n3);
    data=external;
    size_dim0=n0,size_dim1=n1,size_dim2=n2,size_dim3=n3;
    ndim0_ndim1=n0*n1;
    ndim0_ndim1_ndim2=n0*n1*n2;
    element_count=count;
    locally_allocated_data_buffer=false;
  }

  array_4d(const array_4d<T>& other) {
    data=NULL;
    size_dim0=size_dim1=size_dim2=size_dim3=0;
    ndim0_ndim1=ndim0_ndim1_ndim2=0;
    element_count=0;
    locally_allocated_data_buffer=false;
    if (other.size()>0) {
      init(other.size_dim0,other.size_dim1,other.size_dim2,other.size_dim3);
      for (int i=0;i<size();i++) data[i]=other.data[i];
    }
  }

  array_4d(array_4d<T>&& other) noexcept {
    data=other.data;
    size_dim0=other.size_dim0,size_dim1=other.size_dim1;
    size_dim2=other.size_dim2,size_dim3=other.size_dim3;
    ndim0_ndim1=other.ndim0_ndim1;
    ndim0_ndim1_ndim2=other.ndim0_ndim1_ndim2;
    element_count=other.element_count;
    locally_allocated_data_buffer=other.locally_allocated_data_buffer;
    other.data=NULL;
    other.size_dim0=other.size_dim1=other.size_dim2=other.size_dim3=0;
    other.ndim0_ndim1=other.ndim0_ndim1_ndim2=0;
    other.element_count=0;
    other.locally_allocated_data_buffer=false;
  }

  void resize(int n0,int n1,int n2,int n3) {
    if (n0==0 && n1==0 && n2==0 && n3==0) {
      Deallocate();
      return;
    }
    (void)CheckedElementCount(n0,n1,n2,n3);
    if (n0==size_dim0 && n1==size_dim1 && n2==size_dim2 && n3==size_dim3 &&
        locally_allocated_data_buffer) return;
    Deallocate();
    init(n0,n1,n2,n3);
  }

  void find_nan() {
    int i,imax=element_count;

    for (i=0;i<imax;i++) {
      if (isfinite(data[i])==false) {
        exit(__LINE__,__FILE__,"Error: a non-finite number is found");
      }
    }
  }

//===================================================
  T   operator () (int i0,int i1,int i2,int i3) const {
    if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)||(i2<0)||(i2>=size_dim2)||(i3<0)||(i3>=size_dim3)) exit(__LINE__,__FILE__,"Error: out of range");
    return data[i0+size_dim0*i1+ndim0_ndim1*i2+ndim0_ndim1_ndim2*i3]; 
  };

//===================================================
  T & operator () (int i0,int i1,int i2,int i3) {
    if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)||(i2<0)||(i2>=size_dim2)||(i3<0)||(i3>=size_dim3)) exit(__LINE__,__FILE__,"Error: out of range");
    return data[i0+size_dim0*i1+ndim0_ndim1*i2+ndim0_ndim1_ndim2*i3]; 
  };

//===================================================
  array_4d<T>& operator = (const array_4d<T>& v) {
    if (this==&v) return *this;
    if (size_dim0!=v.size_dim0 || size_dim1!=v.size_dim1 ||
        size_dim2!=v.size_dim2 || size_dim3!=v.size_dim3 ||
        !locally_allocated_data_buffer) {
      Deallocate();
      if (v.size()>0) init(v.size_dim0,v.size_dim1,v.size_dim2,v.size_dim3);
    }
    int i,imax=v.size();
    for (i=0;i<imax;i++) data[i]=v.data[i];
 
    return *this;
  };

  array_4d<T>& operator = (array_4d<T>&& v) noexcept {
    if (this==&v) return *this;
    Deallocate();
    data=v.data;
    size_dim0=v.size_dim0,size_dim1=v.size_dim1;
    size_dim2=v.size_dim2,size_dim3=v.size_dim3;
    ndim0_ndim1=v.ndim0_ndim1;
    ndim0_ndim1_ndim2=v.ndim0_ndim1_ndim2;
    element_count=v.element_count;
    locally_allocated_data_buffer=v.locally_allocated_data_buffer;
    v.data=NULL;
    v.size_dim0=v.size_dim1=v.size_dim2=v.size_dim3=0;
    v.ndim0_ndim1=v.ndim0_ndim1_ndim2=0;
    v.element_count=0;
    v.locally_allocated_data_buffer=false;
    return *this;
  }

//===================================================
  array_4d<T>& operator = (T f) {
    int i,imax;

    imax=element_count;
    for (i=0;i<imax;i++) data[i]=f;
 
    return *this;
  };

//===================================================
  friend array_4d<T> operator + (const array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two array_4d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.element_count;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]+v2.data[i]; 
   
    return v3;
  };

//===================================================
  friend array_4d<T> operator - (const array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two array_4d<T> of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.element_count;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]-v2.data[i];

    return v3;
  };

//===================================================
  friend array_4d<T> operator * (const array_4d<T> &v1, const T t) {
    int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.element_count;
    for (i=0;i<imax;i++) v3.data[i]=t*v1.data[i];

    return v3;
  };

//===================================================
  friend array_4d<T> operator / (const array_4d<T> &v1, const T t) {
    if (t == 0) {
      printf("Error: divide vector by 0.\n");
      exit(__LINE__,__FILE__);
    }
    int i,imax;
    array_4d<T> v3(v1.size_dim0,v1.size_dim1,v1.size_dim2,v1.size_dim3);

    imax=v1.element_count;
    for (i=0;i<imax;i++) v3.data[i]=v1.data[i]/t;

    return v3;
  };

//===================================================
  friend array_4d<T>& operator += (array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2) 
    ||(v1.size_dim3!=v2.size_dim3)) {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.element_count;
    for (i=0;i<imax;i++) v1.data[i]+=v2.data[i];

    return v1;
  };

//===================================================
  friend array_4d<T>& operator -= (array_4d<T> &v1,const array_4d<T> &v2) {
    if ((v1.size_dim0!=v2.size_dim0)||
    (v1.size_dim1!=v2.size_dim1)||(v1.size_dim2!=v2.size_dim2)
    ||(v1.size_dim3!=v2.size_dim3))  {
      printf("Error: add two vectors of different length.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.element_count;
    for (i=0;i<imax;i++) v1.data[i]-=v2.data[i];

    return v1;
  };

//===================================================
  friend array_4d<T>& operator *= (array_4d<T> &v1,const T t) {
    int i,imax;

    imax=v1.element_count;
    for (i=0;i<imax;i++) v1.data[i]*=t;

    return v1;
  };

//===================================================
  friend array_4d<T>& operator /= (array_4d<T> &v1,const T t) {
    if (t == 0) {
      printf("Error: divide array_4d<T> by 0.\n");
      exit(__LINE__,__FILE__);
    }

    int i,imax;

    imax=v1.element_count;
    for (i=0;i<imax;i++) v1.data[i]/=t;

    return v1;
  };

  //===================================================
  //get pointer to an element of the array
  T* GetPtr(int i0,int i1,int i2,int i3) {
    if ((i0<0)||(i0>=size_dim0)||(i1<0)||(i1>=size_dim1)||(i2<0)||(i2>=size_dim2)||(i3<0)||(i3>=size_dim3)) exit(__LINE__,__FILE__,"Error: out of range");
    return data+i0+size_dim0*i1+ndim0_ndim1*i2+ndim0_ndim1_ndim2*i3;
  }

  T* GetPtr() {
    return data;
  }

  int size() const {
    return element_count;
  }

  int size(int dimension) const {
    switch (dimension) {
    case 0: return size_dim0;
    case 1: return size_dim1;
    case 2: return size_dim2;
    case 3: return size_dim3;
    default: return 0;
    }
  }

  void reduce(int reduce_rank, MPI_Op operation, MPI_Comm comm) {
    // Get current process rank
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Get total data size (same for all processes)
    int dataSize = size();
    
    // Verify data buffer exists
    if (data == nullptr) {
        exit(__LINE__, __FILE__, "Error: Null data buffer in reduce operation");
    }
    
    // Determine MPI datatype based on template parameter T
    MPI_Datatype mpiType;
    if (std::is_same<T, char>::value) mpiType = MPI_CHAR;
    else if (std::is_same<T, unsigned char>::value) mpiType = MPI_UNSIGNED_CHAR;
    else if (std::is_same<T, short>::value) mpiType = MPI_SHORT;
    else if (std::is_same<T, unsigned short>::value) mpiType = MPI_UNSIGNED_SHORT;
    else if (std::is_same<T, int>::value) mpiType = MPI_INT;
    else if (std::is_same<T, unsigned int>::value) mpiType = MPI_UNSIGNED;
    else if (std::is_same<T, long>::value) mpiType = MPI_LONG;
    else if (std::is_same<T, unsigned long>::value) mpiType = MPI_UNSIGNED_LONG;
    else if (std::is_same<T, float>::value) mpiType = MPI_FLOAT;
    else if (std::is_same<T, double>::value) mpiType = MPI_DOUBLE;
    else if (std::is_same<T, long double>::value) mpiType = MPI_LONG_DOUBLE;
    else {
        exit(__LINE__, __FILE__, "Error: Unsupported data type for MPI reduce");
    }
    
    // Allocate temporary buffer for reduction
    T* recvBuffer = nullptr;
    if (rank == reduce_rank) {
        try {
            recvBuffer = new T[dataSize]();  // Zero-initialize the buffer
        } catch (std::bad_alloc&) {
            exit(__LINE__, __FILE__, "Error: Failed to allocate receive buffer in reduce");
        }
    }

    // Perform reduction
    int error = MPI_Reduce(
        data,           // Send buffer
        recvBuffer,     // Receive buffer
        dataSize,       // Count
        mpiType,        // Datatype
        operation,      // Reduction operation (e.g., MPI_SUM, MPI_MAX)
        reduce_rank,    // Root rank
        comm           // Communicator
    );

    if (error != MPI_SUCCESS) {
        if (rank == reduce_rank && recvBuffer != nullptr) {
            delete[] recvBuffer;
        }
        exit(__LINE__, __FILE__, "Error: MPI_Reduce failed");
    }

    // On reduction rank, copy reduced data back to array
    if (rank == reduce_rank) {
        // Copy received data back to the array
        memcpy(data, recvBuffer, dataSize * sizeof(T));
        
        // Clean up receive buffer
        delete[] recvBuffer;
    }

    // Ensure all processes wait for reduction to complete
    MPI_Barrier(comm);
}

//===================================================
};

#endif
