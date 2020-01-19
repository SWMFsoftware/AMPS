//functionality used for constructing of stencils

/*
 * stencil.h
 *
 *  Created on: Dec 9, 2017
 *      Author: vtenishe
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>


#ifndef _GENERAL_STENCIL_H_
#define _GENERAL_STENCIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <list>

using namespace std;


class cFrac {
public:
  int nominator,denominator;
  
  void Simplify() {
    int i=2;
    
    while (i<=denominator) {
      if ((nominator%i==0)&&(denominator%i==0)) {
        nominator/=i;
        denominator/=i;
      }
      else i++;
    }
  }
  
  void Print(FILE* fout=stdout) {
    if (nominator%denominator==0) {
      fprintf(fout," %i ",nominator/denominator);
    }
    else{
      fprintf(fout," %i/%i ",nominator,denominator);
    }
  }
  
  cFrac(int n,int d) {nominator=n,denominator=d;}
  cFrac() {nominator=0,denominator=0;}
  
  void Set(int n,int d) {nominator=n,denominator=d;}

  int Convert2Int() {
    if (nominator%denominator!=0) {
      printf("Error: cannot be converted (%s@%i\n",__FILE__,__LINE__);
      exit(0);
    }

    return nominator/denominator;
  }
  
  cFrac& operator = (const cFrac& v) {
    nominator=v.nominator;
    denominator=v.denominator;
    
    return *this;
  }

  cFrac& operator = (const int& v) {
    nominator=v;
    denominator=1;

    return *this;
  }
  
  friend cFrac operator + (const cFrac &v1,const cFrac &v2) {
    cFrac v3;
    
    v3.nominator=v1.nominator*v2.denominator+v2.nominator*v1.denominator;
    v3.denominator=v1.denominator*v2.denominator;
    
    v3.Simplify();
    
    return v3;
  }
  
  friend cFrac& operator += (cFrac &v1,const cFrac &v2) {
    cFrac v3;
    
    v3.nominator=v1.nominator*v2.denominator+v2.nominator*v1.denominator;
    v3.denominator=v1.denominator*v2.denominator;
    
    v1=v3;
    v1.Simplify();
    
    return v1;
  }
  
  friend cFrac& operator -= (cFrac &v1,const cFrac &v2) {
    cFrac v3;
    
    v3.nominator=v1.nominator*v2.denominator-v2.nominator*v1.denominator;
    v3.denominator=v1.denominator*v2.denominator;
    
    v1=v3;
    v1.Simplify();
    
    return v1;
  }
  
  friend cFrac operator - (const cFrac &v1,const cFrac &v2) {
    cFrac v3;
    
    v3.nominator=v1.nominator*v2.denominator-v2.nominator*v1.denominator;
    v3.denominator=v1.denominator*v2.denominator;
    
    v3.Simplify();
    
    return v3;
  }
  
  bool operator == (const cFrac& rhs) {
    return (nominator*rhs.denominator==denominator*rhs.nominator);
  }
  
  friend bool operator < (const cFrac& l,const cFrac& r) {
    return (l.nominator*r.denominator<l.denominator*r.nominator);
  }
  
  friend bool operator > (const cFrac& l,const cFrac& r) {
    return (l.nominator*r.denominator>l.denominator*r.nominator);
  }
};


class cStencil {
public:
  cFrac i,j,k;
  char symbol[100];
  
  struct cStencilElement {
    cFrac i,j,k;
    double a; //coefficient used in the stencil
  };

  class cStencilData {
  public:

    class cElementData {
    public:
      int i,j,k;
      double a;

      cElementData() {i=0,j=0,k=0,a=0.0;}
    };

    cElementData *Data;
    int Length;

    void print(const char* msg=NULL) {
      if (msg!=NULL) printf("%s:\n",msg); 
      printf("Stencil length=%i\ni\tj\tk\ta\n",Length);

      for (int it=0;it<Length;it++) {
        printf("%i\t%i\t%i\t%e\n",Data[it].i,Data[it].j,Data[it].k,Data[it].a); 
      }
    }

    cStencilData() {
      Length=0,Data=NULL;
    }

    void remove() {
      if (Data!=NULL) {
        delete [] Data;
        Length=0,Data=NULL;
      }
    }

    ~cStencilData() {
      remove();
    }

    double GetElementValue(int i,int j,int k,bool *SucessFlag=NULL) {
      for (int it=0;it<Length;it++) if ((i==Data[it].i)&&(j==Data[it].j)&&(k==Data[it].k)) {
        if (SucessFlag!=NULL) *SucessFlag=true;
        return Data[it].a;
      }

      if (SucessFlag!=NULL) *SucessFlag=false;

      return 0.0;
    }
 
  };

  list<cStencilElement> StencilData;
  cStencilElement *Stencil;
  int StencilLength,AllocatedStencilLength;

  cStencil() {
    Stencil=NULL;
    StencilLength=0,AllocatedStencilLength=0;
    
    i=0,j=0,k=0;
    sprintf(symbol,"");
  }

  cStencil(const char* s) {
    Stencil=NULL;
    StencilLength=0,AllocatedStencilLength=0;
    
    i=0.0;
    j=0.0;
    k=0.0;

    sprintf(symbol,"%s",s);
  }
  
  cStencil(const char* s,cFrac iIn,cFrac jIn,cFrac kIn) {
    Stencil=NULL;
    StencilLength=0,AllocatedStencilLength=0;
    
    i=iIn,j=jIn,k=kIn;
    sprintf(symbol,"%s",s);
  }
  
  ~cStencil() {
    if (Stencil!=NULL) delete [] Stencil;
    StencilData.clear();
  }

  void SetSymbol(const char* s) {
    sprintf(symbol,"%s",s);
  }
  
  void MoveBase(cFrac iIn,cFrac jIn,cFrac kIn) {
    cFrac di,dj,dk;
    
    di=iIn-i;
    dj=jIn-j;
    dk=kIn-k;
    
    for (list<cStencilElement>::iterator it=StencilData.begin();it!=StencilData.end();it++) {
      it->i-=di;
      it->j-=dj;
      it->k-=dk;
    }
    
    i=iIn,j=jIn,k=kIn;
  }

  void SetBase(cFrac iIn,cFrac jIn,cFrac kIn) {
    i=iIn,j=jIn,k=kIn;
  }
  
  void PrintBase() {
    i.Print();
    j.Print();
    k.Print();
    printf("\n");
  }
  
  
  static bool SortList(const cStencilElement& first, cStencilElement& second) {
    if (first.i<second.i) return true;
    else if (first.i>second.i) return false;
    
    if (first.j<second.j) return true;
    else if (first.j>second.j) return false;
    
    if (first.k<second.k) return true;
    else if (first.k>second.k) return false;
    
    return true;
  }
  
  void Print(const char *fname=NULL) {
    cFrac t;
    int cnt=0;
    FILE* fout=NULL;

    fout=(fname!=NULL) ? fopen(fname,"w") : stdout;

    Simplify();
    
    //sort the list
    StencilData.sort(SortList);
    
    //print the list
    for (list<cStencilElement>::iterator it=StencilData.begin();it!=StencilData.end();it++) {
      fprintf(fout,"%i ", cnt);
      cnt++;
      
      t=i+it->i;
      t.Print(fout);
      
      t=j+it->j;
      t.Print(fout);
      
      t=k+it->k;
      t.Print(fout);
      
      fprintf(fout," %e\n",it->a);
    }

    if (fname!=NULL) fclose(fout);
  }
  
  //remove elements of the list that have the same combination of i,j, and k
  void Simplify() {
    cFrac i,k,j;
    list<cStencilElement>::iterator l0,l1;
    
    for (l0=StencilData.begin();l0!=StencilData.end();l0++) {
      list<list<cStencilElement>::iterator> CleaningList;
      
      i=l0->i,j=l0->j,k=l0->k;
      l1=l0;
      
      for (l1++;l1!=StencilData.end();l1++) {
        if ((i==l1->i)&&(j==l1->j)&&(k==l1->k)) {
          //found another element in the list that point to the point with the same combination of i,j, and k
          //combine l0 and l1
          l0->a+=l1->a;
          
          CleaningList.push_back(l1);
        }
      }
      
      for (list<list<cStencilElement>::iterator>::iterator t=CleaningList.begin();t!=CleaningList.end();t++) StencilData.erase(*t);
      
      CleaningList.clear();
    }
    
    //remove those elements from the list that have zero coefficients
    l0=StencilData.begin();
    
    while (l0!=StencilData.end()) {
      l1=l0;
      l1++;
      
      if (l0->a==0.0) {
        StencilData.erase(l0);
      }
      
      l0=l1;
    }
  }

  //in order to speed up usage of the stencil the data are saved in an aray rather than probided as a list
  cStencilElement* GetStencil() {
    Simplify();
    StencilLength=StencilData.size();

    //allocate the data buffer
    if (Stencil!=NULL) {
      if (AllocatedStencilLength<=StencilLength) {
        delete [] Stencil;
        Stencil=NULL;
      }
    }

    if (Stencil==NULL) {
      AllocatedStencilLength=(int)(1.2*StencilLength);
      Stencil=new cStencilElement[AllocatedStencilLength];
    }

    //copy the contant of the list into the array
    int i;
    list<cStencilElement>::iterator l0;

    for (i=0,l0=StencilData.begin();l0!=StencilData.end();i++,l0++) Stencil[i]=*l0;

    return this->Stencil;
  }

  void ExportStencil(cStencilData* s) {
    int it;
    list<cStencilElement>::iterator l0;

    Simplify();
    StencilLength=StencilData.size();

    s->remove();
    s->Length=StencilLength;
    s->Data=new cStencilData::cElementData[StencilLength];

    for (it=0,l0=StencilData.begin();l0!=StencilData.end();it++,l0++) { 
      s->Data[it].a=l0->a;

      s->Data[it].i=i.Convert2Int()+l0->i.Convert2Int();
      s->Data[it].j=j.Convert2Int()+l0->j.Convert2Int();
      s->Data[it].k=k.Convert2Int()+l0->k.Convert2Int();
    }
  }

  //multiply the entire stencil by a constatnt
  friend cStencil& operator *= (cStencil &v1,const double t) {
    for (list<cStencilElement>::iterator l0=v1.StencilData.begin();l0!=v1.StencilData.end();l0++) l0->a*=t;
    return v1;
  };

  //add a new element to the stencil
  void add(double a,cFrac i,cFrac j,cFrac k) {
    cStencilElement NewElement;

    NewElement.a=a;
    NewElement.i=i;
    NewElement.j=j;
    NewElement.k=k;

    StencilData.push_back(NewElement);
    Simplify();
  }

  void add(double a,int i,int j,int k) {
    cFrac iIn,jIn,kIn;

    iIn=i,jIn=j;kIn=k;
    add(a,iIn,jIn,kIn);
  }

  //shift the entire stencil
  void shift(cFrac di,cFrac dj,cFrac dk) {
    i+=di;
    j+=dj;
    k+=dk;
  }

  void shift(int di,int dj,int dk) {
   cFrac diIn,djIn,dkIn; 

   diIn=di; 
   djIn=dj;
   dkIn=dk;

   shift(diIn,djIn,dkIn);
 }  

  //copy the stencil
  cStencil& operator = (const cStencil& v) {
    i=v.i;
    j=v.j;
    k=v.k;
    
    for (list<cStencilElement>::const_iterator l0=v.StencilData.begin();l0!=v.StencilData.end();l0++) StencilData.push_back(*l0);
    return *this;
  };

  //add and substract another stencil
  friend cStencil& operator += (cStencil &v1,const cStencil &v2) {
    cStencilElement NewElement;

    for (list<cStencilElement>::const_iterator l0=v2.StencilData.begin();l0!=v2.StencilData.end();l0++) {
      NewElement=*l0;
      
      NewElement.i+=v2.i-v1.i;
      NewElement.j+=v2.i-v1.j;
      NewElement.k+=v2.i-v1.k;
  
      v1.StencilData.push_back(NewElement);
    }

    return v1;
  };

  void SwitchAxes(int Axis0,int Axis1) {  
    cFrac t;

    for (list<cStencilElement>::iterator l0=StencilData.begin();l0!=StencilData.end();l0++) {
      switch (3*Axis0+Axis1) {
      case 3*0+1:
        t=l0->i;
        l0->i=l0->j;
        l0->j=t;
        break;

      case 3*0+2:
        t=l0->i;
        l0->i=l0->k;
        l0->k=t;

        break;
      case 3*1+0:
        SwitchAxes(0,1);

        break;
      case 3*1+2:
        t=l0->j;
        l0->j=l0->k;
        l0->k=t;

        break; 
      case 3*2+0:
        SwitchAxes(0,2);

        break;
      case 3*2+1:
        SwitchAxes(1,2);
        break;
      }
    }
  }  


  void AddShifted(cStencil& v,cFrac di,cFrac dj,cFrac dk,double c=1.0) {
    cStencilElement NewElement;

    for (std::list<cStencilElement>::iterator l0=v.StencilData.begin();l0!=v.StencilData.end();l0++) {
      NewElement=*l0;

      NewElement.a*=c;
      
      NewElement.i+=v.i-i+di;
      NewElement.j+=v.j-j+dj;
      NewElement.k+=v.k-k+dk;

      StencilData.push_back(NewElement);
    }
  }

  void AddShifted(cStencil& v,int di,int dj,int dk,double c=1.0) {
    cFrac diIn,djIn,dkIn;

    diIn=di;
    djIn=dj;
    dkIn=dk; 

    AddShifted(v,diIn,djIn,dkIn,c);
  } 

  friend cStencil operator * (cStencil &v1,double t) {
    cStencil v3;

    v3=v1;

    for (std::list<cStencilElement>::iterator it=v3.StencilData.begin();it!=v3.StencilData.end();it++) {
      it->a*=t;
    } 

    return v3;
  }

  friend cStencil operator * (double t,cStencil &v1) {
    return v1*t;
  }

  friend cStencil operator + (cStencil &v1,cStencil &v2) {
    cStencil v3;

    v3=v1;
    v3.AddShifted(v2,0,0,0,1);

    return v3;
  }

  friend cStencil operator - (cStencil &v1,cStencil &v2) {
    cStencil v3;

    v3=v1;
    v3.AddShifted(v2,0,0,0,-1);

    return v3;
  }

  void SubstractShifted(cStencil& v,cFrac di,cFrac dj,cFrac dk,double c=1.0) {
    AddShifted(v,di,dj,dk,-c); 
  }

  void SubstractShifted(cStencil& v,int di,int dj,int dk,double c=1.0) {
    AddShifted(v,di,dj,dk,-c);
  }
};




#endif /* _GENERAL_STENCIL_H_ */
