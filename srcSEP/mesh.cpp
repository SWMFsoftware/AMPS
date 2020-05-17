/*
 * mesh.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */





#include "sep.h"


double** SEP::Mesh::FieldLineTable=NULL;
int SEP::Mesh::FieldLineTableLength=0;


double SEP::Mesh::localSphericalSurfaceResolution(double *x) {
  double res,r,l[3];
  int idim;
  double SubsolarAngle;

  return 0.1*_RADIUS_(_SUN_);
}

double SEP::Mesh::localResolution(double *x) {
  return 10.0*_RADIUS_(_SUN_);
}


void SEP::Mesh::LoadFieldLine(list<SEP::cFieldLine> *field_line) {
  int i;
  list<SEP::cFieldLine>::iterator it;
  FILE *fLine;

  if (FieldLineTableLength!=0) exit(__LINE__,__FILE__,"Error: the filed line is alreadu deined");

  FieldLineTableLength=field_line->size();

  FieldLineTable=new double *[FieldLineTableLength];
  FieldLineTable[0]=new double [3*FieldLineTableLength];

  for (int i=0;i<FieldLineTableLength;i++) FieldLineTable[i]=FieldLineTable[0]+3*i;

  if (PIC::ThisThread==0) {
    fLine=fopen("magnetic-line.dat","w");
    fprintf(fLine,"VARIABLES=\"x\",\"y\",\"z\"\nZONE T=\"MAgnetic Field Line\", I=%i\n",FieldLineTableLength);
  }

  for (i=0,it=field_line->begin();it!=field_line->end();it++,i++) {
    for (int idim=0;idim<3;idim++) FieldLineTable[i][idim]=it->x[idim];

    if (PIC::ThisThread==0) fprintf(fLine,"%e %e %e\n",FieldLineTable[i][0],FieldLineTable[i][1],FieldLineTable[i][2]);
  }

  if (PIC::ThisThread==0) fclose(fLine);
}

bool SEP::Mesh::NodeSplitCriterion(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double *xmin=startNode->xmin;
  double *xmax=startNode->xmax;
  int i;
  double d2;

  bool field_line_intersect_node=false;

  double d2min=pow(2.0*_RADIUS_(_SUN_),2);
  double cell_res2,cell_res; ////=0.1+(0.5-0.1)/(200-1.2)*(sqrt(pow(0.5*(xmin[0]+xmax[0]),2)+pow(0.5*(xmin[1]+xmax[1]),2)+pow(0.5*(xmin[2]+xmax[2]),2))/_RADIUS_(_SUN_)-1.2);

  double r=sqrt(pow(0.5*(xmin[0]+xmax[0]),2)+pow(0.5*(xmin[1]+xmax[1]),2)+pow(0.5*(xmin[2]+xmax[2]),2))/_RADIUS_(_SUN_);

  if (r<1.1) return false;

  cell_res=(r>0.0) ? 0.1*exp(log(0.5/0.1)/log(200.0/1.2)*log(r/1.2)) : 0.1;


  if (cell_res<0.1) cell_res=0.1;
  if (cell_res>0.5) cell_res=0.5;

  cell_res2=pow(cell_res*_RADIUS_(_SUN_),2);

  for (i=0;i<FieldLineTableLength;i++) {
    field_line_intersect_node=true;
    d2=0.0;

    for (int idim=0;idim<3;idim++) {
      if ((FieldLineTable[i][idim]<xmin[idim]) || (FieldLineTable[i][idim]>xmax[idim])) field_line_intersect_node=false;

      double t=0.5*(xmin[idim]+xmax[idim])-FieldLineTable[i][idim];

      d2+=t*t;
    }

    if ((d2<d2min)||(field_line_intersect_node==true)) {
      //the block is close to the field line
      double t0,t1,t2,l2;

      t0=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
      t1=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
      t2=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

      l2=t0*t0+t1*t1+t2*t2;

      if (l2>cell_res2) {
        //the block needs to be refined
        return true;
      }
    }
  }

  return false;
}


double SEP::Mesh::localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize,CharacteristicSpeed=400E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}
