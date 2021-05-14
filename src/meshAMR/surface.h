//the class allows to load a triabgulated surface and store it ourside of AMPS 
//class cNodeDataContainer is a container keping the data at the node locations. cNodeDataContainer should have methods Save() and Read()
//for saving/reading the data to/from a file.

#ifndef _SURFACE_TRIANGULATION_ 
#define _SURFACE_TRIANGULATION_ 

#include "meshAMRcutcell.h"

template <class cNodeDataContainer> 
class cSurfaceTriangulation {
  vector<CutCell::cNASTRANnode> Nodes;
  vector<CutCell::cNASTRANface> Faces;

  vector<cNodeDataContainer> NodeDataVector; 

public:
  CutCell::cNASTRANnode* GetNode(int iNode) {return &Nodes[iNode];}
  cNodeDataContainer* GetNodeData(int iNode) {return &NodeDataVector[iNode];}
  double* GetNodeX(int iNode) {return Nodes[iNode].x;}  

  void SaveData(const char* fname) {
    FILE *fout=fopen(fname,"w");


    for (auto& it : NodeDataVector) {
       it.Save(fout);
    }

    fclose(fout);
  } 

  void ReadData(const char* fname) {
    FILE *fout=fopen(fname,"r");


    for (auto& it : NodeDataVector) {
       it.Read(fout);
    }

    fclose(fout);
  }

  void ReadNastranSurfaceMeshLongFormat(const char* fname,double UnitConversitonFactor=1.0) {
   CiFileOperations ifile;
   char str[10000],dat[10000],*endptr;
   int idim,nnodes=0,nfaces=0,NodeNumberOffset=0,ThisMeshNodeCounter;
   CutCell::cNASTRANnode node;
   CutCell::cNASTRANface face;

   ifile.openfile(fname);
   ifile.GetInputStr(str,sizeof(str));

   //read nodes from the file
   ifile.GetInputStr(str,sizeof(str));
   ifile.CutInputStr(dat,str);

   //load the nodes
   while (strcmp("GRID*",dat)==0) {
     ifile.CutInputStr(dat,str);
     node.id=strtol(dat,&endptr,10);

     ifile.CutInputStr(dat,str);
     ifile.CutInputStr(dat,str);
     node.x[0]=UnitConversitonFactor*strtod(dat,&endptr);

     ifile.CutInputStr(dat,str);
     node.x[1]=UnitConversitonFactor*strtod(dat,&endptr);

     ifile.GetInputStr(str,sizeof(str));
     ifile.CutInputStr(dat,str);

     ifile.CutInputStr(dat,str);
     node.x[2]=UnitConversitonFactor*strtod(dat,&endptr);

     Nodes.push_back(node);

     ifile.GetInputStr(str,sizeof(str));
     ifile.CutInputStr(dat,str);
   }

   //find the beginig for the face information
   while (strcmp("CBAR",dat)==0) {
     ifile.GetInputStr(str,sizeof(str));
     ifile.CutInputStr(dat,str);
   }

   //read the face information
   while (strcmp("CTRIA3",dat)==0) {
     ifile.CutInputStr(dat,str);

     ifile.CutInputStr(dat,str);
     face.faceat=strtol(dat,&endptr,10);

     for (idim=0;idim<3;idim++) {
       ifile.CutInputStr(dat,str);
       face.node[idim]=strtol(dat,&endptr,10)-1;
     }

     Faces.push_back(face);

     ifile.GetInputStr(str,sizeof(str));
     ifile.CutInputStr(dat,str);
   }

   ifile.closefile();
   NodeDataVector.resize(Nodes.size());
  }
};


#endif
