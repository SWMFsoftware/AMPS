//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "Gravity.h"
#include "global.h"
#include "specfunc.h"

#include "ifileopr.h"

#define MAXSTR 1000

using namespace std;

nucleusGravity::cNASTRANnode * nucleusGravity::nodes;
nucleusGravity::cNASTRANtetra * nucleusGravity::tetras;

long int nucleusGravity::nnodes,nucleusGravity::ntetras; ///,ncells;

double nucleusGravity::density;


void nucleusGravity::setDensity(double d) {
  density=d;
}

void nucleusGravity::readMesh_longformat(const char *fname,const char *path) {
  char fullname[500];

  sprintf(fullname,"%s/%s",path,fname);
  readMesh_longformat(fullname);
}

void nucleusGravity::readMesh_longformat(const char *fname) {
  CiFileOperations ifile;
  char str[MAXSTR],dat[MAXSTR],*endptr;
  long int i,j,idim;
  
  cNASTRANnode node;
  cNASTRANtetra tetra;

  //check whether the binary file with the mesh is available 
  FILE *fData=NULL;

  sprintf(str,"%s.bin",fname);
  fData=fopen(str,"r");
 
  if (fData!=NULL) {
    //the binary file exists -> use it directly
    fread(&nnodes,sizeof(long int),1,fData); 
    fread(&ntetras,sizeof(long int),1,fData);

    nodes = new cNASTRANnode [nnodes];
    tetras = new cNASTRANtetra [ntetras];

    fread(nodes,sizeof(cNASTRANnode),nnodes,fData); 
    fread(tetras,sizeof(cNASTRANtetra),ntetras,fData);

    fclose(fData);
    return;
  }
    
      
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
    node.x[0]=strtod(dat,&endptr);
    
    ifile.CutInputStr(dat,str);
    node.x[1]=strtod(dat,&endptr);
    
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);
    
    ifile.CutInputStr(dat,str);
    node.x[2]=strtod(dat,&endptr);
    nnodes++;
    
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);
  }
  nodes = new cNASTRANnode [nnodes];

  rewind(ifile.openfile(fname));
      ifile.GetInputStr(str,sizeof(str));
      
      //read nodes from the file 
      
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(dat,str);


      nnodes=0;
      while (strcmp("GRID*",dat)==0) {
	ifile.CutInputStr(dat,str);
	node.id=strtol(dat,&endptr,10);
	
	ifile.CutInputStr(dat,str);
	ifile.CutInputStr(dat,str);
	node.x[0]=strtod(dat,&endptr);
	//node.x[0]=strtod(dat,&endptr)*5.8;

	ifile.CutInputStr(dat,str);
	node.x[1]=strtod(dat,&endptr);
	//node.x[1]=strtod(dat,&endptr)*5.8;
	
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
	
	ifile.CutInputStr(dat,str);
	node.x[2]=strtod(dat,&endptr);
	//node.x[2]=strtod(dat,&endptr)*5.8;

      //check if coordinates of the node are too close to the nodes of the mesh 
      double dX,dXmin;
      /*
      for (idim=0;idim<3;idim++) {
        dXmin=rootBlock_dX[idim]/(1<<maxLevel);
        dX=(node.x[idim]-rootBlock_Xmin[idim])/dXmin; 

/*
        dX=(dX-(long int)dX)*dXmin;

        if (dX<EPS) node.x[idim]+=2.0*EPS; 
        else if (dXmin-dX<EPS) node.x[idim]-=2.0*EPS; 


        dX-=(long int)dX;
        if (dX<0.1) node.x[idim]+=0.1*dXmin;
        else if (dX>0.9) node.x[idim]-=0.1*dXmin;


      }*/
      
      nodes[nnodes]=node;//nodes.push_back(node);
      nnodes++;

      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(dat,str);
    }
      
      while (strcmp("CBAR",dat)==0) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      } 
      while (strcmp("CTRIA3",dat)==0) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      } 
       
    //read the tetrahedron information
    while (strcmp("CTETRA",dat)==0) {
      //ifile.CutInputStr(dat,str);
      //ifile.CutInputStr(dat,str);
      
      ntetras++;
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(dat,str);

    }
    
      tetras = new cNASTRANtetra [ntetras];
      ntetras=0;

      rewind(ifile.openfile(fname));
      ifile.GetInputStr(str,sizeof(str));
      
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(dat,str);
            
            
      while (strcmp("GRID*",dat)==0) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);	
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      }
      while (strcmp("CBAR",dat)==0) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      } 
      while (strcmp("CTRIA3",dat)==0) {
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      }
      while (strcmp("CTETRA",dat)==0) {
	ifile.CutInputStr(dat,str);
	ifile.CutInputStr(dat,str);
	for (idim=0;idim<4;idim++) {
	  ifile.CutInputStr(dat,str);
	  tetra.node[idim]=strtol(dat,&endptr,10);
	}
	tetras[ntetras]=tetra;
	ntetras++;
	ifile.GetInputStr(str,sizeof(str));
	ifile.CutInputStr(dat,str);
      }
      ifile.closefile();

      register long int nd,nfc,id;
      
      //renumerate nodes
      //determine the maximum valus of node->id and allocate the map for relating node->id with the node number
      int idMax=0;
      int *idMap=NULL;

      for (nd=0;nd<nnodes;nd++) {
        if (nodes[nd].id>idMax) idMax=nodes[nd].id;

        if (nodes[nd].id<0) exit(__LINE__,__FILE__,"Error: The procedure works only when nodes[nd].id>=0; Solution: need to implement processing of nodes[nd].id<0");
      }


      idMap=new int[idMax+1];
      for (nd=0;nd<nnodes;nd++) idMap[nodes[nd].id]=nd;


      for (nfc=0;nfc<ntetras;nfc++) for (idim=0;idim<4;idim++) {
	      id=tetras[nfc].node[idim];
	      tetras[nfc].node[idim]=idMap[id];
	    }

      for (nd=0;nd<nnodes;nd++) nodes[nd].id=nd;

      delete [] idMap;


      //save the mesh binary file
      int rank;

      MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);

      if (rank==0) {
        sprintf(str,"%s.bin",fname);
        fData=fopen(str,"w"); 

        fwrite(&nnodes,sizeof(long int),1,fData);
        fwrite(&ntetras,sizeof(long int),1,fData);

        fwrite(nodes,sizeof(cNASTRANnode),nnodes,fData);
        fwrite(tetras,sizeof(cNASTRANtetra),ntetras,fData);

        fclose(fData);
      }
}

void nucleusGravity::readMesh(const char *fname) {
    FILE *fd;
    long int initPosition,position,nd,i,idim,nfc;
    char str[MAXSTR],dat[MAXSTR],*endptr;

    //open the file 
    errno=0;
    fd=fopen(fname,"r");

    if (errno!=0) {
      printf("cannot open the mesh file %s\n",fname);
      exit(0);
    }

    //cout the number of nodes
    bool foundflag=false;

    while (!feof(fd)) {
      position=ftell(fd);
      fgets(str,MAXSTR,fd);

      if ((str[0]=='G')&&(str[3]=='D')) {
        if (foundflag==false) foundflag=true,initPosition=position;

        ++nnodes; 
      }
      else {
        if (foundflag==true) break;
      }
    } 

    if (foundflag==false) {
      //cout << "ERROR: something wrong with the NASTRAN file -> too short" << endl;
      exit(0);
    }

    //load the nodes
    fseek(fd,initPosition,SEEK_SET);

    nodes = new cNASTRANnode [nnodes];
	  
    for (nd=0;nd<nnodes;nd++) {
      fgets(str,MAXSTR,fd);

      for (i=0;i<12;i++) dat[i]=str[i+4];
      dat[i]='\n';
      nodes[nd].id=strtol(dat,&endptr,10);

      for (idim=0;idim<3;idim++) {
        for (i=0;i<8;i++) dat[i]=str[i+24+8*idim];
        dat[i]='\n';

        //search for an exponential part in the node coordinte 
        for (i=1;i<8;i++) if ((dat[i]=='+')||(dat[i]=='-')) {
          bool foundflag=false;

          //search for symbols 'E' or 'e' in the line 
          for (i=1;i<8;i++) if ((dat[i]=='E')||(dat[i]=='e')) foundflag=true; 
          if (foundflag==true) break;

          //if 'E' or 'e'  is not found in the line -> insert them
          for (i=1;i<8;i++) if ((dat[i]=='+')||(dat[i]=='-')) {
            for (int j=9;j>i;j--) dat[j]=dat[j-1];
            dat[i]='E';
            break;
          }

          break;
        }

        for (i=1;i<8;i++) if (dat[i]==' ') {
          dat[i]='\n';
          break;
        }

        nodes[nd].x[idim]=strtod(dat,&endptr);
      }
    }


    //count the number of tetras 
    foundflag=false;

    while (!feof(fd)) {
      position=ftell(fd);
      fgets(str,MAXSTR,fd);

      if ((str[0]=='C')&&(str[3]=='T')) {
        if (foundflag==false) foundflag=true,initPosition=position;

        ++ntetras;
      }
      else {
        if (foundflag==true) break;
      }
    }

    if (foundflag==false) {
      //cout << "ERROR: something wrong with the NASTRAN file -> too short" << endl;
      exit(0);
    }

    //load the tetras 
    fseek(fd,initPosition,SEEK_SET);
    tetras = new cNASTRANtetra [ntetras];
    
    for (nfc=0;nfc<ntetras;nfc++) {
      fgets(str,MAXSTR,fd);

      for (i=0;i<8;i++) dat[i]=str[i+8];
      dat[i]='\n';
      /// cell.id=strtol(dat,&endptr,10);

      for (i=0;i<8;i++) dat[i]=str[i+16];
      dat[i]='\n';
      //tetras[nfc].faceat=strtol(dat,&endptr,10);

      for (i=0;i<8;i++) dat[i]=str[i+24];
      dat[i]='\n';

      for (idim=0;idim<4;idim++) {
        for (i=0;i<8;i++) dat[i]=str[i+24+8*idim];
        dat[i]='\n';
        tetras[nfc].node[idim]=strtol(dat,&endptr,10);
      }
    }

    //renumerate nodes
    register long int id;

    for (nfc=0;nfc<ntetras;nfc++) for (idim=0;idim<4;idim++) {
      id=tetras[nfc].node[idim];

      for (nd=0;nd<nnodes;nd++) if (nodes[nd].id==id) tetras[nfc].node[idim]=nd;
    }
  }


  
double  nucleusGravity::gravity(double * gravitypermass, double * position) {
      double center[3];
      double a[3],b[3],c[3];
      double forcepermass[3];
      double r_squared=0,r=0;
      double volumetetra=0;
      double volume=0;
      int nfc,idim;

      static double TotalMass=0.0;

      //the static tables used in the gravity calcualtion
      static double **xThetrahedronCenter=NULL;
      static double *ThetrahedronVolume=NULL;

      if (ThetrahedronVolume==NULL) {
        //allocate the static tables
        ThetrahedronVolume=new double [ntetras];
        xThetrahedronCenter=new double* [ntetras];
        xThetrahedronCenter[0]=new double [3*ntetras];

        for (int n=1;n<ntetras;n++) xThetrahedronCenter[n]=xThetrahedronCenter[n-1]+3; 

        //init the tables with the mass and center points of teh tetrahedrons
        for (nfc=0;nfc<ntetras;nfc++){
          r_squared=0;
          for (idim=0;idim<3;idim++) {
            xThetrahedronCenter[nfc][idim]=0;

            for (int n=0;n<4;n++) xThetrahedronCenter[nfc][idim]+=0.25*nodes[tetras[nfc].node[n]].x[idim];
          }

          for (idim=0;idim<3;idim++) {
            a[idim]=nodes[tetras[nfc].node[1]].x[idim]-nodes[tetras[nfc].node[0]].x[idim];
            b[idim]=nodes[tetras[nfc].node[2]].x[idim]-nodes[tetras[nfc].node[0]].x[idim];
            c[idim]=nodes[tetras[nfc].node[3]].x[idim]-nodes[tetras[nfc].node[0]].x[idim];
          }

          ThetrahedronVolume[nfc]=6.67259e-11*density*fabs(c[0]*(a[1]*b[2]-b[1]*a[2])+c[1]*(a[2]*b[0]-b[2]*a[0])+c[2]*(a[0]*b[1]-b[0]*a[1]))/6;
          TotalMass+=density*fabs(c[0]*(a[1]*b[2]-b[1]*a[2])+c[1]*(a[2]*b[0]-b[2]*a[0])+c[2]*(a[0]*b[1]-b[0]*a[1]))/6;
         
        }
      }
  

     //clean the acceleartion vector
     for (idim=0;idim<3;idim++) gravitypermass[idim]=0.0;

     for (nfc=0;nfc<ntetras;nfc++){
        double r_3,r2Center[3],r_squared=0.0;

	for (idim=0;idim<3;idim++) {
	  r_squared+=pow((r2Center[idim]=position[idim]-xThetrahedronCenter[nfc][idim]),2.0);
	}

        r_3=1.0/pow(r_squared,3.0/2.0); 

	for (idim=0;idim<3;idim++) {
	  gravitypermass[idim]-=ThetrahedronVolume[nfc]*r2Center[idim]*r_3;				       
	}
    
     }

     return TotalMass;
}
  
void nucleusGravity::printMeshFile() {
    FILE *mesh;
    mesh = fopen("Mesh.dat","w");
    fprintf(mesh,"ZONE N=%li, E=%li, F=FEPOINT, ET=TETRAHEDRON \n",nnodes,ntetras);
    for (int i=0;i<nnodes;i++) {
      fprintf(mesh,"%e %e %e \n",nodes[i].x[0],nodes[i].x[1],nodes[i].x[2]);
    }
    for (int i=0;i<ntetras;i++) {
      fprintf(mesh,"%li %li %li %li \n",tetras[i].node[0]+1,tetras[i].node[1]+1,tetras[i].node[2]+1,tetras[i].node[3]+1);
    }
    fclose(mesh);
  }
  
void nucleusGravity::clearArrays() {
    free(tetras);
    free(nodes);
  }
