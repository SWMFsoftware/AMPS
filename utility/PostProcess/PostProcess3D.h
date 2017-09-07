//$Id$
//the library for reading and post-processing of AMPS' output files

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <algorithm>
#include <boost/iostreams/device/mapped_file.hpp> // for mmap 

#include "meshAMRcutcell.h"
#include "specfunc.h"
#include "global.h"


#ifndef _POST_PROCESS_3D_
#define _POST_PROCESS_3D_

#define _POST_PROCESS_MAP_SIZE_ (1<<8)

class cPostProcess3D {
public:

  //methods for reading and reconstruction of the output file
  class cCell {
  public:
    int n[8];
    vector<int> TrajectoryPoints,IndividualTrajectories;
    int LastTrajectoryProcessed;

    cCell() {
      LastTrajectoryProcessed=-1;
    }
  };

  class cBlock {
  public:
    cCell ***cell;
    double xmin[3],xmax[3],dx[3];
    int nCellX,nCellY,nCellZ,id;

    void Init(int nx,int ny,int nz);
  };


  cCell *ConnectivityList;
  cBlock *Block;
  vector<cBlock*> Mesh[_POST_PROCESS_MAP_SIZE_][_POST_PROCESS_MAP_SIZE_][_POST_PROCESS_MAP_SIZE_]; //contains pointers to the blocks
  double xmin[3],xmax[3],dxBlockMin[3];
  int nNodes,nCells,nBlocks;
  int nBlockCellX,nBlockCellY,nBlockCellZ;

  //data structure
  class sDataStructure {
  public:
    int nNodes,nVariables;
    std::vector<std::string> VariableList;
    double** data;

    void CopyVariable(int iTarget,sDataStructure* Source,int iSource);
    void InitDataBuffer(int nnodes,int nvars);
  };

  sDataStructure X,data;

  //interpolation stencil
  struct cStencil {
    int Node[8];
    double Weight[8];
  };

  //surface data
  class cSurfaceData {
  public:
    cPostProcess3D *PostProcess3D;

    //data for the data file
    std::vector<std::string> VariableList;
    int nNodes,nCells,nVariables;
    double **data;
    int **ConnectivityList;
    std::vector<int>*  NodeBall;

    void LoadDataFile(const char* fname,const char* path);
    void PrintVariableList();
    void PrintBall(int node);

    cSurfaceData() {
      PostProcess3D=NULL;
      nNodes=0,nCells=0,nVariables=0;
      data=NULL,ConnectivityList=NULL,NodeBall=NULL;
    }
  };

  cSurfaceData SurfaceData;

  
  map<int, int> *CellTrajMap;  //<iTrajecotry,nTrajPnt>
  map<int, int> *TrajCellMap;  //<CellId, nTrajPnts>
 
  void ComputeDtIntegration(double & dtIntegration);
  // for CellTrajMap
  void ReadTrajCellData(int iBuffer);
  void GetCellTrajInfo();
  int  sumTrajPoints(int iCell);
  void AddOneCellTrajPointToTraj(std::map<int, int>  &TrajCellMap, double* x);



  bool initTrajCellMapFlag;
  //trajectories of the infividual particles
  class cParticleTrajectory {
  public:
    cPostProcess3D *PostProcess3D;

    //data for all trajectories
    std::vector<std::string> VariableList;
    std::vector<std::string> BufferNameArray;
    bool LastBufferFull;
    int AvailBufferSize;
   
    bool initBufferFileSizeFlag;
    int  BufferFileSize;
    std::vector<long int *> TrajOffsetInBuffer;
    
    int TrajectoryStartingFaceOffset;
    static int nTrajectoryVariables;
    static int ParticleWeightOverTimeStepOffset;
    int nTotalTrajectories;

    //Individual Particle Trajectory
    class cIndividualTrajectoryData {
    public:
      double **Data;
      int nDataPoints;

      cIndividualTrajectoryData() {
        Data=NULL,nDataPoints=0;
      };
      
      ~cIndividualTrajectoryData(){

	if (Data!=NULL){
	  delete [] Data[0];
	  delete [] Data;
	};
      };

      void AllocateDataBuffer(int n);
    };

    std::vector<cIndividualTrajectoryData> IndividualTrajectories;

    //Load and output trajectory files
    void LoadDataFile(const char* fname,const char* path);
    void AddIndividualTrajectoryData(int& nDataPoints,std::vector<double>& data);

    void PrintDataFileHeader(const char* fname);
    void AddTrajectoryDataFile(cIndividualTrajectoryData* Trajectory,int TrajectoryNumber,const char* fname);

    //plot the surface mesh showing averaged propserties of the trajectories
    void PrintSurfaceData(const char *fname,int (*GetVariableNumber)(),void (*PrintVariableList)(FILE*),void (*GetFaceDataVector)(double*,CutCell::cTriangleFace*,int));

    //print the variable list
    void PrintVariableList();
    
    //functions for using boost and buffer files
    void GetTrajectoryLocation_boost(int & traj_num, std::vector<const char *> & ZoneCharLoc,  boost::iostreams::mapped_file & mmap);
    void WriteBinaryBuffer_boost(const char * bufferName, int head_traj_num, int buffer_traj_size, const char ** ZoneCharLocArr, boost::iostreams::mapped_file & mmap);
    void AppendBinaryBuffer_boost(const char * bufferName, int head_traj_num, int add_buffer_traj_size, const char ** ZoneCharLocArr, boost::iostreams::mapped_file & mmap);
    void WriteIndividualTrajectory_boost(const char * ZoneBegin, const char * ZoneEnd,  boost::iostreams::mapped_file & mmap, FILE* fBinaryOut);
    void InitLoadBufferFile_boost(const char *fname,const char* path, int buffer_traj_num);
    void InitLoadBufferFile_boost(const char *fname,const char* path, int buffer_traj_num, int inputnTraj);
    void LoadBufferHeader(int iBuffer, FILE * & fBinaryIn); 
    void ReadSequentialTrajFromBuffer(int nTraj, FILE * & fBinaryIn);
    
      // functions for using buffer files
    void InitLoadBufferFile(const char *fname,const char* path, int buffer_traj_num, int nTrajInFile);
    void InitLoadBufferFile(const char *fname,const char* path, int buffer_traj_num);
    void ReadBinaryBuffer(const char * bufferName);
    void LoadBufferData(int iBuffer);
    
    void LoadAllBufferOneLine();
    void LoadBufferDataOneLine(int iBuffer, int head_traj_num);


    void ReadIndividualTrajFromBuffer(cIndividualTrajectoryData &traj, const char * bufferName, long int offset);
    void ReadIndividualTrajOneLineFromBuffer(cIndividualTrajectoryData &traj, const char * bufferName, long int offset);
    void LoadIndividualTrajData(cIndividualTrajectoryData &traj, int iTraj);
    void LoadIndividualTrajDataOneLine(cIndividualTrajectoryData &traj, int iTraj);
    void WriteIndividualTrajectory(int LineBegin, int LineEnd, std::ifstream &infile, int& cnt, FILE* fBinaryOut);
    void WriteBinaryBuffer(const char * bufferName, int head_traj_num, int buffer_traj_size, int* ZoneLineNum, std::ifstream &infile, int & currentLineNum);
    void GetIndividualTrajOffset(int iTraj, long int & offset);
    void GetTrajectoryLineNumber(int & traj_num, std::vector<int> & ZoneLineNum, const char* filename);
    void InitTrajectoryVarList(const char* filename);
    void AppendBinaryBuffer(const char * bufferName, int head_traj_num, int add_buffer_traj_size, int* ZoneLineNum, std::ifstream &infile, int & currentLineNum);

    //functions for TrajCellMap
    void ComputeOneTrajsCellInfo(int nTraj, double dtIntegration, std::map<int, int> & TrajCellMap);
    void WriteTrajsCellInfoBuffer(int iBuffer);    
    void WriteTrajsCellInfoTraj(int jTraj, FILE * fBinaryOut, FILE * fBinaryOutOffset);
    void WriteTrajsCellInfoTraj(int jTraj, FILE * fBinaryOut, FILE * fBinaryOutOffset, std::map<int,int> & TrajCellMap);
    void SaveAllTrajsCellInfo();
    void ReadIndividualTrajCellInfo(int nTraj);
 
    cParticleTrajectory() {
      nTotalTrajectories=0,TrajectoryStartingFaceOffset=-1;
      LastBufferFull = true;
    }

  };


  cParticleTrajectory ParticleTrajectory;

  //print particle trajectories into a file
  //nTrajectories -> the number of the trajectories that will be printed
  //OutputMode -> weiting of the trajectories
  #define _OUTPUT_MODE__UNIFORM_ 0
  #define _OUTPUT_MODE__FLUX_    1
  void PrintParticleTrajectory(int nTrajectories,int OutputMode,double (*TrajectoryAcceptableProbability)(int),const char* fname);
  void AssignParticleTrajectoriesToCells();
  void AssignParticleTrajectoriesToCellsOptimize();
  //column integration
  class cColumnIntegral {
  public:
    cPostProcess3D *PostProcess3D;

    //class that contains the user-defined integtated parametes
    struct cColumnIntegrationSet {
      int (*IntegrantVectorLength)();
      void (*IntegrantVector)(double* Data,double *Location);
      void (*PrintVariableList)(FILE* fout);
      void (*PostProcessColumnIntegralVector)(double*);
    };

    class cMap {
    public:
      cColumnIntegral *ColumnIntegral;

      //Types of the column maps integration
      void Circular(double *xObservation,double *PrimaryAxisDirection,double *SecondaryAxisDirection,double HalfMapAngle,int nAzimuthPoints,int nZenithPoints,const char *fname,cColumnIntegrationSet* IntegrationSet);
      void Rectangular(double *xObservation,double *PrimaryAxisDirection,double *SecondaryAxisDirection,double xRange, double yRange, double *ImageOrigin, int nAzimuthPoints,int nZenithPoints,const char *fname,cColumnIntegrationSet* IntegrationSet);

    };

    cMap Map;

  public:
    //define 3 nodes on the surface of a bounding plane; index value: 0 -> xmin component of the coordinate, 1 -> xmax component of the coordinate
    static int x0PlaneNodeIndex[6][3]; //={ {0,0,0},{1,0,0},       {0,0,0},{0,1,0},           {0,0,0},{0,0,1}};
    static int x1PlaneNodeIndex[6][3]; //={ {0,1,0},{1,1,0},       {1,0,0},{1,1,0},           {1,0,0},{1,0,1}};
    static int x2PlaneNodeIndex[6][3]; //={ {0,0,1},{1,0,1},       {0,0,1},{0,1,1},           {0,1,0},{0,1,1}};
    static int PlaneNormal[6][3]; //=     { {1,0,0},{1,0,0},       {0,1,0},{0,1,0},           {0,0,1},{0,0,1}};

    struct cBoundingBoxFace {
      double x[3];
      double e0[3];
      double e1[3];
      double Normal[3];
      double e0Length;
      double e1Length;
    };

    cBoundingBoxFace BoundingBoxFace[6];

    //operations with the column integrals
    void Init(double *xGlobalMin,double *xGlobalMax);
    bool FindIntegrationLimits(double *x0,double *l,double& IntegrationPathLength,double *xStart,double *xFinish,double *xGlobalMin,double *xGlobalMax);

    //calculate the integral
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,double *l,double IntegrationPathLength,void (*Integrand)(double*,double*));
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,double*));


    //constructor
    cColumnIntegral() {
    }
  };

  cColumnIntegral ColumnIntegral;

  //mesh operations
  cBlock* GetBlock(double *x);
  cCell* GetCell(double *x);
  int CellId(double *x);
  void ConvertCellIdToLocalID(int CellId, int & iBlock, int & i, int & j, int & k);
  double CharacteristicCellSize(double *x);
  bool IsInsideDomain(double *x);

  //data file operations
  void LoadDataFile(const char *fname,const char* path);
  void SaveDataFile(const char *fname,sDataStructure &dat);
  void PrintVariableList();

  void SetBlockSize(int nx,int ny,int nz) {nBlockCellX=nx,nBlockCellY=ny,nBlockCellZ=nz;}
  void GetInterpolationStencil(double *x,cStencil* stencil);

  //init MPI
  int size,rank;

  inline void InitMPI() {

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
    MPI_Init(NULL, NULL);

    #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    int provided;
    MPI_Init_thread(NULL,NULL,MPI_THREAD_FUNNELED,&provided);

    #else
    #error Unknown option
    #endif //_COMPILATION_MODE_

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    MPI_GLOBAL_COMMUNICATOR=MPI_COMM_WORLD;

    //seed rnd generator
    rnd_seed(100);
  }

  inline void FinalizeMPI() {
    MPI_Finalize();
  }



  //constructor
  cPostProcess3D() {
    data.nVariables=0,data.nNodes=0;

    ColumnIntegral.PostProcess3D=this;
    ColumnIntegral.Map.ColumnIntegral=&ColumnIntegral;

    ParticleTrajectory.PostProcess3D=this;
    SurfaceData.PostProcess3D=this;
   
    initTrajCellMapFlag=false;
  }
};

#endif //_POST_PROCESS_3D_
