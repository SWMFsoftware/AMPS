/*
 * Regression test of the cell measures (volumes) of a mesh with an internal body
 * (cMeshAMRgeneric::InitCellMeasure(), InitCellMeasureBlock()).
 *
 * Usage: amps sphere|nastran
 *   sphere   the body is an analytic sphere (cInternalSphericalData);
 *   nastran  the body is a triangulated sphere (icosphere, 1280 faces) read from a NASTRAN file.
 *
 * The domain is the box [-1,1]^3 m; the body has radius 0.6 m and is centered off the
 * coordinate planes.  The mesh is refined near the surface, coarser inside the body (blocks
 * entirely inside it) and outside, so it has jumps of resolution next to cut blocks.
 *
 * Checks (one line "PASS|FAIL <check> ..." each in cell-measure.<case>.np<N>.report):
 *   total       sum of the measures of the interior cells = box volume - body volume:
 *               sphere: analytic, to the round-off 1e-9 of the body volume; nastran: the exact polyhedron volume
 *               (divergence theorem), within the sum of the integrator's per-cell bound
 *               1e-2*V_cell over the cut cells;
 *   inside      every cell entirely inside the body (all corners inside) has measure <= 0;
 *   outside     every cell entirely outside the body has the full cell volume;
 *   cut         every other cell agrees with an independent estimate within 1e-2*V_cell: the
 *               measure is integrated along z exactly (chord of the sphere / interval of the
 *               convex polyhedron cut by its face planes) with a 64x64 midpoint rule in x,y;
 *   noaccel     the per-block calculation InitCellMeasureBlock() (single-block path, no
 *               bounding-volume tree of the surface), after the measures of the block are
 *               reset, gives bit-identical measures of its interior cells and the same measures
 *               of its ghost cells within 1e-9 V_cell (a ghost cell can be calculated by any
 *               block that has it): for all blocks with the sphere, for 3 cut and 2 other blocks
 *               of each process with the triangulated sphere (there a point classification
 *               tests every face);
 *   cache       InitCellMeasure() loads the saved cache file and gets the same measures, a
 *               cache file with a wrong header is not loaded, and InitCellMeasure() then
 *               calculates the same measures again.
 * cell-measure.<case>.np<N>.cells holds the measures of all interior cells in a canonical
 * order: the files of runs with different numbers of MPI processes must be identical.
 */
#include "pic.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cNode;

const double kRadius=0.6;
const double kCenter[3]={0.0137,-0.0211,0.0173};

bool CaseNastran=false;
std::string Case;
FILE *Report=NULL;
int nFailed=0;

void Check(bool pass,const char *name,const char *format,...) {
  char message[2000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    printf("%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    fflush(Report);
  }

  if (!pass) nFailed++;
}

//---------------------------------- the body ----------------------------------
//the icosphere: vertices on the sphere, faces with outward normals (convex)
struct cPlane {double n[3],d;};
std::vector<double> IcoVertex;
std::vector<int> IcoFace;
std::vector<cPlane> IcoPlane;
double IcoVolume=0.0;

void BuildIcosphere(int nSubdivisions) {
  const double t=(1.0+std::sqrt(5.0))/2.0;
  double v0[12][3]={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
  int f0[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
                 {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
  std::vector<double> v;
  std::vector<int> f;

  for (int i=0;i<12;i++) {
    double l=std::sqrt(v0[i][0]*v0[i][0]+v0[i][1]*v0[i][1]+v0[i][2]*v0[i][2]);
    for (int d=0;d<3;d++) v.push_back(v0[i][d]/l);
  }
  for (int i=0;i<20;i++) for (int d=0;d<3;d++) f.push_back(f0[i][d]);

  for (int s=0;s<nSubdivisions;s++) {
    std::map<std::pair<int,int>,int> middle;
    std::vector<int> fnew;

    auto Middle=[&](int a,int b) {
      std::pair<int,int> key(std::min(a,b),std::max(a,b));
      std::map<std::pair<int,int>,int>::iterator it=middle.find(key);
      if (it!=middle.end()) return it->second;
      double m[3],l=0.0;
      for (int d=0;d<3;d++) m[d]=0.5*(v[3*a+d]+v[3*b+d]),l+=m[d]*m[d];
      l=std::sqrt(l);
      for (int d=0;d<3;d++) v.push_back(m[d]/l);
      int id=v.size()/3-1;
      middle[key]=id;
      return id;
    };

    for (size_t n=0;n<f.size();n+=3) {
      int a=f[n],b=f[n+1],c=f[n+2];
      int ab=Middle(a,b),bc=Middle(b,c),ca=Middle(c,a);
      int add[12]={a,ab,ca, b,bc,ab, c,ca,bc, ab,bc,ca};
      fnew.insert(fnew.end(),add,add+12);
    }

    f.swap(fnew);
  }

  IcoVertex.resize(v.size());
  for (size_t n=0;n<v.size();n+=3) for (int d=0;d<3;d++) IcoVertex[n+d]=kCenter[d]+kRadius*v[n+d];
  IcoFace=f;

  //outward planes and the exact volume of the polyhedron
  IcoVolume=0.0;
  for (size_t n=0;n<IcoFace.size();n+=3) {
    const double *a=&IcoVertex[3*IcoFace[n]],*b=&IcoVertex[3*IcoFace[n+1]],*c=&IcoVertex[3*IcoFace[n+2]];
    double e1[3],e2[3],nn[3],l=0.0,s=0.0;
    cPlane p;

    for (int d=0;d<3;d++) e1[d]=b[d]-a[d],e2[d]=c[d]-a[d];
    nn[0]=e1[1]*e2[2]-e1[2]*e2[1],nn[1]=e1[2]*e2[0]-e1[0]*e2[2],nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
    for (int d=0;d<3;d++) l+=nn[d]*nn[d],s+=nn[d]*(a[d]-kCenter[d]);
    l=std::sqrt(l);
    if (s<0.0) l=-l;
    for (int d=0;d<3;d++) p.n[d]=nn[d]/l;
    p.d=p.n[0]*a[0]+p.n[1]*a[1]+p.n[2]*a[2];
    IcoPlane.push_back(p);

    //tetrahedron with the center
    double ra[3],rb[3],rc[3];
    for (int d=0;d<3;d++) ra[d]=a[d]-kCenter[d],rb[d]=b[d]-kCenter[d],rc[d]=c[d]-kCenter[d];
    IcoVolume+=std::fabs(ra[0]*(rb[1]*rc[2]-rb[2]*rc[1])-ra[1]*(rb[0]*rc[2]-rb[2]*rc[0])+ra[2]*(rb[0]*rc[1]-rb[1]*rc[0]))/6.0;
  }
}

void WriteNastranFile(const char *fname) {
  FILE *f=fopen(fname,"w");

  fprintf(f,"$ AMPS cell-measure test: icosphere radius %.17e\n",kRadius);
  for (size_t n=0;n<IcoVertex.size();n+=3) {
    fprintf(f,"GRID* %i 0 %.17E %.17E\n* %.17E\n",(int)(n/3+1),IcoVertex[n],IcoVertex[n+1],IcoVertex[n+2]);
  }
  for (size_t n=0;n<IcoFace.size();n+=3) fprintf(f,"CTRIA3 %i 1 %i %i %i\n",(int)(n/3+1),IcoFace[n]+1,IcoFace[n+1]+1,IcoFace[n+2]+1);
  fprintf(f,"ENDDATA\n");
  fclose(f);
}

double BodyVolume() {return CaseNastran ? IcoVolume : 4.0/3.0*Pi*pow(kRadius,3);}

//inside the body (closed set for the polyhedron, open for the sphere: only used with margins)
bool PointInsideBody(const double *x) {
  if (!CaseNastran) {
    double r2=0.0;
    for (int d=0;d<3;d++) r2+=pow(x[d]-kCenter[d],2);
    return r2<kRadius*kRadius;
  }

  for (size_t n=0;n<IcoPlane.size();n++) {
    if (IcoPlane[n].n[0]*x[0]+IcoPlane[n].n[1]*x[1]+IcoPlane[n].n[2]*x[2]>IcoPlane[n].d) return false;
  }
  return true;
}

//classification of a cell: 1 entirely inside the body, 0 entirely outside, 2 otherwise
int ClassifyCell(const double *xmin,const double *xmax,std::vector<int>& planes) {
  planes.clear();

  if (!CaseNastran) {
    double dmin2=0.0,dmax2=0.0;
    for (int d=0;d<3;d++) {
      double a=xmin[d]-kCenter[d],b=xmax[d]-kCenter[d];
      double gap=(a>0.0) ? a : ((b<0.0) ? -b : 0.0);
      dmin2+=gap*gap,dmax2+=std::max(a*a,b*b);
    }
    if (dmax2<kRadius*kRadius*(1.0-1.0E-12)) return 1;
    if (dmin2>kRadius*kRadius*(1.0+1.0E-12)) return 0;
    return 2;
  }

  bool allInside=true;
  for (size_t n=0;n<IcoPlane.size();n++) {
    const cPlane& p=IcoPlane[n];
    double c=0.0,h=0.0;
    for (int d=0;d<3;d++) c+=p.n[d]*0.5*(xmin[d]+xmax[d]),h+=std::fabs(p.n[d])*0.5*(xmax[d]-xmin[d]);
    if (c-h>p.d+1.0E-12) return 0;            //the box is on the outer side of a face plane
    if (c+h>p.d-1.0E-12) planes.push_back(n),allInside=false; //the plane cuts (or touches) the box
  }
  return allInside ? 1 : 2;
}

//independent estimate of the volume of the part of the cell outside the body: exact along z, midpoint rule in x,y
double EstimateCellVolume(const double *xmin,const double *xmax,const std::vector<int>& planes) {
  const int n=64;
  double inside=0.0;
  const double hx=(xmax[0]-xmin[0])/n,hy=(xmax[1]-xmin[1])/n;

  for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
    double x=xmin[0]+(i+0.5)*hx,y=xmin[1]+(j+0.5)*hy,zlo,zhi;

    if (!CaseNastran) {
      double h2=kRadius*kRadius-pow(x-kCenter[0],2)-pow(y-kCenter[1],2);
      if (h2<=0.0) continue;
      zlo=kCenter[2]-std::sqrt(h2),zhi=kCenter[2]+std::sqrt(h2);
    }
    else {
      zlo=-1.0E30,zhi=1.0E30;
      for (size_t m=0;m<planes.size();m++) {
        const cPlane& p=IcoPlane[planes[m]];
        double r=p.d-p.n[0]*x-p.n[1]*y;
        if (p.n[2]>0.0) zhi=std::min(zhi,r/p.n[2]);
        else if (p.n[2]<0.0) zlo=std::max(zlo,r/p.n[2]);
        else if (r<0.0) zhi=zlo;
      }
    }

    double a=std::max(zlo,xmin[2]),b=std::min(zhi,xmax[2]);
    if (b>a) inside+=(b-a)*hx*hy;
  }

  return (xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])-inside;
}

//---------------------------------- the mesh ----------------------------------
//the cell size near the surface: 0.05 m for the sphere, 0.1 m for the triangulated sphere (whose cut cells are integrated numerically)
double SurfaceCellSize=0.06;

double Resolution(double *x) {
  double r=0.0;
  for (int d=0;d<3;d++) r+=pow(x[d]-kCenter[d],2);
  r=std::sqrt(r);

  if (std::fabs(r-kRadius)<2.0*SurfaceCellSize) return SurfaceCellSize;
  if (r<kRadius) return 0.15;                    //coarser inside the body
  return 0.25;
}

double SurfaceResolution(CutCell::cTriangleFace*) {return SurfaceCellSize;}
double UnitLoad(cNode*) {return 1.0;}

void CollectBlocks(cNode* node,std::vector<cNode*>& blocks) {
  if (node->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    for (int n=0;n<8;n++) if (node->downNode[n]!=NULL) CollectBlocks(node->downNode[n],blocks);
  }
  else if (node->block!=NULL) blocks.push_back(node);
}

//the interior cells of the blocks of this process
struct cCell {
  long int key[4];
  double xmin[3],xmax[3],measure;
  bool operator<(const cCell& c) const {return std::lexicographical_compare(key,key+4,c.key,c.key+4);}
};

void GetOwnCells(std::vector<cCell>& cells) {
  cells.clear();

  for (cNode* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    double dx[3];
    for (int d=0;d<3;d++) dx[d]=(node->xmax[d]-node->xmin[d])/((d==0) ? _BLOCK_CELLS_X_ : ((d==1) ? _BLOCK_CELLS_Y_ : _BLOCK_CELLS_Z_));

    for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
      int idx[3]={i,j,k};
      cCell c;

      if (cell==NULL) continue;
      c.key[0]=node->RefinmentLevel;
      for (int d=0;d<3;d++) {
        c.xmin[d]=node->xmin[d]+idx[d]*dx[d],c.xmax[d]=node->xmin[d]+(idx[d]+1)*dx[d];
        c.key[d+1]=std::lround((c.xmin[d]+1.0)/dx[d]);
      }
      c.measure=cell->Measure;
      cells.push_back(c);
    }
  }
}

//the measures of the center nodes of the allocated blocks 'blocks' (interior and ghost cells); a node is
//'interior' when it is an interior cell of one of them; 'volume' is the full volume of the cell
struct cNodeMeasure {double measure,volume; bool interior;};
typedef std::map<PIC::Mesh::cDataCenterNode*,cNodeMeasure> cNodeMeasureMap;

void GetNodeMeasures(const std::vector<cNode*>& blocks,cNodeMeasureMap& m) {
  m.clear();

  for (size_t n=0;n<blocks.size();n++) {
    double v=1.0;
    for (int d=0;d<3;d++) v*=(blocks[n]->xmax[d]-blocks[n]->xmin[d])/((d==0) ? _BLOCK_CELLS_X_ : ((d==1) ? _BLOCK_CELLS_Y_ : _BLOCK_CELLS_Z_));

    for (int k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) for (int j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) for (int i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      PIC::Mesh::cDataCenterNode *cell=blocks[n]->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
      bool interior=(i>=0)&&(i<_BLOCK_CELLS_X_)&&(j>=0)&&(j<_BLOCK_CELLS_Y_)&&(k>=0)&&(k<_BLOCK_CELLS_Z_);

      if (cell==NULL) continue;

      cNodeMeasureMap::iterator it=m.find(cell);
      if (it==m.end()) {cNodeMeasure t={cell->Measure,v,interior}; m[cell]=t;}
      else if (interior) it->second.interior=true;
    }
  }
}

//compare the measures of the nodes in 'a' with their present values: interior cells must be bit-identical,
//the other (ghost) cells agree within 1e-9 of the cell volume (they can be calculated by any block that has them,
//with the coordinates of the cell in that block)
void CompareNodeMeasures(const cNodeMeasureMap& a,long int& nInteriorDiff,long int& nGhostDiff,double& maxGhostDiff) {
  long int lv[2]={0,0},gv[2];
  double maxDiff=0.0;

  for (cNodeMeasureMap::const_iterator it=a.begin();it!=a.end();it++) {
    double m=it->first->Measure;

    if (it->second.interior) {
      if (std::memcmp(&m,&it->second.measure,sizeof(double))!=0) lv[0]++;
    }
    else {
      double diff=std::fabs(std::max(m,0.0)-std::max(it->second.measure,0.0))/it->second.volume;
      if (((m>0.0)!=(it->second.measure>0.0))||(diff>1.0E-9)) lv[1]++;
      maxDiff=std::max(maxDiff,diff);
    }
  }

  MPI_Allreduce(lv,gv,2,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&maxDiff,&maxGhostDiff,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  nInteriorDiff=gv[0],nGhostDiff=gv[1];
}

}  // namespace


int main(int argc,char **argv) {
  PIC::Init_BeforeParser();
  PIC::Mesh::initCellSamplingDataBuffer();

  Case=(argc>1) ? argv[1] : "sphere";
  if ((Case!="sphere")&&(Case!="nastran")) {
    if (PIC::ThisThread==0) printf("usage: %s sphere|nastran\n",argv[0]);
    MPI_Finalize();
    return 2;
  }
  CaseNastran=(Case=="nastran");
  SurfaceCellSize=CaseNastran ? 0.12 : 0.06;

  char fname[400];
  sprintf(fname,"cell-measure.%s.np%i.report",Case.c_str(),PIC::nTotalThreads);
  if (PIC::ThisThread==0) Report=fopen(fname,"w");

  //the body
  if (CaseNastran) {
    BuildIcosphere(3);
    if (PIC::ThisThread==0) WriteNastranFile("cell-measure-icosphere.bdf");
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    CutCell::ReadNastranSurfaceMeshLongFormat("cell-measure-icosphere.bdf");
    PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
    PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;
  }
  else {
    PIC::BC::InternalBoundary::Sphere::Init();
    cInternalBoundaryConditionsDescriptor descriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    cInternalSphericalData *sphere=(cInternalSphericalData*)descriptor.BoundaryElement;
    double x0[3]={kCenter[0],kCenter[1],kCenter[2]};

    sphere->SetSphereGeometricalParameters(x0,kRadius);
    sphere->localResolution=Resolution;
  }

  //the mesh: generated by process 0 and read from the file by the others, as in the applications
  double xmin[3]={-1.0,-1.0,-1.0},xmax[3]={1.0,1.0,1.0};

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Resolution);

  if (PIC::ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile("cell-measure.mesh");
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh->readMeshFile("cell-measure.mesh");
  }

  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  //remove a cache of an earlier run
  unsigned long signature=PIC::Mesh::mesh->getMeshSignature();
  char cache[600];
  PIC::Mesh::mesh->GetCenterNodeMeasureFileName(cache,NULL,signature);
  if (PIC::ThisThread==0) remove(cache);
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //1. the full initialization (root: the bounding-volume tree of the surface is active)
  double t0=MPI_Wtime();
  PIC::Mesh::mesh->InitCellMeasure();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  double tInit=MPI_Wtime()-t0;

  std::vector<cNode*> blocks;
  cNodeMeasureMap reference;
  CollectBlocks(PIC::Mesh::mesh->rootTree,blocks);
  GetNodeMeasures(blocks,reference);

  std::vector<cCell> cells;
  GetOwnCells(cells);

  //2. the checks of the interior cells
  long int nInside=0,nOutside=0,nCut=0,badInside=0,badOutside=0,badCut=0;
  double sum=0.0,cutVolume=0.0,maxCutErr=0.0,sumCutErr=0.0;
  std::vector<int> planes;

  for (size_t n=0;n<cells.size();n++) {
    const cCell& c=cells[n];
    double v=(c.xmax[0]-c.xmin[0])*(c.xmax[1]-c.xmin[1])*(c.xmax[2]-c.xmin[2]);
    int type=ClassifyCell(c.xmin,c.xmax,planes);

    if (c.measure>0.0) sum+=c.measure;

    if (type==1) {
      nInside++;
      if (c.measure>0.0) badInside++;
    }
    else if (type==0) {
      nOutside++;
      if (std::fabs(c.measure-v)>1.0E-12*v) badOutside++;
    }
    else {
      double estimate=EstimateCellVolume(c.xmin,c.xmax,planes),err=std::fabs(std::max(c.measure,0.0)-estimate)/v;
      nCut++,cutVolume+=v,sumCutErr+=err;
      if (err>maxCutErr) maxCutErr=err;
      if (err>1.0E-2) badCut++;
    }
  }

  long int lv[7]={(long int)cells.size(),nInside,nOutside,nCut,badInside,badOutside,badCut},glv[7];
  double dv[3]={sum,cutVolume,sumCutErr},gdv[3],gmax;
  MPI_Allreduce(lv,glv,7,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(dv,gdv,3,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&maxCutErr,&gmax,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);

  const double box=8.0,body=BodyVolume(),exact=box-body;
  //nastran: the bound 1e-2*V_cell of the cut-cell integration summed over the cut cells;
  //sphere: the round-off of the analytic volume of the cut cells (a few 1e-11 of the body volume are observed)
  const double tolerance=CaseNastran ? 1.0E-2*gdv[1] : 1.0E-9*body;

  if (PIC::ThisThread==0) printf("cell-measure %s np=%i: %li cells (%li inside the body, %li outside, %li cut), InitCellMeasure %.3f s\n",
      Case.c_str(),PIC::nTotalThreads,glv[0],glv[1],glv[2],glv[3],tInit);

  Check(std::fabs(gdv[0]-exact)<=tolerance,"total","sum=%.17e box-body=%.17e difference=%.3e (%.3e of the body volume) tolerance=%.3e",
      gdv[0],exact,gdv[0]-exact,(gdv[0]-exact)/body,tolerance);
  Check((glv[1]>0)&&(glv[4]==0),"inside","%li cells entirely inside the body, %li of them with measure>0",glv[1],glv[4]);
  Check((glv[2]>0)&&(glv[5]==0),"outside","%li cells entirely outside the body, %li of them without the full volume",glv[2],glv[5]);
  Check((glv[3]>0)&&(glv[6]==0),"cut","%li other cells: |measure-estimate|/V_cell max=%.3e mean=%.3e, %li above 1e-2",
      glv[3],gmax,(glv[3]>0) ? gdv[2]/glv[3] : 0.0,glv[6]);

  //3. the interior cells in the canonical order (identical for any number of processes)
  {
    std::vector<long int> keys;
    std::vector<double> measures;
    for (size_t n=0;n<cells.size();n++) {keys.insert(keys.end(),cells[n].key,cells[n].key+4);measures.push_back(cells[n].measure);}

    int nLocal=cells.size(),*counts=NULL,*displs=NULL,*counts4=NULL,*displs4=NULL,nTotal=0;
    if (PIC::ThisThread==0) counts=new int[PIC::nTotalThreads],displs=new int[PIC::nTotalThreads],counts4=new int[PIC::nTotalThreads],displs4=new int[PIC::nTotalThreads];
    MPI_Gather(&nLocal,1,MPI_INT,counts,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    std::vector<long int> allKeys;
    std::vector<double> allMeasures;
    if (PIC::ThisThread==0) {
      for (int t=0;t<PIC::nTotalThreads;t++) displs[t]=nTotal,counts4[t]=4*counts[t],displs4[t]=4*nTotal,nTotal+=counts[t];
      allKeys.resize(4*nTotal),allMeasures.resize(nTotal);
    }
    MPI_Gatherv(keys.data(),4*nLocal,MPI_LONG,allKeys.data(),counts4,displs4,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gatherv(measures.data(),nLocal,MPI_DOUBLE,allMeasures.data(),counts,displs,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      std::vector<cCell> all(nTotal);
      for (int n=0;n<nTotal;n++) {for (int d=0;d<4;d++) all[n].key[d]=allKeys[4*n+d]; all[n].measure=allMeasures[n];}
      std::sort(all.begin(),all.end());

      sprintf(fname,"cell-measure.%s.np%i.cells",Case.c_str(),PIC::nTotalThreads);
      FILE *f=fopen(fname,"w");
      for (int n=0;n<nTotal;n++) {fwrite(all[n].key,sizeof(long int),4,f);fwrite(&all[n].measure,sizeof(double),1,f);}
      fclose(f);
      delete [] counts; delete [] displs; delete [] counts4; delete [] displs4;
    }
  }

  //4. the per-block calculation without the bounding-volume tree, for every block (sphere) or, with the triangulated
  //surface (the classification of a point tests all faces), for 3 blocks that the surface cuts and 2 other blocks of
  //each process; the measures of these blocks are reset
  {
    std::vector<cNode*> selected,cutBlocks,otherBlocks;
    cNodeMeasureMap selectedReference;

    if (CaseNastran) {
      for (size_t n=0;n<blocks.size();n++) ((blocks[n]->FirstTriangleCutFace!=NULL) ? cutBlocks : otherBlocks).push_back(blocks[n]);
      for (size_t n=0;n<cutBlocks.size();n+=std::max<size_t>(1,cutBlocks.size()/3)) selected.push_back(cutBlocks[n]);
      for (size_t n=0;n<otherBlocks.size();n+=std::max<size_t>(1,otherBlocks.size()/2)) selected.push_back(otherBlocks[n]);
    }
    else selected=blocks;
    GetNodeMeasures(selected,selectedReference);
    for (cNodeMeasureMap::iterator it=selectedReference.begin();it!=selectedReference.end();it++) it->first->Measure=-1.0;

    t0=MPI_Wtime();
    for (size_t n=0;n<selected.size();n++) PIC::Mesh::mesh->InitCellMeasureBlock(selected[n]);
    double tBlock=MPI_Wtime()-t0;

    long int nInteriorDiff,nGhostDiff;
    double maxGhostDiff;
    CompareNodeMeasures(selectedReference,nInteriorDiff,nGhostDiff,maxGhostDiff);
    Check((nInteriorDiff==0)&&(nGhostDiff==0),"noaccel","InitCellMeasureBlock() for %li of the %li local blocks of process 0 (%.3f s): %li interior cells differ from InitCellMeasure(), %li ghost cells by more than 1e-9 V_cell (max %.3e)",
        (long int)selected.size(),(long int)blocks.size(),tBlock,nInteriorDiff,nGhostDiff,maxGhostDiff);
  }

  //5. the cache file
  {
    bool loaded=PIC::Mesh::mesh->LoadCenterNodeMeasure(NULL,NULL);
    long int nInteriorDiff,nGhostDiff;
    double maxGhostDiff;

    CompareNodeMeasures(reference,nInteriorDiff,nGhostDiff,maxGhostDiff);
    Check(loaded&&(nInteriorDiff==0)&&(nGhostDiff==0),"cache","the saved file %s is loaded: %s; %li interior cells differ, %li ghost cells by more than 1e-9 V_cell (max %.3e)",
        cache,loaded ? "yes" : "no",nInteriorDiff,nGhostDiff,maxGhostDiff);

    //a file with the right name and size but without the header of this version
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) {
      FILE *f=fopen(cache,"r+");
      double minusOne=-1.0;
      for (int n=0;n<8;n++) fwrite(&minusOne,sizeof(double),1,f);
      fclose(f);
    }
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    loaded=PIC::Mesh::mesh->LoadCenterNodeMeasure(NULL,NULL);
    int loadedAny=loaded ? 1 : 0,loadedGlobal;
    MPI_Allreduce(&loadedAny,&loadedGlobal,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    Check(loadedGlobal==0,"cache-header","a cache file with a wrong header is %s",(loadedGlobal==0) ? "not loaded" : "LOADED");

    PIC::Mesh::mesh->InitCellMeasure();
    CompareNodeMeasures(reference,nInteriorDiff,nGhostDiff,maxGhostDiff);
    Check((nInteriorDiff==0)&&(nGhostDiff==0),"recalculate","InitCellMeasure() with the rejected cache file: %li interior cells differ, %li ghost cells by more than 1e-9 V_cell (max %.3e)",
        nInteriorDiff,nGhostDiff,maxGhostDiff);
  }

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s: %i checks failed\n",(nFailed==0) ? "PASS" : "FAIL",Case.c_str(),nFailed);
    fclose(Report);
  }

  MPI_Finalize();
  return (nFailed==0) ? 0 : 1;
}
