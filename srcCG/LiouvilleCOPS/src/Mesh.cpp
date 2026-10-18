#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <array>

#include "Mesh.h"
#include "SpiceUsr.h"

//convert a NASTRAN real field: accepts the "1.5-3" / "1.5+3" short exponent notation
static double NastranReal(std::string s) {
  s.erase(std::remove(s.begin(),s.end(),' '),s.end());

  if (s.empty()) {
    fprintf(stderr,"Error: empty NASTRAN real field\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i=1;i<s.size();i++) if (((s[i]=='-')||(s[i]=='+'))&&(s[i-1]!='E')&&(s[i-1]!='e')&&(s[i-1]!='D')&&(s[i-1]!='d')) {
    s.insert(i,"E");
    break;
  }

  for (auto& ch : s) if ((ch=='D')||(ch=='d')) ch='E';

  char *end;
  double res=strtod(s.c_str(),&end);

  if (*end!=0) {
    fprintf(stderr,"Error: cannot convert NASTRAN real field '%s'\n",s.c_str());
    exit(EXIT_FAILURE);
  }

  return res;
}

static std::string Field(const std::string& line,size_t begin,size_t width) {
  return (line.size()>begin) ? line.substr(begin,width) : std::string();
}

void cMesh::ReadNastran(const std::string& fname) {
  std::ifstream fin(fname);
  std::string line,next;
  std::unordered_map<long,int> NodeIndex;
  std::vector<long> TriNodes;
  long nTetra=0;

  if (!fin) {
    fprintf(stderr,"Error: cannot open the shape model file '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  while (std::getline(fin,line)) {
    if (line.compare(0,5,"GRID*")==0) {
      //large field format: 8-character key, 16-character fields, continuation line holds z
      if (!std::getline(fin,next)) {
        fprintf(stderr,"Error: missing continuation of a GRID* card in '%s'\n",fname.c_str());
        exit(EXIT_FAILURE);
      }

      long id=atol(Field(line,8,16).c_str());
      Vec3 x(NastranReal(Field(line,40,16)),NastranReal(Field(line,56,16)),NastranReal(Field(next,8,16)));

      NodeIndex[id]=(int)Nodes.size();
      Nodes.push_back(x);
    }
    else if ((line.compare(0,4,"GRID")==0)&&(line.compare(0,5,"GRID*")!=0)) {
      //small field format
      long id=atol(Field(line,8,8).c_str());
      Vec3 x(NastranReal(Field(line,24,8)),NastranReal(Field(line,32,8)),NastranReal(Field(line,40,8)));

      NodeIndex[id]=(int)Nodes.size();
      Nodes.push_back(x);
    }
    else if (line.compare(0,6,"CTRIA3")==0) {
      std::istringstream ss(line);
      std::string key;
      long id,pid,n0,n1,n2;

      if (!(ss >> key >> id >> pid >> n0 >> n1 >> n2)) {
        fprintf(stderr,"Error: cannot parse the CTRIA3 card: %s\n",line.c_str());
        exit(EXIT_FAILURE);
      }

      TriNodes.push_back(n0),TriNodes.push_back(n1),TriNodes.push_back(n2);
    }
    else if (line.compare(0,6,"CTETRA")==0) nTetra++;
  }

  //convert the node IDs into the local indexes; keep only the nodes used by the surface triangulation
  std::vector<int> NewIndex(Nodes.size(),-1);
  std::vector<Vec3> SurfaceNodes;

  Faces.resize(TriNodes.size()/3);

  for (size_t i=0;i<TriNodes.size();i++) {
    auto it=NodeIndex.find(TriNodes[i]);

    if (it==NodeIndex.end()) {
      fprintf(stderr,"Error: CTRIA3 references the undefined GRID %li\n",TriNodes[i]);
      exit(EXIT_FAILURE);
    }

    int old=it->second;

    if (NewIndex[old]==-1) {
      NewIndex[old]=(int)SurfaceNodes.size();
      SurfaceNodes.push_back(Nodes[old]);
    }

    Faces[i/3].node[i%3]=NewIndex[old];
  }

  printf("Shape model '%s': %zu GRID nodes, %li CTETRA (ignored), %zu CTRIA3 surface triangles, %zu surface nodes\n",
      fname.c_str(),Nodes.size(),nTetra,Faces.size(),SurfaceNodes.size());

  Nodes.swap(SurfaceNodes);

  if (Faces.empty()) {
    fprintf(stderr,"Error: the shape model has no CTRIA3 surface triangles\n");
    exit(EXIT_FAILURE);
  }
}

void cMesh::ReadDSK(const std::string& fname,const std::string& expectedFrame) {
  SpiceInt handle,nv,np;
  SpiceDLADescr dladsc;
  SpiceDSKDescr dskdsc;
  SpiceBoolean found;
  SpiceChar frname[64];

  dasopr_c(fname.c_str(),&handle);
  dlabfs_c(handle,&dladsc,&found);

  if (failed_c()||!found) {
    fprintf(stderr,"Error: cannot read the DSK file '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  dskgd_c(handle,&dladsc,&dskdsc);
  frmnam_c((SpiceInt)dskdsc.frmcde,sizeof(frname),frname);

  if (dskdsc.dtype!=2) {
    fprintf(stderr,"Error: DSK '%s' is of type %i (only type 2 plate models are supported)\n",fname.c_str(),(int)dskdsc.dtype);
    exit(EXIT_FAILURE);
  }

  if (expectedFrame!=frname) {
    fprintf(stderr,"Error: the DSK '%s' is given in frame '%s' (frame code %i), but BodyFrame is '%s'\n",fname.c_str(),frname,(int)dskdsc.frmcde,expectedFrame.c_str());
    exit(EXIT_FAILURE);
  }

  //only single-segment DSK files are supported (a second segment would be a different surface/level of detail)
  SpiceDLADescr next;
  SpiceBoolean more;
  dlafns_c(handle,&dladsc,&next,&more);

  if (more) {
    fprintf(stderr,"Error: the DSK '%s' has more than one segment\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  dskz02_c(handle,&dladsc,&nv,&np);

  std::vector<double> v(3*(size_t)nv);
  std::vector<int> pl(3*(size_t)np);
  const SpiceInt chunk=100000;
  SpiceInt n;

  for (SpiceInt i=0;i<nv;i+=n) {
    dskv02_c(handle,&dladsc,i+1,std::min(chunk,nv-i),&n,(SpiceDouble(*)[3])(v.data()+3*(size_t)i));
    if (failed_c()||(n<=0)) {fprintf(stderr,"Error: cannot read the DSK vertices\n"); exit(EXIT_FAILURE);}
  }

  for (SpiceInt i=0;i<np;i+=n) {
    dskp02_c(handle,&dladsc,i+1,std::min(chunk,np-i),&n,(SpiceInt(*)[3])(pl.data()+3*(size_t)i));
    if (failed_c()||(n<=0)) {fprintf(stderr,"Error: cannot read the DSK plates\n"); exit(EXIT_FAILURE);}
  }

  dascls_c(handle);

  Nodes.resize(nv);
  for (SpiceInt i=0;i<nv;i++) Nodes[i]=Vec3(1.0E3*v[3*i],1.0E3*v[3*i+1],1.0E3*v[3*i+2]);

  Faces.resize(np);
  for (SpiceInt i=0;i<np;i++) for (int j=0;j<3;j++) {
    int k=pl[3*i+j]-1;   //DSK vertex indices are 1-based

    if ((k<0)||(k>=nv)) {
      fprintf(stderr,"Error: DSK plate %i references the vertex %i (nv=%i)\n",(int)i+1,k+1,(int)nv);
      exit(EXIT_FAILURE);
    }

    Faces[i].node[j]=k;
  }

  printf("Shape model '%s': DSK type 2, frame %s, %i vertices, %i plates\n",fname.c_str(),frname,(int)nv,(int)np);
}

void cMesh::MakeIcosphere(double R,int nRefinement) {
  const double t=(1.0+std::sqrt(5.0))/2.0;
  const double v[12][3]={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
  const int f[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
                      {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};

  Nodes.clear(),Faces.clear();
  for (int i=0;i<12;i++) Nodes.push_back(Normalized(Vec3(v[i][0],v[i][1],v[i][2])));

  std::vector<std::array<int,3> > tri;
  for (int i=0;i<20;i++) tri.push_back({f[i][0],f[i][1],f[i][2]});

  for (int level=0;level<nRefinement;level++) {
    std::map<std::pair<int,int>,int> mid;
    std::vector<std::array<int,3> > newtri;

    auto Midpoint=[&](int a,int b) {
      std::pair<int,int> key(std::min(a,b),std::max(a,b));
      auto it=mid.find(key);
      if (it!=mid.end()) return it->second;

      Nodes.push_back(Normalized(Nodes[a]+Nodes[b]));
      mid[key]=(int)Nodes.size()-1;
      return (int)Nodes.size()-1;
    };

    for (auto& tr : tri) {
      int a=Midpoint(tr[0],tr[1]),b=Midpoint(tr[1],tr[2]),c=Midpoint(tr[2],tr[0]);

      newtri.push_back({tr[0],a,c});
      newtri.push_back({tr[1],b,a});
      newtri.push_back({tr[2],c,b});
      newtri.push_back({a,b,c});
    }

    tri.swap(newtri);
  }

  for (auto& x : Nodes) x*=R;

  Faces.resize(tri.size());
  for (size_t i=0;i<tri.size();i++) for (int j=0;j<3;j++) Faces[i].node[j]=tri[i][j];
}

void cMesh::FinalizeGeometry(bool verbose,bool requireClosed) {
  //1. topology: every edge must be shared by exactly two triangles, traversed in opposite directions
  std::unordered_map<unsigned long long,int> DirectedEdge;
  long nInconsistent=0;

  auto Key=[](int a,int b) {return ((unsigned long long)(unsigned int)a<<32)|(unsigned int)b;};

  DirectedEdge.reserve(3*Faces.size());
  for (auto& f : Faces) for (int i=0;i<3;i++) DirectedEdge[Key(f.node[i],f.node[(i+1)%3])]++;

  nOpenEdges=0;

  for (auto& e : DirectedEdge) {
    int a=(int)(e.first>>32),b=(int)(e.first&0xffffffffULL);

    if (e.second>1) nInconsistent++;

    auto rev=DirectedEdge.find(Key(b,a));
    if ((rev==DirectedEdge.end())||(rev->second+e.second!=2)) nOpenEdges++;
  }

  if (nInconsistent!=0) {
    fprintf(stderr,"Error: the surface triangulation is not consistently oriented (%li directed edges are used more than once)\n",nInconsistent);
    exit(EXIT_FAILURE);
  }

  if (nOpenEdges!=0) {
    if (requireClosed) {
      fprintf(stderr,"Error: the surface is not closed (%li open/non-manifold edges)\n",nOpenEdges);
      exit(EXIT_FAILURE);
    }

    printf("WARNING: the surface is not closed: %li open/non-manifold edges; the inside-nucleus test and the orientation from the signed volume are approximate\n",nOpenEdges);
  }

  //2. orientation: the signed volume must be positive for outward normals; otherwise flip all triangles
  double vol=0.0;

  for (auto& f : Faces) vol+=Dot(Nodes[f.node[0]],Cross(Nodes[f.node[1]],Nodes[f.node[2]]))/6.0;

  bool flipped=false;

  if (vol<0.0) {
    for (auto& f : Faces) std::swap(f.node[1],f.node[2]);
    vol=-vol;
    flipped=true;
  }

  //3. normals, areas, centroids
  TotalArea=0.0;

  for (auto& f : Faces) {
    Vec3 a=Nodes[f.node[0]],b=Nodes[f.node[1]],c=Nodes[f.node[2]];
    Vec3 cr=Cross(b-a,c-a);
    double l=Length(cr);

    if (l<=0.0) {
      fprintf(stderr,"Error: degenerate surface triangle\n");
      exit(EXIT_FAILURE);
    }

    f.Normal=cr/l;
    f.Area=0.5*l;
    f.Centroid=(a+b+c)/3.0;
    f.MaxEdge=std::max(Length(b-a),std::max(Length(c-b),Length(a-c)));
    TotalArea+=f.Area;
  }

  Volume=vol;

  if (verbose) {
    printf("Surface: %s and consistently oriented; triangle node order %s to obtain outward normals\n",(nOpenEdges==0) ? "closed" : "NOT closed",flipped ? "was reversed" : "kept");
    printf("Surface: total area = %.6e m^2, enclosed volume = %.6e m^3\n",TotalArea,Volume);
  }
}

//closest point computation (Ericson, Real-Time Collision Detection, 5.1.5)
double PointTriangleDistance(const Vec3& p,const Vec3& a,const Vec3& b,const Vec3& c) {
  Vec3 ab=b-a,ac=c-a,ap=p-a;
  double d1=Dot(ab,ap),d2=Dot(ac,ap);
  if ((d1<=0.0)&&(d2<=0.0)) return Length(p-a);

  Vec3 bp=p-b;
  double d3=Dot(ab,bp),d4=Dot(ac,bp);
  if ((d3>=0.0)&&(d4<=d3)) return Length(p-b);

  double vc=d1*d4-d3*d2;
  if ((vc<=0.0)&&(d1>=0.0)&&(d3<=0.0)) {
    double v=d1/(d1-d3);
    return Length(p-(a+v*ab));
  }

  Vec3 cp=p-c;
  double d5=Dot(ab,cp),d6=Dot(ac,cp);
  if ((d6>=0.0)&&(d5<=d6)) return Length(p-c);

  double vb=d5*d2-d1*d6;
  if ((vb<=0.0)&&(d2>=0.0)&&(d6<=0.0)) {
    double w=d2/(d2-d6);
    return Length(p-(a+w*ac));
  }

  double va=d3*d6-d5*d4;
  if ((va<=0.0)&&((d4-d3)>=0.0)&&((d5-d6)>=0.0)) {
    double w=(d4-d3)/((d4-d3)+(d5-d6));
    return Length(p-(b+w*(c-b)));
  }

  double denom=1.0/(va+vb+vc);
  double v=vb*denom,w=vc*denom;
  return Length(p-(a+ab*v+ac*w));
}
