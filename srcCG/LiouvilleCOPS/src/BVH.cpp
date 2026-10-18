#include <algorithm>
#include <cmath>
#include <limits>

#include "BVH.h"

static const int LeafSize=4;

void cBVH::Build(const cMesh* mesh) {
  Mesh=mesh;
  Nodes.clear();
  FaceOrder.resize(mesh->Faces.size());

  std::vector<Vec3> centroids(mesh->Faces.size());

  for (size_t i=0;i<mesh->Faces.size();i++) {
    FaceOrder[i]=(int)i;
    centroids[i]=mesh->Faces[i].Centroid;
  }

  Nodes.reserve(2*mesh->Faces.size()/LeafSize+16);
  BuildNode(0,(int)FaceOrder.size(),centroids);
}

int cBVH::BuildNode(int first,int count,std::vector<Vec3>& centroids) {
  cNode node;
  double inf=std::numeric_limits<double>::max();

  node.bmin=Vec3(inf,inf,inf),node.bmax=Vec3(-inf,-inf,-inf);

  for (int i=first;i<first+count;i++) for (int j=0;j<3;j++) {
    Vec3 v=Mesh->Vertex(FaceOrder[i],j);

    node.bmin=Min(node.bmin,v);
    node.bmax=Max(node.bmax,v);
  }

  node.left=-1,node.right=-1,node.first=first,node.count=count;

  int inode=(int)Nodes.size();
  Nodes.push_back(node);

  if (count<=LeafSize) return inode;

  //split at the median of the centroids along the longest extent of the centroid bounding box
  Vec3 cmin(inf,inf,inf),cmax(-inf,-inf,-inf);

  for (int i=first;i<first+count;i++) {
    cmin=Min(cmin,centroids[FaceOrder[i]]);
    cmax=Max(cmax,centroids[FaceOrder[i]]);
  }

  Vec3 ext=cmax-cmin;
  int axis=(ext.x>ext.y) ? ((ext.x>ext.z) ? 0 : 2) : ((ext.y>ext.z) ? 1 : 2);
  int mid=first+count/2;

  std::nth_element(FaceOrder.begin()+first,FaceOrder.begin()+mid,FaceOrder.begin()+first+count,
      [&](int a,int b) {return centroids[a][axis]<centroids[b][axis];});

  int left=BuildNode(first,mid-first,centroids);
  int right=BuildNode(mid,first+count-mid,centroids);

  Nodes[inode].left=left,Nodes[inode].right=right;
  return inode;
}

//Moller-Trumbore ray/triangle intersection
bool cBVH::RayTriangle(int iface,const Vec3& x0,const Vec3& d,double tMin,double tMax,double& t) const {
  const cTriangle& f=Mesh->Faces[iface];
  Vec3 a=Mesh->Nodes[f.node[0]],e1=Mesh->Nodes[f.node[1]]-a,e2=Mesh->Nodes[f.node[2]]-a;
  Vec3 p=Cross(d,e2);
  double det=Dot(e1,p);

  if (det==0.0) return false;

  double inv=1.0/det;
  Vec3 s=x0-a;
  double u=Dot(s,p)*inv;
  if ((u<0.0)||(u>1.0)) return false;

  Vec3 q=Cross(s,e1);
  double v=Dot(d,q)*inv;
  if ((v<0.0)||(u+v>1.0)) return false;

  t=Dot(e2,q)*inv;
  return (t>tMin)&&(t<tMax);
}

static inline bool RayBox(const Vec3& bmin,const Vec3& bmax,const Vec3& x0,const Vec3& invd,double tMin,double tMax) {
  for (int i=0;i<3;i++) {
    double t0=(bmin[i]-x0[i])*invd[i],t1=(bmax[i]-x0[i])*invd[i];

    if (invd[i]<0.0) std::swap(t0,t1);
    tMin=std::fmax(tMin,t0);
    tMax=std::fmin(tMax,t1);
    if (tMax<tMin) return false;
  }

  return true;
}

int cBVH::FirstHit(const Vec3& x0,const Vec3& d,double tMin,double tMax,double& tHit,int ignore) const {
  Vec3 invd(1.0/d.x,1.0/d.y,1.0/d.z);
  int stack[128],nstack=0,res=-1;
  double best=tMax,t;

  stack[nstack++]=0;

  while (nstack>0) {
    const cNode& n=Nodes[stack[--nstack]];

    if (!RayBox(n.bmin,n.bmax,x0,invd,tMin,best)) continue;

    if (n.left<0) {
      for (int i=n.first;i<n.first+n.count;i++) {
        int iface=FaceOrder[i];

        if ((iface!=ignore)&&RayTriangle(iface,x0,d,tMin,best,t)) best=t,res=iface;
      }
    }
    else {
      stack[nstack++]=n.left;
      stack[nstack++]=n.right;
    }
  }

  if (res!=-1) tHit=best;
  return res;
}

bool cBVH::AnyHit(const Vec3& x0,const Vec3& d,double tMin,double tMax,int ignore) const {
  Vec3 invd(1.0/d.x,1.0/d.y,1.0/d.z);
  int stack[128],nstack=0;
  double t;

  stack[nstack++]=0;

  while (nstack>0) {
    const cNode& n=Nodes[stack[--nstack]];

    if (!RayBox(n.bmin,n.bmax,x0,invd,tMin,tMax)) continue;

    if (n.left<0) {
      for (int i=n.first;i<n.first+n.count;i++) {
        int iface=FaceOrder[i];

        if ((iface!=ignore)&&RayTriangle(iface,x0,d,tMin,tMax,t)) return true;
      }
    }
    else {
      stack[nstack++]=n.left;
      stack[nstack++]=n.right;
    }
  }

  return false;
}

int cBVH::FirstHitBruteForce(const Vec3& x0,const Vec3& d,double tMin,double tMax,double& tHit,int ignore) const {
  int res=-1;
  double best=tMax,t;

  for (int i=0;i<(int)Mesh->Faces.size();i++) if ((i!=ignore)&&RayTriangle(i,x0,d,tMin,best,t)) best=t,res=i;

  if (res!=-1) tHit=best;
  return res;
}

static inline double PointBoxDistanceSquared(const Vec3& x,const Vec3& bmin,const Vec3& bmax) {
  double d2=0.0;

  for (int i=0;i<3;i++) {
    double d=(x[i]<bmin[i]) ? bmin[i]-x[i] : ((x[i]>bmax[i]) ? x[i]-bmax[i] : 0.0);
    d2+=d*d;
  }

  return d2;
}

double cBVH::NearestDistance(const Vec3& x) const {
  if ((Mesh==nullptr)||Nodes.empty()) return std::numeric_limits<double>::infinity();

  int stack[128],nstack=0;
  double best2=std::numeric_limits<double>::infinity();
  stack[nstack++]=0;

  while (nstack>0) {
    const cNode& n=Nodes[stack[--nstack]];
    if (PointBoxDistanceSquared(x,n.bmin,n.bmax)>=best2) continue;

    if (n.left<0) {
      for (int i=n.first;i<n.first+n.count;i++) {
        int iface=FaceOrder[i];
        const cTriangle& f=Mesh->Faces[iface];
        double d=PointTriangleDistance(x,Mesh->Nodes[f.node[0]],Mesh->Nodes[f.node[1]],Mesh->Nodes[f.node[2]]);
        best2=std::min(best2,d*d);
      }
    }
    else {
      double dl=PointBoxDistanceSquared(x,Nodes[n.left].bmin,Nodes[n.left].bmax);
      double dr=PointBoxDistanceSquared(x,Nodes[n.right].bmin,Nodes[n.right].bmax);
      int near=(dl<=dr) ? n.left : n.right;
      int far=(dl<=dr) ? n.right : n.left;

      //LIFO order: visit the nearer child first so it tightens the pruning distance.
      if (std::max(dl,dr)<best2) stack[nstack++]=far;
      if (std::min(dl,dr)<best2) stack[nstack++]=near;
    }
  }

  return sqrt(best2);
}
