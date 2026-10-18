//bounding volume hierarchy for ray tracing through the nucleus triangulation
#ifndef BVH_H
#define BVH_H

#include <vector>

#include "Mesh.h"

class cBVH {
public:
  void Build(const cMesh* mesh);

  //nearest intersection of the ray x0+t*d, tMin<t<tMax; returns the face index or -1; tHit is set on success.
  //the face 'ignore' is excluded (use -1 to test all faces)
  int FirstHit(const Vec3& x0,const Vec3& d,double tMin,double tMax,double& tHit,int ignore=-1) const;

  //true when any face intersects the ray on tMin<t<tMax (face 'ignore' excluded)
  bool AnyHit(const Vec3& x0,const Vec3& d,double tMin,double tMax,int ignore=-1) const;

  //brute force version of FirstHit() used to verify the tree
  int FirstHitBruteForce(const Vec3& x0,const Vec3& d,double tMin,double tMax,double& tHit,int ignore=-1) const;

  //exact shortest distance from a point to the triangulated surface, with BVH branch-and-bound pruning
  double NearestDistance(const Vec3& x) const;

private:
  struct cNode {
    Vec3 bmin,bmax;
    int left,right;    //children (internal node) or -1
    int first,count;   //range in FaceOrder (leaf)
  };

  const cMesh* Mesh=nullptr;
  std::vector<cNode> Nodes;
  std::vector<int> FaceOrder;

  int BuildNode(int first,int count,std::vector<Vec3>& centroids);
  bool RayTriangle(int iface,const Vec3& x0,const Vec3& d,double tMin,double tMax,double& t) const;
};

#endif
