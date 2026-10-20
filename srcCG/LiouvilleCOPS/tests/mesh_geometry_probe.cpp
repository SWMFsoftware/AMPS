// Read-only, production-parser check for externally stored high-resolution
// meshes. This deliberately uses the same closed-surface finalization as the
// solver. It is not a gravity admission or a self-intersection test.
#include "Mesh.h"
#include <cstdio>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "Usage: mesh_geometry_probe surface.bdf\n");
    return 2;
  }
  cMesh mesh;
  mesh.ReadNastran(argv[1]);
  mesh.FinalizeGeometry(true, true);
  std::printf("Production mesh validation PASS: %zu nodes, %zu faces\n",
              mesh.Nodes.size(), mesh.Faces.size());
  return 0;
}
