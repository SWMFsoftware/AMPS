#include <gtest/gtest.h>
#include "pic.h"
#include "meshAMRcutcell.h"
#include "marching_cubes.h"


void Orbiter::link_geomentry_test() {
  using namespace MarchingCubes;

      // Define unit cube vertices
    std::vector<Point3D> cubeVertices = {
        Point3D(0,0,0), Point3D(1,0,0), Point3D(1,1,0), Point3D(0,1,0),  // bottom face
        Point3D(0,0,1), Point3D(1,0,1), Point3D(1,1,1), Point3D(0,1,1)   // top face
    };

    std::vector<bool> vertexAttributes = {
        false, false, true, false,   // bottom face vertices (0,1,2,3)
        false, true, true, false     // top face vertices (4,5,6,7)
    };

    std::vector<EdgeIntersection> intersections;

    // Create triangulator
    CubeTriangulator triangulator(cubeVertices, vertexAttributes);

    // Get surface triangulation
    auto triangles = triangulator.triangulateOutsideSurface();
    std::cout << "Number of triangles: " << triangles.size() << std::endl;
    for(const auto& tri : triangles) {
        std::cout << "Triangle: " << tri.v1 << " " << tri.v2 << " " << tri.v3 << std::endl;
    }

    // Get all intersection points
    const auto& points = triangulator.getIntersectionPoints();
    std::cout << "Number of intersection points: " << points.size() << std::endl;
    for(size_t i = 0; i < points.size(); i++) {
        std::cout << "Point " << i << ": " << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
    }

    // Get tetrahedralization
    auto tetrahedra = triangulator.createTetrahedralization();
    std::cout << "Number of tetrahedra: " << tetrahedra.size() << std::endl;
    for(const auto& tet : tetrahedra) {
        std::cout << "Tetrahedron: " << tet.v1 << " " << tet.v2 << " " << tet.v3 << " " << tet.v4 << std::endl;
    }

    triangulator.exportToTecplot("test.dat");
}


class Point3D {
public:
    double x, y, z;

    void GetRandom(double min, double max) {
      x=min+rnd()*(max-min);
      y=min+rnd()*(max-min);
      z=min+rnd()*(max-min);
    } 
};


TEST(SimpleTest, Addition) {
    int a = 2;
    int b = 3;
    EXPECT_EQ(a + b, 5); // Check if the sum is correct
}

// Define the test fixture class
class TriangleIntersectionTest : public ::testing::Test {
public:
    static constexpr double cubeMin = -5.0;
    static constexpr double cubeMax = 5.0;
    static constexpr double tolerance = 1.0E-6;

    Point3D vertices[3];

    void SetUp() override {
        // Generate a random triangle
        for (int i=0;i<3;i++) vertices[i].GetRandom(cubeMin, cubeMax);
    }

    void TearDown() override {
    }

    bool IntervalIntersection(double *x0, double *x1, double *xIntersection) {
        double direction[3];
        for (int i = 0; i < 3; ++i) {
            direction[i] = x1[i] - x0[i];
        }

        // Check if the line segment intersects the triangle plane
        Point3D edge0 = {vertices[1].x - vertices[0].x, vertices[1].y - vertices[0].y, vertices[1].z - vertices[0].z};
        Point3D edge1 = {vertices[2].x - vertices[0].x, vertices[2].y - vertices[0].y, vertices[2].z - vertices[0].z};

        // Compute the normal of the triangle plane
        Point3D normal = {
            edge0.y * edge1.z - edge0.z * edge1.y,
            edge0.z * edge1.x - edge0.x * edge1.z,
            edge0.x * edge1.y - edge0.y * edge1.x
        };

        double dotDirectionNormal = direction[0] * normal.x + direction[1] * normal.y + direction[2] * normal.z;

        // If the dot product is near zero, the line is parallel to the plane
        if (std::abs(dotDirectionNormal) < tolerance) {
            return false;
        }

        // Calculate the intersection point using the plane equation
        double d = normal.x * vertices[0].x + normal.y * vertices[0].y + normal.z * vertices[0].z;
        double t = (d - (normal.x * x0[0] + normal.y * x0[1] + normal.z * x0[2])) / dotDirectionNormal;

        // Ensure the intersection point is within the line segment bounds
        if (t < 0.0 || t > 1.0) {
            return false;
        }

        // Calculate the intersection point
        for (int i = 0; i < 3; ++i) {
            xIntersection[i] = x0[i] + t * direction[i];
        }

        // Check if the intersection point lies inside the triangle using barycentric coordinates
        Point3D intersect = {xIntersection[0], xIntersection[1], xIntersection[2]};
        Point3D v0 = vertices[0], v1 = vertices[1], v2 = vertices[2];

        Point3D edgeToIntersect = {intersect.x - v0.x, intersect.y - v0.y, intersect.z - v0.z};

        double dot00 = edge0.x * edge0.x + edge0.y * edge0.y + edge0.z * edge0.z;
        double dot01 = edge0.x * edge1.x + edge0.y * edge1.y + edge0.z * edge1.z;
        double dot02 = edge0.x * edgeToIntersect.x + edge0.y * edgeToIntersect.y + edge0.z * edgeToIntersect.z;
        double dot11 = edge1.x * edge1.x + edge1.y * edge1.y + edge1.z * edge1.z;
        double dot12 = edge1.x * edgeToIntersect.x + edge1.y * edgeToIntersect.y + edge1.z * edgeToIntersect.z;

        double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= 0) && (v >= 0) && (u + v <= 1);
    }    


bool RayIntersection(double *x0, double *l, double& IntersectionTime) {
    constexpr double tolerance = 1e-10;  // Numerical tolerance for floating-point comparisons
    
    // Calculate triangle edges from vertex[0]
    Point3D edge0 = {
        vertices[1].x - vertices[0].x,
        vertices[1].y - vertices[0].y,
        vertices[1].z - vertices[0].z
    };
    
    Point3D edge1 = {
        vertices[2].x - vertices[0].x,
        vertices[2].y - vertices[0].y,
        vertices[2].z - vertices[0].z
    };

    // Compute the normal of the triangle plane using cross product
    Point3D normal = {
        edge0.y * edge1.z - edge0.z * edge1.y,
        edge0.z * edge1.x - edge0.x * edge1.z,
        edge0.x * edge1.y - edge0.y * edge1.x
    };

    // Calculate dot product between ray direction and triangle normal
    double dotDirectionNormal = l[0] * normal.x + l[1] * normal.y + l[2] * normal.z;

    // If the dot product is near zero, the ray is parallel to the plane
    if (std::abs(dotDirectionNormal) < tolerance) {
        return false;
    }

    // Calculate the intersection time using the plane equation
    // Plane equation: normal · (point - vertex0) = 0
    // Ray equation: point = x0 + t * l
    // Solving for t: normal · (x0 + t * l - vertex0) = 0
    double d = normal.x * vertices[0].x + normal.y * vertices[0].y + normal.z * vertices[0].z;
    IntersectionTime = (d - (normal.x * x0[0] + normal.y * x0[1] + normal.z * x0[2])) / dotDirectionNormal;

    // For a ray, we only want intersections in the positive direction
    if (IntersectionTime < 0.0) {
        return false;
    }

    // Calculate the intersection point
    Point3D intersect = {
        x0[0] + IntersectionTime * l[0],
        x0[1] + IntersectionTime * l[1],
        x0[2] + IntersectionTime * l[2]
    };

    // Check if the intersection point lies inside the triangle using barycentric coordinates
    Point3D edgeToIntersect = {
        intersect.x - vertices[0].x,
        intersect.y - vertices[0].y,
        intersect.z - vertices[0].z
    };

    // Calculate dot products for barycentric coordinate computation
    double dot00 = edge0.x * edge0.x + edge0.y * edge0.y + edge0.z * edge0.z;
    double dot01 = edge0.x * edge1.x + edge0.y * edge1.y + edge0.z * edge1.z;
    double dot02 = edge0.x * edgeToIntersect.x + edge0.y * edgeToIntersect.y + edge0.z * edgeToIntersect.z;
    double dot11 = edge1.x * edge1.x + edge1.y * edge1.y + edge1.z * edge1.z;
    double dot12 = edge1.x * edgeToIntersect.x + edge1.y * edgeToIntersect.y + edge1.z * edgeToIntersect.z;

    // Compute barycentric coordinates
    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Return true if point lies inside the triangle
    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

};



// Then update the test to properly inherit from the fixture
TEST_F(TriangleIntersectionTest, TriangleSegmentIntersection) {
   CutCell::cTriangleFace tr;

   tr.SetFaceNodes(&this->vertices[0].x,&this->vertices[1].x,&this->vertices[2].x); 


    for (int i = 0; i < 1000000; ++i) {
        // Generate random interval endpoints
        Point3D x0,x1;
       
	x0.GetRandom(TriangleIntersectionTest::cubeMin, TriangleIntersectionTest::cubeMax);
        x1.GetRandom(TriangleIntersectionTest::cubeMin, TriangleIntersectionTest::cubeMax);

        // Prepare intersection point storage
        double xIntersection[3];

        // Execute the intersection test
        bool isIntersected = tr.IntervalIntersection(&x0.x, &x1.x, xIntersection, TriangleIntersectionTest::tolerance);



	double xIntersectionTest[3];
	bool flag=this->IntervalIntersection(&x0.x, &x1.x,xIntersectionTest);

	EXPECT_EQ(flag,isIntersected); 

	if (flag!=isIntersected) {
          double a=2;

	  isIntersected = tr.IntervalIntersection(&x0.x, &x1.x, xIntersection, TriangleIntersectionTest::tolerance);
        }

        if (isIntersected) {
            // Verify the intersection point lies within the line segment bounds
            for (int j = 0; j < 3; ++j) {
                ASSERT_GE(xIntersection[j], xIntersectionTest[j] - TriangleIntersectionTest::tolerance);
                ASSERT_LE(xIntersection[j], xIntersectionTest[j] + TriangleIntersectionTest::tolerance);
            }
        }

        double l[3],IntersectionTime,IntersectionTimeTheory;
	Vector3D::Substract(l,&x1.x,&x0.x);
	Vector3D::Normalize(l);

	flag=this->RayIntersection(&x0.x, l, IntersectionTimeTheory);
	isIntersected=tr.RayIntersection(&x0.x,l,IntersectionTime,PIC::Mesh::mesh->EPS);
	EXPECT_EQ(flag,isIntersected);

        if (flag!=isIntersected) {
          double a=2;

	 isIntersected = tr.IntervalIntersection(&x0.x, &x1.x, xIntersection, TriangleIntersectionTest::tolerance);
	  isIntersected=tr.RayIntersection(&x0.x,l,IntersectionTime,PIC::Mesh::mesh->EPS);

	}


	if (isIntersected) {
            // Verify the intersection point lies within the line segment bounds
            ASSERT_LE(fabs(IntersectionTimeTheory-IntersectionTime)/(IntersectionTimeTheory+IntersectionTime),PIC::Mesh::mesh->EPS);
        }
    }
}

// Then update the test to properly inherit from the fixture
class RayIntersectionTest : public TriangleIntersectionTest {}; 
TEST_F(RayIntersectionTest, TriangleSegmentIntersection) {
  double x0[3],l[3];
  int nIntersections,nIntereseactionAnalytic;

  const double R=1737000.0*5.0E-7;

  auto calculateIntersections = [&]() -> int {
    // Calculate coefficients for the quadratic equation
    double a = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
    double b = 2.0 * (x0[0] * l[0] + x0[1] * l[1] + x0[2] * l[2]);
    double c = x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2] - R * R;
    int res;

    // Discriminant of the quadratic equation
    double discriminant = b * b - 4 * a * c;

    // Determine the number of intersections based on the discriminant
    if (discriminant > 0) {
      // Two intersections
         res=2;
      } else if (discriminant == 0) {
         // Tangent to the sphere (one intersection)
         res=1;
      } else {
         // No intersection
         res= 0;
      }

      return res;
   };

   for (int i = 0; i < 1000000; ++i) {
     // Generate random interval endpoints

     do {
       for (int idim=0;idim<3;idim++) {
         x0[idim]=-1.3*R+rnd()*2.6*R;
       } 
     }
     while (Vector3D::DotProduct(x0,x0)<R*R);

     Vector3D::Distribution::Uniform(l);

     double xTarget[3];

     for (int idim=0;idim<3;idim++) xTarget[idim]=x0[idim]+3*R*l[idim];

     //get the analytic number of intersection 
     nIntereseactionAnalytic=calculateIntersections();


     //get numerical number of intersection  
     nIntersections=PIC::RayTracing::CountFaceIntersectionNumber(x0,xTarget,-1,false,NULL);

     EXPECT_EQ(nIntereseactionAnalytic,nIntersections);


     //get numerical number of intersection
     nIntersections=PIC::RayTracing::CountFaceIntersectionNumber(x0,xTarget,-1,false,NULL);
   }
}


