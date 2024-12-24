#include <gtest/gtest.h>
#include "pic.h"
#include "meshAMRcutcell.h"


void Orbiter::link_geomentry_test() {}


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

        if (isIntersected) {
            // Verify the intersection point lies within the line segment bounds
            for (int j = 0; j < 3; ++j) {
                ASSERT_GE(xIntersection[j], xIntersectionTest[j] - TriangleIntersectionTest::tolerance);
                ASSERT_LE(xIntersection[j], xIntersectionTest[j] + TriangleIntersectionTest::tolerance);
            }

            // Additional verification can be added here based on triangle intersection logic
        }
    }

    int a = 2;
    int b = 3;
    EXPECT_EQ(a + b, 5); // Check if the sum is correct
}
