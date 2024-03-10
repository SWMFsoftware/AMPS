//read tophological map derived from LRO/LOLA 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <map>
#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <string>
#include <algorithm>

#include "../../src/general/constants.PlanetaryData.h"
#include "../../src/pic/pic.h"

using namespace std;

string fNameLON="ldem_4_float_LON.txt";
string fNameLAT="ldem_4_float_LAT.txt";
string fNameALT="ldem_4_float.txt";

char fNameSphere[100]="sphere-coarce.mail";

vector<double> vLON,vLAT;
vector<vector<double> >vALT;

void ReadFileLON_LAT(string fName,vector<double>& vRes) {
  ifstream f(fName);
  string line;
  double  number;
  
  getline(f,line);
  std::istringstream iss(line);

  // Extract each float from the line
  while (iss >> number) { 
    vRes.push_back(number);
  }
 
  f.close();
}

void ReadFileALT(string fName,vector<vector<double> >& vRes) {
  ifstream f(fName);
  string line;
  double  number;

  for (int i=0;i<vRes.size();i++) {
    getline(f,line);
    std::istringstream iss(line);

    // Extract each float from the line
    while (iss >> number) {
      vRes[i].push_back(number);
    }
  }

  f.close();
}



int main () {

  PIC::InitMPI();

cout << "sdsfsfasdfsadfsd\n";
  

  //read the sphere 
  PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(fNameSphere,_MOON__RADIUS_);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("imported-sphere.dat");

  //read the map	
  ReadFileLON_LAT(fNameLON,vLON);
  ReadFileLON_LAT(fNameLAT,vLAT);


  vALT.resize(vLAT.size());
  ReadFileALT(fNameALT,vALT);

  return 0;
}
