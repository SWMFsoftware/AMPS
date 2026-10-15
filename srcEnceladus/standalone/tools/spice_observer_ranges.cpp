// Compute Cassini-to-tangent distances for UVIS rows whose time is NAIF ET.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "SpiceUsr.h"

static void fail(const char* where) {
  if (!failed_c()) return;
  SpiceChar message[1841]; getmsg_c("LONG", sizeof(message), message);
  std::fprintf(stderr, "SPICE ERROR during %s: %s\n", where, message);
  std::exit(2);
}

static std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> result; std::stringstream stream(line); std::string item;
  while (std::getline(stream, item, ',')) result.push_back(item);
  return result;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s META INPUT.csv OUTPUT.csv\n", argv[0]); return 2;
  }
  SpiceChar action[]="RETURN", device[]="NULL";
  erract_c("SET",0,action); errdev_c("SET",0,device); furnsh_c(argv[1]); fail("kernel loading");
  std::ifstream input(argv[2]); if (!input) { std::perror(argv[2]); return 2; }
  std::ofstream output(argv[3], std::ios::out|std::ios::trunc);
  if (!output) { std::perror(argv[3]); return 2; }
  std::string line; if (!std::getline(input,line)) return 2;
  const auto header=split(line); int it=-1, ix=-1, iy=-1, iz=-1;
  for (int i=0;i<(int)header.size();++i) {
    if (header[i]=="time_source_s") it=i;
    if (header[i]=="tangent_x_m") ix=i;
    if (header[i]=="tangent_y_m") iy=i;
    if (header[i]=="tangent_z_m") iz=i;
  }
  if (it<0||ix<0||iy<0||iz<0) { std::fprintf(stderr,"missing input column\n"); return 2; }
  output << "time_source_s,spacecraft_radius_m,tangent_radius_m,"
            "observer_coordinate_from_tangent_m,ray_ux_table,ray_uy_table,ray_uz_table\n";
  output.precision(15);
  while (std::getline(input,line)) {
    if (line.empty()) continue;
    const auto row=split(line);
    const double et=std::stod(row[it]);
    const double tx=std::stod(row[ix]),ty=std::stod(row[iy]),tz=std::stod(row[iz]);
    const double tangent_radius=std::sqrt(tx*tx+ty*ty+tz*tz);
    SpiceDouble state[6],lt; spkezr_c("CASSINI",et,"IAU_ENCELADUS","NONE","ENCELADUS",state,&lt);
    fail("Cassini state");
    const double spacecraft_radius=1000.0*std::sqrt(state[0]*state[0]+state[1]*state[1]+state[2]*state[2]);
    SpiceDouble sun[3]; spkpos_c("SUN",et,"IAU_ENCELADUS","LT+S","CASSINI",sun,&lt);
    fail("Sun direction"); vhat_c(sun,sun);
    const double observer_coordinate=1000.0*(state[0]*sun[0]+state[1]*sun[1]+state[2]*sun[2]);
    // The SwRI release negates the IAU-frame x and y coordinates before
    // constructing its UVIS tangent tracks; retain that documented convention.
    output << et << ',' << spacecraft_radius << ',' << tangent_radius << ','
           << observer_coordinate << ',' << -sun[0] << ',' << -sun[1] << ',' << sun[2] << '\n';
  }
  if (!input.eof() || !output) return 2;
  kclear_c(); return 0;
}
