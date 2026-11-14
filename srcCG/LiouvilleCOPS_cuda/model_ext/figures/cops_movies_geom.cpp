// ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cops_movies_geom.cpp (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cops_movies_geom.cpp)
// SHA256 of the run copy = SHA256 of this file without these header lines: 60324a2a795d7cf43faafce949569a894bc7c37be04b7c321e5fc47992a2b2fb
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: SPICE helper of the movies (CSPICE, read only): geom, shield, fov/fovinfo, scgeom; 2026-11-13 hgageom (HGA
//   dish centre, ROS_HGA rotation, gimbal, solar-array frames and gimbals; HGA held after its C-kernel ends, solar
//   arrays after theirs).
// RUN: built by jobs/cops_movies_prep.sbatch: g++ -O2 -std=c++17 -I<cspice>/include cops_movies_geom.cpp
//   <cspice>/lib/cspice.a -lm
// ==== end of header ====
// cops_movies_geom.cpp -- observation geometry for the COPS movies (user request 2026-11-12), read only.
//
//   cops_movies_geom geom <kernel list file> <gx> <gy> <gz>  < utc_list  > geometry.txt
//   cops_movies_geom shield <Rmax> <Hmax> <Rtube> <nBins> <nSamples>     > shield_table.txt
//
// geom: for every UTC line on stdin one output line
//   utc et gx gy gz vx vy vz sx sy sz rh ngx ngy ngz rgx rgy rgz
// exactly as LiouvilleCOPS src/Geometry.cpp cSpiceGeometry::Evaluate (unchanged reference code, copied expressions):
//   spacecraft state: spkezr_c(ROSETTA, et, 67P/C-G_CK, NONE, CHURYUMOV-GERASIMENKO); gauge position = spacecraft position
//   + pxform(ROS_SPACECRAFT -> 67P/C-G_CK) (GaugeOffsetSC = gx gy gz [m]); Sun: spkpos_c(SUN, et, 67P/C-G_CK, LT+S,
//   CHURYUMOV-GERASIMENKO), unit vector and distance rh [m]; gauge boresights: the IK boresights of ROS_ROSINA_COPS_G0
//   (NG, -226331: (0,-1,0)) and _G1 (RG, -226332: (0,0,1)) read with getfov_c and rotated with pxform(ROS_ROSINA_COPS ->
//   67P/C-G_CK). Positions in m, velocity in m/s (relative to the body-fixed frame).
// shield: the NG shield fraction table of src/Gauges.cpp cGaugeModel::BuildShieldTable, code copied verbatim (same RNG,
//   seed and draw order), so the movie's NG weight equals the fitted model's.
// Added 2026-11-13 (observable-part masks of the movies; the modes above are unchanged):
//   cops_movies_geom fovinfo <kernel list>                   > fovinfo.txt
//     per gauge (NG -226331, RG -226332) the getfov_c result in the instrument frame: id shape frame n, boresight,
//     n bound vectors (as LiouvilleCOPS src/Geometry.cpp cSpiceGeometry::ReadFOV, normalized)
//   cops_movies_geom fov <kernel list>  < utc_list           > fov.txt
//     per UTC line: utc et, then for NG and RG: the boresight and the n bound vectors rotated into 67P/C-G_CK with
//     pxform(<FOV frame> -> 67P/C-G_CK, et) and normalized (as Geometry.cpp Evaluate / Rotate), then the 3x3 rotation
//     matrix (row-major) itself
//   cops_movies_geom scgeom <kernel list> <SA hold UTC>  < utc_list  > sc.txt        (spacecraft shadow, movies only)
//     per UTC line: utc et et_sa held(0/1) | ROSETTA position in 67P/C-G_CK [m] (spkpos ROSETTA, NONE, observer
//     CHURYUMOV-GERASIMENKO) | rotation ROS_SPACECRAFT -> 67P/C-G_CK (9, row-major) | bus corners ROS_SPACECRAFT_C1..C8
//     (-226001..-226008) in ROS_SPACECRAFT relative to its origin ROS_SPACECRAFT (-226000) [m] | +Y wing ROS_SA+Y_C1..C4
//     (-226011..-226014) and -Y wing ROS_SA-Y_C1..C4 (-226021..-226024) in ROS_SPACECRAFT relative to -226000 at
//     et_sa = min(et, SA hold time) (the solar-array orientation is held after the SA C-kernel ends; user decision
//     2026-11-13) | check of the rotation: the bus corners at et and (when et <= the hold time) the wing corners at et
//     directly in 67P/C-G_CK relative to -226000 [m] (nan when held) | Sun direction in ROS_SPACECRAFT (spkpos SUN LT+S
//     from ROSETTA) | r_h [m]. NOTE (ROS_V38.TF): the origin of ROS_SPACECRAFT (-226000) is the launch-vehicle interface
//     point, the origin of ROSETTA (-226) the centre of gravity; ROS_STRUCT_V8.BSP gives the structure relative to
//     -226000 and no loaded kernel links -226000 to -226: the caller places -226000 relative to ROSETTA (option).
//   cops_movies_geom hgageom <kernel list> <extra CK> <HGA hold UTC> <SA hold UTC>  < utc_list  > hga.txt
//     (2026-11-13, detach-200 m movies: spacecraft model) per UTC line: utc et et_hga held_hga(0/1) et_sa | HGA dish rim
//     centre ROS_HGA (-226075) relative to ROS_SPACECRAFT (-226000) in ROS_SPACECRAFT [m] at et_hga | rotation ROS_HGA ->
//     ROS_SPACECRAFT (9, row-major) at et_hga | HGA gimbal ROS_HGA_GIMBAL (-226070) [m] | rotations ROS_SA+Y and ROS_SA-Y ->
//     ROS_SPACECRAFT (9 + 9) at et_sa | SA gimbals ROS_SA+Y_GIMBAL (-226010), ROS_SA-Y_GIMBAL (-226020) [m]; et_hga =
//     min(et, HGA hold) (the HGA C-kernel ends there; orientation held, user decision 2026-11-13), et_sa = min(et, SA hold)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include "SpiceUsr.h"

static void check(const char* where) {
  if (failed_c()) {
    SpiceChar msg[1841];
    getmsg_c("LONG", sizeof(msg), msg);
    fprintf(stderr, "SPICE error in %s: %s\n", where, msg);
    exit(EXIT_FAILURE);
  }
}

static int shield(double Rmax, double Hmax, double Rtube, int nBins, int nSamples) {
  // ---- copy of cGaugeModel::BuildShieldTable (src/Gauges.cpp) ----
  std::mt19937_64 rng(20160930);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<double> ShieldTable(nBins);
  for (int ib = 0; ib < nBins; ib++) {
    double mu = (double)ib / (nBins - 1);
    double sinth = sqrt(std::max(0.0, 1.0 - mu * mu));
    double sx = sinth, sz = -mu;
    int cnt = 0;
    for (int i = 0; i < nSamples; i++) {
      double r = sqrt(U(rng)) * Rmax, phi = 2.0 * M_PI * U(rng);
      double px = r * cos(phi), py = r * sin(phi), pz = U(rng) * Hmax;
      if (sz >= 0.0) {cnt++; continue;}
      double t = -pz / sz;
      double cx = px + sx * t, cy = py;
      if (cx * cx + cy * cy > Rtube * Rtube) cnt++;
    }
    ShieldTable[ib] = (double)cnt / nSamples;
  }
  // ---- end of copy ----
  for (int ib = 0; ib < nBins; ib++) printf("%.17g\n", ShieldTable[ib]);
  return 0;
}

static void load_kernels(const char* list) {
  FILE* k = fopen(list, "r");
  if (!k) {fprintf(stderr, "cannot open %s\n", list); exit(2);}
  char line[1024];
  while (fgets(line, sizeof line, k)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    furnsh_c(line);
    check(line);
  }
  fclose(k);
}

struct FovIK {
  int id; std::string shape, frame; double bs[3]; int n; double b[16][3];
};

static FovIK read_fov(int id) {
  // ---- as cSpiceGeometry::ReadFOV (src/Geometry.cpp): getfov_c, boresight and bounds normalized ----
  FovIK f; f.id = id;
  SpiceChar shape[64], frame[64];
  SpiceDouble bs[3], bounds[16][3];
  SpiceInt n;
  getfov_c(id, 16, sizeof(shape), sizeof(frame), shape, frame, bs, &n, bounds);
  check("getfov_c");
  f.shape = shape; f.frame = frame; f.n = (int)n;
  double l = sqrt(bs[0] * bs[0] + bs[1] * bs[1] + bs[2] * bs[2]);
  for (int i = 0; i < 3; i++) f.bs[i] = bs[i] / l;
  for (int j = 0; j < f.n; j++) {
    double m = sqrt(bounds[j][0] * bounds[j][0] + bounds[j][1] * bounds[j][1] + bounds[j][2] * bounds[j][2]);
    for (int i = 0; i < 3; i++) f.b[j][i] = bounds[j][i] / m;
  }
  return f;
}

static void rot_norm(double r[3][3], const double* v, double* o) {
  // ---- as Geometry.cpp Rotate: Normalized(R(v)) ----
  double t[3];
  for (int i = 0; i < 3; i++) t[i] = r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2];
  double l = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
  for (int i = 0; i < 3; i++) o[i] = t[i] / l;
}

static int fov_modes(int argc, char** argv) {
  erract_c("SET", 0, (SpiceChar*)"RETURN");
  errprt_c("SET", 0, (SpiceChar*)"NONE");
  load_kernels(argv[2]);
  FovIK g[2] = {read_fov(-226331), read_fov(-226332)};
  if (!strcmp(argv[1], "fovinfo")) {
    for (int k = 0; k < 2; k++) {
      printf("%d %s %s %d %.17g %.17g %.17g", g[k].id, g[k].shape.c_str(), g[k].frame.c_str(), g[k].n, g[k].bs[0],
             g[k].bs[1], g[k].bs[2]);
      for (int j = 0; j < g[k].n; j++) printf(" %.17g %.17g %.17g", g[k].b[j][0], g[k].b[j][1], g[k].b[j][2]);
      printf("\n");
    }
    return 0;
  }
  char line[1024];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    SpiceDouble et, r[3][3];
    utc2et_c(line, &et); check("utc2et_c");
    printf("%s %.6f", line, et);
    for (int k = 0; k < 2; k++) {
      pxform_c(g[k].frame.c_str(), "67P/C-G_CK", et, r); check("pxform_c(fov)");
      double o[3];
      rot_norm(r, g[k].bs, o);
      printf(" %d %.17g %.17g %.17g", g[k].n, o[0], o[1], o[2]);
      for (int j = 0; j < g[k].n; j++) {
        rot_norm(r, g[k].b[j], o);
        printf(" %.17g %.17g %.17g", o[0], o[1], o[2]);
      }
      for (int i = 0; i < 3; i++) printf(" %.17g %.17g %.17g", r[i][0], r[i][1], r[i][2]);
    }
    printf("\n");
  }
  return 0;
}

static void spk_rel(const char* target, double et, const char* frame, const char* obs, double* out) {
  SpiceDouble pos[3], lt;
  spkpos_c(target, et, frame, "NONE", obs, pos, &lt);
  check(target);
  for (int i = 0; i < 3; i++) out[i] = 1.0E3 * pos[i];
}

static int scgeom(int argc, char** argv) {
  erract_c("SET", 0, (SpiceChar*)"RETURN");
  errprt_c("SET", 0, (SpiceChar*)"NONE");
  load_kernels(argv[2]);
  SpiceDouble et_hold;
  utc2et_c(argv[3], &et_hold); check("utc2et_c(hold)");
  const char* bus[8] = {"-226001", "-226002", "-226003", "-226004", "-226005", "-226006", "-226007", "-226008"};
  const char* wing[8] = {"-226011", "-226012", "-226013", "-226014", "-226021", "-226022", "-226023", "-226024"};
  char line[1024];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    SpiceDouble et, r[3][3];
    utc2et_c(line, &et); check("utc2et_c");
    double et_sa = et <= et_hold ? et : et_hold;
    printf("%s %.6f %.6f %d", line, et, et_sa, et > et_hold ? 1 : 0);
    double x[3];
    spk_rel("ROSETTA", et, "67P/C-G_CK", "CHURYUMOV-GERASIMENKO", x);
    printf(" %.9f %.9f %.9f", x[0], x[1], x[2]);
    pxform_c("ROS_SPACECRAFT", "67P/C-G_CK", et, r); check("pxform_c(ROS_SPACECRAFT)");
    for (int i = 0; i < 3; i++) printf(" %.17g %.17g %.17g", r[i][0], r[i][1], r[i][2]);
    double v[3];
    for (int k = 0; k < 8; k++) {
      spk_rel(bus[k], et, "ROS_SPACECRAFT", "-226000", v);
      printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    }
    for (int k = 0; k < 8; k++) {
      spk_rel(wing[k], et_sa, "ROS_SPACECRAFT", "-226000", v);
      printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    }
    for (int k = 0; k < 8; k++) {
      spk_rel(bus[k], et, "67P/C-G_CK", "-226000", v);
      printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    }
    for (int k = 0; k < 8; k++) {
      if (et <= et_hold) {
        spk_rel(wing[k], et, "67P/C-G_CK", "-226000", v);
        printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
      } else {
        printf(" nan nan nan");
      }
    }
    SpiceDouble sp[3], lt;
    spkpos_c("SUN", et, "ROS_SPACECRAFT", "LT+S", "ROSETTA", sp, &lt); check("spkpos_c(SUN, ROS_SPACECRAFT)");
    double rh = sqrt(sp[0] * sp[0] + sp[1] * sp[1] + sp[2] * sp[2]);
    printf(" %.17g %.17g %.17g %.6f\n", sp[0] / rh, sp[1] / rh, sp[2] / rh, 1.0E3 * rh);
  }
  return 0;
}

static int hgageom(int argc, char** argv) {
  erract_c("SET", 0, (SpiceChar*)"RETURN");
  errprt_c("SET", 0, (SpiceChar*)"NONE");
  load_kernels(argv[2]);
  furnsh_c(argv[3]); check(argv[3]);
  SpiceDouble et_hh, et_sh;
  utc2et_c(argv[4], &et_hh); check("utc2et_c(hga hold)");
  utc2et_c(argv[5], &et_sh); check("utc2et_c(sa hold)");
  char line[1024];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    SpiceDouble et, r[3][3];
    utc2et_c(line, &et); check("utc2et_c");
    double et_h = et <= et_hh ? et : et_hh, et_s = et <= et_sh ? et : et_sh;
    printf("%s %.6f %.6f %d %.6f", line, et, et_h, et > et_hh ? 1 : 0, et_s);
    double v[3];
    spk_rel("-226075", et_h, "ROS_SPACECRAFT", "-226000", v);
    printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    pxform_c("ROS_HGA", "ROS_SPACECRAFT", et_h, r); check("pxform_c(ROS_HGA)");
    for (int i = 0; i < 3; i++) printf(" %.17g %.17g %.17g", r[i][0], r[i][1], r[i][2]);
    spk_rel("-226070", et_h, "ROS_SPACECRAFT", "-226000", v);
    printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    const char* saf[2] = {"ROS_SA+Y", "ROS_SA-Y"};
    for (int k = 0; k < 2; k++) {
      pxform_c(saf[k], "ROS_SPACECRAFT", et_s, r); check(saf[k]);
      for (int i = 0; i < 3; i++) printf(" %.17g %.17g %.17g", r[i][0], r[i][1], r[i][2]);
    }
    spk_rel("-226010", et_s, "ROS_SPACECRAFT", "-226000", v);
    printf(" %.9f %.9f %.9f", v[0], v[1], v[2]);
    spk_rel("-226020", et_s, "ROS_SPACECRAFT", "-226000", v);
    printf(" %.9f %.9f %.9f\n", v[0], v[1], v[2]);
  }
  return 0;
}

int main(int argc, char** argv) {
  if (argc >= 6 && !strcmp(argv[1], "hgageom"))
    return hgageom(argc, argv);
  if (argc >= 4 && !strcmp(argv[1], "scgeom"))
    return scgeom(argc, argv);
  if (argc >= 7 && !strcmp(argv[1], "shield"))
    return shield(atof(argv[2]), atof(argv[3]), atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
  if (argc >= 3 && (!strcmp(argv[1], "fov") || !strcmp(argv[1], "fovinfo")))
    return fov_modes(argc, argv);
  if (argc < 6 || strcmp(argv[1], "geom")) {
    fprintf(stderr, "usage: %s geom <kernel list> gx gy gz < utc | shield Rmax Hmax Rtube nBins nSamples | "
            "fovinfo <kernel list> | fov <kernel list> < utc\n", argv[0]);
    return 2;
  }
  erract_c("SET", 0, (SpiceChar*)"RETURN");
  errprt_c("SET", 0, (SpiceChar*)"NONE");
  FILE* k = fopen(argv[2], "r");
  if (!k) {fprintf(stderr, "cannot open %s\n", argv[2]); return 2;}
  char line[1024];
  while (fgets(line, sizeof line, k)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    furnsh_c(line);
    check(line);
  }
  fclose(k);
  double off[3] = {atof(argv[3]), atof(argv[4]), atof(argv[5])};
  SpiceChar shape[64], frame[64];
  SpiceDouble bs[2][3], bounds[16][3];
  SpiceInt n;
  int ids[2] = {-226331, -226332};
  std::string frames[2];
  for (int g = 0; g < 2; g++) {
    getfov_c(ids[g], 16, sizeof(shape), sizeof(frame), shape, frame, bs[g], &n, bounds);
    check("getfov_c");
    frames[g] = frame;
    double l = sqrt(bs[g][0] * bs[g][0] + bs[g][1] * bs[g][1] + bs[g][2] * bs[g][2]);
    for (int i = 0; i < 3; i++) bs[g][i] /= l;
  }
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\r\n")] = 0;
    if (!line[0]) continue;
    SpiceDouble et, st[6], lt, m[3][3], pos[3], r[3][3], bb[2][3];
    utc2et_c(line, &et); check("utc2et_c");
    spkezr_c("ROSETTA", et, "67P/C-G_CK", "NONE", "CHURYUMOV-GERASIMENKO", st, &lt); check("spkezr_c");
    pxform_c("ROS_SPACECRAFT", "67P/C-G_CK", et, m); check("pxform_c(spacecraft)");
    double x[3];
    for (int i = 0; i < 3; i++) x[i] = 1.0E3 * st[i] + m[i][0] * off[0] + m[i][1] * off[1] + m[i][2] * off[2];
    spkpos_c("SUN", et, "67P/C-G_CK", "LT+S", "CHURYUMOV-GERASIMENKO", pos, &lt); check("spkpos_c(SUN)");
    double xs[3] = {1.0E3 * pos[0], 1.0E3 * pos[1], 1.0E3 * pos[2]};
    double rh = sqrt(xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2]);
    for (int g = 0; g < 2; g++) {
      pxform_c(frames[g].c_str(), "67P/C-G_CK", et, r); check("pxform_c(gauge)");
      mxv_c(r, bs[g], bb[g]);
      double l = sqrt(bb[g][0] * bb[g][0] + bb[g][1] * bb[g][1] + bb[g][2] * bb[g][2]);
      for (int i = 0; i < 3; i++) bb[g][i] /= l;
    }
    printf("%s %.6f %.9f %.9f %.9f %.9f %.9f %.9f %.17g %.17g %.17g %.9f %.17g %.17g %.17g %.17g %.17g %.17g\n", line, et,
           x[0], x[1], x[2], 1.0E3 * st[3], 1.0E3 * st[4], 1.0E3 * st[5], xs[0] / rh, xs[1] / rh, xs[2] / rh, rh,
           bb[0][0], bb[0][1], bb[0][2], bb[1][0], bb[1][1], bb[1][2]);
  }
  return 0;
}
