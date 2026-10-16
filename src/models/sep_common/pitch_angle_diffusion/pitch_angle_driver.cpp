// ============================================================================
// Command driver for the Python reference checks (reference_checks.py).
//
// Reads one command per line on stdin and prints one result line on stdout,
// with values at 17 significant digits, so an independent mpmath
// implementation and the companion bundles' JSON fixtures can be compared
// with the library without transcribing numbers by hand.  Commands:
//
//   kinematics MASS_KG CHARGE_C MOMENTUM      -> gamma speed rigidity_V T_J
//   qform Q H                                 -> I_closed I_quadrature
//   phi EPS                                   -> phi_closed phi_quadrature
//   shape MODEL D0 MU Q H EPS                 -> D dD/dmu (nan if absent)
//   lambda MODEL Q H EPS LAMBDA_M MASS CHARGE P
//          -> lambda_Eq23 D0 nu_max  (constant law: a = b = 0, refs 1, r = 1 AU)
//   transmissibility MODEL Q H EPS D0 E0 E1 E2  -> T between the centres of
//          the cells [E0,E1] and [E1,E2] on a grid -1..E0..E2..1
//
// MODEL is an inventory ID (PA-ISO, PA-QFORM, PA-EPS, PA-EPREM).  An error is
// printed as "error <StatusCode> <detail>".
// ============================================================================

#include "pitch_angle_diffusion.h"
#include "pitch_angle_kernels.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace PA = SEP::PitchAngleDiffusion;
namespace K = SEP::PitchAngleDiffusion::Kernels;

namespace {

void PrintValues(const std::vector<double>& values) {
  for (std::size_t i = 0; i < values.size(); ++i)
    std::printf(i ? " %.17g" : "%.17g", values[i]);
  std::printf("\n");
}

void PrintError(const PA::Status& status) {
  std::printf("error %s %s\n", PA::Name(status.code), status.detail.c_str());
}

}  // namespace

int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string command;
    if (!(in >> command)) continue;
    if (command == "kinematics") {
      PA::ParticleState p;
      in >> p.massKg >> p.chargeC >> p.momentumKgMPerS;
      PA::Kinematics k;
      const PA::Status s = PA::ComputeKinematics(p, &k);
      if (!s.ok()) { PrintError(s); continue; }
      PrintValues({k.gamma, k.speedMPerS, k.rigidityV, k.kineticEnergyJ});
    } else if (command == "qform") {
      double q = 0, h = 0, closed = 0, quadrature = 0;
      in >> q >> h;
      std::string method;
      PA::Status s = PA::QFormIntegralClosedForm(q, h, &closed, &method);
      if (!s.ok()) { PrintError(s); continue; }
      s = PA::QFormIntegralQuadrature(q, h, &quadrature);
      if (!s.ok()) { PrintError(s); continue; }
      PrintValues({closed, quadrature});
    } else if (command == "phi") {
      double eps = 0, closed = 0, quadrature = 0;
      in >> eps;
      PA::Status s = PA::EpsilonPhiClosedForm(eps, &closed);
      if (!s.ok()) { PrintError(s); continue; }
      s = PA::EpsilonPhiQuadrature(eps, &quadrature);
      if (!s.ok()) { PrintError(s); continue; }
      PrintValues({closed, quadrature});
    } else if (command == "shape" || command == "lambda" ||
               command == "transmissibility") {
      std::string id;
      in >> id;
      PA::ModelId model;
      const PA::Status parsed = PA::ParseModelId(id, &model);
      if (!parsed.ok()) { PrintError(parsed); continue; }
      if (command == "shape") {
        double d0 = 0, mu = 0;
        PA::ShapeParameters shape;
        in >> d0 >> mu >> shape.q >> shape.gapH >> shape.epsilon;
        const PA::CoefficientResult r = PA::EvaluateShape(mu, d0, model, shape);
        if (!r.status.ok()) { PrintError(r.status); continue; }
        PrintValues({r.dMuMuPerS, r.dDmuMuDmuPerS
                                      ? *r.dDmuMuDmuPerS
                                      : std::numeric_limits<double>::quiet_NaN()});
      } else if (command == "lambda") {
        PA::ModelConfiguration c;
        c.model = model;
        double lambda = 0;
        PA::ParticleState p;
        in >> c.shape.q >> c.shape.gapH >> c.shape.epsilon >> lambda >> p.massKg >>
            p.chargeC >> p.momentumKgMPerS;
        c.lambdaLaw.lambdaReferenceM = lambda;
        c.lambdaLaw.referenceRigidityV = 1.0;
        c.lambdaLaw.referenceRadiusM = 1.0;
        c.step.rule = PA::StepRule::RateFraction;
        c.step.rateFraction = 0.001;
        PA::LocalState local;
        local.heliocentricPositionM = {{149597870700.0, 0.0, 0.0}};
        const PA::CoefficientResult r = PA::Evaluate(0.25, p, local, c);
        if (!r.status.ok()) { PrintError(r.status); continue; }
        PrintValues({r.lambdaParallelUnfocusedM, r.amplitudeD0PerS,
                     r.maximumScatteringFrequencyPerS});
      } else {
        double d0 = 0, e0 = 0, e1 = 0, e2 = 0;
        PA::ShapeParameters shape;
        in >> shape.q >> shape.gapH >> shape.epsilon >> d0 >> e0 >> e1 >> e2;
        std::vector<double> edges = {-1.0};
        if (e0 > -1.0) edges.push_back(e0);
        edges.push_back(e1);
        if (e2 < 1.0) edges.push_back(e2);
        edges.push_back(1.0);
        K::FvGenerator g;
        const PA::Status s = K::BuildFvGenerator(model, shape, d0, edges, &g);
        if (!s.ok()) { PrintError(s); continue; }
        // The face between the cells [e0,e1] and [e1,e2].
        const std::size_t face = e0 > -1.0 ? 1 : 0;
        PrintValues({g.conductance[face]});
      }
    } else {
      std::printf("error InvalidInput unknown command %s\n", command.c_str());
    }
    std::fflush(stdout);
  }
  return 0;
}
