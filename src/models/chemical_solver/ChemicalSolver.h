/*
 * Dependency-light stiff chemistry integrator.
 *
 * Chemical state is represented by non-negative number densities y_i
 * [m^-3]. A model callback supplies a non-negative production rate P_i
 * [m^-3 s^-1] and loss frequency L_i [s^-1], so the integrated ODE is
 *
 *   dy_i/dt = P_i(y,t) - L_i(y,t) y_i .
 *
 * The production/loss split is the natural AMPS chemistry contract and lets
 * the solver distinguish an invalid negative rate from a physically negative
 * net derivative. The implementation uses fully implicit Euler substeps, a
 * damped Newton solve, step-doubling error control, and Richardson
 * extrapolation.  The accepted state is second-order accurate when the
 * extrapolated state remains in the non-negative physical domain. It is not
 * an equilibrium or reaction-network generator. Elemental, charge, and mass
 * conservation follow only when the callback describes a conservative
 * network; the solver never invents a correction that could conceal a bad
 * reaction table.
 */

#ifndef AMPS_MODELS_CHEMICAL_SOLVER_H
#define AMPS_MODELS_CHEMICAL_SOLVER_H

namespace ChemicalSolver {

enum class Status {
  Ok = 0,
  InvalidArgument,
  InvalidState,
  InvalidNetwork,
  DisconnectedNetwork,
  CallbackFailure,
  NegativeRate,
  NonFiniteRate,
  SingularJacobian,
  NewtonFailure,
  StepUnderflow,
  MaximumStepsExceeded
};

const char* StatusName(Status status);

/* Dependency-light reaction topology used by the strict preflight checker.
 * coefficient is the positive, dimensionless stoichiometric multiplier; a
 * species index is in [0,species_count). Repeated species on one side are
 * rejected rather than silently summed because duplicate parser entries most
 * often indicate a malformed input deck. The same species may occur on both
 * sides (a catalyst).
 *
 * ValidateNetwork intentionally requires one connected reaction graph. An
 * application with independent chemical subsystems must validate and solve
 * those subsystems separately. The callback-only Integrate API cannot infer
 * topology from a rate evaluated at one state, so validation is an explicit
 * precondition rather than an unreliable runtime guess. */
struct ReactionTerm {
  int species;
  double coefficient;
};

struct ReactionDefinition {
  const ReactionTerm* reactants;
  int reactant_count;
  const ReactionTerm* products;
  int product_count;
};

Status ValidateNetwork(int species_count, const ReactionDefinition* reactions,
                       int reaction_count);

/* The callback must fill all n entries of production and loss_frequency.
 * Returning a non-Ok status aborts the integration without converting the
 * model failure into a numerical result. State is immutable and time is in
 * seconds. user_data may be null when the model needs no parameters. */
using RateFunction = Status (*)(const double* state, int n, double time,
                                double* production,
                                double* loss_frequency, void* user_data);

struct Options {
  double relative_tolerance = 1.0e-7;
  double absolute_tolerance = 1.0e-30;  // density units [m^-3]
  double initial_step = 0.0;            // s; 0 selects the full interval
  double minimum_step = 0.0;            // s; 0 selects a roundoff limit
  double maximum_step = 0.0;            // s; 0 selects the full interval
  int maximum_steps = 200000;
  int maximum_newton_iterations = 16;
};

struct Diagnostics {
  int accepted_steps = 0;
  int rejected_steps = 0;
  int rate_evaluations = 0;
  int newton_iterations = 0;
  double final_time = 0.0;               // s
  double minimum_accepted_step = 0.0;    // s; zero for a zero-time call
  double maximum_accepted_step = 0.0;    // s
};

/* Advance state in place from start_time through start_time+duration.
 *
 * On success every returned density is finite and non-negative. On failure
 * state is restored to the last accepted time and diagnostics.final_time
 * identifies that time. Invalid initial state is rejected before invoking
 * the callback. The routine owns no static storage and is reentrant.
 */
Status Integrate(double* state, int n, double start_time, double duration,
                 RateFunction rates, void* user_data,
                 const Options& options = Options(),
                 Diagnostics* diagnostics = nullptr);

}  // namespace ChemicalSolver

#endif  // AMPS_MODELS_CHEMICAL_SOLVER_H
