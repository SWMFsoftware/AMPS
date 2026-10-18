#ifndef SHIELDSIM_RESPONSE_MATRIX_HH
#define SHIELDSIM_RESPONSE_MATRIX_HH

/*
 * Versioned Stage-3 response-matrix data model.
 *
 * A response cell is the mean raw transport response to one primary at an
 * exact species/energy/direction node.  Matrix assembly is allowed only when
 * every cell has the same response key.  Folding consumes explicit normalized
 * node weights; it never guesses an interpolation, silently extrapolates, or
 * folds a nonlinear ratio.  A spectrum adapter or quadrature driver must
 * construct and archive those weights separately.
 */

#include "StrictJson.hh"

#include <map>
#include <string>
#include <vector>

namespace shield_sim {

struct ResponseScore {
  std::string id;
  std::string unit;
  double meanPerPrimary = 0.0;
  double standardError = 0.0;
  double histories = 0.0;
};

struct ResponseCell {
  std::string incidentId;
  std::string particle;
  int atomicNumber = 0;
  int massNumber = 0;
  double totalKineticEnergyMeV = 0.0;
  double mu = 1.0;
  double phiRad = 0.0;
  std::map<std::string,ResponseScore> scores;
  // Covariance of the estimated per-primary means. Keys are canonical
  // lexicographic upper-triangle score-ID pairs.
  std::map<std::pair<std::string,std::string>,double> meanCovariance;
};

class ResponseMatrix {
public:
  static ResponseMatrix Assemble(const std::vector<std::string>& responseCellPaths);
  // Assemble a matrix containing only the explicitly requested primitive
  // scores.  A Geant4 response cell records many diagnostic scores and the
  // complete score covariance triangle can dominate the matrix size.  Stage-7
  // qualification needs only its registered dose and grid-monitor scores, so
  // this operation validates the complete cells first and then removes all
  // other scores and covariance pairs.  It fails if any requested score is
  // absent from any cell; selection can therefore never create a ragged
  // response matrix or silently substitute a different observable.
  static ResponseMatrix AssembleSelected(
      const std::vector<std::string>& responseCellPaths,
      const std::vector<std::string>& scoreIds);
  static ResponseMatrix Load(const std::string& manifestPath);

  // Writes the canonical manifest and a long-form CSV table beside it.  JSON
  // remains authoritative; CSV is a review/analysis derivative whose SHA-256
  // is recorded in the manifest.
  void Write(const std::string& manifestPath) const;

  // Fold exact response nodes using a shieldSim.response_weights.v1 document.
  // Optional source-weight covariance is propagated as r^T C r independently
  // for every primitive score.  Matrix-cell Monte Carlo variance is the sum of
  // squared weights times independent cell standard-error variances.
  Json Fold(const Json& weightsDocument) const;

  const std::string& ResponseKey() const { return responseKey_; }
  const std::string& MatrixHash() const { return matrixHash_; }
  const std::vector<ResponseCell>& Cells() const { return cells_; }
  Json CanonicalDocument() const;

private:
  std::string responseKey_;
  std::string matrixHash_;
  Json responseKeyDocument_;
  std::vector<ResponseCell> cells_;
};

} // namespace shield_sim

#endif
